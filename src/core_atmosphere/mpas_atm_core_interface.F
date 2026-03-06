! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_core_interface

   use mpas_attlist
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units

   contains


   !***********************************************************************
   !
   !  routine atm_setup_core
   !
   !> \brief   Atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use atm_core, only : atm_core_init, atm_core_run, atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => atm_core_init
      core % core_run => atm_core_run
      core % core_finalize => atm_core_finalize
      core % define_packages => atm_define_packages
      core % setup_packages => atm_setup_packages
      core % setup_decompositions => atm_setup_decompositions
      core % setup_clock => atm_setup_clock
      core % setup_log => atm_setup_log
      core % get_mesh_stream => atm_get_mesh_stream
      core % setup_immutable_streams => atm_setup_immutable_streams
      core % setup_derived_dimensions => atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => atm_setup_decomposed_dimensions
      core % setup_block => atm_setup_block
      core % setup_namelist => atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "core_variables.inc"

   end subroutine atm_setup_core


   !***********************************************************************
   !
   !  routine atm_setup_domain
   !
   !> \brief   Atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "domain_variables.inc"

   end subroutine atm_setup_domain


   !***********************************************************************
   !
   !  function atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function atm_setup_packages(configs, streamInfo, packages, iocontext) result(ierr)

      use mpas_dmpar
      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

#ifdef DO_PHYSICS
      use mpas_atmphys_control
      use mpas_atmphys_packages
#endif

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      logical, pointer :: iauActive
      character(len=StrKIND), pointer :: config_iau_option
      logical, pointer :: limited_areaActive
      logical, pointer :: config_apply_lbcs
      logical, pointer :: config_jedi_da, jedi_daActive
      logical, pointer :: no_invariant_streamActive
      logical, pointer :: ugwp_orog_streamActive
      logical, pointer :: ugwp_ngw_streamActive
      logical, pointer :: ugwp_diags_streamActive
      character(len=StrKIND), pointer :: config_gwdo_scheme
      logical, pointer :: config_ngw_scheme
      logical, pointer :: config_ugwp_diags
      character(len=StrKIND) :: attvalue
      integer :: local_ierr


      ierr = atm_setup_packages_when(configs, packages)
      if (ierr /= 0) then
         return
      end if

      !
      ! Incremental analysis update
      !
      nullify(config_iau_option)
      call mpas_pool_get_config(configs, 'config_IAU_option', config_iau_option)

      nullify(iauActive)
      call mpas_pool_get_package(packages, 'iauActive', iauActive)

      if (trim(config_iau_option) /= 'off') then
         iauActive = .true.
      else
         iauActive = .false.
      end if

      !
      ! Limited-area
      !
      nullify(config_apply_lbcs)
      call mpas_pool_get_config(configs, 'config_apply_lbcs', config_apply_lbcs)

      nullify(limited_areaActive)
      call mpas_pool_get_package(packages, 'limited_areaActive', limited_areaActive)

      if (config_apply_lbcs) then
         limited_areaActive = .true.
      else
         limited_areaActive = .false.
      end if

      !
      ! JEDI data assimilation
      !
      nullify(config_jedi_da)
      call mpas_pool_get_config(configs, 'config_jedi_da', config_jedi_da)

      nullify(jedi_daActive)
      call mpas_pool_get_package(packages, 'jedi_daActive', jedi_daActive)

      if (associated(config_jedi_da) .and. associated(jedi_daActive)) then
         jedi_daActive = config_jedi_da
      else
         ierr = ierr + 1
         call mpas_log_write('Package setup failed for ''jedi_da''. '// &
              'Either ''jedi_da'' is not a package, or ''config_jedi_da'' is not a namelist option.', &
              messageType=MPAS_LOG_ERR)
      end if

      !
      ! Separate time-invariant stream
      !
      nullify(no_invariant_streamActive)
      call mpas_pool_get_package(packages, 'no_invariant_streamActive', no_invariant_streamActive)

      if (associated(no_invariant_streamActive)) then
         no_invariant_streamActive = .true.
         if (streamInfo % query('invariant', attname='input_interval', attvalue=attvalue)) then
            if (trim(attvalue) == 'initial_only') then
               no_invariant_streamActive = .false.
            end if
         end if
      else
         ierr = ierr + 1
         call mpas_log_write("Package setup failed for 'no_invariant_stream'. 'no_invariant_stream' is not a package.", &
                             messageType=MPAS_LOG_ERR)
      end if

#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at
      !least one physics parameterization is called (using the logical moist_physics):
      call physics_namelist_check(configs)

      local_ierr = atmphys_setup_packages(configs, packages, iocontext)
      if (local_ierr /= 0) then
         ierr = ierr + 1
         call mpas_log_write('Package setup failed for atmphys in core_atmosphere', messageType=MPAS_LOG_ERR)
      end if

      !
      ! Optional gravity wave drag parameterization streams
      !
      call mpas_pool_get_config(configs, 'config_gwdo_scheme', config_gwdo_scheme)
      nullify(ugwp_orog_streamActive)
      call mpas_pool_get_package(packages, 'ugwp_orog_streamActive', ugwp_orog_streamActive)
      if ( associated(config_gwdo_scheme) .and. associated(ugwp_orog_streamActive) ) then
         if (trim(config_gwdo_scheme) == "bl_ugwp_gwdo") then
            ugwp_orog_streamActive = .true.
         else
            ugwp_orog_streamActive = .false.
         endif
      else
         ierr = ierr + 1
         call mpas_log_write("Package setup failed for 'ugwp_orog_stream'. 'ugwp_orog_stream' is not a package.", &
                             messageType=MPAS_LOG_ERR)
      end if

      call mpas_pool_get_config(configs, 'config_ngw_scheme', config_ngw_scheme)
      nullify(ugwp_ngw_streamActive)
      call mpas_pool_get_package(packages, 'ugwp_ngw_streamActive', ugwp_ngw_streamActive)
      if ( associated(config_ngw_scheme) .and. associated(ugwp_ngw_streamActive) ) then
            ugwp_ngw_streamActive = config_ngw_scheme
      else
         ierr = ierr + 1
         call mpas_log_write("Package setup failed for 'ugwp_ngw_stream'. 'ugwp_ngw_stream' is not a package.", &
                             messageType=MPAS_LOG_ERR)
      end if

      call mpas_pool_get_config(configs, 'config_ugwp_diags', config_ugwp_diags)
      nullify(ugwp_diags_streamActive)
      call mpas_pool_get_package(packages, 'ugwp_diags_streamActive', ugwp_diags_streamActive)
      if ( associated(config_ugwp_diags) .and. associated(ugwp_diags_streamActive) ) then
            ugwp_diags_streamActive = config_ugwp_diags
      else
         ierr = ierr + 1
         call mpas_log_write("Package setup failed for 'ugwp_diags_stream'. 'ugwp_diags_stream' is not a package.", &
                             messageType=MPAS_LOG_ERR)
      end if

#endif

   end function atm_setup_packages


   !***********************************************************************
   !
   !  function atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_derived_types, only : MPAS_Clock_type, mpas_pool_type
      use atm_core, only : atm_simulation_clock_init

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      ierr = 0

      call atm_simulation_clock_init(core_clock, configs, ierr)

   end function atm_setup_clock


   !***********************************************************************
   !
   !  function atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function atm_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types, only : mpas_log_type, domain_type
      use mpas_log, only : mpas_log_init, mpas_log_open
      use mpas_framework, only : mpas_framework_report_settings

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo    !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain          !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call mpas_log_write('')
      call mpas_log_write('MPAS-Atmosphere Version '//trim(domain % core % modelVersion))
      call mpas_log_write('')

      call mpas_framework_report_settings(domain)

   end function atm_setup_log!}}}


   !***********************************************************************
   !
   !  function atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function atm_get_mesh_stream(configs, streamInfo, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type, MPAS_streamInfo_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (MPAS_streamInfo_type), intent(inout) :: streamInfo
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart
      character(len=StrKIND) :: attvalue

      ierr = 0

      !
      ! If the 'invariant' stream is defined in the streams XML file with an
      ! input_interval of 'initial_only', then use the 'invariant' stream to
      ! get mesh information
      !
      if (streamInfo % query('invariant', attname='input_interval', attvalue=attvalue)) then
         if (trim(attvalue) == 'initial_only') then
            write(stream,'(a)') 'invariant'
            return
         end if
      end if

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         call mpas_log_write('config_do_restart was not found when defining mesh stream.', messageType=MPAS_LOG_ERR)
         ierr = 1
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end function atm_get_mesh_stream


   !***********************************************************************
   !
   !  function atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function atm_setup_decompositions


   !***********************************************************************
   !
   !  function atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be
   !>  an atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type
      use mpas_pool_routines, only : mpas_pool_get_config
      use mpas_log, only : mpas_log_write

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      integer, pointer :: cam_pcnst
      integer :: err_level

      ierr = 0

      call atm_generate_structs(block, block % structs, block % dimensions, block % packages)

      !
      ! When MPAS-A is operating as a dycore in CAM, the scalars/scalars_tend var_arrays are
      ! allocated by the call to atm_allocate_scalars, below. The CAM-MPAS interface layer
      ! should have added a config, cam_pcnst, to the configs pool to indicate how many scalars
      ! are to be allocated.
      !
      nullify(cam_pcnst)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, 'cam_pcnst', cam_pcnst)
      call mpas_pool_set_error_level(err_level)
      if (associated(cam_pcnst)) then
         call mpas_log_write('')
         call mpas_log_write('** Config ''cam_pcnst'' is defined with a value of $i', intArgs=[cam_pcnst])
         call mpas_log_write('   Scalars will be allocated separately from Registry-defined variables')
         call mpas_log_write('')
         ierr = atm_allocate_scalars(block, cam_pcnst)
      end if

#ifdef MPAS_USE_MUSICA
      ierr = atm_extend_scalars_for_chemistry(block)
      if (ierr /= 0) then
         call mpas_log_write('** Error extending scalars for chemistry', messageType=MPAS_LOG_ERR)
         return
      end if
#endif

   end function atm_setup_block


   !***********************************************************************
   !
   !  function atm_allocate_scalars
   !
   !> \brief   Allocate scalars and scalars_tend var_arrays
   !> \author  Michael G. Duda
   !> \date    20 May 2020
   !> \details 
   !>  When MPAS-A is operating as a dycore for CAM, the scalars and
   !>  scalars_tend var_arrays are allocated separately from other Registry-
   !>  defined variables, since the set of scalars to be handled by the dycore
   !>  is not known until runtime. This routine allocates these var_arrays,
   !>  but it does not define which constituent is at which position in
   !>  var_arrays; this is defined later in the CAM-MPAS interface layer.
   !
   !-----------------------------------------------------------------------
   function atm_allocate_scalars(block, num_scalars) result(ierr)

      use mpas_derived_types, only : block_type

      use mpas_derived_types, only : mpas_pool_type, field3dReal, MPAS_LOG_ERR
      use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_add_dimension, mpas_pool_add_field
      use mpas_log, only : mpas_log_write

      implicit none

      ! Arguments
      type (block_type), pointer :: block
      integer, intent(in) :: num_scalars

      ! Return value
      integer :: ierr

      ! Local variables
      integer :: i, j, timeLevs
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tendPool
      type (field3dReal), dimension(:), pointer :: scalarsField


      ierr = 0

      !
      ! Allocate scalars var_array
      !
      nullify(statePool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)

      if (.not. associated(statePool)) then
         call mpas_log_write('No pool named ''state'' was found in atm_allocate_scalars', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      timeLevs = 2

      call mpas_pool_add_dimension(statePool, 'num_scalars', num_scalars)

      allocate(scalarsField(timeLevs))

      do i = 1, timeLevs
         scalarsField(i) % block => block
         scalarsField(i) % fieldName = 'scalars'
         scalarsField(i) % dimNames(1) = 'num_scalars'
         scalarsField(i) % dimNames(2) = 'nVertLevels'
         scalarsField(i) % dimNames(3) = 'nCells'
         scalarsField(i) % defaultValue = 0.0
         scalarsField(i) % missingValue = -1.0
         scalarsField(i) % isDecomposed = .true.
         scalarsField(i) % hasTimeDimension = .true.
         scalarsField(i) % isActive = .true.
         scalarsField(i) % isVarArray = .true.
         scalarsField(i) % isPersistent = .true.

         allocate(scalarsField(i) % constituentNames(num_scalars))

         allocate(scalarsField(i) % attLists(num_scalars))
         do j = 1, num_scalars
            allocate(scalarsField(i) % attLists(j) % attList)
         end do

      end do

      call mpas_pool_add_field(statePool, 'scalars', scalarsField)
      call mpas_pool_add_field(block % allFields, 'scalars', scalarsField)


      !
      ! Allocate scalars_tend var_array
      !
      nullify(tendPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      if (.not. associated(tendPool)) then
         call mpas_log_write('No pool named ''tend'' was found in atm_allocate_scalars', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      timeLevs = 1

      call mpas_pool_add_dimension(tendPool, 'num_scalars_tend', num_scalars)

      allocate(scalarsField(timeLevs))

      do i = 1, timeLevs
         scalarsField(i) % block => block
         scalarsField(i) % fieldName = 'scalars_tend'
         scalarsField(i) % dimNames(1) = 'num_scalars_tend'
         scalarsField(i) % dimNames(2) = 'nVertLevels'
         scalarsField(i) % dimNames(3) = 'nCells'
         scalarsField(i) % defaultValue = 0.0
         scalarsField(i) % missingValue = -1.0
         scalarsField(i) % isDecomposed = .true.
         scalarsField(i) % hasTimeDimension = .true.
         scalarsField(i) % isActive = .true.
         scalarsField(i) % isVarArray = .true.
         scalarsField(i) % isPersistent = .true.

         allocate(scalarsField(i) % constituentNames(num_scalars))

         allocate(scalarsField(i) % attLists(num_scalars))
         do j = 1, num_scalars
            allocate(scalarsField(i) % attLists(j) % attList)
         end do

      end do

      call mpas_pool_add_field(tendPool, 'scalars_tend', scalarsField)
      call mpas_pool_add_field(block % allFields, 'scalars_tend', scalarsField)

   end function atm_allocate_scalars


#ifdef MPAS_USE_MUSICA
   !-----------------------------------------------------------------------
   !
   !  function atm_extend_scalars_for_chemistry
   !
   !> \brief   Extend scalars/scalars_tend with chemistry tracers from MICM
   !> \author  CheMPAS-A Developers
   !> \date    March 2026
   !> \details
   !>  Queries MICM config for species names, then extends the existing
   !>  scalars and scalars_tend var_arrays by appending chemistry tracers.
   !>  Called from atm_setup_block after atm_generate_structs.
   !
   !-----------------------------------------------------------------------
   function atm_extend_scalars_for_chemistry(block) result(ierr)

      use mpas_atm_chemistry, only : chemistry_query_species
      use mpas_kind_types,    only : StrKIND, RKIND
      use mpas_log,           only : mpas_log_write

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ! Local variables
      character(len=StrKIND), pointer :: config_micm_file
      integer :: num_chem_species
      character(len=64), allocatable :: chem_names(:)
      integer :: error_code
      character(len=:), allocatable :: error_message

      type (mpas_pool_type), pointer :: statePool, tendPool
      type (field3dReal), pointer :: fieldPtr

      real (kind=RKIND), dimension(:,:,:), pointer :: old_array, new_array
      character(len=StrKIND), dimension(:), pointer :: old_names, new_names
      type (att_lists_type), dimension(:), pointer :: old_atts, new_atts

      integer :: old_count, new_count, t, i, j
      integer :: nVertLevels, nCells
      integer, pointer :: nVertLevels_ptr, nCells_ptr
      integer :: err_level

      ierr = 0

      ! Get MICM config path from namelist
      nullify(config_micm_file)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, 'config_micm_file', config_micm_file)
      call mpas_pool_set_error_level(err_level)

      if (.not. associated(config_micm_file)) then
         call mpas_log_write('[Chemistry] No config_micm_file found; skipping dynamic tracer allocation')
         return
      end if

      if (len_trim(config_micm_file) == 0) then
         call mpas_log_write('[Chemistry] config_micm_file is empty; skipping dynamic tracer allocation')
         return
      end if

      ! Query MICM for species names
      call chemistry_query_species(config_micm_file, num_chem_species, chem_names, error_code, error_message)
      if (error_code /= 0) then
         call mpas_log_write('[Chemistry] Error querying MICM species: ' // trim(error_message), &
                             messageType=MPAS_LOG_ERR)
         ierr = error_code
         return
      end if

      if (num_chem_species == 0) then
         call mpas_log_write('[Chemistry] MICM config has no species; no tracers to add')
         return
      end if

      call mpas_log_write('[Chemistry] Extending scalars with $i chemistry tracers', &
                          intArgs=[num_chem_species])

      ! Get dimensions
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells_ptr)
      nVertLevels = nVertLevels_ptr
      nCells = nCells_ptr

      ! --- Extend scalars in state pool ---
      nullify(statePool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      if (.not. associated(statePool)) then
         call mpas_log_write('[Chemistry] No state pool found', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      ! Get time level 1 to determine current size
      nullify(fieldPtr)
      call mpas_pool_get_field(statePool, 'scalars', fieldPtr, 1)
      if (.not. associated(fieldPtr)) then
         call mpas_log_write('[Chemistry] No scalars field found in state pool', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      old_count = size(fieldPtr % array, 1)
      new_count = old_count + num_chem_species

      call mpas_log_write('[Chemistry] scalars: old_count=$i new_count=$i', &
                          intArgs=[old_count, new_count])

      ! Update num_scalars dimension
      call mpas_pool_add_dimension(statePool, 'num_scalars', new_count)

      ! Extend each time level (scalars has 2 time levels)
      do t = 1, 2
         nullify(fieldPtr)
         call mpas_pool_get_field(statePool, 'scalars', fieldPtr, t)

         call extend_field3d(fieldPtr, old_count, new_count, num_chem_species, &
                             chem_names, nVertLevels, nCells)
      end do

      ! Add index dimensions for each chemistry tracer
      do i = 1, num_chem_species
         call mpas_pool_add_dimension(statePool, 'index_' // trim(chem_names(i)), old_count + i)
         call mpas_log_write('[Chemistry] Added index_' // trim(chem_names(i)) // ' = $i', &
                             intArgs=[old_count + i])
      end do

      ! --- Extend scalars_tend in tend pool ---
      nullify(tendPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      if (.not. associated(tendPool)) then
         call mpas_log_write('[Chemistry] No tend pool found', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      nullify(fieldPtr)
      call mpas_pool_get_field(tendPool, 'scalars_tend', fieldPtr, 1)
      if (.not. associated(fieldPtr)) then
         call mpas_log_write('[Chemistry] No scalars_tend field found', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      old_count = size(fieldPtr % array, 1)
      new_count = old_count + num_chem_species

      call mpas_pool_add_dimension(tendPool, 'num_scalars_tend', new_count)

      ! scalars_tend has 1 time level
      call extend_field3d(fieldPtr, old_count, new_count, num_chem_species, &
                          chem_names, nVertLevels, nCells)

      ! Add index dimensions for tend pool too
      do i = 1, num_chem_species
         call mpas_pool_add_dimension(tendPool, 'index_' // trim(chem_names(i)), old_count + i)
      end do

      if (allocated(chem_names)) deallocate(chem_names)

      call mpas_log_write('[Chemistry] Dynamic tracer allocation complete')

      contains

      subroutine extend_field3d(fld, old_n, new_n, n_chem, names, nVert, nCell)

         type (field3dReal), pointer, intent(inout) :: fld
         integer, intent(in) :: old_n, new_n, n_chem, nVert, nCell
         character(len=64), intent(in) :: names(:)

         real (kind=RKIND), dimension(:,:,:), pointer :: old_arr, new_arr
         character(len=StrKIND), dimension(:), pointer :: onames, nnames
         type (att_lists_type), dimension(:), pointer :: oatts, natts
         integer :: ii, jj

         ! Extend data array
         old_arr => fld % array
         allocate(new_arr(new_n, nVert, nCell))
         new_arr = 0.0_RKIND
         new_arr(1:old_n, :, :) = old_arr(1:old_n, :, :)
         fld % array => new_arr
         deallocate(old_arr)
         fld % dimSizes(1) = new_n

         ! Extend constituentNames
         onames => fld % constituentNames
         allocate(nnames(new_n))
         nnames(1:old_n) = onames(1:old_n)
         do ii = 1, n_chem
            nnames(old_n + ii) = trim(names(ii))
         end do
         fld % constituentNames => nnames
         deallocate(onames)

         ! Extend attLists
         oatts => fld % attLists
         allocate(natts(new_n))
         do jj = 1, old_n
            natts(jj) = oatts(jj)
         end do
         do jj = old_n + 1, new_n
            allocate(natts(jj) % attList)
         end do
         fld % attLists => natts
         ! Shallow-copied elements; only deallocate the container
         deallocate(oatts)

      end subroutine extend_field3d

   end function atm_extend_scalars_for_chemistry
#endif


#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "setup_packages.inc"

#include "structs_and_variables.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

end module atm_core_interface
