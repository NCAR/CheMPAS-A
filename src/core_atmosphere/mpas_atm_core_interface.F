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
   use mpas_kind_types, only : StrKIND, ShortStrKIND

   implicit none

#ifdef MPAS_USE_MUSICA
   integer :: n_runtime_chem_species = 0
   integer :: n_runtime_photo_rates = 0
   integer :: n_runtime_emission_species = 0
   integer :: n_runtime_emission_sectors = 0
   integer :: n_runtime_emission_categories = 0
   character(len=64), allocatable :: runtime_chem_species(:)
   character(len=64), allocatable :: runtime_photo_rates(:)
   character(len=64), allocatable :: runtime_emission_species(:)
   character(len=64), allocatable :: runtime_emission_sectors(:)
   character(len=64), allocatable :: runtime_emission_sector_labels(:)
   integer, allocatable :: runtime_emission_categories(:)
   logical, save :: runtime_emission_layered = .false.
   integer, save :: runtime_emission_max_fields = 256
   character(len=StrKIND), save :: runtime_emission_sector_text = ''
   character(len=StrKIND), save :: runtime_emission_category_text = ''
   logical, save :: runtime_vars_prepared = .false.
#endif

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
      core % setup_derived_dimensions => atm_setup_derived_dimensions_with_runtime_chemistry
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
   function atm_setup_packages(configs, streamInfo, packages, iocontext, dminfo) result(ierr)

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
      type (dm_info), intent(in) :: dminfo
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


#ifdef MPAS_USE_MUSICA
   !-----------------------------------------------------------------------
   !
   !  function atm_prepare_runtime_chemistry_vars
   !
   !> \brief   Cache chemistry tracer names for runtime-var generation
   !
   !-----------------------------------------------------------------------
   function atm_prepare_runtime_chemistry_vars(block) result(ierr)

      use iso_fortran_env, only : real64
      use mpas_atm_chemistry, only : chemistry_query_species, chemistry_query_photo_rates, &
                                     chemistry_query_emission_species, &
                                     chemistry_prepare_miem_diagnostic_selection, &
                                     chemistry_prepare_tracer_outputs, chemistry_clear_tracer_outputs
      use mpas_derived_types, only : block_type, MPAS_LOG_CRIT, MPAS_LOG_ERR
      use mpas_kind_types,    only : StrKIND, RKIND
      use mpas_log,           only : mpas_log_write

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      character(len=StrKIND), pointer :: config_micm_file
      character(len=StrKIND), pointer :: config_miem_file
      character(len=StrKIND), pointer :: config_miem_diagnostic_sectors
      character(len=StrKIND), pointer :: config_miem_diagnostic_categories
      character(len=StrKIND), pointer :: config_tuvx_config_file
      character(len=StrKIND), pointer :: config_chem_tracer_output_unit
      character(len=StrKIND), pointer :: config_chem_tracer_output_overrides
      logical, pointer :: config_apply_lbcs
      logical, pointer :: config_miem_layered_diagnostics
      integer, pointer :: config_miem_max_diagnostic_fields
      character(len=64), allocatable :: emission_names(:)
      character(len=64), allocatable :: emission_sectors(:)
      character(len=64), allocatable :: emission_sector_labels(:)
      integer, allocatable :: emission_categories(:)
      character(len=StrKIND) :: emission_sector_text, emission_category_text
      logical :: emission_layered
      integer :: emission_max_fields
      character(len=64), allocatable :: micm_names(:)
      character(len=64), allocatable :: mpas_names(:)
      character(len=64), allocatable :: photo_names(:)
      real(real64), allocatable :: molar_masses(:)
      real(real64), allocatable :: fixed_vmrs(:)
      logical, allocatable :: host_bound(:)
      real(kind=RKIND), pointer :: config_j_no2_max
      integer :: num_micm_species
      integer :: num_photo_rates
      integer :: num_emission_species
      integer :: error_code
      integer :: err_level
      integer :: i, runtime_idx
      character(len=:), allocatable :: error_message
      character(len=StrKIND) :: output_unit, output_overrides

      ierr = 0
      if (runtime_vars_prepared) then
         ierr = 0
         return   ! already queried (multi-block/second call); reuse the cache
      end if

      n_runtime_chem_species = 0
      n_runtime_photo_rates = 0
      n_runtime_emission_species = 0
      n_runtime_emission_sectors = 0
      n_runtime_emission_categories = 0
      runtime_emission_layered = .false.
      runtime_emission_max_fields = 256
      runtime_emission_sector_text = ''
      runtime_emission_category_text = ''
      if (allocated(runtime_chem_species)) deallocate(runtime_chem_species)
      if (allocated(runtime_photo_rates)) deallocate(runtime_photo_rates)
      if (allocated(runtime_emission_species)) deallocate(runtime_emission_species)
      if (allocated(runtime_emission_sectors)) deallocate(runtime_emission_sectors)
      if (allocated(runtime_emission_sector_labels)) deallocate(runtime_emission_sector_labels)
      if (allocated(runtime_emission_categories)) deallocate(runtime_emission_categories)
      call chemistry_clear_tracer_outputs()

      nullify(config_micm_file)
      nullify(config_miem_file)
      nullify(config_miem_diagnostic_sectors)
      nullify(config_miem_diagnostic_categories)
      nullify(config_miem_layered_diagnostics)
      nullify(config_miem_max_diagnostic_fields)
      nullify(config_tuvx_config_file)
      nullify(config_chem_tracer_output_unit)
      nullify(config_chem_tracer_output_overrides)
      nullify(config_j_no2_max)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, 'config_micm_file', config_micm_file)
      call mpas_pool_get_config(block % domain % configs, 'config_miem_file', config_miem_file)
      call mpas_pool_get_config(block % domain % configs, 'config_miem_diagnostic_sectors', &
                                config_miem_diagnostic_sectors)
      call mpas_pool_get_config(block % domain % configs, 'config_miem_diagnostic_categories', &
                                config_miem_diagnostic_categories)
      call mpas_pool_get_config(block % domain % configs, 'config_miem_layered_diagnostics', &
                                config_miem_layered_diagnostics)
      call mpas_pool_get_config(block % domain % configs, 'config_miem_max_diagnostic_fields', &
                                config_miem_max_diagnostic_fields)
      call mpas_pool_get_config(block % domain % configs, 'config_tuvx_config_file', config_tuvx_config_file)
      call mpas_pool_get_config(block % domain % configs, 'config_chem_tracer_output_unit', &
                                config_chem_tracer_output_unit)
      call mpas_pool_get_config(block % domain % configs, 'config_chem_tracer_output_overrides', &
                                config_chem_tracer_output_overrides)
      call mpas_pool_get_config(block % domain % configs, 'config_j_no2_max', config_j_no2_max)
      call mpas_pool_set_error_level(err_level)

      output_unit = 'mmr'
      output_overrides = ''
      if (associated(config_chem_tracer_output_unit)) output_unit = config_chem_tracer_output_unit
      if (associated(config_chem_tracer_output_overrides)) output_overrides = config_chem_tracer_output_overrides

      if (.not. associated(config_micm_file)) then
         call mpas_log_write('[Chemistry] No config_micm_file found; skipping runtime chemistry tracers')
      else if (len_trim(config_micm_file) == 0) then
         call mpas_log_write('[Chemistry] config_micm_file is empty; skipping runtime chemistry tracers')
      else
         call chemistry_query_species(config_micm_file, num_micm_species, micm_names, mpas_names, &
                                      molar_masses, host_bound, fixed_vmrs, error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Error querying MICM species: ' // trim(error_message), &
                                messageType=MPAS_LOG_ERR)
            ierr = error_code
            return
         end if

         call chemistry_prepare_tracer_outputs(output_unit, output_overrides, num_micm_species, &
                                               micm_names, mpas_names, molar_masses, fixed_vmrs, &
                                               error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Error configuring MICM tracer output units: ' // &
                                trim(error_message), messageType=MPAS_LOG_ERR)
            ierr = error_code
            return
         end if

         n_runtime_chem_species = count(.not. host_bound)
         if (n_runtime_chem_species == 0) then
            call mpas_log_write('[Chemistry] MICM config has no runtime chemistry tracers to add')
         else
            config_apply_lbcs => null()
            call mpas_pool_get_config(block % domain % configs, 'config_apply_lbcs', config_apply_lbcs)
            if (associated(config_apply_lbcs)) then
               if (config_apply_lbcs) then
                  call mpas_log_write('Runtime chemistry tracers are incompatible with config_apply_lbcs=true. ' // &
                                      'lbc_scalars is statically sized and cannot accommodate dynamic tracers.', &
                                      messageType=MPAS_LOG_CRIT)
                  ierr = 1
                  return
               end if
            end if

            allocate(runtime_chem_species(n_runtime_chem_species))
            runtime_idx = 0
            do i = 1, num_micm_species
               if (host_bound(i)) cycle
               runtime_idx = runtime_idx + 1
               runtime_chem_species(runtime_idx) = mpas_names(i)
            end do

            call mpas_log_write('[Chemistry] Prepared $i runtime chemistry tracers', &
                                intArgs=[n_runtime_chem_species])
         end if

         if (allocated(micm_names)) deallocate(micm_names)
         if (allocated(mpas_names)) deallocate(mpas_names)
         if (allocated(molar_masses)) deallocate(molar_masses)
         if (allocated(fixed_vmrs)) deallocate(fixed_vmrs)
         if (allocated(host_bound)) deallocate(host_bound)
      end if

      ! config_miem_file is only the enable switch. Diagnostic names are a
      ! MICM contract and must be discovered from config_micm_file, even when
      ! the MIEM file has an intentionally empty mechanism section.
      num_emission_species = 0
      if (associated(config_miem_file) .and. len_trim(config_miem_file) > 0) then
         if (.not. associated(config_micm_file) .or. len_trim(config_micm_file) == 0) then
            call mpas_log_write('[Chemistry] MIEM diagnostics require a nonempty config_micm_file', &
                                messageType=MPAS_LOG_CRIT)
            ierr = 1
            return
         end if

         call chemistry_query_emission_species(config_micm_file, num_emission_species, &
                                               emission_names, error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Error querying MICM emission diagnostics: ' // &
                                trim(error_message), messageType=MPAS_LOG_ERR)
            ierr = error_code
            return
         end if
         if (num_emission_species < 1) then
            call mpas_log_write('[Chemistry] MIEM is enabled but config_micm_file contains no EMIS.* parameters', &
                                messageType=MPAS_LOG_CRIT)
            ierr = 1
            if (allocated(emission_names)) deallocate(emission_names)
            return
         end if

         n_runtime_emission_species = num_emission_species
         allocate(runtime_emission_species(n_runtime_emission_species))
         runtime_emission_species(:) = emission_names(:)
         if (allocated(emission_names)) deallocate(emission_names)

         emission_sector_text = ''
         emission_category_text = ''
         emission_layered = .false.
         emission_max_fields = 256
         if (associated(config_miem_diagnostic_sectors)) then
            emission_sector_text = config_miem_diagnostic_sectors
         end if
         if (associated(config_miem_diagnostic_categories)) then
            emission_category_text = config_miem_diagnostic_categories
         end if
         if (associated(config_miem_layered_diagnostics)) then
            emission_layered = config_miem_layered_diagnostics
         end if
         if (associated(config_miem_max_diagnostic_fields)) then
            emission_max_fields = config_miem_max_diagnostic_fields
         end if
         ! nVertLevels is read from the mesh after block descriptors are built.
         ! Validate names and a one-level lower bound here; the derived-dimension
         ! callback performs the definitive cap check before field allocation.
         call chemistry_prepare_miem_diagnostic_selection(emission_sector_text, &
                                                          emission_category_text, 1, &
                                                          n_runtime_emission_species, &
                                                          emission_layered, emission_max_fields, &
                                                          emission_sectors, &
                                                          emission_sector_labels, &
                                                          emission_categories, &
                                                          error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Invalid MIEM diagnostic selection: ' // &
                                trim(error_message), messageType=MPAS_LOG_CRIT)
            ierr = error_code
            return
         end if
         n_runtime_emission_sectors = size(emission_sectors)
         n_runtime_emission_categories = size(emission_categories)
         runtime_emission_layered = emission_layered
         runtime_emission_max_fields = emission_max_fields
         runtime_emission_sector_text = emission_sector_text
         runtime_emission_category_text = emission_category_text
         call move_alloc(emission_sectors, runtime_emission_sectors)
         call move_alloc(emission_sector_labels, runtime_emission_sector_labels)
         call move_alloc(emission_categories, runtime_emission_categories)
         call atm_validate_runtime_emission_diag_names(error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Invalid MIEM diagnostic field names: ' // &
                                trim(error_message), messageType=MPAS_LOG_CRIT)
            ierr = error_code
            return
         end if
         call mpas_log_write('[Chemistry] Prepared $i runtime emissions diagnostics', &
                             intArgs=[n_runtime_emission_species])
         call mpas_log_write('[Chemistry] Prepared $i sector and $i category emissions groups', &
                             intArgs=[n_runtime_emission_sectors, n_runtime_emission_categories])
      else
         call mpas_log_write('[Chemistry] MIEM disabled; no runtime emissions diagnostics to add')
      end if

      if (associated(config_tuvx_config_file) .and. len_trim(config_tuvx_config_file) > 0) then
         call chemistry_query_photo_rates(config_tuvx_config_file, num_photo_rates, photo_names, &
                                          error_code, error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Error querying TUV-x photolysis rates: ' // &
                                trim(error_message), messageType=MPAS_LOG_ERR)
            ierr = error_code
            return
         end if
      else if (associated(config_j_no2_max) .and. config_j_no2_max > 0.0_RKIND) then
         num_photo_rates = 1
         allocate(photo_names(1))
         photo_names(1) = 'jNO2'
      else
         num_photo_rates = 0
      end if

      if (num_photo_rates > 0) then
         n_runtime_photo_rates = num_photo_rates
         allocate(runtime_photo_rates(n_runtime_photo_rates))
         runtime_photo_rates(:) = photo_names(:)
         if (allocated(photo_names)) deallocate(photo_names)
         call mpas_log_write('[Chemistry] Prepared $i runtime photolysis diagnostics', &
                             intArgs=[n_runtime_photo_rates])
      else
         call mpas_log_write('[Chemistry] No runtime photolysis diagnostics to add')
         if (allocated(photo_names)) deallocate(photo_names)
      end if

      runtime_vars_prepared = .true.

   end function atm_prepare_runtime_chemistry_vars


   !-----------------------------------------------------------------------
   !
   !  subroutine atm_validate_runtime_emission_diag_names
   !
   !> \brief Reject overlength or colliding generated MPAS field names
   !
   !-----------------------------------------------------------------------
   subroutine atm_validate_runtime_emission_diag_names(error_code, error_message)

      implicit none

      integer, intent(out) :: error_code
      character(len=:), allocatable, intent(out) :: error_message

      character(len=ShortStrKIND), allocatable :: field_names(:)
      character(len=StrKIND) :: candidate
      character(len=32) :: category_text
      integer :: fields_per_species, i, j, name_count

      error_code = 0
      error_message = ''
      fields_per_species = 1 + n_runtime_emission_sectors + n_runtime_emission_categories
      if (runtime_emission_layered) fields_per_species = 2 * fields_per_species
      allocate(field_names(n_runtime_emission_species * fields_per_species))
      field_names = ''
      name_count = 0

      do i = 1, n_runtime_emission_species
         candidate = 'emis_' // trim(runtime_emission_species(i))
         call append_name(candidate)
         if (error_code /= 0) return
         if (runtime_emission_layered) then
            candidate = 'emis_' // trim(runtime_emission_species(i)) // '_layer'
            call append_name(candidate)
            if (error_code /= 0) return
         end if
         do j = 1, n_runtime_emission_sectors
            candidate = 'emis_' // trim(runtime_emission_species(i)) // &
                        '__sector_' // trim(runtime_emission_sector_labels(j))
            call append_name(candidate)
            if (error_code /= 0) return
            if (runtime_emission_layered) then
               candidate = 'emis_' // trim(runtime_emission_species(i)) // &
                           '__sector_' // trim(runtime_emission_sector_labels(j)) // '_layer'
               call append_name(candidate)
               if (error_code /= 0) return
            end if
         end do
         do j = 1, n_runtime_emission_categories
            write(category_text, '(i0)') runtime_emission_categories(j)
            candidate = 'emis_' // trim(runtime_emission_species(i)) // &
                        '__category_' // trim(category_text)
            call append_name(candidate)
            if (error_code /= 0) return
            if (runtime_emission_layered) then
               candidate = 'emis_' // trim(runtime_emission_species(i)) // &
                           '__category_' // trim(category_text) // '_layer'
               call append_name(candidate)
               if (error_code /= 0) return
            end if
         end do
      end do

   contains

      subroutine append_name(name)

         character(len=*), intent(in) :: name
         integer :: previous

         if (len_trim(name) > ShortStrKIND) then
            error_code = 1
            error_message = '[MIEM] generated field name exceeds the MPAS 64-character limit: ' // &
                            trim(name)
            return
         end if
         do previous = 1, name_count
            if (trim(field_names(previous)) == trim(name)) then
               error_code = 1
               error_message = '[MIEM] generated emissions diagnostic field name collides: ' // &
                               trim(name)
               return
            end if
         end do
         name_count = name_count + 1
         field_names(name_count) = trim(name)

      end subroutine append_name

   end subroutine atm_validate_runtime_emission_diag_names


   !-----------------------------------------------------------------------
   !
   !  subroutine atm_add_runtime_chemistry_vars
   !
   !> \brief   Runtime-var callback that adds chemistry tracers
   !
   !-----------------------------------------------------------------------
   subroutine atm_add_runtime_chemistry_vars(block, structPool, packagePool, newVarArray, poolName, varArrayName)

      use mpas_atm_chemistry, only : chemistry_add_tracer_output_vars
      use mpas_derived_types, only : block_type, mpas_pool_type, MPAS_LOG_CRIT
      use mpas_runtime_vars,  only : MPAS_var_array, MPAS_var, MPAS_VAR_REAL, MPAS_VAR_OK
      use mpas_log,           only : mpas_log_write

      implicit none

      type (block_type), target, intent(inout) :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(in) :: packagePool
      type (MPAS_var_array) :: newVarArray
      character(len=*), intent(in) :: poolName
      character(len=*), intent(in) :: varArrayName

      type (MPAS_var) :: newVar
      character(len=64) :: tracer_name
      character(len=128) :: diag_name
      character(len=32) :: category_text
      integer :: i, j, addErr

      if (trim(poolName) == 'state' .and. trim(varArrayName) == 'scalars') then
         if (n_runtime_chem_species <= 0) return
         do i = 1, n_runtime_chem_species
            tracer_name = atm_runtime_chem_tracer_name(runtime_chem_species(i))

            newVar % name = trim(tracer_name)
            newVar % hasTimeDimension = .true.
            newVar % vartype = MPAS_VAR_REAL
            newVar % dimensions = 'nVertLevels nCells'
            newVar % units = 'kg kg^{-1}'
            newVar % description = 'Chemistry tracer ' // trim(tracer_name) // ' mixing ratio'
            newVar % isPersistent = .true.
            call newVarArray % add_constituent(newVar, 'chemistry', ierr=addErr)
            if (addErr /= MPAS_VAR_OK) then
               call mpas_log_write('Failed to add runtime tracer ' // trim(tracer_name) // &
                                   ' to state/scalars', messageType=MPAS_LOG_CRIT)
               return
            end if

            call mpas_log_write('[Chemistry] Added runtime tracer ' // trim(tracer_name) // &
                                ' to state/scalars')
            call newVar % reset()
         end do
      else if (trim(poolName) == 'tend' .and. trim(varArrayName) == 'scalars_tend') then
         if (n_runtime_chem_species <= 0) return
         do i = 1, n_runtime_chem_species
            tracer_name = atm_runtime_chem_tracer_name(runtime_chem_species(i))

            newVar % name = 'tend_' // trim(tracer_name)
            newVar % name_in_code = trim(tracer_name)
            newVar % hasTimeDimension = .true.
            newVar % vartype = MPAS_VAR_REAL
            newVar % dimensions = 'nVertLevels nCells'
            newVar % units = 'kg m^{-3} s^{-1}'
            newVar % description = 'Tendency of chemistry tracer ' // trim(tracer_name)
            newVar % isPersistent = .true.
            call newVarArray % add_constituent(newVar, 'chemistry', ierr=addErr)
            if (addErr /= MPAS_VAR_OK) then
               call mpas_log_write('Failed to add runtime tracer tendency tend_' // &
                                   trim(tracer_name) // ' to tend/scalars_tend', &
                                   messageType=MPAS_LOG_CRIT)
               return
            end if

            call mpas_log_write('[Chemistry] Added runtime tracer tendency tend_' // &
                                trim(tracer_name) // ' to tend/scalars_tend')
            call newVar % reset()
         end do
      else if (trim(poolName) == 'diag' .and. len_trim(varArrayName) == 0) then
         call chemistry_add_tracer_output_vars(block, structPool, packagePool)
         if (n_runtime_photo_rates <= 0 .and. n_runtime_emission_species <= 0) return
         do i = 1, n_runtime_photo_rates
            diag_name = 'j_' // trim(runtime_photo_rates(i))

            newVar % name = trim(diag_name)
            newVar % hasTimeDimension = .true.
            newVar % vartype = MPAS_VAR_REAL
            newVar % dimensions = 'nVertLevels nCells'
            newVar % units = 's^{-1}'
            newVar % description = 'Photolysis rate ' // trim(runtime_photo_rates(i)) // &
                                   ' used by chemistry'
            newVar % isPersistent = .true.
            newVar % nTimeLevels = 1
            call newVar % add_to_pool(block, structPool, packagePool, useNameInCode=.true., ierr=addErr)
            if (addErr /= MPAS_VAR_OK) then
               call mpas_log_write('Failed to add runtime photolysis diagnostic ' // &
                                   trim(diag_name) // ' to diag', messageType=MPAS_LOG_CRIT)
               return
            end if
            call newVar % add_to_pool(block, block % allFields, packagePool, ierr=addErr)
            if (addErr /= MPAS_VAR_OK) then
               call mpas_log_write('Failed to add runtime photolysis diagnostic ' // &
                                   trim(diag_name) // ' to allFields', messageType=MPAS_LOG_CRIT)
               return
            end if

            call mpas_log_write('[Chemistry] Added runtime photolysis diagnostic ' // &
                                trim(diag_name) // ' to diag')
            call newVar % reset()
         end do

         do i = 1, n_runtime_emission_species
            diag_name = 'emis_' // trim(runtime_emission_species(i))
            call add_runtime_emission_diag(trim(diag_name), 'nCells', &
                 'MIEM column flux applied to chemistry for ' // &
                 trim(runtime_emission_species(i)), addErr)
            if (addErr /= MPAS_VAR_OK) return

            if (runtime_emission_layered) then
               diag_name = 'emis_' // trim(runtime_emission_species(i)) // '_layer'
               call add_runtime_emission_diag(trim(diag_name), 'nVertLevels nCells', &
                    'MIEM per-layer flux applied to chemistry for ' // &
                    trim(runtime_emission_species(i)), addErr)
               if (addErr /= MPAS_VAR_OK) return
            end if

            do j = 1, n_runtime_emission_sectors
               diag_name = 'emis_' // trim(runtime_emission_species(i)) // &
                           '__sector_' // trim(runtime_emission_sector_labels(j))
               call add_runtime_emission_diag(trim(diag_name), 'nCells', &
                    'MIEM column flux for sector ' // trim(runtime_emission_sectors(j)) // &
                    ', species ' // trim(runtime_emission_species(i)), addErr)
               if (addErr /= MPAS_VAR_OK) return
               if (runtime_emission_layered) then
                  diag_name = 'emis_' // trim(runtime_emission_species(i)) // &
                              '__sector_' // trim(runtime_emission_sector_labels(j)) // '_layer'
                  call add_runtime_emission_diag(trim(diag_name), 'nVertLevels nCells', &
                       'MIEM per-layer flux for sector ' // trim(runtime_emission_sectors(j)) // &
                       ', species ' // trim(runtime_emission_species(i)), addErr)
                  if (addErr /= MPAS_VAR_OK) return
               end if
            end do

            do j = 1, n_runtime_emission_categories
               write(category_text, '(i0)') runtime_emission_categories(j)
               diag_name = 'emis_' // trim(runtime_emission_species(i)) // &
                           '__category_' // trim(category_text)
               call add_runtime_emission_diag(trim(diag_name), 'nCells', &
                    'MIEM column flux for category ' // trim(category_text) // &
                    ', species ' // trim(runtime_emission_species(i)), addErr)
               if (addErr /= MPAS_VAR_OK) return
               if (runtime_emission_layered) then
                  diag_name = 'emis_' // trim(runtime_emission_species(i)) // &
                              '__category_' // trim(category_text) // '_layer'
                  call add_runtime_emission_diag(trim(diag_name), 'nVertLevels nCells', &
                       'MIEM per-layer flux for category ' // trim(category_text) // &
                       ', species ' // trim(runtime_emission_species(i)), addErr)
                  if (addErr /= MPAS_VAR_OK) return
               end if
            end do
         end do
      end if

   contains

      subroutine add_runtime_emission_diag(field_name, field_dimensions, field_description, status)

         character(len=*), intent(in) :: field_name, field_dimensions, field_description
         integer, intent(out) :: status

         newVar % name = trim(field_name)
         newVar % hasTimeDimension = .true.
         newVar % vartype = MPAS_VAR_REAL
         newVar % dimensions = trim(field_dimensions)
         newVar % units = 'kg m^{-2} s^{-1}'
         newVar % description = trim(field_description)
         newVar % isPersistent = .true.
         newVar % nTimeLevels = 1
         call newVar % add_to_pool(block, structPool, packagePool, useNameInCode=.true., ierr=status)
         if (status /= MPAS_VAR_OK) then
            call mpas_log_write('Failed to add runtime emissions diagnostic ' // &
                                trim(field_name) // ' to diag', messageType=MPAS_LOG_CRIT)
            call newVar % reset()
            return
         end if
         call newVar % add_to_pool(block, block % allFields, packagePool, ierr=status)
         if (status /= MPAS_VAR_OK) then
            call mpas_log_write('Failed to add runtime emissions diagnostic ' // &
                                trim(field_name) // ' to allFields', messageType=MPAS_LOG_CRIT)
            call newVar % reset()
            return
         end if
         call mpas_log_write('[Chemistry] Added runtime emissions diagnostic ' // &
                             trim(field_name) // ' to diag')
         call newVar % reset()

      end subroutine add_runtime_emission_diag

   end subroutine atm_add_runtime_chemistry_vars


   function atm_runtime_chem_tracer_name(species_name) result(tracer_name)

      implicit none

      character(len=*), intent(in) :: species_name
      character(len=64) :: tracer_name

      tracer_name = trim(species_name)

   end function atm_runtime_chem_tracer_name
#endif


   !***********************************************************************
   !
   !  function atm_setup_derived_dimensions_with_runtime_chemistry
   !
   !> \brief Assign derived dimensions and validate dimension-dependent
   !>        runtime chemistry diagnostics before field allocation.
   !
   !-----------------------------------------------------------------------
   function atm_setup_derived_dimensions_with_runtime_chemistry(readDimensions, dimensionPool, &
                                                                 configPool) result(ierr)

      use mpas_derived_types, only : mpas_pool_type, MPAS_LOG_ERR
      use mpas_pool_routines, only : mpas_pool_get_dimension
      use mpas_log, only : mpas_log_write
#ifdef MPAS_USE_MUSICA
      use mpas_atm_chemistry, only : chemistry_prepare_miem_diagnostic_selection
#endif

      implicit none

      type (mpas_pool_type), intent(inout) :: readDimensions
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(inout) :: configPool
      integer :: ierr

#ifdef MPAS_USE_MUSICA
      integer, pointer :: nVertLevels
      character(len=64), allocatable :: sectors(:), sector_labels(:)
      integer, allocatable :: categories(:)
      integer :: error_code
      character(len=:), allocatable :: error_message
#endif

      ierr = atm_setup_derived_dimensions(readDimensions, dimensionPool, configPool)
      if (ierr /= 0) return

#ifdef MPAS_USE_MUSICA
      if (runtime_vars_prepared .and. n_runtime_emission_species > 0) then
         nullify(nVertLevels)
         call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
         if (.not. associated(nVertLevels)) then
            call mpas_log_write('[Chemistry] MIEM diagnostic cap validation requires nVertLevels', &
                                messageType=MPAS_LOG_ERR)
            ierr = 1
            return
         end if

         call chemistry_prepare_miem_diagnostic_selection(runtime_emission_sector_text, &
                                                          runtime_emission_category_text, &
                                                          nVertLevels, &
                                                          n_runtime_emission_species, &
                                                          runtime_emission_layered, &
                                                          runtime_emission_max_fields, sectors, &
                                                          sector_labels, categories, error_code, &
                                                          error_message)
         if (error_code /= 0) then
            call mpas_log_write('[Chemistry] Invalid MIEM diagnostic selection for mesh dimensions: ' // &
                                trim(error_message), messageType=MPAS_LOG_ERR)
            ierr = error_code
            return
         end if
         if (size(sectors) /= n_runtime_emission_sectors .or. &
             size(categories) /= n_runtime_emission_categories) then
            call mpas_log_write('[Chemistry] MIEM diagnostic selection changed during block setup', &
                                messageType=MPAS_LOG_ERR)
            ierr = 1
            return
         end if
         if (allocated(sectors)) deallocate(sectors)
         if (allocated(sector_labels)) deallocate(sector_labels)
         if (allocated(categories)) deallocate(categories)
      end if
#endif

   end function atm_setup_derived_dimensions_with_runtime_chemistry


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
#ifdef MPAS_USE_MUSICA
      use mpas_runtime_vars, only : MPAS_var_add_callback, MPAS_var_clear_callback
#endif

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      integer, pointer :: cam_pcnst
      integer :: err_level

      ierr = 0

#ifdef MPAS_USE_MUSICA
      call MPAS_var_clear_callback()
      ierr = atm_prepare_runtime_chemistry_vars(block)
      if (ierr /= 0) then
         return
      end if
      call MPAS_var_add_callback(atm_add_runtime_chemistry_vars)
#endif

      call atm_generate_pools(block, block % structs, block % packages)

#ifdef MPAS_USE_MUSICA
      call MPAS_var_clear_callback()
#endif

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



#include "setup_immutable_streams.inc"

#include "block_dimension_routines.inc"

#include "define_packages.inc"

#include "setup_packages.inc"

#include "namelist_call.inc"

#include "namelist_defines.inc"

#include "pools_and_vars.inc"

end module atm_core_interface
