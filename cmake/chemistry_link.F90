! Link-only preflight: exercise both the selected chemistry API and MUSICA MIEM.
program chemistry_link
#ifdef MPAS_USE_TMT
  use tmt_micm, only: get_micm_version
  use tmt_tuvx, only: get_tuvx_version
  use tmt_util, only: string_t
#else
  use musica_micm, only: get_micm_version
  use musica_tuvx, only: get_tuvx_version
  use musica_util, only: string_t
#endif
  use musica_util, only: error_t
  use musica_emissions, only: mechanism_t, emissions_t
  implicit none
  type(string_t) :: version
  type(error_t) :: error
  type(mechanism_t), pointer :: mechanism
  type(emissions_t), pointer :: emissions

  version = get_micm_version()
  print *, version%value_
  version = get_tuvx_version()
  print *, version%value_
  mechanism => mechanism_t('link-probe-does-not-run.yaml', error)
  if (associated(mechanism)) then
    emissions => emissions_t(mechanism, 1, 1, error)
    if (associated(emissions)) deallocate(emissions)
    deallocate(mechanism)
  end if
end program chemistry_link
