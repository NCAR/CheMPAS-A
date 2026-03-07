module mpas_solar_geometry

    use mpas_kind_types, only : RKIND

    implicit none
    private
    public :: solar_cos_sza

    real(kind=RKIND), parameter :: PI = 3.14159265358979323846_RKIND
    real(kind=RKIND), parameter :: DEG2RAD = PI / 180.0_RKIND

contains

    !> Compute cos(SZA) from date/time and location.
    !>
    !> Uses Spencer (1971) formula for solar declination.
    !> Returns cos(SZA); negative values mean sun is below horizon.
    !>
    !> @param[in] DoY        Day of year (1-366)
    !> @param[in] hour_utc   Hour of day in UTC (fractional, 0-24)
    !> @param[in] lat_deg    Latitude in degrees North
    !> @param[in] lon_deg    Longitude in degrees East
    !> @return cos_sza       Cosine of solar zenith angle
    function solar_cos_sza(DoY, hour_utc, lat_deg, lon_deg) result(cos_sza)

        integer, intent(in)           :: DoY
        real(kind=RKIND), intent(in)  :: hour_utc
        real(kind=RKIND), intent(in)  :: lat_deg
        real(kind=RKIND), intent(in)  :: lon_deg
        real(kind=RKIND)              :: cos_sza

        real(kind=RKIND) :: gamma, decl, eqtime, hour_angle, lat_rad

        ! Fractional year [radians] — Spencer (1971)
        gamma = 2.0_RKIND * PI * (real(DoY - 1, RKIND) + (hour_utc - 12.0_RKIND) / 24.0_RKIND) / 365.0_RKIND

        ! Solar declination [radians]
        decl = 0.006918_RKIND &
             - 0.399912_RKIND * cos(gamma) &
             + 0.070257_RKIND * sin(gamma) &
             - 0.006758_RKIND * cos(2.0_RKIND * gamma) &
             + 0.000907_RKIND * sin(2.0_RKIND * gamma) &
             - 0.002697_RKIND * cos(3.0_RKIND * gamma) &
             + 0.001480_RKIND * sin(3.0_RKIND * gamma)

        ! Equation of time [minutes]
        eqtime = 229.18_RKIND * ( &
                   0.000075_RKIND &
                 + 0.001868_RKIND * cos(gamma) &
                 - 0.032077_RKIND * sin(gamma) &
                 - 0.014615_RKIND * cos(2.0_RKIND * gamma) &
                 - 0.040849_RKIND * sin(2.0_RKIND * gamma) )

        ! Hour angle [radians]
        ! Solar noon occurs when hour_angle = 0
        hour_angle = DEG2RAD * ((hour_utc * 60.0_RKIND + eqtime + 4.0_RKIND * lon_deg) / 4.0_RKIND - 180.0_RKIND)

        ! Latitude in radians
        lat_rad = lat_deg * DEG2RAD

        ! Cosine of solar zenith angle
        cos_sza = sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(hour_angle)

    end function solar_cos_sza

end module mpas_solar_geometry
