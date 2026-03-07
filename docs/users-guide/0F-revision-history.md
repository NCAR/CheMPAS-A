# Appendix F: Revision History

### 2 June 2025

- Add new namelist options and fields for MPAS v8.3.0.
- Update physics options in Chapter 6

### 27 June 2024

- Add new namelist options and fields for MPAS v8.2.0.
- Update physics options in Chapter 6

### 18 April 2024

- Add new Section 8.3 for invariant stream capability.
- Add new namelist options and fields for MPAS v8.1.0.

### 6 July 2023

- Update units and description for several fields following the MPAS-A v8.0.1 release.
- Add missing diagnostic fields to Appendix D.

### 16 June 2023

- Updates for MPAS v8.0.0, including new namelist options.
- Updated physics scheme versions in Table 6.3

### 8 July 2019

- Removed references to the `double_to_float_grid` utility from Section 3.4.

### 8 June 2019

- Significant updates to reflect MPAS v7.0, principally, new namelist options, model fields, and a description of the limited-area simulation capability
- Moved description of SST and sea-ice updates to a new Chapter 8
- Updated physics scheme versions in Table 6.3

### 11 May 2019

- Update the default value for `config_o3climatology` to match the v6.3 release: the default is now 'true'.

### 17 April 2018

- Minor updates for MPAS v6.0, including new `config_topo_data` option and new `initial_time` field.

### 22 March 2018

- Update Appendix D with correct units for `scalars_tend`.

### 1 August 2017

- Update documentation in Chapter 3 for building with PIO 2.x versions
- Correct other minor typographical errors and unclear wording

### 12 May 2017

- Add documentation for new `config_extrap_airtemp` initialization namelist option
- Remove documentation for `config_smdiv_p_forward` model option, which was removed in v5.1
- Add note that `config_len_disp` is also used by 3-d divergence damping in v5.1
- Change default value of `config_smdiv` model option to 0.1, reflecting changes in v5.1
- Other minor changes to match the MPAS-Atmosphere v5.1 release

### 23 December 2016

- Added description of the new 'convection_permitting' physics suite, as well as a list of all possible physics parameterizations, to Chapter 6.
- Updated the chapter on building MPAS to mention specific versions of PIO that are known to work with MPAS, information on the method used to build single-precision executables in v5.0, and a mention of the *experimental* OpenMP capability in v5.0.
- Appendices A, B, and D are now automatically generated based on XML Registry files.
- Many other small changes to reflect the state of the v5.0 release.
- Minor wording changes and corrections of typographical errors throughout document.

### 19 May 2015

- Added chapter describing the MPAS-Atmosphere physics suites introduced in MPAS v4.0.
- Added a section on building the model with single-precision reals.
- Updated the I/O chapter to include the new `io_type` option introduced in MPAS v4.0.
- Updated a few namelist options to match MPAS v4.0 code.
- Minor wording changes and corrections of typographical errors throughout document.

### 18 November 2014

- Added chapter describing the MPAS runtime I/O system available in MPAS v3.0.
- Updated the chapter on running MPAS-Atmosphere to match MPAS v3.0 code.
- Updated a few namelist options to match MPAS v3.0 code.
- Minor wording changes throughout document.

### 13 June 2013

Initial version.
