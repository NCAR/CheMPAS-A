# MICM Mechanism Configurations

Chemistry mechanism configs for use with `config_micm_file` in the `&musica`
namelist group. Runtime tracer discovery reads species from these files —
switching mechanisms requires no Fortran or Registry changes.

## Configs

| File | Mechanism | Species | Notes |
|------|-----------|---------|-------|
| `abba.yaml` | ABBA | A, B, AB | Synthetic reversible reaction for coupling validation |
| `chapman.yaml` | Chapman (simple) | O2, O, O3 | 2 photolysis + 2 Arrhenius, O2 constant |
| `chapman_full.yaml` | Chapman (full) | M, O2, O, O1D, O3 | jO2 + jO3_O + jO3_O1D; pair with `tuvx_chapman.json` |
| `chapman.json` | Chapman (full) | M, O2, O, O1D, O3 | JSON version of `chapman_full.yaml` |
| `chapman_nox.yaml` | Chapman + NOx | +NO, NO2 | Chapman cycle with NOx catalytic destruction; pair with `tuvx_chapman_nox.json` |
| `lnox_o3.yaml` | LNOx-O3 | NO, NO2, O3 | Tropospheric photochemistry; pair with `tuvx_no2.json` |

## Usage

Copy the desired config to your run directory and set:

```
&musica
    config_micm_file = 'abba.yaml'
/
```
