#!/usr/bin/env python3
"""Standalone MUSICA-Python box model for LNOx + O3 photochemistry.

Mirrors section 2.6 of the CheMPAS-A tutorial without MPAS and without
the lightning-NOx source (which is a CheMPAS-side operator-split
injection in mpas_lightning_nox.F, not part of the MICM mechanism).
Single grid cell at mid-tropospheric conditions, hardcoded
jNO2 = 0.01 s^-1 (matches CheMPAS-A's config_lnox_j_no2), 2 h
integration. Output: lnox_box.nc, lnox_box.png.

Run:
    ~/miniconda3/envs/mpas/bin/python scripts/musica_python/lnox_box.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import musica
from musica.constants import GAS_CONSTANT
from musica.mechanism_configuration import Parser
from musica.micm.solver_result import SolverState

import style

REPO_ROOT = Path(__file__).parent.parent.parent
LNOX_YAML = REPO_ROOT / "micm_configs" / "lnox_o3.yaml"
OUT_DIR = Path(__file__).parent

# Mid-troposphere conditions (~5-8 km), where the CheMPAS-A LNOx source
# is active.
T_REF = 240.0   # K
P_REF = 5e4     # Pa

# Photolysis rate matches CheMPAS-A's config_lnox_j_no2 setting.
J_NO2 = 0.01    # s^-1


def ppb_to_mol_m3(ppb: float, T: float, P: float) -> float:
    """Volume mixing ratio in ppb to concentration in mol m^-3."""
    return ppb * 1e-9 * P / (GAS_CONSTANT * T)


def main() -> None:
    style.apply_ncar_style()

    parser = Parser()
    mechanism = parser.parse(str(LNOX_YAML))
    solver = musica.MICM(
        mechanism=mechanism,
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(number_of_grid_cells=1)
    state.set_conditions(temperatures=T_REF, pressures=P_REF)

    # 1 ppb total NOx, 50/50 NO/NO2; 50 ppb O3 background.
    nox_each = ppb_to_mol_m3(0.5, T_REF, P_REF)
    o3 = ppb_to_mol_m3(50.0, T_REF, P_REF)
    state.set_concentrations({
        "NO":  [nox_each],
        "NO2": [nox_each],
        "O3":  [o3],
    })

    # PHOTOLYSIS reaction rate (constant for the whole run).
    user = state.get_user_defined_rate_parameters()
    user["PHOTO.jNO2"] = [J_NO2]
    state.set_user_defined_rate_parameters(user)

    dt_out = 60.0
    t_end = 7200.0
    times = [0.0]
    history = [{k: float(v[0]) for k, v in state.get_concentrations().items()}]

    while times[-1] < t_end:
        elapsed = 0.0
        while elapsed < dt_out:
            r = solver.solve(state, dt_out - elapsed)
            elapsed += r.stats.final_time
            if r.state != SolverState.Converged:
                print(f"  solver state: {r.state} at t={times[-1] + elapsed:.1f}")
        times.append(times[-1] + dt_out)
        history.append({k: float(v[0]) for k, v in state.get_concentrations().items()})

    minutes = np.array(times) / 60.0
    ds = xr.Dataset(
        {sp: ("time", np.array([h[sp] for h in history])) for sp in ("NO", "NO2", "O3")},
        coords={"time": minutes},
        attrs={
            "mechanism": "lnox_o3.yaml",
            "T_K": T_REF,
            "P_Pa": P_REF,
            "jNO2_s-1": J_NO2,
            "note": "Lightning-NOx source omitted (CheMPAS operator-split, not in MICM).",
        },
    )
    ds["time"].attrs["units"] = "minutes"
    for sp in ("NO", "NO2", "O3"):
        ds[sp].attrs["units"] = "mol m-3"
    ds.to_netcdf(OUT_DIR / "lnox_box.nc", engine="scipy")

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    palette = style.get_palette(3)
    for ax, sp, color in zip(axes, ("NO", "NO2", "O3"), palette):
        ax.plot(minutes, ds[sp].values, color=color)
        ax.set_ylabel(f"[{style.species_label(sp)}] [mol m$^{{-3}}$]")
        ax.grid(True, alpha=0.4)
    axes[-1].set_xlabel("Time [min]")
    axes[0].set_title(style.format_title("Standalone LNOx + O3 box model"))
    fig.tight_layout()
    fig.savefig(OUT_DIR / "lnox_box.png", dpi=150)
    print(f"Wrote {OUT_DIR / 'lnox_box.nc'} and {OUT_DIR / 'lnox_box.png'}.")


if __name__ == "__main__":
    main()
