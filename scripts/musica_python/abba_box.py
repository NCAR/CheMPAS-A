#!/usr/bin/env python3
"""Standalone MUSICA-Python box model for the ABBA toy mechanism.

Mirrors section 2.5 of the CheMPAS-A tutorial (the supercell + ABBA
run) without MPAS in the loop. Single grid cell, slow two-way reaction
qAB <-> qA + qB, 2 h integration. Output: abba_box.nc, abba_box.png.

Run:
    ~/miniconda3/envs/mpas/bin/python scripts/musica_python/abba_box.py
"""
from __future__ import annotations

import sys
from pathlib import Path

# scripts/style.py lives one level up.
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import musica
from musica.mechanism_configuration import Parser
from musica.micm.solver_result import SolverState

import style

REPO_ROOT = Path(__file__).parent.parent.parent
ABBA_YAML = REPO_ROOT / "micm_configs" / "abba.yaml"
OUT_DIR = Path(__file__).parent

# Reference conditions for the box.
T_REF = 273.0     # K
P_REF = 101325.0  # Pa


def main() -> None:
    style.apply_ncar_style()

    parser = Parser()
    mechanism = parser.parse(str(ABBA_YAML))
    solver = musica.MICM(
        mechanism=mechanism,
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(number_of_grid_cells=1)
    state.set_conditions(temperatures=T_REF, pressures=P_REF)
    state.set_concentrations({"A": [0.0], "B": [0.0], "AB": [1.0]})

    # USER_DEFINED reaction rates (scaling factors from abba.yaml).
    user = state.get_user_defined_rate_parameters()
    user["USER.forward_AB_to_A_B"] = [2.0e-3]
    user["USER.reverse_A_B_to_AB"] = [1.0e-3]
    state.set_user_defined_rate_parameters(user)

    dt_out = 60.0       # output cadence (s)
    t_end  = 7200.0     # 2 h
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
        {sp: ("time", np.array([h[sp] for h in history])) for sp in ("A", "B", "AB")},
        coords={"time": minutes},
        attrs={"mechanism": "abba.yaml", "T_K": T_REF, "P_Pa": P_REF},
    )
    ds["time"].attrs["units"] = "minutes"
    for sp in ("A", "B", "AB"):
        ds[sp].attrs["units"] = "mol m-3"
    ds.to_netcdf(OUT_DIR / "abba_box.nc", engine="scipy")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    palette = style.get_palette(3)
    for sp, color in zip(("AB", "A", "B"), palette):
        ax.plot(minutes, ds[sp].values, color=color, label=style.species_label(sp))
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Concentration [mol m$^{-3}$]")
    ax.set_title(style.format_title("Standalone ABBA box model"))
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abba_box.png", dpi=150)
    print(f"Wrote {OUT_DIR / 'abba_box.nc'} and {OUT_DIR / 'abba_box.png'}.")


if __name__ == "__main__":
    main()
