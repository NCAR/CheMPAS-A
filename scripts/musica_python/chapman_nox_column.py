#!/usr/bin/env python3
"""Standalone MUSICA-Python column model for Chapman + NOx PSS.

Mirrors Chapter 3 of the CheMPAS-A tutorial without MPAS in the loop.
Same chapman_nox.yaml MICM mechanism; the column grid is independent
of the MPAS mesh (TUV-x dictates it via the bundled vTS1 calculator).

TUV-x photolysis sourcing
-------------------------
The CheMPAS-A project ships micm_configs/tuvx_chapman_nox.json paired
with chapman_nox.yaml. Loading that JSON directly via the
musica.tuvx.TUVX path-driven constructor would require also building
GridMap / ProfileMap / RadiatorMap (the constructor takes those as
arguments alongside the config path) — essentially reimplementing
mpas_tuvx.F's profile-construction with USSA76 + AFGL data. For
pedagogical simplicity, this script instead reuses MUSICA's bundled
vTS1 calculator (which already provides jO2, jO3->O(3P), jO3->O(1D),
jNO2) and translates TS1's reaction labels to chapman_nox.yaml's
PHOTO.* parameter names via the TS1_LABEL_MAP table below. The
Leighton PSS demonstration is unchanged.

Output: chapman_nox_column.nc, chapman_nox_column.png.

Run:
    ~/miniconda3/envs/mpas/bin/python scripts/musica_python/chapman_nox_column.py
"""
from __future__ import annotations

import sys
from datetime import datetime, time, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

# scripts/ lives one level up; both style.py and init_chapman.py live there.
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import ephem
import ussa1976

import musica
from musica.constants import GAS_CONSTANT
from musica.mechanism_configuration import Parser
from musica.micm.solver_result import SolverState
from musica.tuvx import vTS1

import style
from init_chapman import afgl_qo3_profile, nox_vmr_profile

REPO_ROOT = Path(__file__).parent.parent.parent
MECH_YAML = REPO_ROOT / "micm_configs" / "chapman_nox.yaml"
OUT_DIR = Path(__file__).parent

# Match the supercell case's nominal location.
LAT, LON = 35.86, -97.93
TZ = ZoneInfo("America/Chicago")  # Norman, OK

# JPL kinetics for NO + O3 -> NO2 + O2 (matches §3.7's Leighton expression).
A_NO_O3, EA_R_NO_O3 = 1.7e-12, 1310.0  # cm^3 molec^-1 s^-1, K

# vTS1 reaction labels for chapman_nox.yaml's photolysis reactions.
# Verified against ~/EarthSystem/MUSICA/configs/tuvx/ts1_tsmlt.json
# (the __CAM options.aliasing.pairs table). Note the j*_a / j*_b swap
# relative to chapman_nox.yaml's naming:
#   jo3_a is the O(1D) pathway (matches jO3_O1D)
#   jo3_b is the O(3P) pathway (matches jO3_O)
# Run `print(list(vTS1.get_tuvx_calculator().run(0.0, 1.0).coords['reaction'].values))`
# to inspect TS1's actual reaction labels if these need adjusting.
TS1_LABEL_MAP = {
    "jO2":     "jo2_b",   # O2 + hv -> O + O   (Chapman ground-state)
    "jO3_O":   "jo3_b",   # O3 + hv -> O2 + O(3P)
    "jO3_O1D": "jo3_a",   # O3 + hv -> O2 + O(1D)
    "jNO2":    "jno2",    # NO2 + hv -> NO + O(3P)
}

# Molar masses (match scripts/init_chapman.py).
M_AIR = 0.02897
M_O3 = 0.048

# Daytime NOx partitioning (matches scripts/init_chapman.py).
NO_FRACTION = 0.30


def vmr_to_mol_m3(vmr: np.ndarray, T: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Volume mixing ratio (mol/mol) -> concentration (mol/m^3)."""
    return vmr * P / (GAS_CONSTANT * T)


def kgkg_to_mol_m3(kgkg: np.ndarray, M_species: float,
                   T: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Mass mixing ratio (kg/kg) -> mol/m^3 at given T (K), P (Pa)."""
    return (kgkg * P / (GAS_CONSTANT * T)) / (M_species / M_AIR)


def get_photolysis(tuvx, utc_time, n_cells, start_idx):
    """Compute jO2, jO3_O, jO3_O1D, jNO2 at each column level."""
    # Solar zenith via ephem. Strings for lat/lon are interpreted as
    # degrees; floats would be radians. ephem.Date wants a naive UTC
    # datetime, so strip the tz-info from the already-UTC sim_time.
    obs = ephem.Observer()
    obs.lat, obs.lon = str(LAT), str(LON)
    obs.date = ephem.Date(utc_time.replace(tzinfo=None))
    sun = ephem.Sun()
    sun.compute(obs)
    sza_rad = np.pi / 2 - float(sun.alt)
    sza_deg = float(np.rad2deg(sza_rad))
    rates = tuvx.run(sza=sza_rad, earth_sun_distance=1.0)
    end = start_idx + n_cells
    out = {}
    for our_name, ts1_label in TS1_LABEL_MAP.items():
        out[f"PHOTO.{our_name}"] = (
            rates.sel(reaction=ts1_label).photolysis_rate_constants.values[start_idx:end]
        )
    return out, sza_deg


def main():
    style.apply_ncar_style()

    # 1. TUV-x: vTS1 dictates the column grid. The bundled chapman.py
    # example treats edge altitudes as cell labels (TUV-x photolysis
    # rates are dimensioned on vertical_edge, not vertical_midpoint),
    # so we follow that convention here.
    tuvx = vTS1.get_tuvx_calculator()
    grids = tuvx.get_grid_map()
    z_edges_km = np.asarray(grids["height", "km"].edges)

    # Use a slice of the TS1 grid: skip the surface cell (start=1), span up to 60 km.
    start = 1
    n_cells = int(np.searchsorted(z_edges_km, 60.0)) - start
    z_cells_km = z_edges_km[start:start + n_cells]

    # 2. USSA76 T, P at column cell altitudes.
    env = ussa1976.compute(z=z_cells_km * 1000.0, variables=["t", "p"])
    T_K = env["t"].values
    P_Pa = env["p"].values

    # 3. Initial profile from scripts/init_chapman.py helpers.
    qo3_kgkg = afgl_qo3_profile(z_cells_km)
    nox_vmr = nox_vmr_profile(z_cells_km)
    initial_concs = {
        "O2":  vmr_to_mol_m3(np.full(n_cells, 0.20946), T_K, P_Pa),
        "O":   np.zeros(n_cells),
        "O1D": np.zeros(n_cells),
        "O3":  kgkg_to_mol_m3(qo3_kgkg, M_O3, T_K, P_Pa),
        "NO":  vmr_to_mol_m3(NO_FRACTION * nox_vmr, T_K, P_Pa),
        "NO2": vmr_to_mol_m3((1.0 - NO_FRACTION) * nox_vmr, T_K, P_Pa),
    }

    # 4. MICM solver setup.
    parser = Parser()
    mechanism = parser.parse(str(MECH_YAML))
    solver = musica.MICM(
        mechanism=mechanism,
        solver_type=musica.SolverType.rosenbrock_standard_order,
    )
    state = solver.create_state(number_of_grid_cells=n_cells)
    state.set_conditions(T_K, P_Pa)
    state.set_concentrations({sp: list(v) for sp, v in initial_concs.items()})

    # 5. 12-hour run starting at 06:00 local, 30-min cadence.
    today = datetime.now(TZ).date()
    sim_time = datetime.combine(today, time(6, 0), tzinfo=TZ).astimezone(ZoneInfo("UTC"))
    dt_out = 30 * 60
    n_steps = 12 * 2  # 12 hours, 30-min cadence

    photo, sza = get_photolysis(tuvx, sim_time, n_cells, start)
    user = state.get_user_defined_rate_parameters()
    for k, v in photo.items():
        user[k] = list(v)
    state.set_user_defined_rate_parameters(user)

    times = [sim_time]
    species_hist = [state.get_concentrations()]
    photo_hist = [photo]
    sza_hist = [sza]

    for _ in range(n_steps):
        elapsed = 0.0
        while elapsed < dt_out:
            r = solver.solve(state, dt_out - elapsed)
            elapsed += r.stats.final_time
            if r.state != SolverState.Converged:
                print(f"  solver state: {r.state} at sim_time={sim_time}")
        sim_time += timedelta(seconds=dt_out)
        photo, sza = get_photolysis(tuvx, sim_time, n_cells, start)
        user = state.get_user_defined_rate_parameters()
        for k, v in photo.items():
            user[k] = list(v)
        state.set_user_defined_rate_parameters(user)
        times.append(sim_time)
        species_hist.append(state.get_concentrations())
        photo_hist.append(photo)
        sza_hist.append(sza)

    # 6. Build xarray Dataset.
    species_names = ("O2", "O", "O1D", "O3", "NO", "NO2")
    photo_names = ("jO2", "jO3_O", "jO3_O1D", "jNO2")
    species_arr = {
        sp: np.array([[h[sp][i] for i in range(n_cells)] for h in species_hist])
        for sp in species_names
    }
    photo_arr = {
        pn: np.array([p[f"PHOTO.{pn}"] for p in photo_hist])
        for pn in photo_names
    }
    data_vars = {
        sp: (["time", "height"], species_arr[sp], {"units": "mol m-3"})
        for sp in species_names
    }
    for pn in photo_names:
        data_vars[pn] = (["time", "height"], photo_arr[pn], {"units": "s-1"})
    data_vars["sza_deg"] = (["time"], np.array(sza_hist), {"units": "degrees"})

    ds = xr.Dataset(
        data_vars,
        coords={
            # Drop tz-info (convert to naive UTC) so np.array doesn't deprecate-warn.
            "time": np.array(
                [t.astimezone(ZoneInfo("UTC")).replace(tzinfo=None) for t in times],
                dtype="datetime64[ns]",
            ),
            "height": z_cells_km,
        },
        attrs={
            "mechanism": "chapman_nox.yaml",
            "tuvx": "vTS1 (custom alias map for chapman_nox.yaml reactions)",
            "lat_deg": LAT, "lon_deg": LON, "tz": str(TZ),
        },
    )
    ds["height"].attrs["units"] = "km"
    ds.to_netcdf(OUT_DIR / "chapman_nox_column.nc", engine="scipy")

    # 7. 2x2 plot (uses scripts/style.py).
    palette = style.get_palette(4)
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2)

    # Solar noon = min SZA in the run.
    noon_idx = int(np.argmin(sza_hist))

    # Top-left: O3 profile at noon (ppb).
    ax_o3 = fig.add_subplot(gs[0, 0])
    o3_ppb_noon = ds["O3"].isel(time=noon_idx).values * GAS_CONSTANT * T_K / P_Pa * 1e9
    ax_o3.plot(o3_ppb_noon, z_cells_km, color=palette[0])
    ax_o3.set_xlabel(f"[{style.species_label('O3')}] [ppb]")
    ax_o3.set_ylabel("Height [km]")
    ax_o3.set_title(style.format_title("Solar-noon O3 profile"))
    ax_o3.grid(True, alpha=0.4)

    # Top-right: NO and NO2 profiles at noon.
    ax_no = fig.add_subplot(gs[0, 1])
    no_ppb = ds["NO"].isel(time=noon_idx).values * GAS_CONSTANT * T_K / P_Pa * 1e9
    no2_ppb = ds["NO2"].isel(time=noon_idx).values * GAS_CONSTANT * T_K / P_Pa * 1e9
    ax_no.plot(no_ppb, z_cells_km, color=palette[1], label=style.species_label("NO"))
    ax_no.plot(no2_ppb, z_cells_km, color=palette[2], label=style.species_label("NO2"))
    ax_no.set_xlabel("Concentration [ppb]")
    ax_no.set_ylabel("Height [km]")
    ax_no.set_title(style.format_title("Solar-noon NOx profile"))
    ax_no.legend()
    ax_no.grid(True, alpha=0.4)

    # Bottom-left: simulated vs analytical Leighton ratio at noon.
    ax_lt = fig.add_subplot(gs[1, 0])
    sim_ratio = (ds["NO"].isel(time=noon_idx).values
                 / np.maximum(ds["NO2"].isel(time=noon_idx).values, 1e-30))
    # Convert [O3] from mol/m^3 to molec/cm^3 so units match the JPL k.
    AVOGADRO = 6.022e23
    o3_molec_cm3 = ds["O3"].isel(time=noon_idx).values * AVOGADRO * 1e-6
    k_NO_O3 = A_NO_O3 * np.exp(-EA_R_NO_O3 / T_K)
    leighton = (ds["jNO2"].isel(time=noon_idx).values
                / np.maximum(k_NO_O3 * o3_molec_cm3, 1e-30))
    ax_lt.plot(sim_ratio, z_cells_km, color=palette[0], label="Simulated")
    ax_lt.plot(leighton, z_cells_km, color=palette[2], linestyle="--", label="Leighton")
    ax_lt.set_xlabel(f"[{style.species_label('NO')}]/[{style.species_label('NO2')}]")
    ax_lt.set_ylabel("Height [km]")
    ax_lt.set_title(style.format_title("Solar-noon Leighton check"))
    ax_lt.set_xscale("log")
    ax_lt.legend()
    ax_lt.grid(True, alpha=0.4)

    # Bottom-right: O3 time series at three altitudes.
    ax_ts = fig.add_subplot(gs[1, 1])
    target_alts_km = [10.0, 30.0, 45.0]
    times_arr = ds["time"].values
    hours_since_start = ((times_arr - times_arr[0]).astype("timedelta64[s]")
                         .astype(float)) / 3600.0
    for alt, color in zip(target_alts_km, palette):
        idx = int(np.argmin(np.abs(z_cells_km - alt)))
        o3_at_alt_ppb = (ds["O3"].isel(height=idx).values
                         * GAS_CONSTANT * T_K[idx] / P_Pa[idx] * 1e9)
        ax_ts.plot(hours_since_start, o3_at_alt_ppb,
                   color=color, label=f"{z_cells_km[idx]:.1f} km")
    ax_ts.set_xlabel("Hours since start")
    ax_ts.set_ylabel(f"[{style.species_label('O3')}] [ppb]")
    ax_ts.set_title(style.format_title("O3 time series"))
    ax_ts.legend()
    ax_ts.grid(True, alpha=0.4)

    fig.suptitle(style.format_title("Chapman + NOx column model"))
    fig.tight_layout()
    fig.savefig(OUT_DIR / "chapman_nox_column.png", dpi=150)
    print(f"Wrote {OUT_DIR / 'chapman_nox_column.nc'} and "
          f"{OUT_DIR / 'chapman_nox_column.png'}.")


if __name__ == "__main__":
    main()
