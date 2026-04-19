#!/usr/bin/env python3
"""Initialize Chapman-cycle tracers in an MPAS init file.

Seeds the chapman_full.yaml and chapman_nox.yaml mechanism tracers with
realistic vertical profiles spanning the full MPAS domain (surface to
~50 km), continuous with the TUV-x upper-atmosphere climatology at the
MPAS top.

Profiles:
  qO2  uniform 0.2313 kg/kg  (air * 0.2095 by volume converted to mass)
  qO3  AFGL mid-latitude-summer profile, interpolated to MPAS grid
  qO   0   (fast radical; chemistry spins up)
  qO1D 0
  qNO  realistic NOx profile with daytime partitioning (~30% NO)
  qNO2 realistic NOx profile (~70% NO2 at 18 UTC)

NOx profile (total NOx = NO + NO2, approximate mid-latitude daytime):
  surface - 10 km:  0.05 ppb   (remote tropospheric background)
  10 - 20 km:       log-linear rise
  20 - 30 km:       10 ppb     (stratospheric N2O-derived peak plateau)
  30 - 40 km:       10 ppb
  40 - 50 km:       log-linear drop to 0.5 ppb at the MPAS top

NO/NO2 daytime split is 30%/70% (close to the Leighton steady state
at ~30 km with jNO2 ≈ 0.01 s^-1, [O3] ≈ 5 ppm); chemistry will
partition within seconds once the solver runs.

Usage:
    python init_chapman.py                          # Default: supercell_init.nc
    python init_chapman.py --nox-peak 15            # Scale the NOx peak
    python init_chapman.py --zero-nox               # Don't seed NOx (Chapman only)
"""
import argparse
import math
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# Molar masses [kg/mol]
M_AIR = 0.02897
M_O2  = 0.032
M_O3  = 0.048
M_O   = 0.016
M_O1D = 0.016
M_NO  = 0.030
M_NO2 = 0.046

# Uniform O2 mass mixing ratio corresponding to 0.2095 by volume in dry air.
QO2_UNIFORM = 0.2095 * (M_O2 / M_AIR)


# AFGL mid-latitude-summer O3 profile (molec/cm^3) sampled at key altitudes.
# Source: Anderson et al. 1986 AFGL atmospheric constituent profiles.
AFGL_MLS_O3 = [
    ( 0.0, 6.80e+11),
    ( 5.0, 5.10e+11),
    (10.0, 1.20e+12),
    (15.0, 2.80e+12),
    (20.0, 5.20e+12),
    (25.0, 4.10e+12),
    (30.0, 2.00e+12),
    (35.0, 8.50e+11),
    (40.0, 3.20e+11),
    (45.0, 1.30e+11),
    (50.0, 6.45e+11),   # matches the extension CSV at 50 km
]


# Matching air density [molec/cm^3] at the same altitudes.
AFGL_MLS_AIR = [
    ( 0.0, 2.55e+19),
    ( 5.0, 1.53e+19),
    (10.0, 8.60e+18),
    (15.0, 4.04e+18),
    (20.0, 1.86e+18),
    (25.0, 8.33e+17),
    (30.0, 3.83e+17),
    (35.0, 1.76e+17),
    (40.0, 8.31e+16),
    (45.0, 4.09e+16),
    (50.0, 2.13e+16),
]


# NOx total volume mixing ratio [ppb] as a function of altitude [km].
# Approximate mid-latitude daytime climatology: low tropospheric
# background, rise through the tropopause, broad stratospheric peak at
# 25-35 km (~10 ppb), decay toward the MPAS top. Values are scaled by
# --nox-peak on the command line; defaults match figure 4.10 of Brasseur
# & Solomon (2005) "Aeronomy of the Middle Atmosphere".
NOX_PROFILE_PPB = [
    ( 0.0, 0.05),
    (10.0, 0.05),
    (15.0, 0.5),
    (20.0, 3.0),
    (25.0, 10.0),
    (30.0, 10.0),
    (35.0, 10.0),
    (40.0, 5.0),
    (45.0, 2.0),
    (50.0, 0.5),   # sets continuity with the (near-zero NOx) extension
]

# Daytime NO fraction of NOx (steady-state Leighton approximation).
NO_FRACTION = 0.30


def loglin(z, table):
    """Log-linear interpolate a strictly-positive quantity."""
    if z <= table[0][0]:
        return table[0][1]
    if z >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        z0, y0 = table[i]
        z1, y1 = table[i + 1]
        if z0 <= z <= z1:
            t = (z - z0) / (z1 - z0)
            return math.exp(math.log(y0) + t * (math.log(y1) - math.log(y0)))
    raise RuntimeError("interpolation out of range")


def afgl_qo3_profile(z_mid_km):
    """Return qO3 (kg/kg) at the given MPAS midpoint altitudes."""
    out = np.zeros_like(z_mid_km, dtype=float)
    for i, z in enumerate(z_mid_km):
        n_o3 = loglin(z, AFGL_MLS_O3)
        n_air = loglin(z, AFGL_MLS_AIR)
        out[i] = (n_o3 / n_air) * (M_O3 / M_AIR)
    return out


def nox_vmr_profile(z_mid_km, peak_scale=1.0):
    """Return total NOx volume mixing ratio (dimensionless, e.g. 1e-9 = 1 ppb)
    at the given MPAS midpoint altitudes. Scaled by peak_scale."""
    out = np.zeros_like(z_mid_km, dtype=float)
    for i, z in enumerate(z_mid_km):
        out[i] = loglin(z, NOX_PROFILE_PPB) * peak_scale * 1.0e-9
    return out


def ensure_tracer(ds, name, values, long_name=None):
    nTimes = ds.dimensions["Time"].size
    nCells = ds.dimensions["nCells"].size
    nVertLevels = ds.dimensions["nVertLevels"].size

    if name not in ds.variables:
        dtype = ds.variables["qv"].dtype if "qv" in ds.variables else "f8"
        var = ds.createVariable(name, dtype, ("Time", "nCells", "nVertLevels"))
        var.units = "kg kg^{-1}"
        var.long_name = long_name or name
        print(f"  Created {name}")
    else:
        var = ds.variables[name]
        print(f"  Found existing {name}")

    arr = np.asarray(values)
    if arr.ndim == 0:
        var[:] = np.full((nTimes, nCells, nVertLevels), float(arr))
        print(f"  Set {name} = {float(arr):.6e} kg/kg (uniform)")
    elif arr.ndim == 1 and arr.size == nVertLevels:
        var[:] = np.broadcast_to(arr, (nTimes, nCells, nVertLevels))
        print(f"  Set {name} from 1-D vertical profile "
              f"(min={arr.min():.3e}, max={arr.max():.3e})")
    else:
        raise ValueError(
            f"{name}: expected scalar or length-{nVertLevels} array, got {arr.shape}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="supercell_init.nc")
    ap.add_argument("--nox-peak", type=float, default=1.0,
                    help="Scale factor on the stratospheric NOx peak (default 1.0 = ~10 ppb)")
    ap.add_argument("--zero-nox", action="store_true",
                    help="Seed qNO = qNO2 = 0 (for chapman_full.yaml without NOx)")
    args = ap.parse_args()

    path = Path(args.input)
    if not path.exists():
        raise SystemExit(f"File not found: {path}")

    with Dataset(path, "r+") as ds:
        zgrid = ds.variables["zgrid"][:]
        if zgrid.ndim == 2:
            edges_m = zgrid[0, :]
        else:
            raise RuntimeError(f"Unexpected zgrid shape {zgrid.shape}")
        z_mid_km = 0.5 * (edges_m[:-1] + edges_m[1:]) / 1000.0

        print(f"Initializing Chapman tracers in {path}")
        print(f"  MPAS top = {edges_m[-1]/1000:.1f} km, nVertLevels = {len(z_mid_km)}")

        ensure_tracer(ds, "qO2", QO2_UNIFORM, "molecular_oxygen")
        ensure_tracer(ds, "qO3", afgl_qo3_profile(z_mid_km), "ozone")
        # Seed [O] with a small non-zero floor (instead of exactly 0.0) to
        # break the degeneracy that lets MICM's implicit solver land on a
        # non-physical fixed point when photolysis is on from t=0. 1e-12
        # kg/kg corresponds to [O] ≈ 7e-11 mol/m³ at surface and ≈ 3e-10
        # mol/m³ at the stratopause — many orders above the MICM absolute
        # tolerance (1e-10 mol/m³), well below any chemical steady state,
        # and negligible for every rate we care about. See
        # docs/plans/2026-04-18-chapman-nox-chem-box-issue.md.
        ensure_tracer(ds, "qO",  1.0e-12, "atomic_oxygen_3P")
        ensure_tracer(ds, "qO1D", 0.0, "atomic_oxygen_1D")

        if args.zero_nox:
            ensure_tracer(ds, "qNO",  0.0, "nitric_oxide")
            ensure_tracer(ds, "qNO2", 0.0, "nitrogen_dioxide")
        else:
            vmr_nox = nox_vmr_profile(z_mid_km, peak_scale=args.nox_peak)
            ensure_tracer(ds, "qNO",  vmr_nox * NO_FRACTION * (M_NO / M_AIR),
                          "nitric_oxide")
            ensure_tracer(ds, "qNO2", vmr_nox * (1.0 - NO_FRACTION) * (M_NO2 / M_AIR),
                          "nitrogen_dioxide")
            print(f"  NOx peak scale: {args.nox_peak:.2f}x "
                  f"(peak total NOx = {10.0 * args.nox_peak:.1f} ppb at 25-35 km)")

        ds.sync()

    print("Done.")


if __name__ == "__main__":
    main()
