#!/usr/bin/env python3
"""Initialize Chapman-cycle tracers in an MPAS init file.

Seeds the five species the chapman_full.yaml mechanism tracks:
  qO2  uniform 0.2313 kg/kg  (air * 0.2095 by volume converted to mass)
  qO3  AFGL mid-latitude-summer profile interpolated to MPAS vertical grid
  qO   0
  qO1D 0

Optional --seed-nox flag (used for the chapman_nox.yaml follow-on) adds
a small stratospheric NOx layer (1 ppbv uniform between 20 and 50 km).

Usage:
    python init_chapman.py                          # Default: supercell_init.nc
    python init_chapman.py --seed-nox               # Plus 1 ppbv stratospheric NOx
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
# Values are approximate — a single log-linear interpolation is adequate
# for an idealized spin-up seed profile.
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


# Matching air density [molec/cm^3] at the same altitudes — needed to
# convert O3 number density back to mass mixing ratio.
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


def loglin(z, table):
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
        # qO3 = (n_o3 / n_air) * (M_O3 / M_AIR)
        out[i] = (n_o3 / n_air) * (M_O3 / M_AIR)
    return out


def stratospheric_nox_profile(z_mid_km, z_min=20.0, z_max=50.0, ppbv=1.0):
    """1 ppbv uniform between z_min and z_max, zero elsewhere (mass ratio)."""
    out = np.zeros_like(z_mid_km, dtype=float)
    # 1 ppbv is a volume mixing ratio; convert to mass for each species.
    # Mass mixing ratio = ppbv * 1e-9 * (M_species / M_AIR)
    in_layer = (z_mid_km >= z_min) & (z_mid_km <= z_max)
    out[in_layer] = ppbv * 1.0e-9   # caller multiplies by M_species/M_AIR
    return out


def ensure_tracer(ds, name, values, long_name=None):
    """Create or overwrite a tracer with either a scalar or a (Time,nCells,nVertLevels) array."""
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
        raise ValueError(f"{name}: expected scalar or length-{nVertLevels} array, got {arr.shape}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="supercell_init.nc")
    ap.add_argument("--seed-nox", action="store_true",
                    help="Seed 1 ppbv NOx between 20 and 50 km (for chapman_nox)")
    ap.add_argument("--nox-ppbv", type=float, default=1.0)
    args = ap.parse_args()

    path = Path(args.input)
    if not path.exists():
        raise SystemExit(f"File not found: {path}")

    with Dataset(path, "r+") as ds:
        # MPAS stores zgrid at edges, shape (nCells, nVertLevels+1). Use the
        # first cell's profile since the supercell case is flat-terrain.
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
        ensure_tracer(ds, "qO",  0.0, "atomic_oxygen_3P")
        ensure_tracer(ds, "qO1D", 0.0, "atomic_oxygen_1D")

        if args.seed_nox:
            vmr = stratospheric_nox_profile(z_mid_km, ppbv=args.nox_ppbv)
            ensure_tracer(ds, "qNO",  vmr * (M_NO  / M_AIR), "nitric_oxide")
            ensure_tracer(ds, "qNO2", vmr * (M_NO2 / M_AIR), "nitrogen_dioxide")

        ds.sync()

    print("Done.")


if __name__ == "__main__":
    main()
