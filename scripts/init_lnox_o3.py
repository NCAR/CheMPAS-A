#!/usr/bin/env python3
"""
Initialize LNOx-O3 tracers (qNO, qNO2, qO3) in an MPAS init file.

Creates the tracer variables if they don't exist, sets:
  qNO  = 0        (no initial NO)
  qNO2 = 0        (no initial NO2)
  qO3  = 50 ppbv  (background ozone, converted to kg/kg)

Usage:
    python init_lnox_o3.py                           # Default: supercell_init.nc
    python init_lnox_o3.py -i supercell_init.nc      # Explicit path
    python init_lnox_o3.py --o3-ppbv 30              # Custom O3 background
"""

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# Molar masses [kg/mol]
M_NO = 0.030
M_NO2 = 0.046
M_O3 = 0.048
M_AIR = 0.029


def ppbv_to_kgkg(ppbv, M_species):
    """Convert ppbv to mass mixing ratio (kg/kg)."""
    return ppbv * 1.0e-9 * (M_species / M_AIR)


def ensure_tracer(ds, name, value, long_name=None):
    """Create or overwrite a tracer variable with a uniform value."""
    nTimes = ds.dimensions['Time'].size
    nCells = ds.dimensions['nCells'].size
    nVertLevels = ds.dimensions['nVertLevels'].size

    if name not in ds.variables:
        dtype = ds.variables['qv'].dtype if 'qv' in ds.variables else 'f8'
        var = ds.createVariable(name, dtype, ('Time', 'nCells', 'nVertLevels'))
        var.units = 'kg kg^{-1}'
        var.long_name = long_name or name
        print(f"  Created {name}")
    else:
        var = ds.variables[name]
        print(f"  Found existing {name}")

    var[:] = np.full((nTimes, nCells, nVertLevels), value)
    print(f"  Set {name} = {value:.6e} kg/kg")
    return var


def main():
    parser = argparse.ArgumentParser(
        description='Initialize LNOx-O3 tracers in MPAS init file')
    parser.add_argument('-i', '--input', default='supercell_init.nc',
                        help='Path to init file (default: supercell_init.nc)')
    parser.add_argument('--o3-ppbv', type=float, default=50.0,
                        help='Background O3 in ppbv (default: 50)')
    args = parser.parse_args()

    filepath = Path(args.input)
    if not filepath.exists():
        raise SystemExit(f"File not found: {filepath}")

    o3_kgkg = ppbv_to_kgkg(args.o3_ppbv, M_O3)

    print(f"Initializing LNOx-O3 tracers in {filepath}")
    print(f"  O3 background: {args.o3_ppbv} ppbv = {o3_kgkg:.6e} kg/kg")

    with Dataset(filepath, 'r+') as ds:
        ensure_tracer(ds, 'qNO', 0.0, 'nitric_oxide')
        ensure_tracer(ds, 'qNO2', 0.0, 'nitrogen_dioxide')
        ensure_tracer(ds, 'qO3', o3_kgkg, 'ozone')
        ds.sync()

    print("Done.")


if __name__ == '__main__':
    main()
