#!/usr/bin/env python3
"""Vertical profile plots for the TUV-x column extension work.

Three panels, all on a shared altitude axis:
  1. Temperature profile [K]
  2. Air number density [molec/cm^3], log scale
  3. O3 number density [molec/cm^3], log scale

Each panel stacks:
  - MPAS state from an output.nc file (0 to MPAS top, usually 50 km),
  - the TUV-x extension climatology CSV above (usually 50 to 100 km).

By default uses a clear-sky-ish column from the run directory to
minimise storm contamination in the profile.

Usage:
    plot_extension_profiles.py -i ~/Data/CheMPAS/supercell/output.nc \
        --csv micm_configs/tuvx_upper_atm.csv \
        -o extension_profiles.png
"""
import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

NA = 6.02214076e23
M_AIR = 0.02897
M_O3 = 0.04800


def read_extension_csv(path: str):
    z_km, T_K, n_air, n_o3 = [], [], [], []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            z_km.append(float(row["z_km"]))
            T_K.append(float(row["T_K"]))
            n_air.append(float(row["n_air_molec_cm3"]))
            n_o3.append(float(row["n_O3_molec_cm3"]))
    return (np.array(z_km), np.array(T_K), np.array(n_air), np.array(n_o3))


def mpas_column(nc_path: str, cell_idx: int, time_idx: int):
    """Return z_mid_km, T, n_air, n_o3 for a single MPAS column."""
    with Dataset(nc_path) as ds:
        zgrid = ds.variables["zgrid"][cell_idx, :]  # nVertLevels+1 edges [m]
        rho = ds.variables["rho"][time_idx, cell_idx, :]  # kg/m^3
        theta = ds.variables["theta"][time_idx, cell_idx, :]  # K (potential)
        pressure = ds.variables["pressure"][time_idx, cell_idx, :]  # Pa
        qO3 = ds.variables["qO3"][time_idx, cell_idx, :]  # kg/kg
    # Temperature from potential temperature and pressure:  T = theta * (p/p0)^(R/cp)
    T = theta * (pressure / 1.0e5) ** (287.0 / 1004.0)
    z_mid_km = 0.5 * (zgrid[:-1] + zgrid[1:]) / 1000.0
    n_air = rho * (NA / M_AIR) * 1.0e-6          # molec/cm^3
    n_o3 = qO3 * rho * (NA / M_O3) * 1.0e-6      # molec/cm^3
    return z_mid_km, T, n_air, n_o3


def ext_midpoints(z_edges, y_edges):
    """Midpoint-average the CSV edges the same way mpas_tuvx.F does."""
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])
    y_mid = 0.5 * (y_edges[:-1] + y_edges[1:])
    return z_mid, y_mid


def pick_clear_cell(nc_path: str, time_idx: int) -> int:
    """Pick a cell with minimum total cloud/rain water at the given time."""
    with Dataset(nc_path) as ds:
        qc = ds.variables["qc"][time_idx, :, :] if "qc" in ds.variables else None
        qr = ds.variables["qr"][time_idx, :, :] if "qr" in ds.variables else None
        if qc is None and qr is None:
            return 0
        total = np.zeros(ds.dimensions["nCells"].size)
        if qc is not None:
            total += qc.sum(axis=1)
        if qr is not None:
            total += qr.sum(axis=1)
    return int(np.argmin(total))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default=str(Path.home() / "Data/CheMPAS/supercell/output.nc"))
    ap.add_argument("--csv", default="micm_configs/tuvx_upper_atm.csv")
    ap.add_argument("-o", "--output", default="extension_profiles.png")
    ap.add_argument("-t", "--time", type=int, default=-1, help="time index (default: last)")
    ap.add_argument("-c", "--cell", type=int, default=None,
                    help="cell index; default = clearest-sky cell at --time")
    args = ap.parse_args()

    # Resolve time index
    with Dataset(args.input) as ds:
        nt = ds.dimensions["Time"].size
    tidx = args.time if args.time >= 0 else nt + args.time

    cell = args.cell if args.cell is not None else pick_clear_cell(args.input, tidx)

    z_mpas, T_mpas, nair_mpas, no3_mpas = mpas_column(args.input, cell, tidx)
    z_ext_edges, T_ext_e, nair_ext_e, no3_ext_e = read_extension_csv(args.csv)

    # Extension midpoints (the values TUV-x uses internally)
    z_ext_mid, T_ext = ext_midpoints(z_ext_edges, T_ext_e)
    _, nair_ext = ext_midpoints(z_ext_edges, nair_ext_e)
    _, no3_ext = ext_midpoints(z_ext_edges, no3_ext_e)

    fig, axes = plt.subplots(1, 3, figsize=(13, 7), sharey=True)

    titles = ["Temperature [K]", "Air number density [molec cm$^{-3}$]",
              "O$_3$ number density [molec cm$^{-3}$]"]
    mpas_data = [T_mpas, nair_mpas, no3_mpas]
    ext_mid_data = [T_ext, nair_ext, no3_ext]
    ext_edge_data = [T_ext_e, nair_ext_e, no3_ext_e]
    logx = [False, True, True]

    for ax, title, mp, emid, eedge, lx in zip(axes, titles, mpas_data,
                                               ext_mid_data, ext_edge_data, logx):
        ax.plot(mp, z_mpas, "o-", color="C0", ms=3, lw=1, label="MPAS")
        ax.plot(emid, z_ext_mid, "s-", color="C3", ms=4, lw=1,
                label="Extension midpoints")
        ax.plot(eedge, z_ext_edges, "x", color="C3", ms=6, alpha=0.4,
                label="Extension CSV edges")
        ax.axhline(z_mpas[-1] + 0.5 * (z_mpas[-1] - z_mpas[-2]), ls="--",
                   color="gray", alpha=0.6, label="MPAS top")
        ax.set_xlabel(title)
        ax.grid(True, alpha=0.3)
        if lx:
            ax.set_xscale("log")

    axes[0].set_ylabel("Altitude [km]")
    axes[0].legend(loc="best", fontsize=9)
    axes[0].set_ylim(0, 105)

    fig.suptitle(
        f"Composite column seen by TUV-x (cell {cell}, time idx {tidx})",
        fontsize=12,
    )
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    pdf = Path(args.output).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {args.output} and {pdf}")
    print(f"cell = {cell}, time idx = {tidx}")
    print(f"MPAS top midpoint = {z_mpas[-1]:.2f} km;  extension spans "
          f"{z_ext_edges[0]:.1f}-{z_ext_edges[-1]:.1f} km ({len(z_ext_mid)} layers)")


if __name__ == "__main__":
    main()
