#!/usr/bin/env python3
"""Vertical profile plots for the TUV-x column extension work.

Three panels, all on a shared altitude axis:
  1. Temperature [K]
  2. Air number density [molec cm^-3], log scale
  3. O3 number density [molec cm^-3], log scale

Each panel shows the composite column TUV-x actually sees, stitched
the same way mpas_tuvx.F does it: MPAS midpoint values for the lower
slice, extension-CSV midpoint values above, with edge-blending at the
MPAS/extension boundary so the profile is continuous. Colored segments
distinguish the MPAS region from the extension. The raw CSV edge
values are overplotted as markers for reference.

Usage:
    plot_extension_profiles.py                              # run-dir defaults
    plot_extension_profiles.py -i output.nc -c 12345        # specific cell
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Reuse the project-wide NCAR style module
sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402


NA = 6.02214076e23
M_AIR = 0.02897
M_O3 = 0.04800


def read_extension_csv(path):
    z_km, T_K, n_air, n_o3 = [], [], [], []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            z_km.append(float(row["z_km"]))
            T_K.append(float(row["T_K"]))
            n_air.append(float(row["n_air_molec_cm3"]))
            n_o3.append(float(row["n_O3_molec_cm3"]))
    return (np.array(z_km), np.array(T_K), np.array(n_air), np.array(n_o3))


def mpas_column(nc_path, cell_idx, time_idx):
    """MPAS midpoint profile for a single column."""
    with Dataset(nc_path) as ds:
        edges_m = ds.variables["zgrid"][cell_idx, :]
        rho = ds.variables["rho"][time_idx, cell_idx, :]
        theta = ds.variables["theta"][time_idx, cell_idx, :]
        pressure = ds.variables["pressure"][time_idx, cell_idx, :]
        qO3 = ds.variables["qO3"][time_idx, cell_idx, :]
    T = theta * (pressure / 1.0e5) ** (287.0 / 1004.0)
    z_mid_km = 0.5 * (edges_m[:-1] + edges_m[1:]) / 1000.0
    n_air = rho * (NA / M_AIR) * 1.0e-6
    n_o3 = qO3 * rho * (NA / M_O3) * 1.0e-6
    return z_mid_km, edges_m / 1000.0, T, n_air, n_o3


def composite_column(z_mpas_mid_km, mpas_values, z_ext_edges, ext_edges_values):
    """Return the (z_mid, value) composite column TUV-x actually sees.

    Matches the stitching mpas_tuvx.F performs in tuvx_compute_photolysis:
    MPAS midpoints on the lower slice, extension midpoints (averages of
    adjacent CSV edges) on the upper slice, with a single shared
    altitude axis. Returns one continuous line so the plot shows exactly
    what the radiative transfer sees.
    """
    z_ext_mid = 0.5 * (z_ext_edges[:-1] + z_ext_edges[1:])
    ext_mid_values = 0.5 * (ext_edges_values[:-1] + ext_edges_values[1:])
    return (
        np.concatenate([z_mpas_mid_km, z_ext_mid]),
        np.concatenate([mpas_values, ext_mid_values]),
    )


def pick_clear_cell(nc_path, time_idx):
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
    ap.add_argument("-o", "--output", default=None,
                    help="output PNG path; default <input_dir>/plots/extension_profiles.png")
    ap.add_argument("-t", "--time", type=int, default=-1)
    ap.add_argument("-c", "--cell", type=int, default=None,
                    help="cell index; default = clearest-sky cell at --time")
    args = ap.parse_args()

    if args.output is None:
        plots_dir = Path(args.input).parent / "plots"
        plots_dir.mkdir(exist_ok=True)
        args.output = str(plots_dir / "extension_profiles.png")

    style.setup()

    with Dataset(args.input) as ds:
        nt = ds.dimensions["Time"].size
    tidx = args.time if args.time >= 0 else nt + args.time
    cell = args.cell if args.cell is not None else pick_clear_cell(args.input, tidx)

    z_mpas_mid, z_mpas_edges, T_mpas, nair_mpas, no3_mpas = mpas_column(
        args.input, cell, tidx)
    z_ext_edges, T_ext_e, nair_ext_e, no3_ext_e = read_extension_csv(args.csv)

    # T-anchor: mpas_tuvx.F shifts the extension T profile so its first
    # midpoint matches the MPAS top-midpoint T. Replicate the same shift
    # here so the plot shows what TUV-x actually sees (otherwise we
    # display an ~30 K step that is not in the RT column).
    ext_first_mid_T_raw = 0.5 * (T_ext_e[0] + T_ext_e[1])
    T_anchor_offset = T_mpas[-1] - ext_first_mid_T_raw
    T_ext_anchored = T_ext_e + T_anchor_offset

    fig, axes = plt.subplots(1, 3, figsize=(14, 7), sharey=True)

    mpas_color = style.NCAR_COLORS["ncar_blue"]
    ext_color = style.NCAR_COLORS["red"]
    edge_color = style.NCAR_COLORS["gray"]

    mpas_top_km = z_mpas_edges[-1]

    panels = [
        ("Temperature [K]",                       T_mpas,     T_ext_anchored,  False),
        (r"Air number density [molec cm$^{-3}$]", nair_mpas,  nair_ext_e,      True),
        (r"O$_3$ number density [molec cm$^{-3}$]", no3_mpas, no3_ext_e,       True),
    ]

    for ax, (xlabel, mp_mid, ext_edges, logx) in zip(axes, panels):
        z_full, v_full = composite_column(z_mpas_mid, mp_mid, z_ext_edges, ext_edges)

        # MPAS segment (lower)
        m = z_full <= mpas_top_km
        ax.plot(v_full[m], z_full[m], "o-", color=mpas_color, ms=3, lw=1.5,
                label="MPAS")
        # Extension segment (upper)
        e = z_full >= mpas_top_km
        ax.plot(v_full[e], z_full[e], "s-", color=ext_color, ms=4, lw=1.5,
                label="Extension (CSV)")
        # Connector: last MPAS midpoint -> first extension midpoint.
        # Uses the edge-blending join at MPAS top so TUV-x's continuous
        # view is visible.
        last_mpas = np.where(m)[0][-1]
        first_ext = np.where(e)[0][0]
        ax.plot(v_full[[last_mpas, first_ext]], z_full[[last_mpas, first_ext]],
                "-", color=edge_color, lw=1.0, alpha=0.7)
        # Raw CSV edge markers for reference
        ax.plot(ext_edges, z_ext_edges, "x", color=ext_color, ms=7, alpha=0.4,
                label="CSV edges")

        ax.axhline(mpas_top_km, ls="--", color=edge_color, alpha=0.6)

        ax.set_xlabel(xlabel)
        if logx:
            ax.set_xscale("log")

    axes[0].set_ylabel("Altitude [km]")
    axes[0].legend(loc="best", fontsize=style.FONT_SIZES_DEFAULT.legend_small)
    axes[0].set_ylim(0, 105)
    axes[0].text(0.02, mpas_top_km / 105 + 0.01, "MPAS top",
                 transform=axes[0].transAxes, color=edge_color,
                 fontsize=style.FONT_SIZES_DEFAULT.annotation_small)

    fig.suptitle(
        style.format_title(
            f"Composite column seen by TUV-x  |  cell {cell}, time idx {tidx}"
        )
    )
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    pdf = Path(args.output).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {args.output} and {pdf}")
    print(f"cell = {cell}, time idx = {tidx}")
    print(f"MPAS top edge = {mpas_top_km:.2f} km;  extension spans "
          f"{z_ext_edges[0]:.1f}-{z_ext_edges[-1]:.1f} km "
          f"({len(z_ext_edges) - 1} layers)")


if __name__ == "__main__":
    main()
