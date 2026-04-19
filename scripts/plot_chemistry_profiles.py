#!/usr/bin/env python3
"""Horizontal-mean vertical profiles of chemistry output.

Plots one panel per variable (species tracer or photolysis rate) with
altitude on the y-axis and the horizontal-mean value on the x-axis.
Each output time gets its own coloured curve so we can see evolution.
A shaded band shows the min/max envelope across cells at each level.

Defaults: all common Chapman / Chapman+NOx species and photolysis rates
found in the file. Missing variables are silently skipped.

Usage:
    plot_chemistry_profiles.py -i ~/Data/CheMPAS/supercell/output.nc \
        -o chemistry_profiles.png
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# Default set of species / rates to look for. Ordered for a tidy panel layout.
DEFAULT_SPECIES = ["qO3", "qO", "qO1D", "qNO", "qNO2", "qO2"]
DEFAULT_JVARS   = ["j_jO2", "j_jO3_O", "j_jO3_O1D", "j_jNO2"]


def horizontal_stats(arr, axis_cells):
    """Return (mean, min, max) over the cells axis at every level."""
    return (arr.mean(axis=axis_cells),
            arr.min(axis=axis_cells),
            arr.max(axis=axis_cells))


def read_var(ds, name):
    if name not in ds.variables:
        return None
    return ds.variables[name][:]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default=str(Path.home() / "Data/CheMPAS/supercell/output.nc"))
    ap.add_argument("-o", "--output", default=None,
                    help="output PNG path; default <input_dir>/plots/chemistry_profiles.png")
    ap.add_argument("--logx", choices=["auto", "always", "never"], default="auto",
                    help="x-axis scale: 'auto' = log when dynamic range > 2 orders")
    args = ap.parse_args()

    if args.output is None:
        input_dir = Path(args.input).parent
        plots_dir = input_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        args.output = str(plots_dir / "chemistry_profiles.png")

    ds = Dataset(args.input, "r")
    zgrid = ds.variables["zgrid"][0, :]  # first-cell edges (flat supercell)
    z_mid_km = 0.5 * (zgrid[:-1] + zgrid[1:]) / 1000.0
    nTimes = ds.dimensions["Time"].size

    # Decide which variables are actually in the file.
    species = [v for v in DEFAULT_SPECIES if v in ds.variables]
    jvars   = [v for v in DEFAULT_JVARS   if v in ds.variables]
    variables = species + jvars
    if not variables:
        raise SystemExit("No known species or j_* diagnostic variables found.")

    ncols = 4
    nrows = int(np.ceil(len(variables) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows),
                              sharey=True)
    axes = np.atleast_2d(axes)

    times = np.arange(nTimes)
    cmap = plt.cm.viridis
    colors = [cmap(i / max(nTimes - 1, 1)) for i in range(nTimes)]

    for idx, name in enumerate(variables):
        ax = axes[idx // ncols, idx % ncols]
        arr = read_var(ds, name)  # (Time, nCells, nVertLevels)

        for t in range(nTimes):
            slab = arr[t, :, :]        # (nCells, nVertLevels)
            mean = slab.mean(axis=0)
            lo = slab.min(axis=0)
            hi = slab.max(axis=0)
            label = f"t={t}" if idx == 0 else None
            ax.plot(mean, z_mid_km, color=colors[t], lw=1.5, label=label)
            ax.fill_betweenx(z_mid_km, lo, hi, color=colors[t], alpha=0.15,
                              linewidth=0)

        ax.set_xlabel(name)
        if idx % ncols == 0:
            ax.set_ylabel("Altitude [km]")
        ax.grid(True, alpha=0.3)

        use_log = False
        if args.logx == "always":
            use_log = (arr[-1] > 0).any()
        elif args.logx == "auto":
            # log when the non-zero dynamic range spans >= 2 decades
            positive = arr[-1][arr[-1] > 0]
            if positive.size > 10:
                lo_val, hi_val = positive.min(), positive.max()
                if hi_val / max(lo_val, 1e-300) > 100.0:
                    use_log = True
        if use_log:
            ax.set_xscale("log")
        ax.set_title(name)

    # Hide any unused panels
    for j in range(len(variables), nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    axes[0, 0].legend(loc="best", fontsize=8, title="time idx")
    fig.suptitle(f"Horizontal-mean profiles from {Path(args.input).name} (shaded = min/max envelope)",
                 fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    pdf = Path(args.output).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {args.output} and {pdf}")
    print(f"variables plotted: {', '.join(variables)}")
    print(f"nTimes: {nTimes}, nVertLevels: {len(z_mid_km)}, "
          f"z range: {z_mid_km[0]:.2f} - {z_mid_km[-1]:.2f} km")

    ds.close()


if __name__ == "__main__":
    main()
