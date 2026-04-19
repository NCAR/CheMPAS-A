#!/usr/bin/env python3
"""Horizontal-mean vertical profiles of chemistry output.

Plots one panel per variable (species tracer or photolysis rate) with
altitude on the y-axis and the horizontal-mean value on the x-axis.
Each output time gets its own coloured curve so we can see evolution.
A shaded band shows the min/max envelope across cells at each level.

Conventions (scripts/style.py):
- NCAR brand palette via style.setup()
- Species-specific colors via style.species_color()
- LaTeX species labels via style.species_label() / style.format_title()
- Chemistry tracers are shown in ppb (converted from MPAS kg/kg using
  per-species molar masses); photolysis rates in s^{-1}
- Units are always in axis labels.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402


# Molar masses [kg/mol] for mass->ppbv conversion
M_AIR = 0.02897
MOLAR_MASS = {
    "qO":   0.016,
    "qO1D": 0.016,
    "qO2":  0.032,
    "qO3":  0.048,
    "qNO":  0.030,
    "qNO2": 0.046,
    "qCO":  0.028,
    "qCH4": 0.016,
    "qOH":  0.017,
    "qHO2": 0.033,
}

DEFAULT_SPECIES = ["qO3", "qO", "qO1D", "qNO", "qNO2", "qO2"]
DEFAULT_JVARS   = ["j_jO2", "j_jO3_O", "j_jO3_O1D", "j_jNO2"]


# Fixed x-axis ranges sized to the physically expected values (not
# idealized spin-up transients or floating-point denormals). Species in
# ppb, photolysis rates in s^-1.
AXIS_SPECS = {
    # Species [ppb]
    "qO3":   {"range": (1.0, 1.0e5),     "log": True},   # 50 ppb tropo → 10 ppm stratospheric peak
    "qO":    {"range": (0.0, 5.0e7),     "log": False},  # idealized spin-up transient; steady-state tiny
    "qO1D":  {"range": (0.0, 1.0),       "log": False},  # idealized spin-up transient
    "qNO":   {"range": (0.0, 5.0),       "log": False},  # stratospheric NO peak ~3 ppb
    "qNO2":  {"range": (0.0, 10.0),      "log": False},  # stratospheric NO2 peak ~8 ppb
    "qO2":   {"range": (1.8e8, 2.3e8),   "log": False},  # ~0.2095 VMR -> ~2.09e8 ppb
    # Photolysis rates [s^-1]
    "j_jO2":     {"range": (0.0, 5.0e-10), "log": False},  # Schumann-Runge weak
    "j_jO3_O":   {"range": (0.0, 1.0e-3),  "log": False},  # Hartley band, O(3P) channel
    "j_jO3_O1D": {"range": (0.0, 1.0e-3),  "log": False},  # Hartley band, O(1D) channel
    "j_jNO2":    {"range": (0.0, 1.5e-2),  "log": False},  # midday NO2 photolysis
}


def to_ppbv(q_kgkg, M_species):
    """Convert mass mixing ratio (kg/kg) to ppbv."""
    return q_kgkg * (M_AIR / M_species) * 1.0e9


def rate_label(name):
    """Return a LaTeX display label for a j_<rxn> diagnostic variable."""
    if not name.startswith("j_"):
        return name
    rxn = name[2:]
    # Known rates: jNO2, jO2, jO3_O, jO3_O1D
    pretty = {
        "jNO2":    r"$j_{NO_2}$",
        "jO2":     r"$j_{O_2}$",
        "jO3_O":   r"$j_{O_3 \to O(^3P)}$",
        "jO3_O1D": r"$j_{O_3 \to O(^1D)}$",
    }
    return pretty.get(rxn, name)


def species_panel_label(name):
    return f"{style.species_label(name)} [ppb]"


def rate_panel_label(name):
    return f"{rate_label(name)} [s$^{{-1}}$]"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    default=str(Path.home() / "Data/CheMPAS/supercell/output.nc"))
    ap.add_argument("-o", "--output", default=None,
                    help="output PNG path; default <input_dir>/plots/chemistry_profiles.png")
    args = ap.parse_args()

    if args.output is None:
        plots_dir = Path(args.input).parent / "plots"
        plots_dir.mkdir(exist_ok=True)
        args.output = str(plots_dir / "chemistry_profiles.png")

    style.setup()

    ds = Dataset(args.input, "r")
    zgrid = ds.variables["zgrid"][0, :]
    z_mid_km = 0.5 * (zgrid[:-1] + zgrid[1:]) / 1000.0
    nTimes = ds.dimensions["Time"].size

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

    # Time-series color ramp: viridis, one shade per output frame
    cmap = plt.cm.viridis
    time_colors = [cmap(i / max(nTimes - 1, 1)) for i in range(nTimes)]

    for idx, name in enumerate(variables):
        ax = axes[idx // ncols, idx % ncols]
        arr = ds.variables[name][:]  # (Time, nCells, nVertLevels)

        is_species = name in species
        if is_species:
            Mw = MOLAR_MASS.get(name)
            if Mw is None:
                # Unknown species — plot raw kg/kg
                arr_plot = arr
                xlabel = f"{style.species_label(name)} [kg kg$^{{-1}}$]"
            else:
                arr_plot = to_ppbv(arr, Mw)
                xlabel = species_panel_label(name)
            title = style.species_label(name)
        else:
            arr_plot = arr
            xlabel = rate_panel_label(name)
            title = rate_label(name)

        for t in range(nTimes):
            slab = arr_plot[t, :, :]
            mean = slab.mean(axis=0)
            lo = slab.min(axis=0)
            hi = slab.max(axis=0)
            label = f"t={t}" if idx == 0 else None
            ax.plot(mean, z_mid_km, color=time_colors[t], lw=1.5, label=label)
            ax.fill_betweenx(z_mid_km, lo, hi, color=time_colors[t], alpha=0.15,
                              linewidth=0)

        ax.set_xlabel(xlabel)
        ax.set_title(title)
        if idx % ncols == 0:
            ax.set_ylabel("Altitude [km]")

        # Fixed x-axis range per the physically expected values in AXIS_SPECS.
        # Data outside the range is clipped visually, which is the intent —
        # numerical transients and denormals are not interesting here.
        spec = AXIS_SPECS.get(name)
        if spec is not None:
            lo, hi = spec["range"]
            ax.set_xlim(lo, hi)
            if spec["log"]:
                ax.set_xscale("log")

    # Hide unused panels
    for j in range(len(variables), nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    axes[0, 0].legend(loc="best",
                       fontsize=style.FONT_SIZES_DEFAULT.legend_small,
                       title="time idx")

    fig.suptitle(
        f"Horizontal-mean profiles from {Path(args.input).name}  |  "
        "shaded = min/max envelope across cells"
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
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
