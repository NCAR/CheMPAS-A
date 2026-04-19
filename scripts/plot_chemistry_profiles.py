#!/usr/bin/env python3
"""Horizontal-mean vertical profiles of chemistry output.

One panel per variable (species or photolysis rate), altitude on y,
horizontal-mean on x. Each output frame gets its own curve labelled
with the simulation time in seconds; a shaded band is the min/max
envelope across cells at each level.

Conventions (match scripts/style.py and scripts/plot_lnox_o3.py):
- style.setup() for NCAR brand palette / fonts
- Simulation time in seconds parsed from xtime (reference = t=0 of
  the init file)
- Tracer mass mixing ratios shown in ppb via per-species molar mass
- NCAR palette colours for the time series
- Photolysis rate tick labels trimmed to avoid overlap

Panels:
- O3, NO, NO2       (main chemistry tracers in ppb)
- O = q_O + q_O1D   (atomic oxygen summed; O1D is rapidly quenched
                      to O(3P) by M so the sum is the useful view)
- O2 is not plotted (constant-ish at 0.2095 VMR by construction)
- jO2, jO3 -> O(3P), jO3 -> O(1D), jNO2  (photolysis rates in s^-1)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402


# Molar masses [kg/mol] for mass->ppbv conversion
M_AIR = 0.02897
M_O2  = 0.032
M_O3  = 0.048
MOLAR_MASS = {
    "qO":   0.016,
    "qO1D": 0.016,
    "qO2":  0.032,
    "qO3":  0.048,
    "qNO":  0.030,
    "qNO2": 0.046,
}

NA = 6.02214076e23  # Avogadro

# Species shown (O = qO + qO1D summed; O2 omitted because near-constant)
SPECIES_PANELS = ["qO3", "qO_total", "qNO", "qNO2"]
# Photolysis: combine jO3 channels since O(1D) is rapidly quenched to
# O(3P) by M, so both produce atomic oxygen on the same timescale.
JVAR_PANELS    = ["j_jO2", "j_jO3_total", "j_jNO2"]

AXIS_SPECS = {
    # Species [ppb]
    "qO3":      {"range": (1.0, 1.0e5),     "log": True},   # AFGL MLS profile
    "qO_total": {"range": (1.0e-4, 1.0e8),  "log": True},   # wide: spin-up (~1e7) + Chapman SS (~1e2)
    "qNO":      {"range": (0.0, 5.0),       "log": False},
    "qNO2":     {"range": (0.0, 10.0),      "log": False},
    # Photolysis rates [s^-1]
    "j_jO2":       {"range": (0.0, 5.0e-10), "log": False},
    "j_jO3_total": {"range": (0.0, 2.0e-3),  "log": False},  # jO3_O + jO3_O1D combined
    "j_jNO2":      {"range": (0.0, 1.5e-2),  "log": False},
}


def to_ppbv(q_kgkg, M_species):
    """Convert mass mixing ratio (kg/kg) to ppbv."""
    return q_kgkg * (M_AIR / M_species) * 1.0e9


def xtime_to_seconds(xt_frame):
    """Decode an xtime char array to absolute seconds since 0000-01-01_00:00:00."""
    s = "".join([c.decode() if isinstance(c, bytes) else c for c in xt_frame]).strip()
    parts = s.split("_")
    dparts = parts[0].split("-")
    hms = parts[-1].split(":")
    days = int(dparts[0]) * 365 + (int(dparts[1]) - 1) * 30 + int(dparts[2]) - 1
    return days * 86400 + int(hms[0]) * 3600 + int(hms[1]) * 60 + int(hms[2])


def read_times_seconds(filename):
    """Return simulation time in seconds since the init-file t=0."""
    ds = Dataset(filename, "r")
    xt = ds.variables["xtime"][:]
    t0 = xtime_to_seconds(xt[0])
    init_file = Path(filename).parent / "supercell_init.nc"
    if init_file.exists():
        with Dataset(init_file, "r") as dsi:
            t0 = xtime_to_seconds(dsi.variables["xtime"][0])
    times = np.array([xtime_to_seconds(xt[t]) - t0 for t in range(len(xt))])
    ds.close()
    return times


def rate_label(name):
    rxn = name[2:] if name.startswith("j_") else name
    pretty = {
        "jNO2":       r"$j_{NO_2}$",
        "jO2":        r"$j_{O_2}$",
        "jO3_O":      r"$j_{O_3 \to O(^3P)}$",
        "jO3_O1D":    r"$j_{O_3 \to O(^1D)}$",
        "jO3_total":  r"$j_{O_3 \to O}$",   # jO3_O + jO3_O1D combined
    }
    return pretty.get(rxn, name)


def panel_title(name):
    if name == "qO_total":
        return "O"   # sum of qO + qO1D, labelled simply as O
    if name.startswith("j_"):
        return rate_label(name)
    return style.species_label(name)


def panel_xlabel(name):
    if name == "qO_total":
        return "O [ppb]"
    if name.startswith("j_"):
        return f"{rate_label(name)} [s$^{{-1}}$]"
    return f"{style.species_label(name)} [ppb]"


def chapman_steady_state_O_ppb(ds, time_idx):
    """Return horizontal-mean Chapman steady-state [O] in ppb on the MPAS grid.

    Solves P_O = L_O at each (cell, level):
      P_O = 2 * jO2 * [O2] + (jO3_O + jO3_O1D) * [O3]
      L_O = k(O+O2+M) * [M] * [O2] * [O] + k(O+O3) * [O3] * [O]
    with JPL-style rate constants:
      k(O+O2+M) = 6.1e-34 * (T/300)^-2.4   cm^6 molec^-2 s^-1
      k(O+O3)   = 8.0e-12 * exp(-2060/T)    cm^3 molec^-1 s^-1

    The product [M]*[O2]*[O] from the 3-body reaction dominates below
    ~60 km. The horizontal mean is taken after the per-cell solve
    because the formula is non-linear in T.
    """
    rho    = ds.variables["rho"][time_idx, :, :]
    theta  = ds.variables["theta"][time_idx, :, :]
    press  = ds.variables["pressure"][time_idx, :, :]
    qO2    = ds.variables["qO2"][time_idx, :, :]
    qO3    = ds.variables["qO3"][time_idx, :, :]
    jO2_a  = ds.variables["j_jO2"][time_idx, :, :]
    jO3a_a = ds.variables["j_jO3_O"][time_idx, :, :]
    jO3b_a = ds.variables["j_jO3_O1D"][time_idx, :, :]
    jO3    = jO3a_a + jO3b_a

    T = theta * (press / 1.0e5) ** (287.0 / 1004.0)

    # Number densities [molec / cm^3]
    n_M  = rho / M_AIR * NA * 1.0e-6
    n_O2 = qO2 * rho / M_O2 * NA * 1.0e-6
    n_O3 = qO3 * rho / M_O3 * NA * 1.0e-6

    k_O_O2_M = 6.1e-34 * (T / 300.0) ** (-2.4)           # cm^6/molec^2/s
    k_O_O3   = 8.0e-12 * np.exp(-2060.0 / T)             # cm^3/molec/s

    prod = 2.0 * jO2_a * n_O2 + jO3 * n_O3
    loss_coeff = k_O_O2_M * n_M * n_O2 + k_O_O3 * n_O3
    with np.errstate(divide="ignore", invalid="ignore"):
        n_O_ss = np.where(loss_coeff > 0, prod / loss_coeff, np.nan)
        ppb_ss = n_O_ss / n_M * 1.0e9

    return np.nanmean(ppb_ss, axis=0)


def load_panel_data(ds, name):
    """Return an (nTimes, nCells, nVertLevels) array in the panel's display units."""
    if name == "qO_total":
        # Sum atomic oxygen channels; both have M = 0.016 so convert together.
        qO  = ds.variables["qO"][:] if "qO" in ds.variables else 0.0
        qO1D = ds.variables["qO1D"][:] if "qO1D" in ds.variables else 0.0
        return to_ppbv(qO + qO1D, MOLAR_MASS["qO"])
    if name == "j_jO3_total":
        # O(1D) is rapidly quenched to O(3P), so both channels produce O.
        jO  = ds.variables["j_jO3_O"][:] if "j_jO3_O" in ds.variables else 0.0
        j1D = ds.variables["j_jO3_O1D"][:] if "j_jO3_O1D" in ds.variables else 0.0
        return jO + j1D
    if name.startswith("j_"):
        return ds.variables[name][:]
    Mw = MOLAR_MASS.get(name)
    arr = ds.variables[name][:]
    return to_ppbv(arr, Mw) if Mw else arr


def _format_seconds(s):
    return f"t = {int(round(float(s)))} s"


def _photolysis_tick_formatter():
    """Compact scientific-notation formatter for photolysis rate ticks."""
    def fmt(x, _pos):
        if x == 0:
            return "0"
        return f"{x:.1e}"
    return mticker.FuncFormatter(fmt)


def apply_tick_style(ax, name, spec):
    """Limit tick density and formatting per panel type."""
    if spec is not None and spec["log"]:
        # Log axes: one tick per decade keeps labels from overlapping.
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=6))
        ax.xaxis.set_minor_locator(mticker.LogLocator(
            base=10.0, subs=np.arange(2, 10) * 0.1, numticks=6))
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        return
    if name.startswith("j_"):
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=3))
        ax.xaxis.set_major_formatter(_photolysis_tick_formatter())
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5))


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
    time_s = read_times_seconds(args.input)
    nTimes = len(time_s)

    # Drop panels whose inputs are missing
    species = [v for v in SPECIES_PANELS
               if v == "qO_total" or v in ds.variables]
    if "qO_total" in species and "qO" not in ds.variables and "qO1D" not in ds.variables:
        species.remove("qO_total")
    jvars = []
    for v in JVAR_PANELS:
        if v == "j_jO3_total":
            # Synthesized: needs both components
            if "j_jO3_O" in ds.variables and "j_jO3_O1D" in ds.variables:
                jvars.append(v)
        elif v in ds.variables:
            jvars.append(v)
    variables = species + jvars
    if not variables:
        raise SystemExit("No known species or j_* diagnostic variables found.")

    ncols = 4
    nrows = int(np.ceil(len(variables) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows),
                              sharey=True)
    axes = np.atleast_2d(axes)

    # Time-series colours from the NCAR palette so they match the rest of
    # the project. If there are more frames than palette colours we cycle.
    time_colors = style.get_palette(nTimes)

    have_ss_inputs = all(v in ds.variables for v in
                         ("rho", "theta", "pressure", "qO2", "qO3",
                          "j_jO2", "j_jO3_O", "j_jO3_O1D"))

    for idx, name in enumerate(variables):
        ax = axes[idx // ncols, idx % ncols]
        arr_plot = load_panel_data(ds, name)

        for t in range(nTimes):
            slab = arr_plot[t, :, :]
            mean = slab.mean(axis=0)
            lo = slab.min(axis=0)
            hi = slab.max(axis=0)
            label = _format_seconds(time_s[t]) if idx == 0 else None
            ax.plot(mean, z_mid_km, color=time_colors[t], lw=1.8, label=label)
            ax.fill_betweenx(z_mid_km, lo, hi, color=time_colors[t], alpha=0.12,
                              linewidth=0)

        # Chapman analytic steady-state [O] reference for the O panel
        if name == "qO_total" and have_ss_inputs:
            ppb_ss = chapman_steady_state_O_ppb(ds, nTimes - 1)
            ax.plot(ppb_ss, z_mid_km, color=style.NCAR_COLORS["gray"],
                    lw=1.5, ls="--", label="Chapman SS")

        ax.set_xlabel(panel_xlabel(name))
        ax.set_title(panel_title(name))
        if idx % ncols == 0:
            ax.set_ylabel("Altitude [km]")

        spec = AXIS_SPECS.get(name)
        if spec is not None:
            lo, hi = spec["range"]
            ax.set_xlim(lo, hi)
            if spec["log"]:
                ax.set_xscale("log")

        apply_tick_style(ax, name, spec)

    for j in range(len(variables), nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    axes[0, 0].legend(loc="best",
                       fontsize=style.FONT_SIZES_DEFAULT.legend_small,
                       title="sim time")

    # O panel gets its own single-entry legend for the steady-state line
    if "qO_total" in variables:
        o_idx = variables.index("qO_total")
        ax_o = axes[o_idx // ncols, o_idx % ncols]
        handles, labels = ax_o.get_legend_handles_labels()
        ss = [(h, l) for h, l in zip(handles, labels) if l == "Chapman SS"]
        if ss:
            ax_o.legend([ss[0][0]], [ss[0][1]], loc="best",
                        fontsize=style.FONT_SIZES_DEFAULT.legend_small)

    fig.suptitle(
        f"Horizontal-mean profiles from {Path(args.input).name}  |  "
        "shaded = min/max envelope across cells"
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    pdf = Path(args.output).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {args.output} and {pdf}")
    print(f"panels: {', '.join(variables)}")
    print(f"times [s]: {time_s.tolist()}")

    ds.close()


if __name__ == "__main__":
    main()
