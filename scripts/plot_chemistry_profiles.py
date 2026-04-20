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
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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
SPECIES_PANELS = ["qO3", "qO_total", "qNO2", "qNO"]
# Photolysis: combine jO3 channels since O(1D) is rapidly quenched to
# O(3P) by M, so both produce atomic oxygen on the same timescale.
JVAR_PANELS    = ["j_jO2", "j_jO3_total", "j_jNO2"]

AXIS_SPECS = {
    # Species [ppb]
    "qO3":      {"range": (10.0, 1.0e5),    "log": True},   # AFGL MLS profile: ~27 to ~2.4e4 ppb
    "qO_total": {"range": (1.0e-4, 1.0e3),  "log": True},   # daytime Chapman-SS peak ~500 ppb; lower bound at numerical-floor threshold
    "qNO":      {"range": "auto",           "log": False},
    "qNO2":     {"range": "auto",           "log": False},
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


def xtime_to_datetime(xt_frame):
    """Decode an xtime char array to a Python datetime; None if year<1."""
    s = "".join([c.decode() if isinstance(c, bytes) else c for c in xt_frame]).strip()
    try:
        return datetime.strptime(s, "%Y-%m-%d_%H:%M:%S")
    except ValueError:
        return None


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
    if name.startswith("j_"):
        return "s$^{-1}$"
    return "ppb"


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
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    default=str(Path.home() / "Data/CheMPAS/supercell/output.nc"))
    ap.add_argument("-o", "--output", default=None,
                    help="output PNG path; default <input_dir>/plots/chemistry_profiles.png")
    ap.add_argument("--endpoints", action="store_true",
                    help="profile plot shows only the first and last frames "
                         "(skipping frame 0 if its j-rates are all zero, which "
                         "happens because MPAS writes the initial snapshot before "
                         "the first chemistry call)")
    ap.add_argument("--frames", default=None,
                    help="comma-separated frame indices for the profile plot "
                         "(supports negatives), e.g. --frames 1,20,40,-1")
    ap.add_argument("--time-series", default=None,
                    help="comma-separated altitudes [km] for a companion time-series "
                         "plot, e.g. --time-series 1,5,10,20,30,40")
    ap.add_argument("--title", default=None,
                    help="custom figure suptitle; overrides the default.")
    ap.add_argument("--tz-offset-hours", type=float, default=None,
                    help="hours offset from UTC. When set, the time-series "
                         "x-axis shows HH:MM local time (e.g. -6 for MDT).")
    ap.add_argument("--tz-label", default="local",
                    help="timezone label in the time-series x-axis (default: local)")
    ap.add_argument("--subtitle-extra", default=None,
                    help="extra text appended to the '{nCells}-hex domain mean' "
                         "subtitle line, separated by ', '")
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
    nCells = len(ds.dimensions["nCells"])
    domain_subtitle = f"{nCells}-hex domain mean"
    if args.subtitle_extra:
        domain_subtitle = f"{domain_subtitle}  {args.subtitle_extra}"

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

    # Frame selection: --frames > --endpoints > all frames.
    if args.frames:
        raw = [int(x) for x in args.frames.split(",") if x.strip()]
        frame_idx = [i if i >= 0 else nTimes + i for i in raw]
    elif args.endpoints and nTimes >= 2:
        # Skip frame 0 if its j-rates are all zero (pre-first-chem-call snapshot).
        first = 0
        jvars_present = [v for v in ("j_jO2", "j_jO3_O", "j_jO3_O1D", "j_jNO2")
                         if v in ds.variables]
        if jvars_present and all(
            float(np.asarray(ds.variables[v][0]).max()) == 0.0
            for v in jvars_present
        ):
            first = 1
        frame_idx = [first, nTimes - 1]
    else:
        frame_idx = list(range(nTimes))
    n_frames = len(frame_idx)

    # Day/night-aware palette: warm (YlOrRd) for frames with nonzero
    # photolysis, cool (Blues) for frames at night. Falls back to the NCAR
    # palette if jO2 is not in the output.
    if "j_jO2" in ds.variables:
        is_day = [float(np.asarray(ds.variables["j_jO2"][t]).max()) > 0.0
                  for t in frame_idx]
    else:
        is_day = [True] * n_frames
    day_count = sum(is_day)
    night_count = n_frames - day_count
    day_pos = (np.linspace(0.35, 0.90, day_count)
               if day_count > 1 else np.array([0.75]))
    night_pos = (np.linspace(0.95, 0.55, night_count)
                 if night_count > 1 else np.array([0.80]))
    time_colors = []
    d_i = n_i = 0
    for day in is_day:
        if day:
            time_colors.append(plt.cm.YlOrRd(day_pos[d_i])); d_i += 1
        else:
            time_colors.append(plt.cm.Blues(night_pos[n_i])); n_i += 1

    have_ss_inputs = all(v in ds.variables for v in
                         ("rho", "theta", "pressure", "qO2", "qO3",
                          "j_jO2", "j_jO3_O", "j_jO3_O1D"))

    t0_utc = xtime_to_datetime(ds.variables["xtime"][0])
    use_local_labels = args.tz_offset_hours is not None and t0_utc is not None
    t0_local = (t0_utc + timedelta(hours=args.tz_offset_hours)
                if use_local_labels else None)

    for idx, name in enumerate(variables):
        ax = axes[idx // ncols, idx % ncols]
        arr_plot = load_panel_data(ds, name)

        for k, t in enumerate(frame_idx):
            slab = arr_plot[t, :, :]
            mean = slab.mean(axis=0)
            lo = slab.min(axis=0)
            hi = slab.max(axis=0)
            if idx == 0:
                if use_local_labels:
                    label = (t0_local + timedelta(
                        seconds=float(time_s[t]))).strftime("%H:%M")
                else:
                    label = _format_seconds(time_s[t])
            else:
                label = None
            ax.plot(mean, z_mid_km, color=time_colors[k], lw=1.8, label=label)
            ax.fill_betweenx(z_mid_km, lo, hi, color=time_colors[k], alpha=0.12,
                              linewidth=0)

        # Chapman analytic steady-state [O] reference for the O panel
        if name == "qO_total" and have_ss_inputs:
            ppb_ss = chapman_steady_state_O_ppb(ds, nTimes - 1)
            ax.plot(ppb_ss, z_mid_km, color=style.NCAR_COLORS["gray"],
                    lw=1.5, ls="--", label="Steady State")

        ax.set_xlabel(panel_xlabel(name))
        ax.set_title(panel_title(name))
        if idx % ncols == 0:
            ax.set_ylabel("Altitude [km]")

        spec = AXIS_SPECS.get(name)
        if spec is not None:
            rng = spec["range"]
            if rng == "auto":
                data = np.asarray(arr_plot[frame_idx, :, :])
                lo_x = float(np.floor(float(np.min(data))))
                hi_x = float(np.ceil(float(np.max(data))))
                if hi_x <= lo_x:
                    hi_x = lo_x + 1.0
            else:
                lo_x, hi_x = rng
            ax.set_xlim(lo_x, hi_x)
            if spec["log"]:
                ax.set_xscale("log")

        apply_tick_style(ax, name, spec)

    for j in range(len(variables), nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    legend_title = f"Time ({args.tz_label})" if use_local_labels else "sim time"
    handles, labels = axes[0, 0].get_legend_handles_labels()

    empty_idx = len(variables)
    if empty_idx < nrows * ncols:
        legend_ax = axes[empty_idx // ncols, empty_idx % ncols]
        leg = legend_ax.legend(handles, labels, loc="center", title=legend_title,
                               fontsize=16, title_fontsize=18,
                               ncol=2, columnspacing=1.5,
                               handlelength=2.5, labelspacing=0.7,
                               borderpad=0.8, frameon=False)
        leg._legend_box.align = "center"
    else:
        axes[0, 0].legend(handles, labels, loc="best", title=legend_title,
                          fontsize=style.FONT_SIZES_DEFAULT.legend_small)

    # Steady-state reference gets its own small legend on the qO_total panel.
    if "qO_total" in variables:
        o_idx = variables.index("qO_total")
        ax_o = axes[o_idx // ncols, o_idx % ncols]
        for h, l in zip(*ax_o.get_legend_handles_labels()):
            if l == "Steady State":
                ax_o.legend([h], [l], loc="best",
                            fontsize=style.FONT_SIZES_DEFAULT.legend_small)
                break

    if args.title:
        fig.suptitle(f"{args.title}\n{domain_subtitle}")
    else:
        subtitle = "shaded = min/max envelope across cells"
        if args.endpoints and nTimes > 2:
            subtitle = f"endpoints only (t=0, t={int(round(time_s[-1]))} s); " + subtitle
        fig.suptitle(
            f"Horizontal-mean profiles from {Path(args.input).name}  |  {subtitle}"
            f"\n{domain_subtitle}"
        )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    pdf = Path(args.output).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {args.output} and {pdf}")
    print(f"panels: {', '.join(variables)}")
    print(f"times [s]: {time_s.tolist()}")

    if args.time_series:
        ts_alts = [float(x) for x in args.time_series.split(",") if x.strip()]
        ts_path = Path(args.output).with_name(
            Path(args.output).stem + "_timeseries.png"
        )
        _plot_time_series_at_levels(
            ds, variables, time_s, z_mid_km, ts_alts, str(ts_path),
            title=args.title, t0_utc=t0_utc,
            tz_offset_hours=args.tz_offset_hours, tz_label=args.tz_label,
            subtitle=domain_subtitle,
        )

    ds.close()


def _plot_time_series_at_levels(ds, variables, time_s, z_mid_km,
                                 target_alts_km, output_path,
                                 title=None, t0_utc=None,
                                 tz_offset_hours=None, tz_label="local",
                                 subtitle=None):
    """Companion plot: domain-mean time series at a set of altitudes."""
    level_idx = [int(np.argmin(np.abs(z_mid_km - a))) for a in target_alts_km]
    actual_alts = [z_mid_km[i] for i in level_idx]

    ncols = 4
    nrows = int(np.ceil(len(variables) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.2 * nrows),
                              sharex=True)
    axes = np.atleast_2d(axes)

    level_colors = style.get_palette(len(level_idx))

    use_local = tz_offset_hours is not None and t0_utc is not None
    if use_local:
        t0_local = t0_utc + timedelta(hours=tz_offset_hours)
        t_plot = [t0_local + timedelta(seconds=float(s)) for s in time_s]
        xlabel = f"Time ({tz_label})"
    else:
        t_plot = time_s / 60.0
        xlabel = "Time [min]"

    for idx, name in enumerate(variables):
        ax = axes[idx // ncols, idx % ncols]
        arr_plot = load_panel_data(ds, name)
        for k, (li, za) in enumerate(zip(level_idx, actual_alts)):
            series = arr_plot[:, :, li].mean(axis=1)
            label = f"{za:.1f} km" if idx == 0 else None
            ax.plot(t_plot, series, color=level_colors[k], lw=1.8, label=label)

        ax.set_title(panel_title(name))
        ax.set_xlabel(xlabel)
        # Each panel carries its own units (ppb vs s^-1), so label every one.
        ax.set_ylabel(panel_xlabel(name))
        spec = AXIS_SPECS.get(name)
        if spec is not None and spec["log"]:
            # Log scale is useful for wide column-range plots (profiles), but in
            # a time-series at a single level the range is bounded. Matplotlib's
            # log autoscaler latches onto tiny nighttime zeros and stretches the
            # axis over ~120 decades, so clamp the lower bound to the smallest
            # physically meaningful value across the plotted panels.
            positives = arr_plot[:, :, level_idx].ravel()
            positives = positives[positives > 0]
            if positives.size > 0:
                ymin = max(positives.min(), 10 ** np.floor(
                    np.log10(positives.max()) - 4))
                ax.set_ylim(bottom=ymin)
            ax.set_yscale("log")
        ax.grid(True, alpha=0.3)

    for j in range(len(variables), nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles, labels = list(reversed(handles)), list(reversed(labels))
    empty_idx = len(variables)
    if empty_idx < nrows * ncols:
        legend_ax = axes[empty_idx // ncols, empty_idx % ncols]
        leg = legend_ax.legend(handles, labels, loc="center", title="Altitude",
                               fontsize=16, title_fontsize=18,
                               handlelength=2.5, labelspacing=0.7,
                               borderpad=0.8, frameon=False)
        leg._legend_box.align = "center"
    else:
        axes[0, 0].legend(handles, labels, loc="best", title="Altitude",
                          fontsize=style.FONT_SIZES_DEFAULT.legend_small)

    for ax in axes.flat:
        ax.set_xlim(t_plot[0], t_plot[-1])
    if use_local:
        for ax in axes.flat:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        fig.autofmt_xdate(rotation=30, ha="right")

    primary = title or f"Domain-mean time series from {Path(ds.filepath()).name}"
    fig.suptitle(f"{primary}\n{subtitle}" if subtitle else primary)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    pdf = Path(output_path).with_suffix(".pdf")
    plt.savefig(pdf, bbox_inches="tight")
    print(f"saved {output_path} and {pdf}")
    print(f"levels picked [km]: {[round(a, 2) for a in actual_alts]}")


if __name__ == "__main__":
    main()
