#!/usr/bin/env python3
"""
Fill an MPAS-A init file tracer (e.g., qAB or qX) with a horizontal sine wave
perturbation around a baseline value (default baseline = 1).

By default the script edits supercell_init.nc in place and overwrites the
requested tracer for all times and vertical levels.
"""

import argparse
import shutil
from pathlib import Path
from typing import Optional

import numpy as np
from netCDF4 import Dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overwrite a tracer in an MPAS init file with a sine wave pattern."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="supercell_init.nc",
        help="Path to the input init file (default: supercell_init.nc).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Optional output path. If omitted, the input file is edited in place.",
    )
    parser.add_argument(
        "-t",
        "--tracer",
        default="qAB",
        help="Tracer variable name to overwrite (default: qAB).",
    )
    parser.add_argument(
        "--amplitude",
        type=float,
        default=0.5,
        help="Amplitude of the sine perturbation about the baseline (default: 0.5).",
    )
    parser.add_argument(
        "--offset",
        type=float,
        default=1.0,
        help="Baseline (mean) value added to the sine pattern (default: 1).",
    )
    parser.add_argument(
        "--waves-x",
        type=float,
        default=1.0,
        help="Number of sine cycles across the domain in x (default: 1).",
    )
    parser.add_argument(
        "--waves-y",
        type=float,
        default=1.0,
        help="Number of sine cycles across the domain in y (default: 1).",
    )
    parser.add_argument(
        "--phase-x",
        type=float,
        default=0.0,
        help="Phase shift in radians applied to the x sine.",
    )
    parser.add_argument(
        "--phase-y",
        type=float,
        default=0.0,
        help="Phase shift in radians applied to the y sine.",
    )
    parser.add_argument(
        "--create",
        action="store_true",
        help="Create the tracer variable if it does not exist.",
    )
    parser.add_argument(
        "--units",
        default="kg kg^{-1}",
        help="Units to assign when creating a new tracer (default: kg kg^{-1}).",
    )
    parser.add_argument(
        "--long-name",
        dest="long_name",
        help="Long name to assign when creating a new tracer. Defaults to the tracer name.",
    )
    parser.add_argument(
        "--spherical",
        action="store_true",
        help="Use lonCell/latCell (radians) instead of xCell/yCell — needed for global meshes.",
    )
    return parser.parse_args()


def copy_if_needed(src: Path, dest: Path) -> None:
    if src.resolve() != dest.resolve():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)


def ensure_tracer_var(
    ds: Dataset, name: str, create: bool, units: str, long_name: Optional[str]
):
    if name in ds.variables:
        return ds.variables[name]
    if not create:
        raise SystemExit(f"Tracer '{name}' not found in {ds.filepath()} (use --create to add it).")
    dtype = ds.variables["qv"].dtype if "qv" in ds.variables else "f4"
    tracer = ds.createVariable(name, dtype, ("Time", "nCells", "nVertLevels"))
    tracer.units = units
    tracer.long_name = long_name or name
    return tracer


def build_sine_pattern(
    x: np.ndarray,
    y: np.ndarray,
    waves_x: float,
    waves_y: float,
    phase_x: float,
    phase_y: float,
    amplitude: float,
    offset: float,
    spherical: bool = False,
) -> np.ndarray:
    if spherical:
        # x=lon in [0, 2π), y=lat in [-π/2, π/2]. Wrap longitudinally with
        # waves_x full cycles around the globe; cos(lat)^waves_y tapers toward
        # the poles so the pattern is smooth there.
        phase_x_total = waves_x * x + phase_x
        horizontal = np.sin(phase_x_total) * np.cos(y) ** max(waves_y, 1.0)
        return offset + amplitude * horizontal

    x_span = float(np.ptp(x))
    y_span = float(np.ptp(y))
    if x_span == 0 or y_span == 0:
        raise SystemExit("Domain spans in x or y are zero; cannot build sine pattern.")

    phase_x_total = 2.0 * np.pi * waves_x * (x - x.min()) / x_span + phase_x
    phase_y_total = 2.0 * np.pi * waves_y * (y - y.min()) / y_span + phase_y

    horizontal = np.sin(phase_x_total)
    horizontal *= np.sin(phase_y_total)

    return offset + amplitude * horizontal


def main() -> None:
    args = parse_args()
    src = Path(args.input)
    dest = Path(args.output) if args.output else src

    if not src.exists():
        raise SystemExit(f"Input file does not exist: {src}")

    copy_if_needed(src, dest)

    with Dataset(dest, "r+") as ds:
        tracer = ensure_tracer_var(ds, args.tracer, args.create, args.units, args.long_name)

        if args.spherical:
            x = ds.variables["lonCell"][:]
            y = ds.variables["latCell"][:]
        else:
            x = ds.variables["xCell"][:]
            y = ds.variables["yCell"][:]
        n_cells = ds.dimensions["nCells"].size
        n_vert = ds.dimensions["nVertLevels"].size

        if tracer.shape[1] != n_cells or tracer.shape[2] != n_vert:
            raise SystemExit(
                f"Tracer '{args.tracer}' has unexpected shape {tracer.shape}; expected (Time, {n_cells}, {n_vert})."
            )

        horizontal_field = build_sine_pattern(
            x=x,
            y=y,
            waves_x=args.waves_x,
            waves_y=args.waves_y,
            phase_x=args.phase_x,
            phase_y=args.phase_y,
            amplitude=args.amplitude,
            offset=args.offset,
            spherical=args.spherical,
        )

        field_3d = np.broadcast_to(horizontal_field[:, None], (n_cells, n_vert))
        field_with_time = np.broadcast_to(field_3d, (tracer.shape[0], n_cells, n_vert))
        tracer[:] = field_with_time
        ds.sync()

    print(f"Updated tracer '{args.tracer}' in {dest}")


if __name__ == "__main__":
    main()
