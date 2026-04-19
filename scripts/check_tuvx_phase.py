#!/usr/bin/env python3
"""
Phase-gate checks for CheMPAS-A LNOx-O3 / TUV-x integration.

Imported from the ancestor MPAS-Model-ACOM-dev tooling and adapted for the
CheMPAS-A phase matrix. Each subcommand prints a compact report and exits
non-zero on failure.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


EPS = 1.0e-30


@dataclass
class CheckResult:
    passed: bool
    messages: list[str]


def parse_csv_list(raw: str) -> list[str]:
    return [v.strip() for v in raw.split(",") if v.strip()]


def open_dataset(path: str):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing NetCDF file: {p}")
    try:
        from netCDF4 import Dataset  # type: ignore
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            "Python module 'netCDF4' is required for check_tuvx_phase.py. "
            "Install with: pip install netCDF4"
        ) from exc
    return Dataset(str(p), "r")


def get_var(ds, name: str) -> np.ndarray:
    if name not in ds.variables:
        raise KeyError(f"Variable '{name}' not found in {ds.filepath()}")
    return np.asarray(ds.variables[name][:])


def to_tcl(arr: np.ndarray, n_cells: int, n_levels: int, name: str) -> np.ndarray:
    """
    Normalize 3D arrays to [time, cell, level].
    Supports [time, cell, level] and [time, level, cell].
    """
    if arr.ndim != 3:
        raise ValueError(f"{name}: expected 3D array, got shape {arr.shape}")
    if arr.shape[1] == n_cells and arr.shape[2] == n_levels:
        return arr
    if arr.shape[1] == n_levels and arr.shape[2] == n_cells:
        return np.transpose(arr, (0, 2, 1))
    raise ValueError(
        f"{name}: could not infer [time,cell,level] ordering from shape {arr.shape}, "
        f"nCells={n_cells}, nVertLevels={n_levels}"
    )


def to_tc(arr: np.ndarray, n_cells: int, name: str) -> np.ndarray:
    """
    Normalize 2D arrays to [time, cell].
    Supports [time, cell] and [cell, time].
    """
    if arr.ndim != 2:
        raise ValueError(f"{name}: expected 2D array, got shape {arr.shape}")
    if arr.shape[1] == n_cells:
        return arr
    if arr.shape[0] == n_cells:
        return np.transpose(arr, (1, 0))
    raise ValueError(
        f"{name}: could not infer [time,cell] ordering from shape {arr.shape}, nCells={n_cells}"
    )


def to_tcl_compare(arr: np.ndarray, n_cells: int, n_levels: int, name: str) -> np.ndarray:
    """
    Normalize fields for comparison to [time, cell, level_like].
    Supports 3D tracer fields and 2D diagnostic fields.
    """
    if arr.ndim == 3:
        return to_tcl(arr, n_cells, n_levels, name)
    if arr.ndim == 2:
        return to_tc(arr, n_cells, name)[:, :, None]
    raise ValueError(f"{name}: expected 2D/3D array, got shape {arr.shape}")


def zgrid_to_cl(arr: np.ndarray, n_cells: int, n_levels_p1: int, name: str) -> np.ndarray:
    """
    Normalize zgrid-like arrays to [cell, level+1].
    Supports [level+1, cell], [cell, level+1], and time-leading 3D variants.
    """
    if arr.ndim == 3:
        arr = arr[0, :, :] if arr.shape[0] != n_cells else arr[:, :, 0]
    if arr.ndim != 2:
        raise ValueError(f"{name}: expected 2D/3D array, got shape {arr.shape}")
    if arr.shape[0] == n_levels_p1 and arr.shape[1] == n_cells:
        return np.transpose(arr, (1, 0))
    if arr.shape[0] == n_cells and arr.shape[1] == n_levels_p1:
        return arr
    raise ValueError(
        f"{name}: could not infer [cell,level+1] ordering from shape {arr.shape}, "
        f"nCells={n_cells}, nVertLevelsP1={n_levels_p1}"
    )


def check_nonnegative(ds, vars_csv: str, tol: float) -> CheckResult:
    vars_ = parse_csv_list(vars_csv)
    messages: list[str] = []
    failed = False
    for name in vars_:
        arr = get_var(ds, name)
        min_val = float(np.nanmin(arr))
        max_val = float(np.nanmax(arr))
        ok = min_val >= -tol
        state = "PASS" if ok else "FAIL"
        messages.append(
            f"[{state}] {name}: min={min_val:.6e} max={max_val:.6e} tol={tol:.3e}"
        )
        failed = failed or (not ok)
    return CheckResult(passed=not failed, messages=messages)


def check_oxygen_budget(
    ds,
    qo: str,
    qo2: str,
    qo3: str,
    rho_var: str,
    zgrid_var: str,
    max_domain_drift: float,
    max_column_p95_drift: float,
) -> CheckResult:
    n_cells = len(ds.dimensions["nCells"])
    n_levels = len(ds.dimensions["nVertLevels"])
    n_levels_p1 = len(ds.dimensions["nVertLevelsP1"])

    q_o = to_tcl(get_var(ds, qo), n_cells, n_levels, qo)
    q_o2 = to_tcl(get_var(ds, qo2), n_cells, n_levels, qo2)
    q_o3 = to_tcl(get_var(ds, qo3), n_cells, n_levels, qo3)
    oxygen_mass = q_o + q_o2 + q_o3

    # Prefer mass-weighted column burden if rho and zgrid are present.
    have_rho = rho_var in ds.variables
    have_zgrid = zgrid_var in ds.variables
    if have_rho and have_zgrid:
        rho = to_tcl(get_var(ds, rho_var), n_cells, n_levels, rho_var)
        z_clp1 = zgrid_to_cl(get_var(ds, zgrid_var), n_cells, n_levels_p1, zgrid_var)
        dz_cl = np.maximum(z_clp1[:, 1:] - z_clp1[:, :-1], 0.0)
        weights = rho * dz_cl[None, :, :]
    else:
        weights = np.ones_like(oxygen_mass)

    weighted = oxygen_mass * weights
    domain_total = np.sum(weighted, axis=(1, 2))
    column_total = np.sum(weighted, axis=2)

    base_domain = float(domain_total[0])
    base_columns = column_total[0, :]
    domain_rel = np.abs(domain_total - base_domain) / max(abs(base_domain), EPS)

    # Avoid exploding drift for columns whose baseline is effectively zero.
    valid = np.abs(base_columns) > EPS
    if np.any(valid):
        col_rel = np.abs(column_total[:, valid] - base_columns[valid][None, :]) / np.abs(
            base_columns[valid][None, :]
        )
        per_col_max = np.max(col_rel, axis=0)
        col_p95 = float(np.percentile(per_col_max, 95.0))
    else:
        col_p95 = 0.0

    max_domain = float(np.max(domain_rel))
    pass_domain = max_domain <= max_domain_drift
    pass_column = col_p95 <= max_column_p95_drift

    mode = "mass-weighted (rho*dz)" if have_rho and have_zgrid else "unweighted"
    messages = [
        f"[INFO] oxygen budget mode: {mode}",
        (
            f"[{'PASS' if pass_domain else 'FAIL'}] domain drift: "
            f"max={max_domain:.6e} threshold={max_domain_drift:.6e}"
        ),
        (
            f"[{'PASS' if pass_column else 'FAIL'}] column drift p95: "
            f"{col_p95:.6e} threshold={max_column_p95_drift:.6e}"
        ),
    ]
    return CheckResult(passed=(pass_domain and pass_column), messages=messages)


def check_night_jzero(
    ds,
    coszr_var: str,
    j_vars_csv: str,
    night_threshold: float,
    abs_tol: float,
    allow_missing: bool,
) -> CheckResult:
    n_cells = len(ds.dimensions["nCells"])
    coszr = to_tc(get_var(ds, coszr_var), n_cells, coszr_var)
    night_mask = coszr <= night_threshold

    messages: list[str] = [
        (
            f"[INFO] night mask points: {int(np.count_nonzero(night_mask))} / "
            f"{night_mask.size} (threshold coszr<={night_threshold:.3f})"
        )
    ]

    failed = False
    for j_name in parse_csv_list(j_vars_csv):
        if j_name not in ds.variables:
            msg = f"[WARN] missing variable: {j_name}"
            if allow_missing:
                messages.append(msg)
                continue
            messages.append(msg.replace("[WARN]", "[FAIL]"))
            failed = True
            continue

        arr = np.asarray(ds.variables[j_name][:])
        if arr.ndim == 2:
            j_tc = to_tc(arr, n_cells, j_name)
            max_night = float(np.max(np.abs(j_tc[night_mask]))) if np.any(night_mask) else 0.0
        elif arr.ndim == 3:
            j_tcl = to_tcl(arr, n_cells, len(ds.dimensions["nVertLevels"]), j_name)
            max_night = (
                float(np.max(np.abs(j_tcl[night_mask, :])))
                if np.any(night_mask)
                else 0.0
            )
        else:
            messages.append(f"[FAIL] {j_name}: unsupported rank {arr.ndim}")
            failed = True
            continue

        ok = max_night <= abs_tol
        messages.append(
            f"[{'PASS' if ok else 'FAIL'}] {j_name}: night max |j|={max_night:.6e} tol={abs_tol:.3e}"
        )
        failed = failed or (not ok)

    return CheckResult(passed=not failed, messages=messages)


def parse_species_floors(raw: str) -> dict[str, float]:
    floors: dict[str, float] = {}
    for token in parse_csv_list(raw):
        if ":" not in token:
            raise ValueError(
                f"Invalid floor spec '{token}'. Expected format like 'qO:1e-18,qO3:1e-14'."
            )
        key, value = token.split(":", maxsplit=1)
        floors[key.strip()] = float(value.strip())
    return floors


def check_top_growth(
    ds,
    vars_csv: str,
    coszr_var: str,
    top_levels: int,
    night_threshold: float,
    min_night_samples: int,
    floors_csv: str,
    monotonic_fraction_threshold: float,
    max_domain_growth: float,
    max_p95_growth: float,
    max_mono_frac: float,
) -> CheckResult:
    n_cells = len(ds.dimensions["nCells"])
    n_levels = len(ds.dimensions["nVertLevels"])
    top_levels_use = min(max(top_levels, 1), n_levels)
    level_start = n_levels - top_levels_use
    floors = parse_species_floors(floors_csv)

    coszr = to_tc(get_var(ds, coszr_var), n_cells, coszr_var)
    night_mask = coszr <= night_threshold

    messages: list[str] = [
        (
            f"[INFO] top-growth setup: top_levels={top_levels_use} "
            f"night_points={int(np.count_nonzero(night_mask))}/{night_mask.size} "
            f"night_threshold={night_threshold:.3f}"
        )
    ]
    failed = False

    if not np.any(night_mask):
        return CheckResult(
            passed=False,
            messages=messages + ["[FAIL] no nighttime samples available for top-growth check."],
        )

    for var_name in parse_csv_list(vars_csv):
        arr = to_tcl(get_var(ds, var_name), n_cells, n_levels, var_name)
        top_arr = arr[:, :, level_start:]
        floor_val = max(floors.get(var_name, EPS), EPS)

        domain_series: list[float] = []
        for t in range(top_arr.shape[0]):
            night_cells = night_mask[t, :]
            if np.any(night_cells):
                domain_series.append(float(np.mean(top_arr[t, night_cells, :])))

        if len(domain_series) < 2:
            messages.append(
                f"[FAIL] {var_name}: insufficient nighttime domain samples ({len(domain_series)} < 2)."
            )
            failed = True
            continue

        domain_series_np = np.asarray(domain_series, dtype=float)
        domain_baseline = float(np.median(domain_series_np))
        domain_growth = float(
            (domain_series_np[-1] - domain_series_np[0]) / max(abs(domain_baseline), floor_val)
        )

        growth_list: list[float] = []
        mono_list: list[float] = []
        valid_count = 0

        for i_cell in range(n_cells):
            night_idx = np.flatnonzero(night_mask[:, i_cell])
            if night_idx.size < min_night_samples:
                continue
            for k_top in range(top_levels_use):
                series = top_arr[night_idx, i_cell, k_top]
                baseline = float(np.median(series))
                if baseline < floor_val:
                    continue
                valid_count += 1
                growth = float((series[-1] - series[0]) / max(baseline, floor_val))
                growth_list.append(max(growth, 0.0))
                if series.size <= 1:
                    mono = 0.0
                else:
                    mono = float(np.mean(np.diff(series) > 0.0))
                mono_list.append(mono)

        if valid_count == 0:
            messages.append(
                f"[FAIL] {var_name}: no valid top-layer points after floor/sample filters."
            )
            failed = True
            continue

        growth_np = np.asarray(growth_list, dtype=float)
        mono_np = np.asarray(mono_list, dtype=float)

        p95_growth = float(np.percentile(growth_np, 95.0))
        mono_frac = float(np.mean(mono_np >= monotonic_fraction_threshold))

        pass_domain = domain_growth <= max_domain_growth
        pass_p95 = p95_growth <= max_p95_growth
        pass_mono = mono_frac <= max_mono_frac
        failed = failed or (not (pass_domain and pass_p95 and pass_mono))

        messages.append(
            (
                f"[{'PASS' if pass_domain else 'FAIL'}] {var_name}: domain_growth={domain_growth:.6e} "
                f"threshold={max_domain_growth:.6e}"
            )
        )
        messages.append(
            (
                f"[{'PASS' if pass_p95 else 'FAIL'}] {var_name}: p95_growth={p95_growth:.6e} "
                f"threshold={max_p95_growth:.6e} valid_points={valid_count}"
            )
        )
        messages.append(
            (
                f"[{'PASS' if pass_mono else 'FAIL'}] {var_name}: monotonic_frac={mono_frac:.6e} "
                f"threshold={max_mono_frac:.6e} monotonic_score_threshold={monotonic_fraction_threshold:.3f}"
            )
        )

    return CheckResult(passed=not failed, messages=messages)


def check_transition_smooth(
    ds,
    coszr_var: str,
    j_vars_csv: str,
    transition_coszr: float,
    day_coszr: float,
    dt_seconds: float,
    j_floor: float,
    max_p99_jump: float,
    max_max_jump: float,
    max_p99_curvature: float,
    skip_curvature: bool,
    allow_missing: bool,
) -> CheckResult:
    n_cells = len(ds.dimensions["nCells"])
    n_levels = len(ds.dimensions["nVertLevels"])

    if dt_seconds <= 0.0:
        raise ValueError("--dt-seconds must be positive.")

    coszr = to_tc(get_var(ds, coszr_var), n_cells, coszr_var)
    transition_mask = np.abs(coszr) <= transition_coszr
    day_mask = coszr >= day_coszr

    messages: list[str] = [
        (
            f"[INFO] transition-smooth setup: transition |coszr|<={transition_coszr:.3f}, "
            f"day coszr>={day_coszr:.3f}, dt={dt_seconds:.3f}s"
        ),
        (
            f"[INFO] transition mask points: {int(np.count_nonzero(transition_mask))}/"
            f"{transition_mask.size}"
        ),
    ]
    failed = False

    jump_scale = 3.0 / dt_seconds
    curvature_scale = jump_scale * jump_scale

    for j_name in parse_csv_list(j_vars_csv):
        if j_name not in ds.variables:
            msg = f"[WARN] missing variable: {j_name}"
            if allow_missing:
                messages.append(msg)
                continue
            messages.append(msg.replace("[WARN]", "[FAIL]"))
            failed = True
            continue

        arr = np.asarray(ds.variables[j_name][:])
        if arr.ndim == 2:
            j_tc = to_tc(arr, n_cells, j_name)
            j_tcl = j_tc[:, :, None]
        elif arr.ndim == 3:
            j_tcl = to_tcl(arr, n_cells, n_levels, j_name)
        else:
            messages.append(f"[FAIL] {j_name}: unsupported rank {arr.ndim}")
            failed = True
            continue

        nt, n_cells_j, n_levels_j = j_tcl.shape
        if nt < 3:
            messages.append(
                f"[FAIL] {j_name}: need at least 3 timesteps for transition smoothness (got {nt})."
            )
            failed = True
            continue

        j_ref = np.full((n_cells_j, n_levels_j), max(j_floor, EPS), dtype=float)
        for i_cell in range(n_cells_j):
            day_idx = np.flatnonzero(day_mask[:, i_cell])
            if day_idx.size == 0:
                continue
            day_vals = j_tcl[day_idx, i_cell, :]
            ref_vals = np.percentile(day_vals, 99.0, axis=0)
            j_ref[i_cell, :] = np.maximum(ref_vals, max(j_floor, EPS))

        jump_values: list[float] = []
        curvature_values: list[float] = []
        transition_pairs = 0
        transition_triples = 0

        for i_cell in range(n_cells_j):
            trans_col = transition_mask[:, i_cell]
            pair_idx = np.flatnonzero(trans_col[:-1] & trans_col[1:])
            transition_pairs += int(pair_idx.size)
            for t in pair_idx:
                delta = np.abs(j_tcl[t + 1, i_cell, :] - j_tcl[t, i_cell, :]) / j_ref[i_cell, :]
                jump_values.extend((delta * jump_scale).tolist())

            if not skip_curvature:
                tri_idx = np.flatnonzero(trans_col[:-2] & trans_col[1:-1] & trans_col[2:])
                transition_triples += int(tri_idx.size)
                for t in tri_idx:
                    curv = np.abs(
                        j_tcl[t + 2, i_cell, :] - 2.0 * j_tcl[t + 1, i_cell, :] + j_tcl[t, i_cell, :]
                    ) / j_ref[i_cell, :]
                    curvature_values.extend((curv * curvature_scale).tolist())

        jump_np = np.asarray(jump_values, dtype=float)
        jump_np = jump_np[np.isfinite(jump_np)]
        if jump_np.size == 0:
            messages.append(
                f"[FAIL] {j_name}: no valid transition jump samples (pairs={transition_pairs})."
            )
            failed = True
            continue

        p99_jump = float(np.percentile(jump_np, 99.0))
        max_jump = float(np.max(jump_np))
        pass_p99_jump = p99_jump <= max_p99_jump
        pass_max_jump = max_jump <= max_max_jump
        failed = failed or (not (pass_p99_jump and pass_max_jump))

        messages.append(
            (
                f"[{'PASS' if pass_p99_jump else 'FAIL'}] {j_name}: p99_jump={p99_jump:.6e} "
                f"threshold={max_p99_jump:.6e}"
            )
        )
        messages.append(
            (
                f"[{'PASS' if pass_max_jump else 'FAIL'}] {j_name}: max_jump={max_jump:.6e} "
                f"threshold={max_max_jump:.6e}"
            )
        )

        if skip_curvature:
            messages.append(f"[INFO] {j_name}: curvature check skipped.")
            continue

        curv_np = np.asarray(curvature_values, dtype=float)
        curv_np = curv_np[np.isfinite(curv_np)]
        if curv_np.size == 0:
            messages.append(
                f"[FAIL] {j_name}: no valid transition curvature samples (triples={transition_triples})."
            )
            failed = True
            continue

        p99_curv = float(np.percentile(curv_np, 99.0))
        pass_curv = p99_curv <= max_p99_curvature
        failed = failed or (not pass_curv)
        messages.append(
            (
                f"[{'PASS' if pass_curv else 'FAIL'}] {j_name}: p99_curvature={p99_curv:.6e} "
                f"threshold={max_p99_curvature:.6e}"
            )
        )

    return CheckResult(passed=not failed, messages=messages)


def rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a.ravel()) + EPS
    return float(np.linalg.norm((a - b).ravel()) / denom)


def compare_history_fields(
    ds_a,
    ds_b,
    vars_csv: str,
    max_rel_l2: float,
    max_final_mean_rel: float,
    max_abs_diff: float | None = None,
) -> CheckResult:
    n_cells = min(len(ds_a.dimensions["nCells"]), len(ds_b.dimensions["nCells"]))
    n_levels = min(len(ds_a.dimensions["nVertLevels"]), len(ds_b.dimensions["nVertLevels"]))

    failed = False
    messages: list[str] = []

    for name in parse_csv_list(vars_csv):
        a = to_tcl_compare(
            get_var(ds_a, name),
            len(ds_a.dimensions["nCells"]),
            len(ds_a.dimensions["nVertLevels"]),
            f"{name}@A",
        )
        b = to_tcl_compare(
            get_var(ds_b, name),
            len(ds_b.dimensions["nCells"]),
            len(ds_b.dimensions["nVertLevels"]),
            f"{name}@B",
        )

        nt = min(a.shape[0], b.shape[0])
        n_levels_use = min(a.shape[2], b.shape[2], n_levels)
        a = a[:nt, :n_cells, :n_levels_use]
        b = b[:nt, :n_cells, :n_levels_use]

        l2 = rel_l2(a, b)
        max_abs = float(np.max(np.abs(a - b)))
        mean_a = np.mean(a, axis=(1, 2))
        mean_b = np.mean(b, axis=(1, 2))
        denom = max(abs(float(mean_a[-1])), abs(float(mean_b[-1])), EPS)
        final_rel = float(abs(mean_a[-1] - mean_b[-1]) / denom)

        ok_l2 = l2 <= max_rel_l2
        ok_final = final_rel <= max_final_mean_rel
        failed = failed or (not (ok_l2 and ok_final))

        messages.append(
            f"[{'PASS' if ok_l2 else 'FAIL'}] {name}: rel_l2={l2:.6e} threshold={max_rel_l2:.6e}"
        )
        messages.append(
            "[{}] {}: final mean rel diff={:.6e} threshold={:.6e}".format(
                "PASS" if ok_final else "FAIL",
                name,
                final_rel,
                max_final_mean_rel,
            )
        )

        if max_abs_diff is not None:
            ok_abs = max_abs <= max_abs_diff
            failed = failed or (not ok_abs)
            messages.append(
                f"[{'PASS' if ok_abs else 'FAIL'}] {name}: max abs diff={max_abs:.6e} "
                f"threshold={max_abs_diff:.6e}"
            )

    return CheckResult(passed=not failed, messages=messages)


def check_decomp_compare(
    ds_a,
    ds_b,
    vars_csv: str,
    max_rel_l2: float,
    max_final_mean_rel: float,
) -> CheckResult:
    return compare_history_fields(
        ds_a=ds_a,
        ds_b=ds_b,
        vars_csv=vars_csv,
        max_rel_l2=max_rel_l2,
        max_final_mean_rel=max_final_mean_rel,
    )


def check_fallback_compare(
    ds_a,
    ds_b,
    vars_csv: str,
    max_rel_l2: float,
    max_final_mean_rel: float,
    max_abs_diff: float,
) -> CheckResult:
    return compare_history_fields(
        ds_a=ds_a,
        ds_b=ds_b,
        vars_csv=vars_csv,
        max_rel_l2=max_rel_l2,
        max_final_mean_rel=max_final_mean_rel,
        max_abs_diff=max_abs_diff,
    )


def run_check(result: CheckResult) -> int:
    for line in result.messages:
        print(line)
    return 0 if result.passed else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run CheMPAS-A TUV-x phase-gate checks.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_nonneg = sub.add_parser("nonnegative", help="Check tracers are non-negative.")
    p_nonneg.add_argument("-i", "--input", required=True, help="History file.")
    p_nonneg.add_argument("--vars", default="qNO,qNO2,qO3", help="Comma-separated tracer vars.")
    p_nonneg.add_argument("--tol", type=float, default=0.0, help="Allowed negative tolerance.")

    p_budget = sub.add_parser("oxygen-budget", help="Legacy Chapman-style oxygen budget drift.")
    p_budget.add_argument("-i", "--input", required=True, help="History file.")
    p_budget.add_argument("--qO", default="qO")
    p_budget.add_argument("--qO2", default="qO2")
    p_budget.add_argument("--qO3", default="qO3")
    p_budget.add_argument("--rho-var", default="rho", help="Density variable name.")
    p_budget.add_argument("--zgrid-var", default="zgrid", help="Height-edge variable name.")
    p_budget.add_argument(
        "--max-domain-drift",
        type=float,
        default=5.0e-3,
        help="Maximum relative drift in domain-integrated oxygen budget.",
    )
    p_budget.add_argument(
        "--max-column-p95-drift",
        type=float,
        default=1.0e-2,
        help="Maximum p95 of per-column relative drift.",
    )

    p_night = sub.add_parser("night-jzero", help="Check nighttime photolysis is ~0.")
    p_night.add_argument("-i", "--input", required=True, help="History/diag file.")
    p_night.add_argument("--coszr-var", default="coszr")
    p_night.add_argument("--j-vars", default="j_no2")
    p_night.add_argument("--night-threshold", type=float, default=0.0)
    p_night.add_argument("--abs-tol", type=float, default=1.0e-20)
    p_night.add_argument(
        "--allow-missing",
        action="store_true",
        help="Do not fail when requested j variables are missing.",
    )

    p_top = sub.add_parser("top-growth", help="Legacy Chapman top-layer nighttime growth artifacts.")
    p_top.add_argument("-i", "--input", required=True, help="History file.")
    p_top.add_argument("--vars", default="qO,qO3", help="Comma-separated tracer vars.")
    p_top.add_argument("--coszr-var", default="coszr")
    p_top.add_argument("--top-levels", type=int, default=3)
    p_top.add_argument("--night-threshold", type=float, default=0.0)
    p_top.add_argument("--min-night-samples", type=int, default=8)
    p_top.add_argument("--floors", default="qO:1e-18,qO3:1e-14")
    p_top.add_argument("--monotonic-fraction-threshold", type=float, default=0.8)
    p_top.add_argument("--max-domain-growth", type=float, default=5.0e-2)
    p_top.add_argument("--max-p95-growth", type=float, default=2.0e-1)
    p_top.add_argument("--max-mono-frac", type=float, default=2.0e-1)

    p_trans = sub.add_parser(
        "transition-smooth", help="Check dawn/dusk photolysis transition smoothness."
    )
    p_trans.add_argument("-i", "--input", required=True, help="History/diag file.")
    p_trans.add_argument("--coszr-var", default="coszr")
    p_trans.add_argument("--j-vars", default="j_no2")
    p_trans.add_argument("--transition-coszr", type=float, default=0.08)
    p_trans.add_argument("--day-coszr", type=float, default=0.20)
    p_trans.add_argument("--dt-seconds", type=float, default=3.0)
    p_trans.add_argument("--j-floor", type=float, default=1.0e-12)
    p_trans.add_argument("--max-p99-jump", type=float, default=0.35)
    p_trans.add_argument("--max-max-jump", type=float, default=0.75)
    p_trans.add_argument("--max-p99-curvature", type=float, default=0.50)
    p_trans.add_argument(
        "--skip-curvature",
        action="store_true",
        help="Skip the second-difference curvature criterion.",
    )
    p_trans.add_argument(
        "--allow-missing",
        action="store_true",
        help="Do not fail when requested j variables are missing.",
    )

    p_cmp = sub.add_parser("decomp-compare", help="Compare two decompositions.")
    p_cmp.add_argument("-a", "--input-a", required=True, help="First history file.")
    p_cmp.add_argument("-b", "--input-b", required=True, help="Second history file.")
    p_cmp.add_argument("--vars", default="qNO,qNO2,qO3")
    p_cmp.add_argument("--max-rel-l2", type=float, default=1.0e-3)
    p_cmp.add_argument("--max-final-mean-rel", type=float, default=1.0e-3)

    p_fallback = sub.add_parser(
        "fallback-compare",
        help="Compare a fallback/Phase-1 run against a TUV-x-disabled reference run.",
    )
    p_fallback.add_argument("-a", "--input-a", required=True, help="First history file.")
    p_fallback.add_argument("-b", "--input-b", required=True, help="Second history file.")
    p_fallback.add_argument("--vars", default="qNO,qNO2,qO3")
    p_fallback.add_argument("--max-rel-l2", type=float, default=1.0e-12)
    p_fallback.add_argument("--max-final-mean-rel", type=float, default=1.0e-12)
    p_fallback.add_argument("--max-abs-diff", type=float, default=1.0e-18)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "nonnegative":
            with open_dataset(args.input) as ds:
                return run_check(check_nonnegative(ds, args.vars, args.tol))

        if args.command == "oxygen-budget":
            with open_dataset(args.input) as ds:
                return run_check(
                    check_oxygen_budget(
                        ds=ds,
                        qo=args.qO,
                        qo2=args.qO2,
                        qo3=args.qO3,
                        rho_var=args.rho_var,
                        zgrid_var=args.zgrid_var,
                        max_domain_drift=args.max_domain_drift,
                        max_column_p95_drift=args.max_column_p95_drift,
                    )
                )

        if args.command == "night-jzero":
            with open_dataset(args.input) as ds:
                return run_check(
                    check_night_jzero(
                        ds=ds,
                        coszr_var=args.coszr_var,
                        j_vars_csv=args.j_vars,
                        night_threshold=args.night_threshold,
                        abs_tol=args.abs_tol,
                        allow_missing=args.allow_missing,
                    )
                )

        if args.command == "decomp-compare":
            with open_dataset(args.input_a) as ds_a, open_dataset(args.input_b) as ds_b:
                return run_check(
                    check_decomp_compare(
                        ds_a=ds_a,
                        ds_b=ds_b,
                        vars_csv=args.vars,
                        max_rel_l2=args.max_rel_l2,
                        max_final_mean_rel=args.max_final_mean_rel,
                    )
                )

        if args.command == "fallback-compare":
            with open_dataset(args.input_a) as ds_a, open_dataset(args.input_b) as ds_b:
                return run_check(
                    check_fallback_compare(
                        ds_a=ds_a,
                        ds_b=ds_b,
                        vars_csv=args.vars,
                        max_rel_l2=args.max_rel_l2,
                        max_final_mean_rel=args.max_final_mean_rel,
                        max_abs_diff=args.max_abs_diff,
                    )
                )

        if args.command == "top-growth":
            with open_dataset(args.input) as ds:
                return run_check(
                    check_top_growth(
                        ds=ds,
                        vars_csv=args.vars,
                        coszr_var=args.coszr_var,
                        top_levels=args.top_levels,
                        night_threshold=args.night_threshold,
                        min_night_samples=args.min_night_samples,
                        floors_csv=args.floors,
                        monotonic_fraction_threshold=args.monotonic_fraction_threshold,
                        max_domain_growth=args.max_domain_growth,
                        max_p95_growth=args.max_p95_growth,
                        max_mono_frac=args.max_mono_frac,
                    )
                )

        if args.command == "transition-smooth":
            with open_dataset(args.input) as ds:
                return run_check(
                    check_transition_smooth(
                        ds=ds,
                        coszr_var=args.coszr_var,
                        j_vars_csv=args.j_vars,
                        transition_coszr=args.transition_coszr,
                        day_coszr=args.day_coszr,
                        dt_seconds=args.dt_seconds,
                        j_floor=args.j_floor,
                        max_p99_jump=args.max_p99_jump,
                        max_max_jump=args.max_max_jump,
                        max_p99_curvature=args.max_p99_curvature,
                        skip_curvature=args.skip_curvature,
                        allow_missing=args.allow_missing,
                    )
                )

        parser.error(f"Unhandled command: {args.command}")
        return 2
    except Exception as exc:  # pylint: disable=broad-except
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
