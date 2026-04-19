"""Plot the three idealized test cases from their output.nc files."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def cell_height(d):
    zg = d.variables["zgrid"][:]
    return 0.5 * (zg[:, :-1] + zg[:, 1:])


def xtime_strings(d):
    xt = d.variables["xtime"][:]
    return [b"".join(row).decode(errors="ignore").strip() for row in xt]


def plot_supercell(path: Path, outdir: Path) -> None:
    d = nc.Dataset(path)
    xC = np.asarray(d.variables["xCell"][:] / 1e3)
    yC = np.asarray(d.variables["yCell"][:] / 1e3)
    w = d.variables["w"][:]                   # (T, nC, nVL+1)
    qAB = d.variables["qAB"][:]               # (T, nC, nVL)
    qA = d.variables["qA"][:]
    theta = d.variables["theta"][:]
    qc = d.variables["qc"][:] if "qc" in d.variables else None
    qr = d.variables["qr"][:] if "qr" in d.variables else None
    zc = cell_height(d)
    zg = d.variables["zgrid"][:]
    times = xtime_strings(d)
    d.close()

    nT = qAB.shape[0]
    z_col = np.asarray(zc.mean(axis=0) / 1e3)          # km
    z_w = np.asarray(zg.mean(axis=0) / 1e3)

    # Panel 1: time evolution of max w and storm diagnostics
    t_min = np.arange(nT) * 2.0   # 2-min output interval
    max_w = np.max(w, axis=(1, 2))
    min_w = np.min(w, axis=(1, 2))
    qAB_min = np.min(qAB, axis=(1, 2))
    qAB_max = np.max(qAB, axis=(1, 2))
    qAB_mean = qAB.mean(axis=(1, 2))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    ax.plot(t_min, max_w, "r-", label="max w")
    ax.plot(t_min, min_w, "b-", label="min w")
    ax.set_xlabel("time (min)"); ax.set_ylabel("w (m/s)")
    ax.set_title("Vertical velocity extremes")
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(t_min, qAB_min, "b-", label="min qAB")
    ax.plot(t_min, qAB_max, "r-", label="max qAB")
    ax.plot(t_min, qAB_mean, "k--", label="mean qAB")
    ax.set_xlabel("time (min)"); ax.set_ylabel("qAB (kg/kg)")
    ax.set_title("qAB tracer extremes (mixing range)")
    ax.legend(); ax.grid(True, alpha=0.3)

    # Vertical profile: max|w| vs z, several times
    ax = axes[1, 0]
    for ti in np.linspace(0, nT - 1, 5).astype(int):
        ax.plot(np.max(np.abs(w[ti]), axis=0), z_w,
                label=f"t={t_min[ti]:.0f} min")
    ax.set_xlabel("max |w| (m/s)"); ax.set_ylabel("z (km)")
    ax.set_title("Domain-max |w| vertical profile")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Cloud/rain amount vs time
    ax = axes[1, 1]
    if qc is not None:
        ax.plot(t_min, np.max(qc, axis=(1, 2)) * 1e3, "b-", label="max qc")
    if qr is not None:
        ax.plot(t_min, np.max(qr, axis=(1, 2)) * 1e3, "g-", label="max qr")
    ax.plot(t_min, (theta.max(axis=(1, 2)) - theta[0].max()), "r-", label="Δ max θ")
    ax.set_xlabel("time (min)"); ax.set_ylabel("g/kg  /  K")
    ax.set_title("Cloud water, rain, θ change")
    ax.legend(); ax.grid(True, alpha=0.3)

    fig.suptitle("Supercell — full 2 h run", y=1.00)
    fig.tight_layout()
    out = outdir / "supercell_timeseries.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # Horizontal slices of qAB at z≈5 km at 6 times
    k5 = int(np.argmin(np.abs(z_col - 5.0)))
    k1 = int(np.argmin(np.abs(z_col - 1.0)))
    k10 = int(np.argmin(np.abs(z_col - 10.0)))

    show_times = np.linspace(0, nT - 1, 6).astype(int)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    flat_axes = list(axes.flat)
    qmin, qmax = float(qAB.min()), float(qAB.max())
    for ax, ti in zip(flat_axes, show_times):
        sc = ax.scatter(xC, yC, c=qAB[ti, :, k5], s=2, cmap="viridis",
                        vmin=qmin, vmax=qmax)
        ax.set_aspect("equal")
        ax.set_title(f"t = {t_min[ti]:.0f} min")
        ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
    fig.colorbar(sc, ax=flat_axes, label=f"qAB @ z≈{z_col[k5]:.1f} km", shrink=0.8)
    fig.suptitle("Supercell qAB advection (horizontal slice at 5 km)", y=1.00)
    out2 = outdir / "supercell_qAB_slices.png"
    fig.savefig(out2, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out2}")

    # Vertical x-z cross-section of w and qAB at final time through max updraft
    i_peak = int(np.argmax(np.max(w[-1], axis=1)))
    y_band = yC[i_peak]
    y_mask = np.abs(yC - y_band) < 1.0
    idx = np.where(y_mask)[0]
    idx = idx[np.argsort(xC[idx])]
    x_line = np.asarray(xC[idx], dtype=float)
    z_line = np.asarray(zg[idx].mean(axis=0) / 1e3, dtype=float)
    z_line_c = np.asarray(zc[idx].mean(axis=0) / 1e3, dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    W = np.asarray(w[-1, idx, :].T, dtype=float)
    wamp = float(np.max(np.abs(W)))
    pcm = axes[0].pcolormesh(x_line, z_line, W, cmap="RdBu_r",
                             vmin=-wamp, vmax=wamp, shading="auto")
    axes[0].set_ylabel("z (km)"); axes[0].set_title(f"w  —  t = {times[-1]}  (y≈{y_band:.1f} km)")
    plt.colorbar(pcm, ax=axes[0], label="m/s")

    Q = np.asarray(qAB[-1, idx, :].T, dtype=float)
    pcm = axes[1].pcolormesh(x_line, z_line_c, Q, cmap="viridis", shading="auto")
    axes[1].set_ylabel("z (km)"); axes[1].set_xlabel("x (km)")
    axes[1].set_title("qAB (kg/kg)")
    plt.colorbar(pcm, ax=axes[1], label="qAB")
    fig.tight_layout()
    out3 = outdir / "supercell_xsection_final.png"
    fig.savefig(out3, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out3}  | peak w @ final = {W.max():.1f} m/s")

    print(f"[supercell] qA range {qA.min():.3e} .. {qA.max():.3e}  — chemistry {'active' if qA.max() > 0 else 'INACTIVE'}")


def plot_mountain_wave(path: Path, outdir: Path) -> None:
    d = nc.Dataset(path)
    xC = np.asarray(d.variables["xCell"][:] / 1e3)
    yC = np.asarray(d.variables["yCell"][:] / 1e3)
    w = d.variables["w"][:]
    qAB = d.variables["qAB"][:]
    qA = d.variables["qA"][:]
    zg = d.variables["zgrid"][:]
    zc = cell_height(d)
    times = xtime_strings(d)
    d.close()

    nT = w.shape[0]

    # y-centerline cells
    y_mid = 0.5 * (yC.min() + yC.max())
    y_unique = np.unique(np.round(yC, 1))
    y_target = y_unique[np.argmin(np.abs(y_unique - y_mid))]
    idx = np.where(np.isclose(yC, y_target, atol=50.0 / 1e3))[0]
    idx = idx[np.argsort(xC[idx])]
    x_line = np.asarray(xC[idx], dtype=float)
    z_line = np.asarray(zg[idx].mean(axis=0) / 1e3, dtype=float)
    z_line_c = np.asarray(zc[idx].mean(axis=0) / 1e3, dtype=float)

    show_times = np.linspace(0, nT - 1, 6).astype(int)

    # w cross-section evolution
    fig, axes = plt.subplots(len(show_times), 1, figsize=(11, 2.5 * len(show_times)), sharex=True)
    wmax = float(np.max(np.abs(w)))
    for ax, ti in zip(axes, show_times):
        W = np.asarray(w[ti, idx, :].T, dtype=float)
        pcm = ax.pcolormesh(x_line, z_line, W, cmap="RdBu_r",
                            vmin=-wmax, vmax=wmax, shading="auto")
        ax.set_ylabel("z (km)")
        ax.set_title(f"w  —  t = {times[ti]}")
        plt.colorbar(pcm, ax=ax, label="m/s")
    axes[-1].set_xlabel("x (km)")
    fig.tight_layout()
    out = outdir / "mountain_wave_w_evolution.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # qAB cross-section evolution (tracer advection by wave)
    fig, axes = plt.subplots(len(show_times), 1, figsize=(11, 2.5 * len(show_times)), sharex=True)
    qmin, qmax = float(qAB.min()), float(qAB.max())
    for ax, ti in zip(axes, show_times):
        Q = np.asarray(qAB[ti, idx, :].T, dtype=float)
        pcm = ax.pcolormesh(x_line, z_line_c, Q, cmap="viridis",
                            vmin=qmin, vmax=qmax, shading="auto")
        ax.set_ylabel("z (km)")
        ax.set_title(f"qAB  —  t = {times[ti]}")
        plt.colorbar(pcm, ax=ax, label="kg/kg")
    axes[-1].set_xlabel("x (km)")
    fig.tight_layout()
    out2 = outdir / "mountain_wave_qAB_evolution.png"
    fig.savefig(out2, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out2}")

    # Time series
    t_min = np.arange(nT) * 10.0   # 10-min output interval
    max_w = np.max(np.abs(w), axis=(1, 2))
    qAB_spread = qAB.max(axis=(1, 2)) - qAB.min(axis=(1, 2))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(t_min / 60.0, max_w, "k-")
    axes[0].set_xlabel("time (h)"); axes[0].set_ylabel("max |w| (m/s)")
    axes[0].set_title("Domain-max |w|"); axes[0].grid(True, alpha=0.3)
    axes[1].plot(t_min / 60.0, qAB_spread, "b-", label="max − min qAB")
    axes[1].plot(t_min / 60.0, qAB.mean(axis=(1, 2)), "k--", label="mean qAB")
    axes[1].set_xlabel("time (h)"); axes[1].set_ylabel("kg/kg")
    axes[1].set_title("qAB spread & mean"); axes[1].legend(); axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    out3 = outdir / "mountain_wave_timeseries.png"
    fig.savefig(out3, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out3}  | final-time max|w| = {max_w[-1]:.3f} m/s")
    print(f"[mountain_wave] qA range {qA.min():.3e} .. {qA.max():.3e}  — chemistry {'active' if qA.max() > 0 else 'INACTIVE'}")


def plot_jw_baroclinic(path: Path, outdir: Path) -> None:
    d = nc.Dataset(path)
    lat = np.degrees(d.variables["latCell"][:])
    lon = np.degrees(d.variables["lonCell"][:])
    lon = np.where(lon > 180, lon - 360, lon)
    ps = d.variables["surface_pressure"][:]
    theta = d.variables["theta"][:]
    pres = d.variables["pressure"][:]
    qAB = d.variables["qAB"][:]
    qA = d.variables["qA"][:]
    times = xtime_strings(d)
    d.close()

    nT = ps.shape[0]

    def theta_on_pressure(pres_col, theta_col, p_target):
        out = np.empty(pres_col.shape[0])
        for i in range(pres_col.shape[0]):
            p = pres_col[i]; th = theta_col[i]
            order = np.argsort(p)
            out[i] = np.interp(p_target, p[order], th[order])
        return out

    show_times = np.linspace(0, nT - 1, min(nT, 4)).astype(int)
    fig, axes = plt.subplots(len(show_times), 3, figsize=(15, 4 * len(show_times)))
    if len(show_times) == 1:
        axes = axes.reshape(1, 3)
    for row, ti in enumerate(show_times):
        sc0 = axes[row, 0].scatter(lon, lat, c=ps[ti] / 1e2, s=2, cmap="viridis")
        axes[row, 0].set_title(f"ps (hPa)  t = {times[ti]}")
        axes[row, 0].set_xlabel("lon"); axes[row, 0].set_ylabel("lat")
        axes[row, 0].set_xlim(-180, 180); axes[row, 0].set_ylim(-90, 90)
        plt.colorbar(sc0, ax=axes[row, 0])

        th850 = theta_on_pressure(pres[ti], theta[ti], 8.5e4)
        sc1 = axes[row, 1].scatter(lon, lat, c=th850, s=2, cmap="RdYlBu_r")
        axes[row, 1].set_title(f"θ @ 850 hPa  t = {times[ti]}")
        axes[row, 1].set_xlabel("lon")
        axes[row, 1].set_xlim(-180, 180); axes[row, 1].set_ylim(-90, 90)
        plt.colorbar(sc1, ax=axes[row, 1])

        # qAB at ~5 km (mid-level) — use column index via pressure level ~500 hPa
        k_mid = qAB.shape[2] // 2
        sc2 = axes[row, 2].scatter(lon, lat, c=qAB[ti, :, k_mid], s=2, cmap="viridis")
        axes[row, 2].set_title(f"qAB mid-level  t = {times[ti]}")
        axes[row, 2].set_xlabel("lon")
        axes[row, 2].set_xlim(-180, 180); axes[row, 2].set_ylim(-90, 90)
        plt.colorbar(sc2, ax=axes[row, 2])
    fig.tight_layout()
    out = outdir / "jw_baroclinic_evolution.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"[jw_baroclinic] qA range {qA.min():.3e} .. {qA.max():.3e}  — chemistry {'active' if qA.max() > 0 else 'INACTIVE'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/home/fillmore/Data/CheMPAS")
    ap.add_argument("--only", choices=["supercell", "mountain_wave", "jw_baroclinic", "all"],
                    default="all")
    args = ap.parse_args()
    root = Path(args.data_root)

    if args.only in ("supercell", "all"):
        plot_supercell(root / "supercell/output.nc", root / "supercell/plots")
    if args.only in ("mountain_wave", "all"):
        plot_mountain_wave(root / "mountain_wave/output.nc", root / "mountain_wave/plots")
    if args.only in ("jw_baroclinic", "all"):
        plot_jw_baroclinic(root / "jw_baroclinic_wave/output.nc", root / "jw_baroclinic_wave/plots")


if __name__ == "__main__":
    main()
