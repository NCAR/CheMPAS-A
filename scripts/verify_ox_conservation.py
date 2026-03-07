#!/usr/bin/env python3
"""
Verify domain-integrated Ox (O3 + NO2) conservation in MPAS-MUSICA output.

Ox = O3 + NO2 is conserved by the LNOx-O3 chemistry when the lightning source
and NOx sink are both disabled.  Advection redistributes species but conserves
domain totals, so the mass-weighted domain integral of Ox should be constant.

The script computes:
    Ox_mass(t) = sum over cells,levels of (qO3/M_O3 + qNO2/M_NO2) * rho * dV

where dV = areaCell * dz (layer thickness from zgrid).

Usage:
    python verify_ox_conservation.py                    # default output.nc
    python verify_ox_conservation.py -i output.nc       # explicit path
    python verify_ox_conservation.py --plot              # save plot
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# Add scripts dir to path for style module
sys.path.insert(0, str(Path(__file__).parent))

# Molar masses [kg/mol]
M_NO = 0.030
M_NO2 = 0.046
M_O3 = 0.048
M_AIR = 0.029


def compute_domain_ox(ds):
    """Compute domain-integrated Ox [mol] at each timestep."""
    areaCell = ds.variables['areaCell'][:]          # (nCells,)
    zgrid = ds.variables['zgrid'][:]                # (nCells, nVertLevelsP1)
    dz = np.diff(zgrid, axis=1)                     # (nCells, nVertLevels)

    nTimes = ds.dimensions['Time'].size
    ox_total = np.zeros(nTimes)
    no_total = np.zeros(nTimes)
    no2_total = np.zeros(nTimes)
    o3_total = np.zeros(nTimes)

    for t in range(nTimes):
        rho = ds.variables['rho'][t, :, :]          # (nCells, nVertLevels)
        qNO = ds.variables['qNO'][t, :, :]
        qNO2 = ds.variables['qNO2'][t, :, :]
        qO3 = ds.variables['qO3'][t, :, :]

        # dV = areaCell * dz  (broadcast areaCell over levels)
        dV = areaCell[:, np.newaxis] * dz            # (nCells, nVertLevels)

        # Mass in each cell-level [kg of air]
        dm = rho * dV

        # Moles of each species = (mass mixing ratio / molar mass) * air mass [kg]
        # Actually: q [kg species / kg air], so moles = q * dm / M_species
        no_mol = np.sum(qNO * dm / M_NO)
        no2_mol = np.sum(qNO2 * dm / M_NO2)
        o3_mol = np.sum(qO3 * dm / M_O3)

        no_total[t] = no_mol
        no2_total[t] = no2_mol
        o3_total[t] = o3_mol
        ox_total[t] = o3_mol + no2_mol

    return ox_total, no_total, no2_total, o3_total


def main():
    parser = argparse.ArgumentParser(
        description='Verify domain-integrated Ox conservation')
    parser.add_argument('-i', '--input', default='output.nc',
                        help='Path to MPAS output file')
    parser.add_argument('--plot', action='store_true',
                        help='Save a time-series plot')
    parser.add_argument('-o', '--output', default='ox_conservation',
                        help='Output filename stem (default: ox_conservation)')
    parser.add_argument('--max-drift-pct', type=float, default=0.01,
                        help='Maximum allowed absolute Ox drift percentage for PASS')
    parser.add_argument('--warn-drift-pct', type=float, default=0.1,
                        help='Warning threshold used for report text when PASS is not met')
    args = parser.parse_args()

    filepath = Path(args.input)
    if not filepath.exists():
        raise SystemExit(f"File not found: {filepath}")

    try:
        from netCDF4 import Dataset  # type: ignore
    except Exception as exc:
        raise SystemExit(
            "Python module 'netCDF4' is required for verify_ox_conservation.py. "
            "Install with: pip install netCDF4"
        ) from exc

    print(f"Reading {filepath} ...")
    with Dataset(filepath) as ds:
        ox, no, no2, o3 = compute_domain_ox(ds)
        dt_output = 30.0  # seconds between output frames
        nTimes = ds.dimensions['Time'].size

    times_min = np.arange(nTimes) * dt_output / 60.0

    # Report
    ox0 = ox[0]
    ox_rel = (ox - ox0) / ox0 * 100.0
    print(f"\nDomain-integrated Ox (O3 + NO2):")
    print(f"  t=0:   {ox0:.6e} mol")
    print(f"  t=end: {ox[-1]:.6e} mol")
    print(f"  Drift: {ox[-1] - ox0:.4e} mol ({ox_rel[-1]:+.6f}%)")
    print(f"  Max deviation: {np.max(np.abs(ox_rel)):.6f}%")

    # Also report NO total (should decrease if titration converts NO->NO2)
    print(f"\nDomain-integrated NO:")
    print(f"  t=0:   {no[0]:.6e} mol")
    print(f"  t=end: {no[-1]:.6e} mol")

    print(f"\nDomain-integrated NOx (NO + NO2):")
    nox = no + no2
    nox_rel = (nox - nox[0]) / nox[0] * 100.0 if nox[0] > 0 else np.zeros_like(nox)
    print(f"  t=0:   {nox[0]:.6e} mol")
    print(f"  t=end: {nox[-1]:.6e} mol")
    print(f"  Drift: {nox[-1] - nox[0]:.4e} mol ({nox_rel[-1]:+.6f}%)")

    # Pass/fail
    max_drift_pct = np.max(np.abs(ox_rel))
    status = 0
    if max_drift_pct <= args.max_drift_pct:
        print(f"\n*** PASS: Ox conserved to {max_drift_pct:.4f}% ***")
    elif max_drift_pct <= args.warn_drift_pct:
        print(f"\n*** FAIL: Ox drift {max_drift_pct:.4f}% exceeds PASS threshold "
              f"{args.max_drift_pct:.4f}% but is below warning threshold "
              f"{args.warn_drift_pct:.4f}% ***")
        status = 2
    else:
        print(f"\n*** FAIL: Ox drift {max_drift_pct:.4f}% exceeds warning threshold "
              f"{args.warn_drift_pct:.4f}% ***")
        status = 2

    if args.plot:
        try:
            import style
            style.apply()
        except Exception:
            pass
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        # Top: absolute Ox
        ax = axes[0]
        ax.plot(times_min, o3, label='O3', color='tab:blue')
        ax.plot(times_min, no2, label='NO2', color='tab:orange')
        ax.plot(times_min, ox, label='Ox (O3+NO2)', color='black', lw=2)
        ax.plot(times_min, no, label='NO', color='tab:red', ls='--')
        ax.set_ylabel('Domain total [mol]')
        ax.legend(fontsize=8)
        ax.set_title('Domain-Integrated Species')

        # Bottom: relative Ox drift
        ax = axes[1]
        ax.plot(times_min, ox_rel, color='black', lw=2)
        ax.axhline(0, color='gray', ls=':', lw=0.5)
        ax.set_ylabel('Ox drift [%]')
        ax.set_xlabel('Time [min]')
        ax.set_title(f'Ox Conservation (max drift: {max_drift_pct:.4f}%)')

        plt.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(f'{args.output}.{ext}', dpi=150)
            print(f"Saved {args.output}.{ext}")
        plt.close()

    return status


if __name__ == '__main__':
    raise SystemExit(main())
