#!/usr/bin/env python3
"""
Visualize MPAS-A chemistry tracer output (A, B, AB concentrations).

Produces:
- Horizontal slices at specified vertical levels
- Vertical cross-sections through domain center
- Time evolution panels

Usage:
    python plot_chemistry.py                    # Default: output.nc
    python plot_chemistry.py -i output.nc -r output_reference.nc  # With reference
    python plot_chemistry.py --level 20 --time 2  # Specific level/time
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset


def save_figure(output_file: str, dpi: int = 300) -> None:
    """Save figure to both PNG and PDF formats."""
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    pdf_file = output_file.replace('.png', '.pdf')
    plt.savefig(pdf_file, bbox_inches='tight')
    print(f"Saved {output_file} and {pdf_file}")


def get_level_height(data: dict, level: int) -> float:
    """Get approximate height in km for a vertical level."""
    if data.get('zCell') is not None:
        # Use first cell's zgrid as reference
        return float(data['zCell'][0, level])
    else:
        # Fallback: assume 0-20 km over nVertLevels
        return level * 20.0 / data['nVertLevels']


def level_label(data: dict, level: int) -> str:
    """Return a label with both level index and height."""
    height = get_level_height(data, level)
    return f"{height:.1f} km"


def load_mpas_data(filename: str) -> dict:
    """Load MPAS output data."""
    ds = Dataset(filename, 'r')
    data = {
        'latCell': np.degrees(ds.variables['latCell'][:]),
        'lonCell': np.degrees(ds.variables['lonCell'][:]),
        'xCell': ds.variables['xCell'][:] / 1000,  # km
        'yCell': ds.variables['yCell'][:] / 1000,  # km
        'zCell': ds.variables['zgrid'][:] / 1000 if 'zgrid' in ds.variables else None,  # km
        'qAB': ds.variables['qAB'][:],
        'qA': ds.variables['qA'][:],
        'qB': ds.variables['qB'][:],
        'nCells': len(ds.dimensions['nCells']),
        'nVertLevels': len(ds.dimensions['nVertLevels']),
        'nTimes': len(ds.dimensions['Time']),
        'times': np.arange(len(ds.dimensions['Time'])),
    }
    # Load wind fields if available
    if 'uReconstructZonal' in ds.variables:
        data['uZonal'] = ds.variables['uReconstructZonal'][:]
        data['uMeridional'] = ds.variables['uReconstructMeridional'][:]
    if 'w' in ds.variables:
        # Trim w to match nVertLevels (w has nVertLevelsP1)
        data['w'] = ds.variables['w'][:, :, :-1]
    ds.close()
    return data


def plot_horizontal_slice(data: dict, tracer: str, level: int, time_idx: int,
                          ax: plt.Axes, title: str = None, vmin: float = None,
                          vmax: float = None, show_wind: bool = False,
                          wind_skip: int = 50, wind_scale: float = 300) -> None:
    """Plot horizontal slice of tracer at given level and time.

    Args:
        show_wind: If True, overlay wind vectors
        wind_skip: Plot every Nth cell for wind vectors (reduces clutter)
        wind_scale: Scale factor for quiver (larger = shorter arrows)
    """
    x = data['xCell']
    y = data['yCell']
    values = data[tracer][time_idx, :, level]

    # Create triangulation for unstructured mesh
    tri = Triangulation(x, y)

    # Plot
    cf = ax.tricontourf(tri, values, levels=50, cmap='viridis',
                        vmin=vmin, vmax=vmax)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_aspect('equal')

    # Add wind vectors if requested
    if show_wind and 'uZonal' in data:
        u = data['uZonal'][time_idx, :, level]
        v = data['uMeridional'][time_idx, :, level]
        # Subsample for clarity
        idx = np.arange(0, len(x), wind_skip)
        ax.quiver(x[idx], y[idx], u[idx], v[idx],
                  scale=wind_scale, alpha=0.7, width=0.003)

    if title:
        ax.set_title(title)
    plt.colorbar(cf, ax=ax, label=f'{tracer} (kg/kg)')


def plot_vertical_cross_section(data: dict, tracer: str, time_idx: int,
                                 y_slice: float, ax: plt.Axes,
                                 title: str = None) -> None:
    """Plot vertical cross-section at given y position."""
    x = data['xCell']
    y = data['yCell']

    # Find cells near the y-slice
    y_tol = (y.max() - y.min()) / 50
    mask = np.abs(y - y_slice) < y_tol

    if not np.any(mask):
        ax.text(0.5, 0.5, 'No cells in slice', transform=ax.transAxes, ha='center')
        return

    x_slice = x[mask]
    values = data[tracer][time_idx, mask, :]  # (nCells_slice, nVertLevels)

    # Sort by x for cleaner plotting
    sort_idx = np.argsort(x_slice)
    x_slice = x_slice[sort_idx]
    values = values[sort_idx, :]

    # Create height coordinate (approximate)
    nLevels = values.shape[1]
    z = np.linspace(0, 20, nLevels)  # Approximate 0-20 km

    X, Z = np.meshgrid(x_slice, z)

    cf = ax.contourf(X, Z, values.T, levels=50, cmap='viridis')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Height (km)')
    if title:
        ax.set_title(title)
    plt.colorbar(cf, ax=ax, label=f'{tracer} (kg/kg)')


def plot_time_evolution(data: dict, tracer: str, level: int, ax: plt.Axes,
                        color: str = 'b', show_legend: bool = True) -> None:
    """Plot domain-mean tracer evolution over time."""
    times = data['times']
    values = data[tracer][:, :, level]  # (nTimes, nCells)

    mean_val = values.mean(axis=1)
    min_val = values.min(axis=1)
    max_val = values.max(axis=1)

    ax.fill_between(times, min_val, max_val, alpha=0.3, color=color, label='min-max')
    ax.plot(times, mean_val, color=color, linestyle='-', linewidth=2, label='mean')
    ax.set_xlabel('Time step')
    ax.set_ylabel(f'{tracer} (kg/kg)')
    if show_legend:
        ax.legend()
    ax.set_title(f'{tracer} at {level_label(data, level)}')


def plot_time_series_standalone(data: dict, level: int, output_file: str,
                                 time_range: tuple = None, dt_seconds: float = 30.0) -> None:
    """Generate standalone time series plots for all tracers.

    Args:
        time_range: Optional (start, end) tuple to limit time axis (in time steps)
        dt_seconds: Time interval between output steps in seconds
    """
    tracers = ['qAB', 'qA', 'qB']
    colors = ['purple', 'blue', 'red']
    times = data['times']

    # Determine time slice
    if time_range:
        t_start, t_end = time_range
        t_mask = (times >= t_start) & (times <= t_end)
    else:
        t_mask = np.ones(len(times), dtype=bool)

    # Convert to minutes for plotting
    times_minutes = times * dt_seconds / 60.0

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    for ax, tracer, color in zip(axes, tracers, colors):
        values = data[tracer][:, :, level]  # (nTimes, nCells)
        mean_val = values.mean(axis=1)[t_mask]
        min_val = values.min(axis=1)[t_mask]
        max_val = values.max(axis=1)[t_mask]
        t_plot = times_minutes[t_mask]

        ax.fill_between(t_plot, min_val, max_val, alpha=0.3, color=color, label='min-max')
        ax.plot(t_plot, mean_val, color=color, linestyle='-', linewidth=2, label='mean')
        ax.set_xlabel('Time (minutes)')
        ax.set_ylabel(f'{tracer} (kg/kg)')
        ax.legend()
        ax.set_title(f'{tracer} at {level_label(data, level)}')

    t_min = times_minutes[t_mask].min()
    t_max = times_minutes[t_mask].max()
    plt.suptitle(f'Tracer Time Evolution ({level_label(data, level)}, t={t_min:.1f}-{t_max:.1f} min)', fontsize=14)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_single_cell_timeseries(data: dict, levels: list, output_file: str,
                                 cell_idx: int = None, time_range: tuple = None,
                                 dt_seconds: float = 30.0) -> None:
    """Plot time series of AB and A at a single grid cell for multiple levels.

    Args:
        levels: List of vertical levels to plot
        cell_idx: Grid cell index (if None, use cell nearest domain center)
        time_range: Optional (start, end) tuple for time axis
        dt_seconds: Time interval between output steps
    """
    x = data['xCell']
    y = data['yCell']
    times = data['times']

    # Find center cell if not specified
    if cell_idx is None:
        x_center = (x.min() + x.max()) / 2
        y_center = (y.min() + y.max()) / 2
        distances = np.sqrt((x - x_center)**2 + (y - y_center)**2)
        cell_idx = np.argmin(distances)

    # Time range
    if time_range:
        t_start, t_end = time_range
        t_mask = (times >= t_start) & (times <= t_end)
    else:
        t_mask = np.ones(len(times), dtype=bool)

    times_minutes = times[t_mask] * dt_seconds / 60.0

    # Create figure: 2 columns (qAB, qA), rows = levels
    fig, axes = plt.subplots(len(levels), 2, figsize=(10, 3 * len(levels)), sharex=True)
    if len(levels) == 1:
        axes = axes.reshape(1, -1)

    tracers = ['qAB', 'qA']
    colors = ['purple', 'blue']

    for row_idx, level in enumerate(levels):
        for col_idx, (tracer, color) in enumerate(zip(tracers, colors)):
            ax = axes[row_idx, col_idx]
            values = data[tracer][t_mask, cell_idx, level]

            ax.plot(times_minutes, values, color=color, linewidth=2)
            ax.set_ylabel(f'{tracer} (kg/kg)')
            ax.grid(True, alpha=0.3)

            if row_idx == 0:
                ax.set_title(f'{tracer}')
            if row_idx == len(levels) - 1:
                ax.set_xlabel('Time (minutes)')

            # Add height label on left side
            if col_idx == 0:
                height_str = level_label(data, level)
                ax.annotate(height_str, xy=(-0.15, 0.5), xycoords='axes fraction',
                           fontsize=10, ha='right', va='center', rotation=90)

    cell_x, cell_y = x[cell_idx], y[cell_idx]
    plt.suptitle(f'Single Cell Time Series (cell {cell_idx}: x={cell_x:.1f}, y={cell_y:.1f} km)',
                 fontsize=12)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_vertical_structure(data: dict, time_idx: int, output_file: str,
                            y_slice: float = None) -> None:
    """Plot vertical cross-section showing convective structure and tracers.

    Args:
        time_idx: Time index to plot
        y_slice: Y coordinate for cross-section (default: location of max updraft)
    """
    x = data['xCell']
    y = data['yCell']
    nLevels = data['nVertLevels']
    z = np.linspace(0, 20, nLevels)  # Approximate height in km

    # Find y_slice at max updraft if not specified
    if y_slice is None and 'w' in data:
        w_horizontal = data['w'][time_idx, :, nLevels // 4]  # Mid-level
        max_idx = np.argmax(w_horizontal)
        y_slice = y[max_idx]

    if y_slice is None:
        y_slice = (y.min() + y.max()) / 2

    # Extract cross-section - use tight tolerance to get single row of cells
    y_tol = 0.2  # km - tight enough for single row on typical MPAS mesh
    mask = np.abs(y - y_slice) < y_tol

    if not np.any(mask):
        print(f"No cells found near y={y_slice:.1f} km")
        return

    x_slice = x[mask]
    sort_idx = np.argsort(x_slice)
    x_slice = x_slice[sort_idx]

    X, Z = np.meshgrid(x_slice, z)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Vertical velocity
    ax = axes[0, 0]
    if 'w' in data:
        w_slice = data['w'][time_idx, mask, :][sort_idx, :]
        w_max = max(abs(w_slice.min()), abs(w_slice.max()))
        levels_w = np.linspace(-w_max, w_max, 41)
        cf = ax.contourf(X, Z, w_slice.T, levels=levels_w, cmap='RdBu_r')
        ax.set_title('Vertical velocity (m/s)')
        plt.colorbar(cf, ax=ax)
    else:
        ax.text(0.5, 0.5, 'No w data', transform=ax.transAxes, ha='center')
    ax.set_ylabel('Height (km)')

    # Zonal wind
    ax = axes[0, 1]
    if 'uZonal' in data:
        u_slice = data['uZonal'][time_idx, mask, :][sort_idx, :]
        cf = ax.contourf(X, Z, u_slice.T, levels=50, cmap='RdBu_r')
        ax.set_title('Zonal wind (m/s)')
        plt.colorbar(cf, ax=ax)
    else:
        ax.text(0.5, 0.5, 'No wind data', transform=ax.transAxes, ha='center')

    # qAB
    ax = axes[1, 0]
    qAB_slice = data['qAB'][time_idx, mask, :][sort_idx, :]
    cf = ax.contourf(X, Z, qAB_slice.T, levels=50, cmap='viridis')
    ax.set_title('qAB (kg/kg)')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Height (km)')
    plt.colorbar(cf, ax=ax)

    # qA
    ax = axes[1, 1]
    qA_slice = data['qA'][time_idx, mask, :][sort_idx, :]
    cf = ax.contourf(X, Z, qA_slice.T, levels=50, cmap='viridis')
    ax.set_title('qA (kg/kg)')
    ax.set_xlabel('X (km)')
    plt.colorbar(cf, ax=ax)

    time_min = time_idx * 0.5  # Assuming 30s output interval
    plt.suptitle(f'Vertical Cross-Section at y={y_slice:.1f} km (t={time_min:.0f} min)', fontsize=14)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_time_series_combined(data: dict, level: int, output_file: str) -> None:
    """Generate combined time series plot with all tracers on one axis."""
    tracers = ['qAB', 'qA', 'qB']
    colors = ['purple', 'blue', 'red']
    times = data['times']

    fig, ax = plt.subplots(figsize=(10, 6))

    for tracer, color in zip(tracers, colors):
        values = data[tracer][:, :, level]  # (nTimes, nCells)
        mean_val = values.mean(axis=1)
        min_val = values.min(axis=1)
        max_val = values.max(axis=1)

        ax.fill_between(times, min_val, max_val, alpha=0.2, color=color)
        ax.plot(times, mean_val, color=color, linestyle='-', linewidth=2, label=tracer)

    ax.set_xlabel('Time step')
    ax.set_ylabel('Mixing ratio (kg/kg)')
    ax.legend()
    ax.set_title(f'Tracer Time Evolution ({level_label(data, level)})')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_spatial_time_evolution(data: dict, tracer: str, level: int,
                                 output_file: str, n_times: int = 6,
                                 show_wind: bool = False, wind_skip: int = 50,
                                 wind_scale: float = 300) -> None:
    """Plot spatial maps at multiple times to show pattern evolution."""
    x = data['xCell']
    y = data['yCell']
    tri = Triangulation(x, y)

    nTimes = data['nTimes']
    # Select evenly spaced time indices
    if nTimes <= n_times:
        time_indices = list(range(nTimes))
    else:
        time_indices = [int(i * (nTimes - 1) / (n_times - 1)) for i in range(n_times)]

    n_cols = min(3, len(time_indices))
    n_rows = (len(time_indices) + n_cols - 1) // n_cols

    # Global color scale
    vmin = data[tracer].min()
    vmax = data[tracer].max()

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    for idx, t in enumerate(time_indices):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        values = data[tracer][t, :, level]
        cf = ax.tricontourf(tri, values, levels=50, cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_aspect('equal')
        ax.set_title(f't = {t} ({level_label(data, level)})')
        plt.colorbar(cf, ax=ax, label=f'{tracer}')

        # Add wind vectors if requested
        if show_wind and 'uZonal' in data:
            u = data['uZonal'][t, :, level]
            v = data['uMeridional'][t, :, level]
            wind_idx = np.arange(0, len(x), wind_skip)
            ax.quiver(x[wind_idx], y[wind_idx], u[wind_idx], v[wind_idx],
                      scale=wind_scale, alpha=0.7, width=0.003)

    # Hide unused subplots
    for idx in range(len(time_indices), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle(f'{tracer} Time Evolution', fontsize=14)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_multispecies_time_evolution(data: dict, tracers: list, level: int,
                                      output_file: str, n_times: int = 6,
                                      show_wind: bool = False, wind_skip: int = 50,
                                      wind_scale: float = 300) -> None:
    """Plot spatial maps for multiple tracers at multiple times.

    Creates a grid with rows = tracers, columns = time steps.
    """
    x = data['xCell']
    y = data['yCell']
    tri = Triangulation(x, y)

    nTimes = data['nTimes']
    # Select evenly spaced time indices
    if nTimes <= n_times:
        time_indices = list(range(nTimes))
    else:
        time_indices = [int(i * (nTimes - 1) / (n_times - 1)) for i in range(n_times)]

    n_cols = len(time_indices)
    n_rows = len(tracers)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for row_idx, tracer in enumerate(tracers):
        # Global color scale per tracer
        vmin = data[tracer].min()
        vmax = data[tracer].max()

        for col_idx, t in enumerate(time_indices):
            ax = axes[row_idx, col_idx]

            values = data[tracer][t, :, level]
            cf = ax.tricontourf(tri, values, levels=50, cmap='viridis', vmin=vmin, vmax=vmax)
            ax.set_aspect('equal')

            # Add wind vectors if requested
            if show_wind and 'uZonal' in data:
                u = data['uZonal'][t, :, level]
                v = data['uMeridional'][t, :, level]
                wind_idx = np.arange(0, len(x), wind_skip)
                ax.quiver(x[wind_idx], y[wind_idx], u[wind_idx], v[wind_idx],
                          scale=wind_scale, alpha=0.7, width=0.003)

            # Labels
            if row_idx == n_rows - 1:
                ax.set_xlabel('X (km)')
            if col_idx == 0:
                ax.set_ylabel(f'{tracer}\nY (km)')
            if row_idx == 0:
                ax.set_title(f't = {t}')

            # Colorbar on rightmost column
            if col_idx == n_cols - 1:
                plt.colorbar(cf, ax=ax, label=f'{tracer}')

    plt.suptitle(f'Multi-species Time Evolution ({level_label(data, level)})', fontsize=14)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def plot_spatial_diff_evolution(data: dict, tracer: str, level: int,
                                 output_file: str, n_times: int = 6,
                                 diff_from_initial: bool = True) -> None:
    """Plot difference maps to show advection patterns.

    If diff_from_initial=True: shows tracer(t) - tracer(0)
    If diff_from_initial=False: shows tracer(t) - tracer(t-1)
    """
    x = data['xCell']
    y = data['yCell']
    tri = Triangulation(x, y)

    nTimes = data['nTimes']
    # Select evenly spaced time indices (skip t=0 for diff)
    if nTimes <= n_times:
        time_indices = list(range(1, nTimes))
    else:
        time_indices = [int(i * (nTimes - 1) / (n_times - 1)) for i in range(1, n_times)]

    n_cols = min(3, len(time_indices))
    n_rows = (len(time_indices) + n_cols - 1) // n_cols

    # Compute all diffs first to get global color scale
    diffs = []
    for t in time_indices:
        if diff_from_initial:
            diff = data[tracer][t, :, level] - data[tracer][0, :, level]
        else:
            diff = data[tracer][t, :, level] - data[tracer][t-1, :, level]
        diffs.append(diff)

    all_diffs = np.array(diffs)
    vmax = np.abs(all_diffs).max()
    vmin = -vmax  # Symmetric colorscale

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    for idx, (t, diff) in enumerate(zip(time_indices, diffs)):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        cf = ax.tricontourf(tri, diff, levels=50, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_aspect('equal')
        if diff_from_initial:
            ax.set_title(f't={t} - t=0')
        else:
            ax.set_title(f't={t} - t={t-1}')
        plt.colorbar(cf, ax=ax, label=f'Δ{tracer}')

    # Hide unused subplots
    for idx in range(len(time_indices), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    diff_type = "from initial" if diff_from_initial else "consecutive"
    plt.suptitle(f'{tracer} Difference ({diff_type}) - {level_label(data, level)}', fontsize=14)
    plt.tight_layout()
    save_figure(output_file)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize MPAS chemistry output')
    parser.add_argument('-i', '--input', default='output.nc',
                        help='Coupled state output file')
    parser.add_argument('-r', '--reference', default=None,
                        help='Reference state output file (optional)')
    parser.add_argument('-l', '--level', type=int, default=10,
                        help='Vertical level for horizontal slices (default: 10)')
    parser.add_argument('-t', '--time', type=int, default=-1,
                        help='Time index (default: -1 = last)')
    parser.add_argument('-o', '--output', default='chemistry_plots.png',
                        help='Output figure filename')
    parser.add_argument('--show', action='store_true',
                        help='Display plot interactively')
    parser.add_argument('--time-series', action='store_true',
                        help='Generate spatial time evolution plot for qAB')
    parser.add_argument('--temporal', action='store_true',
                        help='Generate standalone time vs mean-value plots')
    parser.add_argument('--multi-species', action='store_true',
                        help='Generate multi-species time evolution (qA and qAB)')
    parser.add_argument('--diff', action='store_true',
                        help='Generate difference plots (t - t0) to show advection')
    parser.add_argument('--diff-consecutive', action='store_true',
                        help='Generate consecutive difference plots (t - t-1)')
    parser.add_argument('--n-times', type=int, default=6,
                        help='Number of time steps for time series plot (default: 6)')
    parser.add_argument('--wind', action='store_true',
                        help='Overlay wind vectors on horizontal plots')
    parser.add_argument('--wind-skip', type=int, default=50,
                        help='Plot every Nth cell for wind vectors (default: 50)')
    parser.add_argument('--wind-scale', type=float, default=300,
                        help='Wind vector scale (larger = shorter arrows, default: 300)')
    parser.add_argument('--time-start', type=int, default=None,
                        help='Start time index for temporal plots')
    parser.add_argument('--time-end', type=int, default=None,
                        help='End time index for temporal plots')
    parser.add_argument('--single-cell', action='store_true',
                        help='Generate single-cell time series at levels 10 and 20')
    parser.add_argument('--cell-idx', type=int, default=None,
                        help='Cell index for single-cell plot (default: domain center)')
    parser.add_argument('--vertical', action='store_true',
                        help='Generate vertical cross-section through updraft')
    parser.add_argument('--y-slice', type=float, default=None,
                        help='Y coordinate for vertical cross-section (default: max updraft)')
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    data = load_mpas_data(args.input)

    time_idx = args.time if args.time >= 0 else data['nTimes'] - 1
    level = args.level

    print(f"Data: {data['nCells']} cells, {data['nVertLevels']} levels, {data['nTimes']} times")
    print(f"Plotting time={time_idx}, level={level}")

    # Generate spatial time evolution if requested
    if args.time_series:
        ts_output = args.output.replace('.png', '_timeseries.png')
        print(f"Generating spatial time series...")
        plot_spatial_time_evolution(data, 'qAB', level, ts_output, args.n_times,
                                    show_wind=args.wind, wind_skip=args.wind_skip,
                                    wind_scale=args.wind_scale)

    # Generate standalone temporal plots if requested
    if args.temporal:
        # Determine time range (default to first 3 min = 6 steps at 30s intervals)
        t_start = args.time_start if args.time_start is not None else 0
        t_end = args.time_end if args.time_end is not None else min(6, data['nTimes'] - 1)
        time_range = (t_start, t_end)

        temporal_output = args.output.replace('.png', '_temporal.png')
        print(f"Generating temporal evolution plots (t={t_start}-{t_end})...")
        plot_time_series_standalone(data, level, temporal_output, time_range=time_range)

    # Generate single-cell time series if requested
    if args.single_cell:
        # Use time range settings
        t_start = args.time_start if args.time_start is not None else 0
        t_end = args.time_end if args.time_end is not None else min(6, data['nTimes'] - 1)
        time_range = (t_start, t_end)

        sc_output = args.output.replace('.png', '_single_cell.png')
        print(f"Generating single-cell time series (levels 10, 20)...")
        plot_single_cell_timeseries(data, [10, 20], sc_output,
                                     cell_idx=args.cell_idx, time_range=time_range)

    # Generate vertical cross-section if requested
    if args.vertical:
        vert_output = args.output.replace('.png', '_vertical.png')
        print(f"Generating vertical cross-section...")
        plot_vertical_structure(data, time_idx, vert_output, y_slice=args.y_slice)

    # Generate multi-species time evolution if requested
    if args.multi_species:
        ms_output = args.output.replace('.png', '_multispecies.png')
        print(f"Generating multi-species time series (qA, qAB)...")
        plot_multispecies_time_evolution(data, ['qA', 'qAB'], level, ms_output, args.n_times,
                                          show_wind=args.wind, wind_skip=args.wind_skip,
                                          wind_scale=args.wind_scale)

    # Generate difference plots if requested
    if args.diff:
        diff_output = args.output.replace('.png', '_diff.png')
        print(f"Generating difference from initial plots...")
        plot_spatial_diff_evolution(data, 'qAB', level, diff_output, args.n_times, diff_from_initial=True)

    if args.diff_consecutive:
        diff_output = args.output.replace('.png', '_diff_consecutive.png')
        print(f"Generating consecutive difference plots...")
        plot_spatial_diff_evolution(data, 'qAB', level, diff_output, args.n_times, diff_from_initial=False)

    # Determine global min/max for consistent color scales
    vmin_ab = min(data['qAB'].min(), 0)
    vmax_ab = data['qAB'].max()
    vmin_a = min(data['qA'].min(), 0)
    vmax_a = max(data['qA'].max(), data['qB'].max())

    # Create figure
    fig = plt.figure(figsize=(16, 12))

    # Row 1: Horizontal slices of AB, A, B
    ax1 = fig.add_subplot(3, 3, 1)
    plot_horizontal_slice(data, 'qAB', level, time_idx, ax1,
                          f'AB (level={level}, t={time_idx})', vmin=vmin_ab, vmax=vmax_ab,
                          show_wind=args.wind, wind_skip=args.wind_skip, wind_scale=args.wind_scale)

    ax2 = fig.add_subplot(3, 3, 2)
    plot_horizontal_slice(data, 'qA', level, time_idx, ax2,
                          f'A (level={level}, t={time_idx})', vmin=vmin_a, vmax=vmax_a,
                          show_wind=args.wind, wind_skip=args.wind_skip, wind_scale=args.wind_scale)

    ax3 = fig.add_subplot(3, 3, 3)
    plot_horizontal_slice(data, 'qB', level, time_idx, ax3,
                          f'B (level={level}, t={time_idx})', vmin=vmin_a, vmax=vmax_a,
                          show_wind=args.wind, wind_skip=args.wind_skip, wind_scale=args.wind_scale)

    # Row 2: Vertical cross-sections
    y_center = (data['yCell'].min() + data['yCell'].max()) / 2

    ax4 = fig.add_subplot(3, 3, 4)
    plot_vertical_cross_section(data, 'qAB', time_idx, y_center, ax4,
                                 f'AB vertical (y={y_center:.0f} km)')

    ax5 = fig.add_subplot(3, 3, 5)
    plot_vertical_cross_section(data, 'qA', time_idx, y_center, ax5,
                                 f'A vertical (y={y_center:.0f} km)')

    ax6 = fig.add_subplot(3, 3, 6)
    plot_vertical_cross_section(data, 'qB', time_idx, y_center, ax6,
                                 f'B vertical (y={y_center:.0f} km)')

    # Row 3: Time evolution
    ax7 = fig.add_subplot(3, 3, 7)
    plot_time_evolution(data, 'qAB', level, ax7)

    ax8 = fig.add_subplot(3, 3, 8)
    plot_time_evolution(data, 'qA', level, ax8)

    ax9 = fig.add_subplot(3, 3, 9)
    plot_time_evolution(data, 'qB', level, ax9)

    plt.suptitle(f'MPAS Chemistry Output: {Path(args.input).name}', fontsize=14)
    plt.tight_layout()

    # Save
    save_figure(args.output)

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
