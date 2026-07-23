#!/usr/bin/env bash
# Run the Jablonowski-Williamson baroclinic wave test case on the 480-km mesh.
#
# Prerequisites:
#   - init_atmosphere_model and atmosphere_model built
#   - 480-km mesh downloaded (scripts/download_data.sh)
#
# Usage: ./quacs/run_jw_test.sh [NPROCS]
#   NPROCS: number of MPI ranks (default: 1)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MPAS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NPROCS="${1:-1}"

# Verify executables exist
for exe in init_atmosphere_model atmosphere_model; do
    if [ ! -x "${MPAS_DIR}/${exe}" ]; then
        echo "ERROR: ${MPAS_DIR}/${exe} not found. Build MPAS first." >&2
        exit 1
    fi
done

# Verify mesh exists
MESH_DIR="${MPAS_DIR}/data"
GRID_FILE="${MESH_DIR}/x1.2562.grid.nc"
if [ ! -f "${GRID_FILE}" ]; then
    echo "ERROR: ${GRID_FILE} not found. Run scripts/download_data.sh first." >&2
    exit 1
fi

# Create working directory
WORK_DIR="${MPAS_DIR}/data/jw_480km"
rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

echo "=== Setting up JW baroclinic wave test (480-km, ${NPROCS} MPI ranks) ==="

# Link executables and mesh
ln -sf "${MPAS_DIR}/init_atmosphere_model" .
ln -sf "${MPAS_DIR}/atmosphere_model" .
ln -sf "${GRID_FILE}" .

# Link graph partition files for this mesh
for f in "${MESH_DIR}"/x1.2562.graph.info.part.*; do
    [ -f "$f" ] && ln -sf "$f" .
done

# ---- Create init namelist ----
cat > namelist.init_atmosphere << 'EOF'
&nhyd_model
    config_start_time = '0000-01-01_00:00:00'
    config_init_case = 2
/

&dimensions
    config_nvertlevels = 26
/

&decomposition
    config_block_decomp_file_prefix = 'x1.2562.graph.info.part.'
/
EOF

# ---- Create init streams ----
cat > streams.init_atmosphere << 'EOF'
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="x1.2562.grid.nc"
                  input_interval="initial_only" />

<immutable_stream name="output"
                  type="output"
                  filename_template="x1.2562.init.nc"
                  packages="initial_conds"
                  output_interval="initial_only" />

<immutable_stream name="surface"
                  type="output"
                  filename_template="sfc_not_needed_for_jw"
                  filename_interval="none"
                  packages="sfc_update"
                  output_interval="86400" />

<immutable_stream name="lbc"
                  type="output"
                  filename_template="lbc_not_needed_for_jw"
                  filename_interval="output_interval"
                  packages="lbcs"
                  output_interval="none" />

</streams>
EOF

# ---- Run init_atmosphere_model ----
echo "=== Creating initial conditions ==="
MPIRUN_OPTS="--oversubscribe -np ${NPROCS}"
# Containers typically run as root; OpenMPI requires an explicit flag for this
if [ "$(id -u)" -eq 0 ]; then
    MPIRUN_OPTS="${MPIRUN_OPTS} --allow-run-as-root"
fi
mpirun ${MPIRUN_OPTS} ./init_atmosphere_model
echo ""

# Check for errors
if grep -q "^    Critical error messages = *[1-9]" log.init_atmosphere.0000.out 2>/dev/null; then
    echo "ERROR: init_atmosphere_model encountered critical errors" >&2
    cat log.init_atmosphere.0000.out
    exit 1
fi

if [ ! -f x1.2562.init.nc ]; then
    echo "ERROR: init_atmosphere_model did not produce x1.2562.init.nc" >&2
    exit 1
fi
echo "  init_atmosphere_model succeeded"

# ---- Create atmosphere namelist ----
# 480-km mesh: dt ~ 1800s (scale from 120-km=450s by ratio 480/120=4)
cat > namelist.atmosphere << 'EOF'
&nhyd_model
    config_dt = 1800.0
    config_start_time = '0000-01-01_00:00:00'
    config_run_duration = '1_00:00:00'
    config_split_dynamics_transport = false
    config_number_of_sub_steps = 6
    config_dynamics_split_steps = 1
    config_h_mom_eddy_visc2    = 0.0
    config_h_mom_eddy_visc4    = 0.0
    config_v_mom_eddy_visc2    = 0.0
    config_h_theta_eddy_visc2  = 0.0
    config_h_theta_eddy_visc4  = 0.0
    config_v_theta_eddy_visc2  = 0.0
    config_horiz_mixing        = '2d_smagorinsky'
    config_len_disp            = 480000.
    config_u_vadv_order      = 3
    config_w_vadv_order      = 3
    config_theta_vadv_order  = 3
    config_scalar_vadv_order = 3
    config_theta_adv_order   = 3
    config_scalar_adv_order  = 3
    config_scalar_advection  = false
    config_positive_definite = false
    config_coef_3rd_order = 1.0
    config_monotonic = false
    config_epssm = 0.1
    config_smdiv = 0.1
/

&damping
    config_zd = 22000.0
    config_xnutr = 0.0
/

&limited_area
    config_apply_lbcs = false
/

&io
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/

&decomposition
    config_block_decomp_file_prefix = 'x1.2562.graph.info.part.'
/

&restart
    config_do_restart = false
/

&printout
    config_print_global_minmax_vel = true
    config_print_global_minmax_sca = false
/

&physics
    config_o3climatology = false
    config_sst_update = false
    config_sstdiurn_update = false
    config_deepsoiltemp_update = false
    config_radtlw_interval = '00:30:00'
    config_radtsw_interval = '00:30:00'
    config_bucket_update = 'none'
    config_physics_suite = 'none'
/
EOF

# ---- Create atmosphere streams ----
cat > streams.atmosphere << 'EOF'
<streams>

<immutable_stream name="input"
                  type="input"
                  filename_template="x1.2562.init.nc"
                  input_interval="initial_only"/>

<immutable_stream name="restart"
                  type="input;output"
                  filename_template="restart.$Y-$M-$D_$h.$m.$s.nc"
                  input_interval="initial_only"
                  output_interval="none"/>

<stream name="output"
        type="output"
        filename_template="output.nc"
        filename_interval="none"
        output_interval="01:00:00">

    <var name="latCell"/>
    <var name="lonCell"/>
    <var name="surface_pressure"/>
    <var name="pressure_base"/>
    <var name="pressure_p"/>
    <var name="theta"/>
    <var name="uReconstructZonal"/>
    <var name="uReconstructMeridional"/>

</stream>

<stream name="surface"
        type="input"
        filename_template="x1.2562.sfc_update.nc"
        filename_interval="none"
        input_interval="none">
</stream>

<immutable_stream name="iau"
                  type="input"
                  filename_template="AmB.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="none"
                  packages="iau"
                  input_interval="initial_only" />

<immutable_stream name="lbc_in"
                  type="input"
                  filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="input_interval"
                  packages="limited_area"
                  input_interval="none" />

</streams>
EOF

# ---- Create empty stream_list file ----
touch stream_list.atmosphere.output
touch stream_list.atmosphere.diagnostics
touch stream_list.atmosphere.diag_ugwp

# ---- Run atmosphere_model ----
echo "=== Running JW baroclinic wave (1 day, 480-km mesh) ==="
mpirun ${MPIRUN_OPTS} ./atmosphere_model
echo ""

# Check for errors
if grep -q "^    Critical error messages = *[1-9]" log.atmosphere.0000.out 2>/dev/null; then
    echo "ERROR: atmosphere_model encountered critical errors" >&2
    tail -20 log.atmosphere.0000.out
    exit 1
fi

if [ ! -f output.nc ]; then
    echo "ERROR: atmosphere_model did not produce output.nc" >&2
    exit 1
fi

echo "  atmosphere_model succeeded"
echo ""
echo "=== JW baroclinic wave test PASSED ==="
echo "  Output: ${WORK_DIR}/output.nc"
ls -lh output.nc
