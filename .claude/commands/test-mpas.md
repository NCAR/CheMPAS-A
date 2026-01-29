# Test MPAS Run

Run an MPAS test case with automatic output archival.

## Arguments

- `$ARGUMENTS` - Optional: test case directory (default: `~/MPAS/supercell`)

## Instructions

Execute the following steps:

### 1. Set Up

Determine the test directory and MPI configuration:

```bash
TEST_DIR="${ARGUMENTS:-$HOME/MPAS/supercell}"
cd "$TEST_DIR"
```

Verify required files exist:
- `namelist.atmosphere`
- `streams.atmosphere`
- Initial conditions file (check `streams.atmosphere` for `input` stream filename)
- Graph partition files (`*.graph.info.part.*`)

Read `namelist.atmosphere` to find the decomposition prefix and determine available partition sizes.

**Important:** Check that `streams.atmosphere` has `io_type="netcdf"` on input and output streams. PnetCDF has compatibility issues on macOS/LLVM builds. If missing, add `io_type="netcdf"` to avoid "Could not open input file" or "PIO error -36" errors.

### 2. Archive Previous Output

Before running, archive any existing output with a timestamp:

```bash
timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
for f in log.atmosphere.*.out; do
    [ -f "$f" ] && mv "$f" "${f%.out}.${timestamp}.out"
done
```

### 3. Run the Model

Run with 8 MPI ranks (optimal for 10-core Mac, avoids oversubscription):

```bash
mpiexec -n 8 ~/EarthSystem/CheMPAS/atmosphere_model 2>&1 | tee run.out
```

**Note:** If user requests different rank count, verify the corresponding partition file exists.

### 4. Verify Results

After completion, check for success:

1. Check `log.atmosphere.0000.out` exists and has content
2. Look for timing statistics near the end of the log
3. Verify `output.nc` was created
4. Report any `CRITICAL ERROR` messages

```bash
# Quick success check
tail -30 log.atmosphere.0000.out
ls -lh output.nc
```

### 5. Report Summary

Provide a summary including:
- Test case run
- Run duration (from log)
- Output file size
- Any warnings or errors encountered

## Available Test Cases

| Directory | Description |
|-----------|-------------|
| `~/MPAS/supercell` | Idealized supercell thunderstorm (2 min run) |

## MPI Rank Recommendations

| Machine Cores | Recommended Ranks |
|---------------|-------------------|
| 8 | 4 or 8 |
| 10 | 8 |
| 12+ | 8 or 12 |

Partition files must exist: `*.graph.info.part.N` where N = rank count.
