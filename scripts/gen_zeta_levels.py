#!/usr/bin/env python3
"""Generate a specified-zeta-levels text file for MPAS idealized init.

Writes nVertLevels+1 edge heights (metres) as a simple stretched power-law
profile: zw(k) = top * ((k-1)/nVertLevels)**stretch for k = 1..nVertLevels+1.

The file format matches what ``read_text_array`` consumes in
``mpas_init_atm_cases.F`` (one value per line, metres).
"""
import argparse


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=float, required=True, help="domain top (m)")
    ap.add_argument("--nlevels", type=int, required=True, help="nVertLevels")
    ap.add_argument("--stretch", type=float, default=1.0,
                    help="power-law exponent (1.0 = uniform)")
    ap.add_argument("--out", required=True, help="output file path")
    args = ap.parse_args()

    nz1 = args.nlevels
    nz = nz1 + 1
    edges = [args.top * (k / nz1) ** args.stretch for k in range(nz)]

    with open(args.out, "w") as fh:
        for z in edges:
            fh.write(f"{z:.6f}\n")

    dz1 = edges[1] - edges[0]
    dz_top = edges[-1] - edges[-2]
    print(f"wrote {args.out}: {nz} edges, dz_surf={dz1:.1f} m, dz_top={dz_top:.1f} m")


if __name__ == "__main__":
    main()
