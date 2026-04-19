#!/usr/bin/env python3
"""Generate the TUV-x upper-atmosphere climatology CSV.

Emits edge-height values for T, air number density, and O3 number density
spanning a requested altitude range at uniform spacing. Defaults to
50-100 km at 5 km spacing, 11 edges (10 layers).

T and air values are interpolated from US Standard Atmosphere 1976 tables.
O3 values are interpolated from AFGL mid-latitude-summer constituent profile
values sampled at the same altitudes.

The produced CSV is consumed by ``load_extension_csv`` in
``src/core_atmosphere/chemistry/mpas_tuvx.F``. Format: header line followed
by one row per edge, columns ``z_km, T_K, n_air_molec_cm3, n_O3_molec_cm3``.
"""
import argparse
import math

K_B = 1.380649e-23  # Boltzmann constant [J/K]


# US Standard Atmosphere 1976 sampled at 5-km spacing over 50-100 km.
# z [km], T [K], P [Pa].
USSA76 = [
    (50.0, 270.65, 7.9779e+01),
    (55.0, 260.77, 4.2522e+01),
    (60.0, 247.02, 2.1958e+01),
    (65.0, 233.29, 1.0930e+01),
    (70.0, 219.58, 5.2200e+00),
    (75.0, 208.40, 2.3880e+00),
    (80.0, 198.64, 1.0520e+00),
    (85.0, 188.89, 4.4570e-01),
    (90.0, 186.87, 1.8360e-01),
    (95.0, 188.42, 7.5970e-02),
    (100.0, 195.08, 3.2010e-02),
]

# AFGL mid-latitude-summer ozone number density [molec/cm^3], sampled
# at matching altitudes. Source: AFGL atmospheric constituent profiles
# (Anderson et al. 1986), interpolated/smoothed to this grid.
O3_AFGL_MLS = [
    (50.0, 6.45e+11),
    (55.0, 2.02e+11),
    (60.0, 7.42e+10),
    (65.0, 3.20e+10),
    (70.0, 9.60e+09),
    (75.0, 3.76e+09),
    (80.0, 1.49e+09),
    (85.0, 6.35e+08),
    (90.0, 2.39e+08),
    (95.0, 8.80e+07),
    (100.0, 2.85e+07),
]


def loglin_interp(z, table):
    """Linear in z for T, log-linear for a positive quantity (P, n_O3)."""
    if z <= table[0][0]:
        return table[0][1:]
    if z >= table[-1][0]:
        return table[-1][1:]
    for i in range(len(table) - 1):
        z0, *y0 = table[i]
        z1, *y1 = table[i + 1]
        if z0 <= z <= z1:
            t = (z - z0) / (z1 - z0)
            out = []
            for a, b in zip(y0, y1):
                if a > 0 and b > 0:
                    out.append(math.exp(math.log(a) + t * (math.log(b) - math.log(a))))
                else:
                    out.append(a + t * (b - a))
            return out
    raise RuntimeError("interpolation out of range")


def n_air_from_PT(P_Pa, T_K):
    """Air number density [molec/cm^3] from pressure [Pa] and T [K]."""
    return P_Pa / (K_B * T_K) * 1.0e-6  # m^-3 -> cm^-3


def sample(z_km):
    T, P = loglin_interp(z_km, USSA76)
    (n_o3,) = loglin_interp(z_km, O3_AFGL_MLS)
    return T, n_air_from_PT(P, T), n_o3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--z-bottom", type=float, default=50.0, help="bottom edge (km)")
    ap.add_argument("--z-top", type=float, default=100.0, help="top edge (km)")
    ap.add_argument("--nlayers", type=int, default=10, help="number of layers (edges = nlayers+1)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    nedges = args.nlayers + 1
    dz = (args.z_top - args.z_bottom) / args.nlayers

    with open(args.out, "w") as fh:
        fh.write("z_km,T_K,n_air_molec_cm3,n_O3_molec_cm3\n")
        for k in range(nedges):
            z = args.z_bottom + k * dz
            T, n_air, n_o3 = sample(z)
            fh.write(f"{z:.2f},{T:.3f},{n_air:.6e},{n_o3:.6e}\n")

    print(f"wrote {args.out}: {nedges} edges ({args.nlayers} layers), {args.z_bottom}-{args.z_top} km")


if __name__ == "__main__":
    main()
