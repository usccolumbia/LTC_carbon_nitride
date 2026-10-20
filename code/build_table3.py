"""Build data/table3_dft.csv (paper Table 3): first-principles kappa at 300 K for the 11 validated structures.

Source: ShengBTE iterative solutions (VASP PBEsol + D3, Pheasy force constants; settings in the paper, Table 1),
from the collaborators' result files. Tensors in W/mK, Cartesian frame of the DFT cell.
Off-diagonal components: computed values for C3N4 I-42m and C6N8 Cmc2_1; for C12N4 C2/c and C6N2 C2/m the
values reported in the paper (Table 3, note a, one decimal); all others vanish by symmetry.
"""
import csv, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data")
# lemat_id: (formula, SG, xx, yy, zz, xy, xz, yz, gamma_bar, IFC cutoff (3rd/4th), supercell, BTE q-grid, off-diagonal source)
DFT = {
 "006232245": ("C12N4", "C2/c",    801.18, 831.99, 655.09, -18.8, 60.2, 28.5, 0.465584548, "3rd NN8 / 4th NN4", "3x3x3", "10x10x10", "paper note a (1 decimal)"),
 "002090448": ("C2N4",  "P-3m1",   882.13, 882.13, 3.698,  0, 0, 0, 0.734237995, "3rd 6 A / 4th 4 A", "6x6x2", "26x26x6", "symmetry"),
 "004092476": ("C11N4", "P-42m",   490.35, 490.35, 403.45, 0, 0, 0, 0.868015042, "3rd NN8 / 4th NN4", "3x3x2", "17x17x9", "symmetry"),
 "000486093": ("C2N4",  "P-42_1m", 599.47, 599.47, 11.638, 0, 0, 0, 0.65229145,  "3rd NN8 / 4th NN4", "4x4x3", "18x18x15", "symmetry"),
 "006423936": ("C6N8",  "Cmc2_1",  296.90, 583.01, 279.96, 81.349, 0, 0, None,   "3rd NN8 / 4th NN4", "", "", "computed"),
 "006305877": ("C6N2",  "C2/m",    248.91, 537.61, 284.03, -121.3, 105.7, 38.5, 0.440643424, "3rd NN8 / 4th NN4", "4x4x2", "17x17x9", "paper note a (1 decimal)"),
 "002016386": ("CN2",   "I-4m2",   427.17, 427.17, 206.48, 0, 0, 0, 0.896231441314249, "3rd NN8 / 4th NN4", "4x4x4", "28x28x28", "symmetry"),
 "002289672": ("C4N4",  "P4_332",  339.87, 339.87, 339.87, 0, 0, 0, -0.082291159947727, "3rd NN8 / 4th NN4", "3x3x3", "20x20x20", "symmetry"),
 "006455880": ("C3N4",  "I-42m",   329.44, 319.45, 309.34, 14.261, 0, 0, None,   "3rd NN8 / 4th NN4", "", "", "computed"),
 "003509749": ("C2N4",  "Cmc2_1",  428.53, 298.11, 186.48, 0, 0, 0, 0.786094914, "3rd NN8 / 4th NN4", "5x4x3", "24x18x13", "symmetry"),
 "006501365": ("C6N8",  "Pmn2_1",  259.38, 264.78, 371.58, 0, 0, 0, 0.783591308, "3rd NN8 / 4th NN4", "3x3x3", "12x14x15", "symmetry"),
}
mlip = {r["lemat_id"]: r for r in csv.DictReader(open(os.path.join(D, "table2_mlip.csv")))}
rows = []
for lid, (f, sg, xx, yy, zz, xy, xz, yz, gam, ifc, sc, qg, src) in DFT.items():
    K = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]); e = np.sort(np.linalg.eigvalsh(K))[::-1]; km = np.trace(K) / 3
    m = mlip[lid]
    rows.append(dict(lemat_id=lid, source_id=m["source_id"], formula=f, space_group=sg,
        kappa_xx=xx, kappa_yy=yy, kappa_zz=zz, kappa_yz=yz, kappa_xz=xz, kappa_xy=xy, offdiag_source=src,
        kappa_mean=round(km, 2), kappa_1=round(e[0], 2), kappa_2=round(e[1], 2), kappa_3=round(e[2], 2), rho=round(e[0] / e[2], 2),
        gamma_bar=gam if gam is not None else "", kappa_mlip=m["kappa_mean"], ratio_dft_mlip=round(km / float(m["kappa_mean"]), 3),
        ifc_cutoff=ifc, supercell=sc, bte_qgrid=qg, temperature_K=300))
rows.sort(key=lambda r: -r["kappa_mean"])
with open(os.path.join(D, "table3_dft.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"table3_dft.csv: {len(rows)} rows")
