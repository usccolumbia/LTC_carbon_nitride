"""Build data/table2_mlip.csv (paper Table 2) from data/mlip_results/*.json."""
import json, glob, os, csv, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data")
from pymatgen.core import Structure
def axes_aligned(stem):
    """Cell axes along x, y, z (90/90/90 or 90/90/120): off-diagonal kappa then vanishes by symmetry for SG >= 16."""
    al, be, ga = Structure.from_file(os.path.join(D, "structures_relaxed", stem + "_nequipL_relaxed.cif")).lattice.angles
    return abs(al - 90) < 1e-3 and abs(be - 90) < 1e-3 and (abs(ga - 90) < 1e-3 or abs(ga - 120) < 1e-3)
# lemat_id: (source ID as in Table 2 = current database entry ID, dHf, E_hull) ; source-database values, not recomputed
T2 = {"004092476":("agm005276040",0.304,0.304),"002016386":("mp-1009818",0.730,0.730),"006232245":("oqmd-2036226",0.354,0.394),
      "000486093":("agm005425724",0.363,0.363),"002090448":("agm005454843",0.587,0.587),"006455880":("oqmd-2032304",0.482,0.505),
      "002728099":("mp-571653",0.498,0.498),"006423936":("oqmd-2032302",0.400,0.423),"006423462":("oqmd-2030885",0.316,0.343),
      "002701325":("mp-1188347",0.265,0.265),"003509749":("mp-1077595",0.707,0.707),"006289931":("oqmd-2053626",0.404,0.427),
      "001490044":("mp-1985",0.336,0.336),"006501365":("oqmd-2032306",0.440,0.462),"006541006":("oqmd-2053624",0.459,0.482),
      "002289672":("agm004449150",0.611,0.611),"004477381":("mp-2852",0.503,0.503),"006243597":("oqmd-2053629",0.545,0.568),
      "006305877":("oqmd-2037964",0.377,0.417),"002375518":("agm005243025",0.559,0.559),"006211723":("oqmd-2053627",0.525,0.548),
      "006290867":("oqmd-1483855",0.735,0.753),"003986562":("mp-9410",0.299,0.299),"006708111":("oqmd-1472226",0.726,0.744),
      "006575614":("oqmd-2053619",0.099,0.122)}
DFT = {"004092476","002016386","006232245","000486093","002090448","006455880","006423936","003509749","006501365","002289672","006305877"}
rows = []
for f in sorted(glob.glob(os.path.join(D, "mlip_results", "T2_*.json"))):
    stem = os.path.basename(f)[:-5]; lid = stem[6:15]; d = json.load(open(f))
    v = {k: d.get(k) for k in ("kappa_xx","kappa_yy","kappa_zz","kappa_yz","kappa_xz","kappa_xy")}
    sgn = int(re.search(r"\((\d+)\)", str(d.get("space_group"))).group(1)); note = ""
    if any(v[k] is None for k in ("kappa_yz", "kappa_xz", "kappa_xy")) and sgn >= 16 and axes_aligned(stem):
        for k in ("kappa_yz", "kappa_xz", "kappa_xy"):
            v[k] = 0.0 if v[k] is None else v[k]
        note = "off-diagonal components not stored by this run; zero by symmetry (orthorhombic, axes along x, y, z)"
    full = all(v[k] is not None for k in v)
    K = np.array([[v["kappa_xx"], v["kappa_xy"] or 0, v["kappa_xz"] or 0], [v["kappa_xy"] or 0, v["kappa_yy"], v["kappa_yz"] or 0],
                  [v["kappa_xz"] or 0, v["kappa_yz"] or 0, v["kappa_zz"]]])
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    sid, dhf, eh = T2[lid]
    rows.append(dict(table2_id=stem[:5], lemat_id=lid, source_id=sid, source_id_lemat=stem.split("_")[3],
        formula=re.sub(r"(?<=[A-Z])1(?![0-9])", "", d.get("formula") or ""), space_group=d.get("space_group"), n_atoms=d.get("n_atoms"),
        settings=d.get("mode"), fc2_supercell=d.get("fc2_supercell"), fc3_supercell=d.get("fc3_supercell"),
        fc3_pair_cutoff_A=d.get("cutoff_pair"), q_mesh=d.get("mesh", [21,21,21]),
        kappa_mean=round(d["kappa_raw"], 2), **{k: (round(x, 3) if x is not None else "") for k, x in v.items()},
        full_tensor=full, tensor_note=note, kappa_1=round(ev[0], 2), kappa_2=round(ev[1], 2), kappa_3=round(ev[2], 2),
        rho=round(ev[0] / ev[2], 2), dynamically_stable=d.get("dynamically_stable"), min_freq_THz=d.get("min_frequency_THz"),
        max_force_after_relax=d.get("max_force_after_relax"), device=d.get("device", "cuda"),
        dHf_eV_atom=dhf, E_hull_eV_atom=eh, dft_validated=lid in DFT))
rows.sort(key=lambda r: -r["kappa_mean"])
with open(os.path.join(D, "table2_mlip.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['table2_id']} {r['formula']:8} {str(r['space_group']):15} {r['kappa_mean']:>8} rho={r['rho']:<7} {r['settings']:17} full={r['full_tensor']}")
