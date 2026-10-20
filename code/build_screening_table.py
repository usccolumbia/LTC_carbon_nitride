"""Build data/screening_initial.csv (all 270 LeMaT C-N structures, initial screen).
Provenance script: reads the raw cluster outputs (first-pass dump + per-structure JSONs of re-run jobs)."""
import csv, json, glob, os, datetime as dt
DUMP = "/work/rchevva/lemat_C_N_screening/all_271_dump.csv"
RETRY = "/work/rchevva/lemat_C_N_screening/results/tier_retry"
import os as _os
REPO = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")
OUT = _os.path.join(REPO, "data", "screening_initial.csv")
# DUMP and RETRY are the raw first-pass outputs on the authors' cluster (not distributed); the CSV is the data product.
MESH_EDIT = dt.datetime(2026, 6, 19, 12, 22)          # worker MESH 11 -> 21
TIER_MODE = {"1": "baseline", "2": "baseline", "3": "baseline", "4": "cutoff_4_fc2_333", "5": "cutoff_4_fc2_333"}
rows = {}
for r in csv.DictReader(open(DUMP)):
    lid = r["lemat_id"].zfill(9)
    rows[lid] = dict(lemat_id=lid, source_id_lemat=r["source_id"], formula=r["formula"],
                     space_group=r["space_group"], n_atoms=r["n_atoms"], mode=TIER_MODE.get(r["tier"], ""),
                     q_mesh=11, kappa_mean=r["kappa_W_mK"], kappa_xx=r["kappa_xx"], kappa_yy=r["kappa_yy"],
                     kappa_zz=r["kappa_zz"], dynamically_stable=r["stable"], n_imaginary=r["n_imag"],
                     min_freq_THz=r["min_freq_THz"], status=r["status_raw"])
for f in glob.glob(f"{RETRY}/*_result.json"):
    lid = os.path.basename(f)[:9]
    try: d = json.load(open(f))
    except Exception: continue                         # unreadable/cancelled job: keep first-pass status
    start = dt.datetime.fromtimestamp(os.path.getmtime(f) - (d.get("runtime_sec") or 0))
    r = rows[lid]
    r.update(formula=d.get("formula") or r["formula"], space_group=d.get("space_group") or r["space_group"],
             n_atoms=d.get("n_atoms") or r["n_atoms"], mode=d.get("mode", r["mode"]), q_mesh=21 if start > MESH_EDIT else 11,
             status=d.get("status"), dynamically_stable=d.get("dynamically_stable"),
             n_imaginary=d.get("n_imaginary_modes"), min_freq_THz=d.get("min_frequency_THz", ""),
             kappa_mean=d.get("kappa_raw"), kappa_xx=d.get("kappa_xx"), kappa_yy=d.get("kappa_yy"),
             kappa_zz=d.get("kappa_zz"))
for r in rows.values():
    try: k = float(r["kappa_mean"])
    except (TypeError, ValueError): k = None
    st = str(r["dynamically_stable"]) == "True"
    r["kappa_physical"] = bool(k is not None and k < 1e3)
    r["stable_physical"] = st and r["kappa_physical"]
    r["above_100"] = r["stable_physical"] and k > 100
cols = ["lemat_id","source_id_lemat","formula","space_group","n_atoms","mode","q_mesh","status","dynamically_stable",
        "n_imaginary","min_freq_THz","kappa_mean","kappa_xx","kappa_yy","kappa_zz","kappa_physical","stable_physical","above_100"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for lid in sorted(rows): w.writerow({c: rows[lid].get(c, "") for c in cols})
R = list(rows.values())
print(f"rows {len(R)} | stable {sum(str(r['dynamically_stable'])=='True' for r in R)} | "
      f"stable+physical {sum(r['stable_physical'] for r in R)} | above 100 {sum(r['above_100'] for r in R)}")
print("q_mesh 21 rows:", [(r['lemat_id'], r['formula'], r['status']) for r in R if r['q_mesh'] == 21])

# Rows with no output file: formula, space group and size from the input CIF; no settings recorded.
import csv
import glob as _glob
from pymatgen.core import Structure as _Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as _SGA
_CIFS = "/work/rchevva/lemat_C_N/lemat_bulk_cifs_elem_C_N"   # 270 LeMaT-Bulk CIFs (not distributed)
_rows = list(csv.DictReader(open(OUT))); _fields = list(_rows[0])
for _r in _rows:
    if _r["status"] in ("EMPTY_FILE", "no_output"):
        _s = _Structure.from_file(_glob.glob(f"{_CIFS}/*{_r['lemat_id']}*.cif")[0]); _g = _SGA(_s, symprec=0.01)
        _r.update(formula=_s.composition.formula.replace(" ", ""), n_atoms=len(_s),
                  space_group=f"{_g.get_space_group_symbol()} ({_g.get_space_group_number()})",
                  mode="", q_mesh="", status="no_output")
with open(OUT, "w", newline="") as _fh:
    _w = csv.DictWriter(_fh, fieldnames=_fields); _w.writeheader(); _w.writerows(_rows)
