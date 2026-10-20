"""Deduplicate the stable, kappa > 100 W/mK structures of the initial screen (paper Table 1 settings):
pymatgen StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, attempt_supercell=True),
compared within groups of equal reduced formula."""
import csv, glob
from collections import defaultdict
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
import os as _os
REPO = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..")
CIFS = _os.path.join(REPO, "data", "structures_above100")
rows = [r for r in csv.DictReader(open(_os.path.join(REPO, "data", "screening_initial.csv"))) if r["above_100"] == "True"]
sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, attempt_supercell=True)
groups = defaultdict(list)
for r in rows:
    s = Structure.from_file(glob.glob(f"{CIFS}/{r['lemat_id']}_*.cif")[0])
    groups[s.composition.reduced_formula].append((r, s))
polys = []
for rf, items in sorted(groups.items()):
    reps = []
    for r, s in items:
        for rep in reps:
            if sm.fit(s, rep[1]): rep[2].append(r); break
        else: reps.append([r, s, [r]])
    polys += [(rf, rep[0]["space_group"], [m["lemat_id"] + ":" + m["source_id_lemat"] for m in rep[2]]) for rep in reps]
print(f"{len(rows)} structures -> {len(polys)} distinct polymorphs")
for rf, sg, mem in polys: print(f"  {rf:7} {sg:16} {len(mem)}x  {', '.join(mem)}")

with open(_os.path.join(REPO, "data", "polymorph_groups.csv"), "w") as fh:
    fh.write("reduced_formula,space_group,n_copies,members_lemat_id:source_id\n")
    for rf, sg, mem in polys: fh.write(f'{rf},"{sg}",{len(mem)},"{" ".join(mem)}"\n')
