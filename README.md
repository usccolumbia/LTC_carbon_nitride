# High lattice thermal conductivity in the bulk carbon-nitride family approaching cubic boron nitride

Code and data for *High lattice thermal conductivity in the bulk carbon-nitride family approaching cubic boron nitride* (R. L. R. Chevva et al.). Citation details will be added on publication.

## Pipeline

1. **Input.** All 270 {C, N} bulk structures in LeMaT-Bulk (Materials Project, OQMD and Alexandria entries).
2. **MLIP screen.** Each structure is relaxed with NequIP-OAM-L (BFGS + ExpCellFilter, forces below 1e-4 eV/Å, at most 300 steps). Second- and third-order force constants come from finite displacements of 0.03 Å, and κ from phono3py in the relaxation-time approximation at 300 K, with three-phonon scattering only. A structure counts as dynamically stable if no mode on the q-mesh lies below −0.05 THz.
3. **Selection.** 270 → 264 completed → 92 dynamically stable → 88 with physical κ (four stable structures gave κ > 10³ W m⁻¹ K⁻¹ and were discarded) → 43 with κ̄ > 100 W m⁻¹ K⁻¹ → 24 distinct polymorphs. Deduplication uses pymatgen StructureMatcher (ltol 0.2, stol 0.3, angle_tol 5°, primitive_cell, attempt_supercell), grouped by reduced formula. C₂₄N₃₂ (Pca2₁, 56 atoms) is added as a large-cell reference, giving 25.
4. **Table 2.** The 25 structures are recomputed at uniform settings: FC2 4×4×4, FC3 2×2×2, no pair cutoff, q-mesh 21×21×21. The exception is C₂₄N₃₂, which uses FC2 3×3×3, FC3 2×2×2 and a 4 Å pair cutoff on the same mesh, because the full calculation is prohibitively expensive for its cell.
5. **DFT validation** of 11 structures (Table 3): VASP PBEsol + D3, Pheasy force constants (third order to the 8th neighbour shell), ShengBTE iterative solution. These settings differ from the MLIP screen; see paper Table 1. The resulting tensors are in `data/table3_dft.csv`. Raw force constants and BTE output are available from the corresponding author.

## Settings modes

| `mode` | FC2 | FC3 | FC3 pair cutoff |
|---|---|---|---|
| `baseline` | 4×4×4 | 2×2×2 | none |
| `cutoff_4_fc2_333` | 3×3×3 | 2×2×2 | 4 Å |

**Initial screen.**
- `baseline` on an 11×11×11 mesh was used for 237 structures (2–29 atoms).
- `cutoff_4_fc2_333` was used for the rest:
  - all cells above 29 atoms;
  - five 28-atom cells and one 24-atom cell;
  - structures whose first calculation did not complete. These were rerun at 21×21×21; all of them except C₂₄N₃₂ are dynamically unstable or failed again.

The `mode` and `q_mesh` columns of `data/screening_initial.csv` give the settings for every structure.

**Incomplete runs.** Six of the 270 have no κ:
- three produced no output (`status = no_output`; formula, space group and size come from the input CIF);
- three failed (`status = failed`).

None of the six enters any count beyond the initial 270.

## Files → paper

| File | Content | Paper |
|---|---|---|
| `data/screening_initial.csv` | Initial screen of all 270: settings, stability, lowest mode, κ̄ and diagonal κ, selection flags | Table 1, Fig. 2, Section 2.1 |
| `data/structures_above100/*.cif` | The 43 LeMaT-Bulk inputs with κ̄ > 100 | Deduplication input |
| `data/polymorph_groups.csv` | 43 → 24 polymorph groups | Section 2.1 |
| `data/structures_input/*.cif` | The 25 LeMaT-Bulk input structures | Table 2 |
| `data/structures_relaxed/*.cif` | The same, relaxed with NequIP-OAM-L | Table 2 |
| `data/mlip_results/*.json` | Raw per-structure output behind Table 2 | Table 2 |
| `data/table2_mlip.csv` | Full MLIP κ tensor, principal values κ₁ ≥ κ₂ ≥ κ₃, ρ = κ₁/κ₃, settings, stability, source data | Table 2, Fig. 3 |
| `data/table3_dft.csv` | First-principles κ tensor (DFT cell frame), principal values, ρ, γ̄, ratio to κ̄_MLIP, DFT settings | Table 3, Figs. 1 and 3 |
| `data/mlip_crosscheck/` | The same structures computed on another device | README only |
| `code/build_screening_table.py` | Builds `screening_initial.csv` from the raw cluster output (not distributed) | |
| `code/dedup_structurematcher.py` | Builds `polymorph_groups.csv` from the repo files | |
| `code/build_table2.py` | Builds `table2_mlip.csv` from `mlip_results/` | |
| `code/build_table3.py` | Builds `table3_dft.csv` from the DFT results (values and their source in the file) | |

**Columns.**
- `source_id`: the current Materials Project, OQMD or Alexandria entry ID, as cited in the paper.
- `source_id_lemat`: the ID stored in the LeMaT-Bulk CIF. For OQMD this is an older snapshot ID; both refer to the same structure.
- `lemat_id`: the LeMaT-Bulk `immutable_id`.
- `kappa_*`: W m⁻¹ K⁻¹ at 300 K, in the Cartesian frame of the corresponding cell.
- `dHf_eV_atom`, `E_hull_eV_atom`: as published by the source database, not recomputed.
- `dft_validated`: one of the 11 structures in Table 3.
- `tensor_note`: the C₂₄N₃₂ run stored only the diagonal. Its cell is orthorhombic with its axes along x, y and z, so the off-diagonal terms are zero by symmetry.
- `offdiag_source` (Table 3): for C₁₂N₄ C2/c and C₆N₂ C2/m, the off-diagonal terms are the values reported in the paper (one decimal). They reproduce the reported principal values.

## Checks

- **Devices:** the same structure computed on CPU and GPU agrees within 0.2% (`data/mlip_crosscheck/`: T2_02 0.014%, T2_12 0.07%, T2_23 0.16%).
- **Reproduction:** A T2_02 rerun from scratch reproduced Table 2 to 0.014%.
- **Reduced force constants:** the two database copies of the 28-atom P31c C₃N₄ polymorph, at full and reduced settings on the same mesh, differ by 1.7%.
- **C₁₆N₂ (Fddd) is marginally stable:** its lowest mode is −0.037 THz, inside the tolerance, and κ₃ = 0.38 W m⁻¹ K⁻¹. It is unstable at the DFT level.

## Rebuilding the tables

The MLIP calculations used NequIP-OAM-L 0.1 (mir-group), https://zenodo.org/records/18775904. The checkpoint's SHA-256 is in `code/nequip_oam_l_metadata.json`. The included scripts rebuild the tables and polymorph groups from the distributed results and structures.

```bash
conda env create -f environment.yml
conda activate nequip-phono3
python code/build_table2.py
python code/build_table3.py
python code/dedup_structurematcher.py
```

`environment.yml` pins the main packages: Python 3.10, phono3py 3.30.1, phonopy 3.5.0, NequIP 0.17.1, PyTorch 2.5.1, ASE 3.28.0, pymatgen 2025.10.7. `environment_full.yml` is the full export of the environment used.

## License

Code: MIT (`LICENSE`). Data in `data/`: CC BY 4.0. The input structures come from LeMaT-Bulk and its source databases; please cite them as in the paper.
