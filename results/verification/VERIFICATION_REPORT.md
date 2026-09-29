# Independent verification of the GBM / Ti3C2O2 manuscript (2026-09-29)

Scripts: `src/verification/verify_docking_identity.py`, `verify_text_numbers.py`
(neither imports the pipeline). The lattice scan and flake tests were verified in the
audit paper (methods-pipeline-2d-nanocarriers/data/VERIFICATION_REPORT.md).

| What | How | Result |
|---|---|---|
| Identity of the 33 structures | InChIKey from SMILES vs table vs PubChem record of the CID (batch) | 33/33 identical; erlotinib SMILES = CID 176870 |
| Vina scores | best REMARK VINA RESULT of each pose file | 32/32 (bortezomib not docked, as the Methods state) |
| Redocking, every mode | RDKit CalcRMS vs AQ4 of 1M17 chain A | max difference 0.0004; top 5.94 / 7.61 Å, crystal-like mode 6 (1.13 Å) and mode 7 (1.44 Å) |
| Mode-score gaps | pose files | 0.21 and 0.27 kcal/mol (text identical); spread 0.275 and 0.451, so the text "within 0.45" became "within 0.46" (value now rounded up) |
| Families, range, medians, rank, Spearman | recomputed | 5/5/19/4; −4.7 to −9.7; best abemaciclib, larotrectinib, palbociclib; medians −6.0 / −8.6; erlotinib −7.3, 25th of 32; ρ = −0.66 (all identical) |
| Residue contacts | re-implemented (4.0 Å heavy; N/O 3.5 Å) | identical to residue_contacts.csv; Leu820, Val702 100%, Met769 75%, Thr766 94% (polar 66%) |
| His964 / Leu977 | distance to erlotinib in 1M17 | 22.2 and 27.3 Å (text 22 and 27) |
| Procarbazine length | lowest of 8 MMFF conformers | 13.8 Å vs cell 12.12 Å (text identical) |
| Lattice constant / metallic | own PBE-D3 primitive cell | 3.022 Å, DOS(E_F) 1.03 /eV/f.u. (replaces an unverified citation) |
| References | Crossref (check_references.py) | 0 problems |
| QSPR | full nested CV + 1,000 permutations re-run | PENDING at time of writing (run started 2026-09-29) |

## Open before submission
- QSPR rerun comparison (see above).
- DFT campaign (Table 1, Fig. 5, S3, S6) not complete.
- Branch `rebuild-2026-09` not pushed; availability statements true only after push.
