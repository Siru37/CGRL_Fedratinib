[README.md](https://github.com/user-attachments/files/32186192/README.md)
# DHRS13–KPNA contact analysis

Calculate residue-level contacts between DHRS13 and KPNA1–7 from AlphaFold 3 server CIF files. Each contact is defined by the minimum heavy-atom distance between a pair of residues.

## Files

- `analyze_dhrs13_kpna.py`: reusable command-line program for contact analysis.
- `requirements.txt`: Python dependency list.
- `original_scripts/`: seven scripts used for the initial analysis, preserved for provenance. They use the original working-directory layout; use the reusable program for new runs.

## Requirements and usage

Requires Python 3.9 or later and NumPy.

```bash
python -m pip install -r requirements.txt
python analyze_dhrs13_kpna.py --input-dir cif_files --output-dir results
```

Place the CIF files in `cif_files`. Each filename must contain `kpna1`, `kpna2`, ..., or `kpna7`, case-insensitively. Include at most one CIF per KPNA in each input directory. A subset of the seven proteins is also supported.

Example filename:

```text
fold_2025_03_28_11_34_dhrs13_iso_1_kpna1_model_0.cif
```

The default assignments are chain A for DHRS13 and chain B for KPNA. Specify alternate labels or a different distance cutoff with:

```bash
python analyze_dhrs13_kpna.py --input-dir cif_files --output-dir results --dhrs13-chain A --kpna-chain B --cutoff 5.0
```

The cutoff is in angstroms. Chain identity must be verified independently from the input sequences; the program does not infer biological identity from coordinates.

## Outputs

For each input model, the program writes:

- `KPNA<N>_interface_contacts_5A.csv`: residue-pair contacts at the default cutoff.
- `KPNA<N>_summary.json`: source filename, sequences, numbering checks, contact count, and statistics for selected DHRS13 regions.

The cutoff value in the CSV filename changes when `--cutoff` is changed. This Python program does not create XLSX files. The accompanying seven-sheet workbook was assembled separately from the contact CSV files.

### Contact CSV columns

| Column | Definition |
|---|---|
| `DHRS13_residue` | DHRS13 residue number from `label_seq_id` |
| `DHRS13_aa` | One-letter amino-acid code |
| `KPNA<N>_residue` | KPNA residue number from `label_seq_id` |
| `KPNA<N>_aa` | One-letter amino-acid code |
| `minimum_heavy_atom_distance_A` | Minimum distance between the residue pair, in angstroms |
| `DHRS13_atom` | DHRS13 atom in the closest atom pair |
| `KPNA<N>_atom` | KPNA atom in the closest atom pair |
| `DHRS13_mean_pLDDT` | Mean per-atom pLDDT within the DHRS13 residue |
| `KPNA<N>_mean_pLDDT` | Mean per-atom pLDDT within the KPNA residue |

## Calculation method

1. Read the `_atom_site` records for the selected chains and exclude hydrogen and deuterium.
2. Compute all heavy-atom Euclidean distances for each interchain residue pair using NumPy float64 arithmetic.
3. Retain a residue pair when its minimum distance is less than or equal to the cutoff. Each retained pair produces one CSV row, not one row per atom pair.
4. Apply the cutoff before rounding. Report distances to three decimal places and pLDDT values to two decimal places. If closest-atom distances tie, retain the first pair in the CIF record order.
5. Read per-atom pLDDT from the AF3 `B_iso_or_equiv` field and average over the heavy atoms within each residue.

The JSON region summaries cover DHRS13 residues 170–190, 309–377, 353–372, and 363–370. Region mean pLDDT is the unweighted mean of residue-level values across all modeled residues in the region, not only contacting residues. These intervals are analysis windows, not automatically detected domains or NLS annotations.

Output numbering uses `label_seq_id`. The JSON reports whether it agrees with `auth_seq_id` for the selected chains.

## Scope and limitations

- The parser targets single-model AF3 server CIF files with standard amino acids and one atom record per line. It is not a general-purpose mmCIF parser.
- The interpretation of `B_iso_or_equiv` as pLDDT is specific to the AF3 prediction files used here. Experimental B factors must not be interpreted as pLDDT.
- The program does not automatically assign ARM repeats or major/minor NLS-binding pockets. Pocket interpretation requires an independent mapping of conserved residues or structural alignment.
- Distance alone does not establish a hydrogen bond, salt bridge, or steric clash. Atom types and geometry require separate assessment.
- Contacts and local pLDDT do not establish binding affinity, physiological interaction, or NLS function.
- The analyzed CIF files do not contain PAE or ipTM. The program does not estimate or generate these metrics.

## Reproducibility

The reusable program was checked against all seven original CIF files. Every CSV header and data row matched the initial results, including residue identities, closest atom pairs, rounded distances, and residue-level pLDDT values.

| Model | Contact residue pairs at 5 Å |
|---|---:|
| KPNA1 | 109 |
| KPNA2 | 109 |
| KPNA3 | 154 |
| KPNA4 | 107 |
| KPNA5 | 140 |
| KPNA6 | 120 |
| KPNA7 | 193 |
| Total | 932 |

These counts apply to the specific input models used for this analysis, not to all possible DHRS13–KPNA predictions.

The archived original scripts expect `work/model.cif` for KPNA6 and `work/kpna1.cif`, etc., for the remaining proteins. The original KPNA6 script also attempts a UniProt download after contact calculation; that step is not required for CSV generation. The reusable program performs no network requests.
