# Changjiang Cup 2026 reference and CSV tools

This repository contains the seven public human UniProt reference regions, seven header-only CSV templates, an offline format validator, and sequence composition reporting. The templates contain zero candidate rows. Candidate design, experimental validation, structure prediction, and ranking are outside this implementation.

## Run locally

Use Python 3.12.16. The CLI uses the standard library and the default checks work offline.

```sh
python cup.py check-targets
python cup.py validate --allow-empty
```

To recreate the seven templates in a new directory:

```sh
python cup.py write-templates --directory ./local-templates
python cup.py validate --allow-empty --directory ./local-templates
```

To validate files you supply, point the validator at a directory with all seven named CSVs. Use `--allow-empty` if some tracks have no rows. Populated rows are checked in both modes. Score method labels have the format `proxy:<method>` or `official-method:<method>`. See the [validation scope](submissions/DESIGN_SUMMARY.md) before using the results.

The FASTA cache comes from UniProtKB release 2026_04. The range is 1-based and inclusive. Restore the cache when online with `python cup.py fetch-targets`. It checks the source accession and human taxonomy, then compares length, endpoint residues, and the SHA-256 of the extracted region against [the manifest](data/targets/tracks.json). A missing accession, network failure, or mismatch exits nonzero. The checked-in FASTAs are sufficient for offline validation.

## Describe sequences

```sh
python cup.py describe reference.fasta another-reference.fasta
```

Pass one or more single-record FASTA files with uppercase standard amino acids. The command writes a JSON array to stdout in input order and works offline without the track data. It reports length, cysteine count, D/E and K/R counts, the K/R minus D/E count, and mean hydropathy. The charge count omits histidine, termini, and pH effects. Mean hydropathy is the arithmetic mean of the published residue values. These descriptors leave the submission CSVs untouched and provide no binding, solubility, or aggregation predictions. [Kyte–Doolittle scale](https://web.expasy.org/protscale/pscale/Hphob.Doolittle.html)

## Tracks

| Track | Target | Type | UniProt region | Reference | Condition |
| --- | --- | --- | --- | --- | --- |
| 1 | EGFR | VHH | P00533-1 Met1–Ser645 | 6ARU A | Affinity at pH 6.0 |
| 2 | TNF-alpha | VHH | P01375 Val77–Leu233 | 1TNF A/B/C | Trimer |
| 3 | GLP-1R | VHH | P43220 Arg24–Ser432 | 6X18 R | Agonist required |
| 4 | TfR1 (CD71) | VHH | P02786 Cys89–Phe760 | 11QC A/B | Dimer |
| 5 | TL1A (TNFSF15) | Minibinder | O95150 Leu72–Leu251 | 2RE9 A/B/C | Trimer |
| 6 | TNFR1 (TNFRSF1A) | VHH | P19438-1 Leu30–Thr211 | 9S3I A | |
| 7 | IL-7Ralpha (CD127) | Minibinder | P16871-1 Glu21–Asp239 | 3DI3 B | |

## Competition rules

Registration and sequence submission both close **November 30, 2026 at 23:59:59**. The seven tracks are independent, with a maximum of 20 distinct de novo sequences per team per track and one molecule type per track. VHHs must be 100–150 residues with a complete variable-domain framework and CDRs. Minibinders must be 40–100 residues in a single chain with one designed binding module. Use the 20 standard amino acids and mature proteins without signal peptides, tags, Fc, or other fusions.

Exclude experimentally validated binders for the same target, including as starting points for affinity optimization. The maximum identity to known same-target binders is 70% for comparisons covering at least 80% of the minibinder sequence or at least 80% of the VHH CDRs, as applicable.

The official open-registration screen uses ESMFold2 and ipSAE with one seed, one sample, no template, and target-only MSA input. Measured KD determines rank from lowest to highest. Only each team's best candidate in each track counts, with six winning entries per track. [Competition site](https://aiprotein.sjtu.edu.cn/)

## Tests

```sh
python -m venv .venv
.venv/bin/python -m pip install --require-hashes -r requirements.txt
.venv/bin/python -m pytest -q
```

The dependency install can use a local wheelhouse with `--no-index --find-links /path/to/wheelhouse`. The tests themselves run offline.
