# Reference data and CSV validation summary

The seven CSVs in this directory are **header-only templates** for the Changjiang Cup 2026. The repository contains zero designs and zero scores.

## Reference methodology

The [track manifest](../data/targets/tracks.json) records the seven human UniProt accessions, inclusive residue ranges, endpoint residues, SHA-256 sequence checksums, structural references, and track-specific conditions. Each FASTA contains the exact region extracted from its named UniProt sequence. Isoform-specific accession suffixes are retained. The cache is from UniProtKB release 2026_04.

## Validation

The offline validator checks the CSV filenames, column order, row size, target metadata, 20-row track limit, sequence ID and sequence uniqueness across the supplied files, standard amino-acid alphabet, and inclusive molecule-specific length limits. It checks that stored lengths match the sequences, row ranks are consecutive, scores are finite numbers, and score methods are labeled. It rejects exact terminal matches to common His, FLAG, HA, Myc, and Strep-II tags.

Passing these checks confirms file and sequence format. Signal peptides, unlisted tags and fusions, domain completeness, CDR boundaries, and same-target novelty need separate evidence. Exact sequence deduplication covers only the supplied CSVs. Score labels are checked for format. The validator does not establish the score's provenance or measure binding, agonism, developability, or affinity.

## Reproducibility

Use Python 3.12.16 and run the following commands from the repository root. Both commands work offline with the checked-in files.

```sh
python cup.py check-targets
python cup.py validate --allow-empty
```

`python cup.py write-templates --directory ./local-templates` recreates the headers without modifying existing CSVs. `python cup.py fetch-targets` restores the reference data over HTTPS, verifying the pinned checksums before replacing cached files. The format tool provides no candidate generator or structure and affinity scoring.
