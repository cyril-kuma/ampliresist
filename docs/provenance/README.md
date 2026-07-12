# Marker provenance

Source material behind the drug-resistance marker definitions the pipeline uses.
**None of these files is read at runtime** — they are the evidence trail for how
`assets/resources/` was derived, and they belong with the pipeline so the marker
calls can be audited.

The pipeline reads exactly five resource files, all in `assets/resources/`:

| File | Used by | For |
|---|---|---|
| `amplicon_gene_positions.csv` | stage 12 | amplicon start/end per gene |
| `DR_variant_info_v2.xlsx` | stage 13 | the SNP panel that is genotyped |
| `genetic_code.csv` | stage 13 | codon → amino-acid translation |
| `k13_seq_annotated.csv` | stage 13 | kelch13 codon annotation |
| `k13_who_artemisinin_partial_resistance_markers.csv` | stage 13 | WHO validated/candidate ART-R markers |

## What's here

- **`MalariaGEN Genomic Epi of ART-R (eLife 2015) k13 mutations table.xlsx`** — the
  published k13 mutation table underpinning the ART-R marker list.
- **`k13_ARTR_interpretation.xlsx`**, `k13_3d7_seq*.txt`, `K13 sequence.docx` —
  kelch13 reference sequence and interpretation working files behind
  `k13_seq_annotated.csv`.
- **`DR_variant_info.xlsx`** — v1 of the SNP panel. The pipeline uses **v2**.
- **`genetic_code.xlsx`** — spreadsheet original of `genetic_code.csv`.
- **`Amplicon target sequences with positions.docx`**,
  `drug_resistance_calls_ref_genomes.xlsx` — amplicon design behind
  `amplicon_gene_positions.csv`.
- **`DR_genotype_calls.xlsx`**, `pf_drug_resistance_calls_.xlsx` — earlier
  hand-curated call sets, kept for comparison.
- **`codon_num_blank.csv`** — a scratch file. It appears in `02_genotype_calls.R`
  only inside a commented-out `write.csv()`; it was an *output*, never an input.
- **Five loose `R script …` files** — exploratory precursors of the numbered
  stages (`R script for handling vcf files.R`, `… - looped.R`,
  `… - non-looped for 8 samples.R`, `R script to call k13 variants.R`,
  `R script to check amplicon coverage (non-looped).R`). Superseded by
  `bin/01_amplicon_coverage.R` … `bin/04_summary_tables.R`.
