# Changelog

## v2.2.0 — 2026-07-11

Pooled-cohort analysis, CSP population genetics, and a self-contained test suite.

### Correctness — sample mis-assignment when runs are pooled

The analysis stages keyed samples on **`ont_barcode` alone**:

```r
ont_barcode <- gsub(".*_barcode", "barcode", vcf_filename)   # -> "barcode01"
merge(df, meta, by = "ont_barcode")
```

That is only safe because the metadata was filtered to a single run. `barcode01`
exists on **every** flowcell, so pooling the runs would have silently multiplied
rows and attached the wrong sample to a genotype. The join key is now
**`(ont_multiplex_group, ont_barcode)`**, with the run recovered from the
filename (`nr_run_from_file()` in `bin/_setup.R`).

Related: `merge()` re-sorts by the join key, but `dataFiles` is in `list.files()`
order and is indexed by row position (`dataFiles[-empty_vcfs]`). Those orderings
happened to agree; the order is now restored explicitly and asserted.

### Analyses run over the pooled cohort

`DRUG_RESISTANCE` and the new `POPULATION_GENETICS` now run **once over every
flowcell in `--input`**, not once per flowcell — allele frequencies and population
structure are properties of the sample set, so per-flowcell denominators were
wrong. New `--cohort_name` labels the pooled set. Per-flowcell QC (`01_qc/`,
`02_coverage/`, `03_variants/`) is still emitted per run.

### New: CSP population genetics (subworkflow 06, modules 17–21)

Integrated from the archived, never-wired analysis. csp is the RTS,S/R21 vaccine
target, so this is a **different question** from drug resistance: parasite
population structure and diversity.

`17 CSP_PREPARE_VCFS` → `18 CSP_MERGE_VCFS` → `19 CSP_POPULATION_MAF` →
`20 CSP_GENEFLOW_PCA` → `21 CSP_REPORT`, producing haplotype diversity, the
TH2R/TH3R haplotype network, per-zone MAF, nucleotide diversity (π), Tajima's D,
pairwise F<sub>ST</sub>, PCA, and an HTML report. Ecological-zone grouping is a
parameter (`--csp_*_regions`), not hard-coded. New container
`containers/Dockerfile.popgen` (bcftools, vcftools, plink, tabix + pegas,
adegenet, vcfR, rmarkdown).

`beagle.29Oct24.c8e.jar` was dropped — no script references it.

### New: self-contained test suite

`assets/test_data/` bundles 2 flowcells × 2 barcodes (subsampled reads, truncated
sequencing summary — the real one is 3 GB), 39 MB total. `-profile test` no longer
depends on the study data, and because it spans **two** runs it exercises the
pooled-cohort path and the (run, barcode) join — a single-run test would have
missed the bug above.

`bin/download_upstream_test_data.py` restores the upstream Sanger fetcher, with a
header stating its scope: it checks the alignment/calling core, but cannot
exercise the resistance or CSP stages (R9.4.1 chemistry, no Pf panel, no sample
sheet).

Run samplesheets now accept **relative** paths (resolved against `projectDir`), so
the bundled test data works from any launch directory.

### Marker provenance moved into the pipeline

`docs/provenance/` — the MalariaGEN eLife-2015 k13 table, ART-R interpretation
sheets, `DR_variant_info.xlsx` v1, amplicon/K13 sequence documents, and the
exploratory R precursors. None is read at runtime, but they are the evidence trail
for `assets/resources/`, so they belong with the pipeline rather than in an
archive. (`codon_num_blank.csv` appears in the R code only inside a commented-out
`write.csv()` — it was an output, never an input.)

### Removed

`04_workflow/.retired/` deleted after verifying every v1 process is ported:

| v1 | now |
|---|---|
| `SAMTOOLS_VIEW_SAM_TO_BAM`, `SAMTOOLS_SORT_AND_INDEX` | merged into the piped `MINIMAP2_ALIGN` |
| `GET_CHROM_SIZES_AND_INDEX` | `SAMTOOLS_FAIDX` |
| `NANOPLOT_QC`, `NORMALISE_FASTAS`, `BGZIP_AND_INDEX_VCF` | `NANOPLOT`, `NORMALISE_FASTA`, `BGZIP_TABIX` |
| `*_VARIANT_CALLING` | `CLAIR3` / `MEDAKA` / `MEDAKA_HAPLOID` / `FREEBAYES` |
| `DOWNSTREAM_ANALYSIS` (monolith) | modules 12–15 |

The upstream `Dockerfile` was **not** rewritten — it was an all-in-one dev image
(minimap2 + sniffles + bedtools + samtools + Java + Nextflow) for running the
pipeline inside one container. Obsolete here: every process pulls its own pinned
biocontainer. It also bundled `sniffles`, which nano-rave never calls.

The upstream `nf-test` harness was deleted: it targeted the removed `main.nf`,
tested medaka, and needed test data that was never present. **Adding a proper
nf-test suite is a recommended follow-up.**

### Known gotcha

**Editing anything in `bin/` invalidates the entire `-resume` cache** — Nextflow
fingerprints the `bin/` directory into every task hash. Observed the hard way: an
R-script fix forced a full recompute of 13,165 tasks.

---

## v2.1.0 — 2026-07-11

Structural consolidation. `06_downstream_analysis/` removed; the resistance
analysis became part of the pipeline rather than a separate phase.

- Modules renamed by execution order (`01_sort_fastqs` … `16_dump_versions`;
  alternative callers share slot `09a`–`09d`). `downstream_` prefixes dropped.
- Subworkflows serialised: `01_prepare_reads` … `05_drug_resistance`.
- R stages moved to `bin/`, made executable and self-locating (Nextflow puts
  `bin/` on PATH). Resources flattened into `assets/resources/`; `nanorave_files/`
  → `assets/references/`.
- Fixed `DUMP_VERSIONS` mis-nesting versions from processes whose script contains
  an unindented heredoc (Nextflow strips only the *common* indentation).

## v2.0.0 — 2026-07-11

Full rebuild as a modular Nextflow DSL2 project.

- **Reference bound to its BAM by a meta map.** Previously each process in the
  alignment chain re-emitted the reference on a *separate* channel and the next
  paired it with the BAM **by channel position** — Clair3's `--ref_fn` was bound to
  its BAM only by ordering.
- **A run where every barcode is filtered out now fails** instead of reporting
  SUCCESS having processed zero samples.
- **Stage 2 no longer crashes when a run has no kelch13 variants**
  (`k13_snps_allsamps_df` was built conditionally but consumed unconditionally).
- Multi-run `--input`; Clair3 model as a staged input rather than a Docker
  bind-mount; `minimap2 → view → sort` collapsed into one piped process;
  nf-schema validation; `versions.yml` everywhere; `stub:` blocks throughout.
- Runs under the strict syntax parser — `NXF_SYNTAX_PARSER=v1` no longer needed.
