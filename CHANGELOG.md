# Changelog

## v2.3.0 — 2026-07-12

Genotype-call correctness. Three defects that each silently produced a *confident
wrong answer* rather than an error, plus a scope cut.

### Fix — `dhps A437G` was reported inverted

The 3D7 reference **carries** 437G (its dhps genotype is `SGKAA`). The frequency
table nevertheless counted "matches reference" as wild-type, so a locus that is
**100% mutant** was reported as 0% — and it contradicted the pipeline's own
haplotype output, which had the polarity hand-encoded correctly.

Polarity is now **derived** from the codon table (`DR_variant_info_v2.xlsx`:
translate `codon_ref` → if it ≠ `aa_ref`, the reference allele *is* the mutant)
and applied in the reporting layer, `bin/03_resistance_frequencies.R`. It is
deliberately **not** applied to the raw genotype matrix, because the dhps
haplotype builder already depends on the "non-reference" semantics. A
`<run>_marker_polarity.csv` audit file is emitted.

### Fix — `AF` was a character column

`vcfR` returns INFO fields as strings. `mean(AF)` therefore returned `NA`
*silently* while `sd(AF)` worked (it coerces internally), and `AF >= 0.8` was a
**lexicographic** string comparison that happened to give the right answer. All
of `AF`/`DP`/`Qual`/`Pos` are now coerced with `as.numeric()` and asserted.

### New — explicit 3-class call rule and a systematic-artefact catalogue

Four positions were called in nearly every specimen, always heterozygous, always
at AF ≈ 0.5 with almost no variance — including k13 **C580Y**, which is
essentially absent from Africa. 89/94 specimens at AF 0.535 ± 0.046 is not
biology; it is a systematic mapping/basecalling artefact.

Every call is now classified `clonal` | `mixed` | `artefact` | `ref`. A position
is an artefact when it is called in ≥40% of field specimens with mean AF in
[0.40, 0.60] and SD < 0.10 (thresholds are env-configurable; the 0.40 recurrence
default was chosen from a 0.80→0.20 sweep on the real cohort, where 0.50–0.20 is
a stable plateau that flags exactly the four known artefacts and never a genuine
variant — real variants sit at AF 0.84–0.95). Artefact positions report **NOT
CALLABLE** (`NA`), not wild-type. Catalogue: `<run>_artefact_catalogue.csv`.

Consequence: `mdr1 Y184F`, `dhps K540N` and `dhps A581G` are now correctly
reported as not callable rather than as confident calls. The SP conclusion is
unaffected — the quintuple/sextuple hinges on **K540E** (codon 1618), a different
position, which is genuinely callable and genuinely 0/413.

### New — `per_call.tsv`

`11b PER_CALL_TABLE` flattens each VCF to one row per (sample, position) with
`GT`, `DP`, `AF` and ref/alt depths, published to `03_variants/<run>/`. The
analysis stages consume a genotype *matrix*; the evidence behind each cell was
previously unauditable.

### Changed — `min_cov` default is 50, not 10

Girgis 2023 (>50×) and Runtuwene 2018 (≥50 reads) converge on 50× independently;
nothing in the curated literature supports 10×. 10× is now the **secondary**
analysis (`-profile sensitivity_10x`), and a conclusion that holds only at 10× is
not a conclusion.

### Removed — CSP population genetics

csp is the RTS,S/R21 vaccine target: a different scientific question. Subworkflow
`06_population_genetics`, modules 17–21, the `csp_*` scripts, and
`containers/Dockerfile.popgen` are deleted, along with `--population_genetics` and
`--csp_*_regions`. (csp remains an amplicon in the panel, so its *coverage* is
still reported as QC.) Recoverable from git history if the question returns.

### Fix — controls collided across flowcells

`Control_KH2`/`NC` are re-sequenced on **every** flowcell, so `sample_id` was not
unique and the six per-gene merges in stage 1 fanned out to 5⁶ = 15,625 rows.
`sample_id` is now disambiguated (`KH2__RUN23`) and both control filters match on
the **base** id — without which the C580Y-carrying control would have leaked into
the field ART-R frequencies.

## v2.2.0 — 2026-07-11

Pooled-cohort analysis and a self-contained test suite.

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

`DRUG_RESISTANCE` and the new `POPULATION_GENETICS` now ran **once over every
flowcell in `--input`**, not once per flowcell — allele frequencies are properties of the sample set, so per-flowcell denominators were
wrong. New `--cohort_name` labels the pooled set. Per-flowcell QC (`01_qc/`,
`02_coverage/`, `03_variants/`) is still emitted per run.

### New: self-contained test suite

`assets/test_data/` bundles 2 flowcells × 2 barcodes (subsampled reads, truncated
sequencing summary — the real one is 3 GB), 39 MB total. `-profile test` no longer
depends on the study data, and because it spans **two** runs it exercises the
pooled-cohort path and the (run, barcode) join — a single-run test would have
missed the bug above.

`bin/download_upstream_test_data.py` restores the upstream Sanger fetcher, with a
header stating its scope: it checks the alignment/calling core, but cannot
exercise the resistance stages (R9.4.1 chemistry, no Pf panel, no sample
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
