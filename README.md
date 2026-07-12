# ampliresist

ONT amplicon variant calling and *Plasmodium falciparum* drug-resistance
genotyping — **one pipeline, one command**.

Adapted from [sanger-pathogens/nano-rave](https://github.com/sanger-pathogens/nano-rave)
(© 2022, 2023 Genome Research Ltd.).

## Overview

Six subworkflows:

```
runs.csv ─> 01 PREPARE_READS ──┐                      ┌─> 05 DRUG_RESISTANCE
            (sort, nanoplot,   │                      │   12 amplicon coverage ─┐
             pycoqc)           │                      │   13 genotype calls ────┼─> 14 frequencies
                               ├─> 03 ALIGN_AND_      │                         └─> 15 summary tables
            02 PREPARE_        │      COVERAGE ───────┤
               REFERENCES ─────┘   (minimap2|sort,    │─> 04 CALL_VARIANTS
            (normalise,             bedtools)         │   (clair3|medaka|freebayes,
             faidx, mmi)                              │    bgzip/tabix, gunzip)
                                                      │
                                                      └─> 06 POPULATION_GENETICS
                                                          17 prepare ─> 18 merge ─> 19 MAF
                                                          ─> 20 pi/Tajima's D/Fst/PCA ─> 21 report
```

Alignment is the fan-out point: every barcode × every reference amplicon
(95 × 7 = 665 alignments per flowcell).

**The two analyses run over the pooled cohort**, not per flowcell. Allele
frequencies and population structure are properties of the *sample set*, so
splitting them per flowcell would compute the wrong denominators. Per-flowcell QC
(`01_qc/`, `02_coverage/`, `03_variants/`) is still emitted per run.

## What you need to supply

This repository contains the **pipeline only**. Study data and large binaries are
deliberately not published here — you provide them:

| Not in the repo | Why | How to supply |
|---|---|---|
| `assets/samplesheet.xlsx` | sample metadata, incl. collection GPS coordinates | build your own from `assets/samplesheet.template.csv`, save as `.xlsx`; pass with `--samplesheet` |
| `assets/test_data/` | real sequencing reads | point `--input` at your own flowcells |
| `assets/runs.csv` | embeds absolute paths to your data | generate with `bin/make_run_samplesheet.sh` (see `assets/runs.example.csv`) |
| `assets/references/clair3_models/` | ~78 MB binary | download the model matching your basecalling chemistry from [ONT Rerio](https://github.com/nanoporetech/rerio); pass with `--clair3_model` |

The reference amplicon FASTAs (`assets/references/`) and the marker resources
(`assets/resources/`) **are** included — the pipeline needs them to run.

`run_name` in `--input` must match `ont_multiplex_group` in the sample sheet.

## Quick start

```bash
cd 04_workflow

# 0. one-off: build the analysis container
docker build -t drag1-downstream:1.0 -f containers/Dockerfile        containers/

# 1. smoke test — self-contained, uses the bundled dataset in assets/test_data/
nextflow run . -profile test,docker

# 2. build the run samplesheet by scanning the data directory
bin/make_run_samplesheet.sh ../02_data > assets/runs.csv

# 3. full run — all flowcells, one command
nextflow run . -profile docker --input assets/runs.csv --outdir ../05_results -resume
```

Only `--input` is required; references, Clair3 model, sample sheet and resources
all default under `assets/`.

## Inputs

**`--input`** — one row per flowcell. Any number of runs; they are pooled into one
cohort for the analyses.

```csv
run_name,sequencing_dir,sequencing_summary
LUC_DRAG1_IMPAVESRUN23_20250920,/…/20250920_1723_X1_FAZ03468_…,/…/sequencing_summary_….txt
```

> `run_name` **must** match `ont_multiplex_group` in `assets/samplesheet.xlsx`.
> Any run missing from the workbook aborts the analysis rather than being silently
> dropped.
>
> Samples are keyed on **(run, barcode)**, not barcode alone — `barcode01` exists
> on every flowcell, so a barcode-only join would mis-assign samples once the runs
> are pooled.

The resistance analysis needs the **full amplicon panel** (crt, dhfr, dhps, mdr1,
k13, csp, msp1). A cut-down `--reference_manifest` is fine for variant calling
only (`--drug_resistance false`).

## Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `--input` | *(required)* | Run samplesheet |
| `--cohort_name` | `DRAG1_cohort` | Names the pooled sample set in the analysis outputs |
| `--variant_caller` | `clair3` | The validated path. `medaka`, `medaka_haploid`, `freebayes` are inherited from upstream and **not validated for this assay**. |
| `--clair3_model` | bundled R10.4.1 model | Must match the basecaller chemistry — the Clair3 image ships only R9.4.1 models. |
| `--clair3_args` | `--no_phasing_for_fa --include_all_ctgs --haploid_precise` | The first two are **required**: Clair3's phasing/contig filters assume human chromosomes and emit no genotypes on small Pf amplicons. `--haploid_precise` because *P. falciparum* is haploid. |
| `--min_cov` | `50` | Amplicon median coverage **and** per-SNP depth threshold. |
| `--min_barcode_dir_size` | `10` | MB. If *every* barcode is filtered out the run fails rather than producing nothing. |
| `--drug_resistance` | `true` | Stages 12–15. |

Parameters are validated against `nextflow_schema.json` before any work starts.

## Results

```
<outdir>/
  01_qc/nanoplot/<run>/ , pycoqc/    per-flowcell QC
  02_coverage/<run>/                 *.bedGraph
  03_variants/<run>/vcf/ , vcf_unzipped/
  04_resistance/                     pooled cohort
    01_coverage_qc/  02_genotype_calls/  03_frequencies/  04_summary/
  pipeline_info/                     timeline, report, trace, DAG, software_versions.yml
```

## Safeguards

A green Nextflow run only means every process exited 0 — it says nothing about
whether the science is right. This pipeline has already shipped, with exit code 0,
a coverage table whose `sample_id` column contained gene names. These guards exist
because each one corresponds to a failure that actually happened:

| Guard | Catches |
|---|---|
| `VERIFY_RUNS_DISTINCT` | the **same flowcell supplied twice** in `--input`. Fingerprints each run by hashing its read UUIDs; identical reads abort the run. Pooling a flowcell twice pseudoreplicates every sample it carries and silently invalidates the cohort frequencies. Override with `-profile dev` (see below). |
| `_setup.R` sample-id uniqueness | controls (KH2, NC) are re-sequenced on **every** flowcell, so `sample_id` is not unique once runs are pooled. Each sequencing instance is disambiguated (`KH2__<run>`); a residual duplicate aborts. |
| stage 12 fan-out check | the coverage table having more rows than the cohort has samples — the signature of a merge fanning out on a non-unique `sample_id` (one control on 5 runs produced **5⁶ = 15,625** rows). |
| stage 12 join check | no coverage row matching any `sample_id` — a broken join, not a result. |
| stage 14 callability check | *every* marker reporting `n_callable = 0`. |

None of these hardcodes a cohort size: they all derive from the sample sheet
filtered to the runs in `--input`, so they scale from 96 samples to 6,000.

Post-run, `tests/validate_cohort.py` re-checks the published results
independently:

```bash
tests/validate_cohort.py --outdir ../05_results/v2 \
    --input assets/runs.csv --samplesheet assets/samplesheet.xlsx
```

> A clean integer scaling factor between a single run and a pooled cohort is **not**
> evidence of correctness — it is exactly what duplicated input data looks like.

### Development mode

While the multi-batch path is being built, the same flowcell is deliberately
copied under several run names so the pooling logic can be exercised before real
second and third batches exist. `VERIFY_RUNS_DISTINCT` would (correctly) refuse
that, so use the `dev` profile:

```bash
nextflow run . -profile docker,dev --input assets/runs.csv --outdir results
```

It downgrades the abort to a warning and prints a banner in the log. **Results
produced this way are pseudoreplicated and are not scientifically usable** —
`dev` is not a switch to reach for on real data.

## Testing

| | |
|---|---|
| `nextflow run . -profile test,docker -stub-run` | validates the whole channel topology in seconds |
| `nextflow run . -profile test,docker` | **self-contained** end-to-end run on `assets/test_data/` (2 flowcells × 3 **disjoint** barcodes — genuinely different reads, so it also passes VERIFY_RUNS_DISTINCT). Spans two runs, so it exercises the pooled-cohort path and the (run, barcode) join. |
| `bin/download_upstream_test_data.py` | fetches the upstream Sanger nano-rave dataset. Checks the alignment/calling core against upstream, but **cannot** exercise the resistance stages (R9.4.1 chemistry, no Pf panel, no sample sheet). |

## Layout

```
main.nf                  entry workflow — thin; composes the six subworkflows
nextflow.config          params, profiles, reporting
nextflow_schema.json     parameter validation (nf-schema)
conf/                    base (resources) · modules (publishDir, ext.args) · test
modules/local/           one process per file, numbered by execution order
  01_sort_fastqs … 08_bedtools_genomecov
  09a_clair3 / 09b_medaka / 09c_medaka_haploid / 09d_freebayes   (alternatives)
  10_bgzip_tabix · 11_gunzip_vcf
  12_amplicon_coverage … 15_summary_tables      (drug resistance)
  16_dump_versions
subworkflows/local/
  01_prepare_reads  02_prepare_references  03_align_and_coverage
  04_call_variants  05_drug_resistance
bin/                     executables, staged onto PATH inside each task
  make_run_samplesheet.sh  download_upstream_test_data.py
  01_amplicon_coverage.R … 04_summary_tables.R   (+ _setup.R, _pipeline_paths.R)
assets/
  samplesheet.xlsx  resources/  references/  test_data/  runs.csv  test_runs.csv
containers/              Dockerfile (R environment for the analysis stages)
docs/                    analysis policy, marker catalogue, provenance/
```

## Notes

- Runs under Nextflow's **strict syntax parser** (tested on 26.04).
- **Editing anything in `bin/` invalidates the entire `-resume` cache** — Nextflow
  fingerprints the `bin/` directory into every task hash. Plan R-script changes
  before a long run, not during one.
- Editing `nextflow.config` also invalidates the cache; use `-c my.config` for
  one-off overrides.
- The two analysis images are built locally and are not in a registry.
- Retired material is in `../07_archive/`.
