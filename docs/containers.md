# Container manifest & strategy

Every process runs in a container. Public tool images are pinned by **immutable
digest** in each module; the two ampliresist-**owned** images are built from
`containers/` and referenced through `--container_registry` (default
`ghcr.io/owner`, a placeholder — replace `owner` with your GitHub org/user).

No process depends on an image that exists only in a local Docker cache: every
image is either a public digest-pinned biocontainer or is reproducibly built from
a version-controlled Dockerfile here.

## Image ↔ process map

| Image | Pin | Processes |
|---|---|---|
| `quay.io/biocontainers/nanoplot` | `@sha256:57cffeb1…` (1.38.0) | NANOPLOT |
| `quay.io/biocontainers/pycoqc` | `@sha256:ea0a0847…` (2.5.2) | PYCOQC |
| `quay.io/biocontainers/biopython` | `@sha256:8bdeb52f…` (1.78) | NORMALISE_FASTA |
| `quay.io/biocontainers/samtools` | `@sha256:c3e0ba2a…` (1.15.1) | SAMTOOLS_FAIDX |
| `quay.io/biocontainers/minimap2` | `@sha256:1f23d5cf…` (2.24) | MINIMAP2_INDEX |
| `quay.io/biocontainers/mulled-v2-6653…` | `@sha256:1a2d2f99…` (minimap2 2.24 + samtools 1.15) | MINIMAP2_ALIGN |
| `quay.io/biocontainers/bedtools` | `@sha256:9199479a…` (2.29.2) | BEDTOOLS_GENOMECOV |
| `docker.io/hkubal/clair3` | `@sha256:57cf5d20…` | CLAIR3 (model supplied separately — see `assets/references/clair3_models/PROVENANCE.md`) |
| `quay.io/biocontainers/tabix` | `@sha256:106e72ca…` (1.11) | BGZIP_TABIX, GUNZIP_VCF |
| `quay.io/biocontainers/python` | `@sha256:f6b44640…` (3.10) | PER_CALL_TABLE, DUMP_VERSIONS |
| **`${container_registry}/ampliresist-r:1.0.0`** | built — `containers/r.Dockerfile` | AMPLICON_COVERAGE, GENOTYPE_CALLS, RESISTANCE_FREQUENCIES, SUMMARY_TABLES, COMPLEXITY_OF_INFECTION, CLEAN_METADATA |
| **`${container_registry}/ampliresist-figures:1.0.0`** | built — `containers/figures.Dockerfile` | PUBLICATION_PLOTS |

The alternative variant callers (`MEDAKA`, `MEDAKA_HAPLOID`, `FREEBAYES`) are
inherited from upstream nano-rave and **not** on the validated (`clair3`) path.

### Provenance of the two owned images

- **`ampliresist-r`** — formerly `drag1-downstream:1.0`. `FROM rocker/tidyverse:4.4.1`
  plus vcfR/data.table/reshape2/cowplot/viridis/yaml. Renamed for ownership; build
  is unchanged, so outputs are identical (verified: R package versions match).
- **`ampliresist-figures`** — **replaces the borrowed `haema-figures:0.4.0`** (a
  HÆMA-project image). `bin/07_build_roadmap_figures.py` imports only
  numpy/pandas/matplotlib (GeoJSON via stdlib `json`; the penalised logistic model
  is hand-rolled numpy IRLS), so the seaborn/geopandas/shapely stack in the old
  image was unused. The lean replacement pins numpy 2.1.3 / pandas 2.2.3 /
  matplotlib 3.9.2 (+ procps for Nextflow metrics). **Verified: figure source CSVs
  are byte-identical between the old and new image.**

## Build, verify, publish (owned images)

```bash
cd 04_workflow
REG=ghcr.io/<your-org-or-user>        # your namespace (lowercase)

# Build
docker build -t $REG/ampliresist-r:1.0.0       -f containers/r.Dockerfile       containers/
docker build -t $REG/ampliresist-figures:1.0.0 -f containers/figures.Dockerfile containers/

# Smoke-test
docker run --rm --entrypoint Rscript $REG/ampliresist-r:1.0.0 \
  -e 'stopifnot(all(c("readxl","dplyr","vcfR") %in% rownames(installed.packages()))); cat("r OK\n")'
docker run --rm --entrypoint python3 $REG/ampliresist-figures:1.0.0 \
  -c 'import numpy,pandas,matplotlib; print("figures OK", numpy.__version__)'

# Publish
docker push $REG/ampliresist-r:1.0.0
docker push $REG/ampliresist-figures:1.0.0

# Run against your registry
nextflow run . -profile docker --container_registry $REG ...
```

## Pull the public images (no local cache required)

```bash
for img in \
  quay.io/biocontainers/nanoplot@sha256:57cffeb19b67d88126a8251069771190a65f73ec694fb2e496e0b4f1a98ef4ab \
  quay.io/biocontainers/pycoqc@sha256:ea0a084751a0b48b5ffe90e9d3adfa8f57473709a1b0a95c9cb38d434ee3a9a2 \
  quay.io/biocontainers/biopython@sha256:8bdeb52fb15b5f61c40292f73d85a3a77cda4bbd95d29e710ddaad7a6bf76720 \
  quay.io/biocontainers/samtools@sha256:c3e0ba2add590177a2e6ea33ae9074dc1f82b99e6913338d3d1c3a70dc78b518 \
  quay.io/biocontainers/minimap2@sha256:1f23d5cfbefb25ef4f9a0ee5b4f78d3b6cb0b3c955028d80e1d8b00bc97e299a \
  quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4@sha256:1a2d2f9958084a835a8897d16d339484202c8a31ad17fded2f14294f26de4d3c \
  quay.io/biocontainers/bedtools@sha256:9199479a4142cb21a7dbd9f215a4f29811d95927cf57b49f7163d3c2d1d70abb \
  docker.io/hkubal/clair3@sha256:57cf5d20f2ee39c1b91493ad1fb5c1b9fa838691efce818c3139caa5e6c6b974 \
  quay.io/biocontainers/tabix@sha256:106e72ca3c7ca98c12b3971ba3d2699f4ec63673976f6037a38ebf1d46727515 \
  quay.io/biocontainers/python@sha256:f6b44640f06e8265ebf5ce85ca12cea53af110c188d6b4acf5f59887c24abb8f ; do
  docker pull "$img"
done
```

## Singularity / Apptainer

The digest-pinned `docker://`/`quay.io` references convert automatically
(`-profile singularity`, `singularity.autoMounts = true`). Build the owned images
to a SIF or push them to a registry Singularity can pull:

```bash
nextflow run . -profile singularity --container_registry docker://$REG ...
```

Images are `linux/amd64`. On arm64 (Apple Silicon) run under emulation or rebuild
the owned images for the local architecture.
