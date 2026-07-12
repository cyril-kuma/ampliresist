# R environment for the drug-resistance analysis stages (modules 12-15).
#
# Build:
#   docker build -t drag1-downstream:1.0 -f containers/Dockerfile containers/
#
# Unrelated to the upstream nano-rave Dockerfile, which was an all-in-one dev
# image (minimap2 + sniffles + bedtools + samtools + Java + Nextflow) for running
# the entire pipeline inside a single container. That approach is obsolete here:
# every other process pulls its own pinned biocontainer, so this image carries
# only R and is used by the analysis stages alone.
#
# rocker/tidyverse ships tidyverse (incl. readxl, dplyr, ggplot2, purrr, readr,
# stringr, tibble, tidyr) and is configured to use Posit Package Manager, so the
# extra packages install as precompiled Ubuntu binaries rather than from source.
FROM rocker/tidyverse:4.4.1

# System libs needed by vcfR (and friends) at build/run time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Packages the downstream R stages load explicitly, beyond tidyverse:
#   vcfR      - stage 2 VCF parsing
#   data.table, reshape2, cowplot, viridis - stages 1-5 wrangling/plots
#   yaml      - run/_common.sh config parsing (Rscript-based)
RUN install2.r --error --skipinstalled --ncpus -1 \
        vcfR \
        data.table \
        reshape2 \
        cowplot \
        viridis \
        yaml

# Fail the build if anything is missing, rather than at pipeline runtime.
RUN Rscript -e 'pkgs <- c("tidyverse","vcfR","readxl","data.table","cowplot","viridis","reshape2","yaml","dplyr","ggplot2","stringr","tibble","tidyr","purrr","readr"); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop("MISSING: ", paste(missing, collapse=", ")); \
    cat("All R packages present\n")'

WORKDIR /analysis
