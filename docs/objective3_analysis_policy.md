# Objective 3 Analysis Policy

This policy locks the downstream interpretation for parasite-positive blood-fed mosquito specimens.

## Cohort

- Field analysis includes parasite-positive, sequenced mosquito specimens only.
- Positive and negative controls are retained for QC and excluded from field estimates.
- Ct is treated as measured Pf qPCR signal, not calibrated parasitaemia.

## Thresholds

- Primary thesis/manuscript callability threshold: `>=50x`.
- Sensitivity/permissive threshold: `>=10x`, reported separately.
- High-depth context threshold: `>=100x`.
- Missing or non-callable loci are never counted as wild type.

## Denominators

- Marker frequencies use marker-callable denominators.
- Gene/haplotype summaries use samples callable at all required constituent codons.
- Cross-gene dhfr/dhps results are genotype patterns, not phased haplotypes.

## Interpretation

- Mosquito-derived genotypes are reported as alleles recoverable from parasite-positive blood-fed mosquitoes.
- Results are not direct human infection prevalence, clinical treatment failure, parasite origin, or selection evidence.
- K13 variants are classified against `k13_who_artemisinin_partial_resistance_markers.csv`; all other non-synonymous variants are variants of unknown or unestablished significance unless independently supported.
