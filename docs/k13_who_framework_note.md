# WHO-aligned PfK13 interpretation framework for Objective 3

This note documents how Objective 3 should interpret *P. falciparum kelch13* (`PfK13`, `kelch13`, `PF3D7_1343700`) variants detected from parasite-positive blood-fed mosquito extracts. It is a human-readable companion to `k13_who_artemisinin_partial_resistance_markers.csv` and should be cited in Methods or supplementary materials when explaining the analysis policy.

## Sources and date checked

The framework is based on the WHO artemisinin partial-resistance Q&A and the WHO compendium of molecular markers for antimalarial drug resistance. Both were checked on 2026-06-26.

- WHO artemisinin partial-resistance guidance: https://www.who.int/news-room/questions-and-answers/item/artemisinin-resistance
- WHO compendium of molecular markers for antimalarial drug resistance: https://www.who.int/tools/compendium-of-molecular-markers-for-antimalarial-drug-resistance

## Interpretation principle

A non-synonymous *kelch13* variant is not automatically an artemisinin partial-resistance marker. Objective 3 should classify each detected *kelch13* amino-acid change into one of four reporting categories:

1. `validated`: a WHO validated PfK13 marker of artemisinin partial resistance.
2. `candidate_or_associated`: a WHO candidate or associated PfK13 marker.
3. `unestablished`: a non-synonymous variant not currently classified by WHO as validated, candidate, or associated.
4. `not_interpretable`: a low-depth, mixed/heterozygous-pattern, low-complexity artefact, control-contaminated, or otherwise non-confident call.

Only confident field calls in categories `validated` or `candidate_or_associated` should contribute to an Objective 3 count of mosquito-derived parasite-positive specimens with a WHO-classified artemisinin partial-resistance marker. The KH2 positive control is used to validate assay behavior, not to estimate field marker frequency.

## WHO-classified PfK13 markers used in this project

The current local marker file lists the following WHO validated markers:

F446I, N458Y, C469Y, M476I, Y493H, R539T, I543T, P553L, R561H, P574L, C580Y, R622I, and A675V.

It lists the following candidate or associated markers:

P441L, G449A, C469F, A481V, R515K, P527H, N537I, N537D, G538V, and V568G.

The marker list should be rechecked against WHO before final thesis submission or manuscript submission because WHO classifications can change as new clinical, in vitro, and genetic evidence accumulates.

## Objective 3 reporting language

Preferred wording:

"A specimen was classified as having a WHO-classified *kelch13* artemisinin partial-resistance marker only when a WHO validated or candidate/associated PfK13 amino-acid change was confidently called at the pre-specified coverage threshold and passed control and artefact review."

Avoid:

"All non-synonymous *kelch13* variants were counted as artemisinin resistance."

"Mosquitoes carried artemisinin-resistant malaria."

"Absence of a *kelch13* marker proves artemisinin resistance is absent from the study sites."

## Quality-control requirements

For Objective 3, *kelch13* calls must satisfy the same conservative calling rules as the rest of the marker panel. A field call should not be promoted to a resistance-marker call if it is sub-threshold, occurs only as a recurrent heterozygous-pattern signal in a low-complexity context, appears in negative controls, or conflicts with expected positive-control behavior. The pilot analysis showed why this is necessary: a true KH2 C580Y positive-control call was recovered as a high-depth high-allele-fraction mutant call, while recurrent field heterozygous-pattern calls at low-complexity contexts required rejection or independent validation.
