# Callability: why three codons are reported NOT CALLABLE

`mdr1 Y184F`, `dhps K540N` and `dhps A581G` are reported as **NOT CALLABLE**
rather than as a frequency. This document is the evidence for that decision.

It matters because the alternative — reporting them — would have been *wrong in
the confident direction*. `mdr1 Y184F` in particular is a real, common West
African mutation (~50–80% in the region), so a spurious call at 65% would have
looked entirely plausible and would have passed review.

Reproduce with:

```bash
bin/05_call_evidence_audit.sh mdr1 551 work/ barcode95 6
```

## The statistical signature

The artefact rule (`bin/02_genotype_calls.R`) keys on **`frac_het`** — the
fraction of specimens in which a position is called heterozygous. *P. falciparum*
is haploid in the human host, so a mono-clonal infection carrying a variant is
`1/1`. Any real polymorphism, at any frequency, therefore produces homozygotes.

On this cohort the separation is absolute:

| position | marker | af_mean | af_sd | **frac_het** | verdict |
|---|---|---|---|---|---|
| `mdr1_551_AT`  | Y184F | 0.655 | 0.030 | **1.000** | artefact |
| `k13_1739_GA`  | C580Y | 0.536 | 0.046 | **1.000** | artefact |
| `dhps_1620_AT` | K540N | 0.541 | 0.062 | **1.000** | artefact |
| `dhps_1742_CG` | A581G | 0.529 | 0.059 | **1.000** | artefact |
| `dhfr_323_GA`  | S108N | 0.886 | 0.069 | 0.079 | real |
| `dhfr_175_TC`  | C59R  | 0.824 | 0.064 | 0.086 | real |
| `csp_902_GA`   | —     | 0.947 | 0.034 | 0.000 | real |

A position called `0/1` in **every** specimen and `1/1` in **none** is not a
polymorphism. No allele frequency and no mixture of clones can produce that.

The k13 result anchors the rule independently: `C580Y` is essentially **absent
from Africa**, so 413 specimens carrying it at AF 0.536 ± 0.046 cannot be real,
whatever the reads say.

## Ruling out the technical explanations

The statistical signature justifies refusing to report a frequency. It does not
say *why*. Four mechanisms could produce a spurious call; the read-level audit
excludes all four for `mdr1 551`.

**1. Insufficient depth — no.** Depth at the position is 5,000–10,000×. The
amplicon is uniformly covered from ~220–600 bp. This was never a coverage
problem, which is why lowering `min_cov` would not have recovered it.

**2. Mismapping or paralogy — no.** Split by the base each read carries, *both*
alleles are equally well supported:

```
field barcode47   depth=10008
    base A  n=3385  (33.8%)  meanMAPQ=60.0  MAPQ0=0  minus-strand=54.8%
    base T  n=6576  (65.7%)  meanMAPQ=60.0  MAPQ0=0  minus-strand=53.6%
```

MAPQ 60, zero MAPQ-0 reads, balanced strands, on both alleles. If the alternate
allele were mismapped noise, its reads would be MAPQ 0 or strand-skewed while the
reference reads stayed clean. They are not. A MAPQ filter would remove nothing.

**3. Basecaller / sequence-context error — no.** Position 551 sits in a TA
dinucleotide repeat (`TTTA·TATAT·TTGGT`), which is exactly the context where ONT
errs. But the **clonal control settles it**: `Control_KH2` is a cultured lab line
— a single genome, which *cannot* be heterozygous.

```
control KH2       depth=5037
    base A  n=5003  (99.3%)  meanMAPQ=59.9  MAPQ0=0  minus-strand=56.9%
```

99.3% reference, 0.5% noise, at the same position, same depth, same context, same
flowcell. The chemistry, the basecaller and the aligner all work correctly here.

**4. Genuine mixed infection — no.** This is the last refuge, and the one the
invariance kills. Mixed-clone infections are normal in mosquito bloodmeals, but
the clone *ratio* is a property of the individual infection and varies freely.
Observed across ~94 unrelated mosquitoes:

```
AF  0.63 ██████████   0.64 █████████████████████████   0.65 ██████████████████████████████████
    0.66 ████████████████████   0.67 ███████   0.68 ███████   0.80 ███
```

~90 of 94 independent infections fall in a 0.05-wide band around 0.65. Unrelated
infections, acquired from different humans at different times, cannot share a
clone ratio. Whatever produces this is a property of the **assay**, not of the
parasites.

## Conclusion

The reads are real, high-quality and uniquely mapped; the clonal control is
clean; the allele fraction is invariant across independent infections. The only
mechanism consistent with all four observations is **co-amplification of a second
template at fixed stoichiometry** in bloodmeal-derived material — present in field
extractions, absent from cultured control DNA.

The correct output is therefore `NA`, not a genotype. Reporting 65% would have
been a fabricated allele frequency for a mutation that is genuinely common in
Ghana, which is precisely the kind of plausible-looking error that survives peer
review.

**Sanger sequencing was considered and is not the right instrument.** It would
tell us the consensus base in ~20 specimens, but it cannot distinguish a 65%
co-amplified template from a 65% minor clone — Sanger reads the same mixed
trace. The discriminating evidence is the clonal control and the cross-specimen
invariance, both of which are already in hand at 10,000× depth. What *would*
resolve the mechanism is a primer/BLAST audit of the `mdr1` amplicon against the
*P. falciparum* and *Anopheles* genomes to identify the second template.

## What this does not affect

The SP conclusion does not depend on any masked codon. The quintuple/sextuple
genotype hinges on **`dhps K540E` at codon 1618** — a *different* position from
the 1620 artefact — which is genuinely callable and genuinely **0/413**.
