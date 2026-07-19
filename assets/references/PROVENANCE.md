# Reference amplicon panel — provenance & loci classification

*Plasmodium falciparum* **3D7** reference coding sequences (CDS) for the targeted
amplicons. FASTA files live in `ref_genomes/target_gene_cds_seqs/`; checksums in
that directory's `SHA256SUMS`. The active panel is defined by
`manifests/reference_manifest.csv` (7 loci).

## Loci — resistance-associated vs. antigen/diversity

The FASTA record IDs self-document the panel: `3D7_DR1_*` = drug-resistance panel,
`3D7_AG1_*` = antigen/diversity panel. **Not all seven amplicons are drug-resistance
genes** — this distinction is load-bearing for interpretation.

| Locus | In manifest | Class | Role |
|---|---|---|---|
| `crt` (pfcrt) | ✓ | **Drug-resistance** | Chloroquine (K76T) |
| `dhfr` (pfdhfr) | ✓ | **Drug-resistance** | Pyrimethamine (N51I, C59R, S108N, I164L) |
| `dhps` (pfdhps) | ✓ | **Drug-resistance** | Sulfadoxine (S436A, A437G, K540E, A581G, A613S) |
| `mdr1` (pfmdr1) | ✓ | **Drug-resistance** | Partner-drug/chloroquine (N86Y, Y184F) |
| `k13` (pfkelch13) | ✓ | **Drug-resistance** | Artemisinin partial resistance (WHO-validated markers) |
| `csp` (pfcsp) | ✓ | Antigen / vaccine | Circumsporozoite protein — sporozoite antigen (RTS,S/R21). **Not a drug-resistance gene.** |
| `msp1` (pfmsp1) | ✓ | Antigen / diversity | Merozoite surface protein 1 — genotyping/diversity & complexity-of-infection. **Not a drug-resistance gene.** |
| `ama1`, `eba175`, `msp2` | ✗ (present, unused) | Antigen / diversity | Retained for reference; **not** in the active manifest. |

Marker-level classification, QC status (accepted / reference_caveat / rejected_artifact /
requires_validation) and interpretation scope are catalogued in
`../../docs/objective3_marker_catalog.csv` and applied via `assets/resources/`
(`DR_variant_info_v2.xlsx`, k13 WHO catalogue). Note `pfdhps` A437G carries a
reference-coding caveat (3D7 itself carries 437G).

## Source, build & orientation

- Organism/strain: *P. falciparum* 3D7 (reference).
- Sequence type: gene CDS (coding sequence), forward/sense orientation as used by
  `minimap2 -ax map-ont --MD` and Clair3; amplicon codon positions are resolved by
  `assets/resources/amplicon_gene_positions.csv` + `genetic_code.csv`.
- `.fai` indexes are generated at runtime (`SAMTOOLS_FAIDX`), not committed.
- **Exact PlasmoDB/GeneDB release to be recorded by the maintainer** — these are
  standard 3D7 CDS references; pin the release used when publishing.

## Integrity & override

Verify with `sha256sum -c SHA256SUMS` in the FASTA directory. To use a different
panel, pass `--reference_manifest <csv>` pointing at your own FASTAs; a cut-down
panel is valid for variant calling only (`--drug_resistance false`) — the pooled
resistance analysis needs the full 7-locus panel.
