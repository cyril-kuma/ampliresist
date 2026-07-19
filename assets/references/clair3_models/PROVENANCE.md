# Clair3 model — provenance

Variant calling (`modules/local/09a_clair3.nf`, `CLAIR3`) uses the model directory
`r1041_e82_400bps_sup_v430`, supplied to the container as a staged input channel
(the pinned `hkubal/clair3` image ships only R9.4.1 models, so the R10.4.1 model
must be provided separately).

| Field | Value |
|---|---|
| Model | `r1041_e82_400bps_sup_v430` |
| Chemistry / pore | R10.4.1, E8.2, 400 bps translocation |
| Basecalling | SUP (super-accurate), model version **v4.3.0** |
| Matches data basecalled with | `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` |
| Files | `full_alignment.data-00000-of-00001`, `full_alignment.index`, `pileup.data-00000-of-00001`, `pileup.index` |
| Size | ~78 MB |
| Integrity | verify against the SHA-256 sums below |
| Clair3 image | `docker.io/hkubal/clair3@sha256:57cf5d20f2ee39c1b91493ad1fb5c1b9fa838691efce818c3139caa5e6c6b974` (digest-pinned) |

## Scientific compatibility

The model **must** match the basecaller chemistry and model version of the input
reads. This model is for R10.4.1/E8.2/400bps SUP v4.3.0 data. Using a mismatched
model (e.g. the R9.4.1 models bundled in the Clair3 image) silently degrades calls.
`main.nf` enforces that `--clair3_model` is set whenever `--variant_caller clair3`.

## Authoritative source & reproducible acquisition

Clair3 R10.4.1 models are published by Oxford Nanopore via **Rerio** and mirrored on
the Clair3 model zoo:

```bash
# Option A — ONT Rerio (authoritative)
git clone https://github.com/nanoporetech/rerio
python3 rerio/download_model.py --clair3 rerio/clair3_models/r1041_e82_400bps_sup_v430_model
# Option B — Clair3 model zoo mirror
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v430.tar.gz
tar -xzf r1041_e82_400bps_sup_v430.tar.gz -C assets/references/clair3_models/
```

After obtaining the model, verify it before use (these sums are tracked here so a
fresh clone can validate a freshly downloaded model — the model binaries
themselves are git-ignored):

```
cad82f0fac02a558e5690da45e46796dcb33a4b2725554404a65f9c94bc8bf0f  full_alignment.data-00000-of-00001
e9b7b0b79976e52970406484b4e88cb21114a57a7dd14ea08309b5eb051d3d32  full_alignment.index
177c921e0dd34e3f3780c76042ddca146ca37b767721d3074fac6c66a1489bc1  pileup.data-00000-of-00001
6d2a0f114a4a4018c25031d190def1d7a5f0371b2d42f349b36ea3ecd526b647  pileup.index
```

Verify with `cd r1041_e82_400bps_sup_v430 && sha256sum -c ../PROVENANCE.sha256`
(or check each hash manually). To use your own model directory (e.g. a different
chemistry), pass `--clair3_model <dir>`; preflight checks the four files exist.

**Redistribution:** the model is bundled here for convenience; confirm the ONT
model licence permits redistribution before publishing the repository, otherwise
ship only this provenance file + the download commands above and keep the model
out of version control.
