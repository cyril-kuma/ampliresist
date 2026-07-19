# Geographic assets — provenance

Boundary polygons for the stage-18 geographic figure
(`bin/07_build_roadmap_figures.py` → `figure1`, DHFR haplotype map).

| File | Content | Features |
|---|---|---|
| `ghana_ADM0.geojson` | Ghana national boundary (ADM0) | 1 |
| `ghana_ADM1.geojson` | Ghana regions (ADM1) | 16 |

- **Source:** geoBoundaries — Open administrative boundaries, country `GHA`, levels ADM0 & ADM1. https://www.geoboundaries.org
- **Citation:** Runfola, D. et al. (2020). geoBoundaries: A global database of political administrative boundaries. *PLoS ONE* 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866
- **Licence:** **CC BY 4.0** — redistribution permitted with attribution (satisfied by this file). https://creativecommons.org/licenses/by/4.0/
- **Redistribution status:** permitted; these files are vendored into the pipeline.
- **CRS:** `urn:ogc:def:crs:OGC:1.3:CRS84` — equivalent to WGS84 / EPSG:4326 (longitude, latitude, decimal degrees).
- **Property schema (geoBoundaries):** `shapeName`, `shapeISO`, `shapeID`, `shapeGroup`, `shapeType`.
- **Preprocessing:** none — files are byte-for-byte as obtained (verify with `SHA256SUMS`).
- **Vendored:** 2026-07-19, from `02_plasmodium_qpcr/data/geo/` (which is *not* a runtime dependency of ampliresist; these copies make ampliresist self-contained).

## Overriding

The default resolves to this directory (`--plot_geo_dir` → `${projectDir}/assets/geo`).
A user-supplied directory may be passed with `--plot_geo_dir <dir>`; it must contain
`ghana_ADM0.geojson` and `ghana_ADM1.geojson` in a CRS84/EPSG:4326 GeoJSON, with the
same polygon geometry structure. Preflight validation checks both files exist and parse.
