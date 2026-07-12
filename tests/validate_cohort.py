#!/usr/bin/env python3
"""Post-run validation of a pooled-cohort result.

A green Nextflow run only means every process exited 0. It says nothing about
whether the science is right: this pipeline has already shipped, with exit code 0,
(a) a coverage table whose sample_id column contained gene names, and (b) a
15,625-row cartesian blow-up caused by a control that is re-sequenced on every
flowcell. These are the checks that would have caught both.

Everything is derived from the inputs -- there are no hardcoded cohort sizes, so
this works for 96 samples or 6,000.

    tests/validate_cohort.py --outdir ../05_results/v2 \
                             --input assets/runs.csv \
                             --samplesheet assets/samplesheet.xlsx \
                             [--baseline <old per-run 04_downstream dir>]
"""
import argparse
import collections
import csv
import pathlib
import sys

CONTROL_IDS = {"Control_HB3", "Control_Dd2", "Control_KH2", "control_KH2",
               "control_Dd2", "KH2", "PC", "NC"}
GENES = {"crt", "dhfr", "dhps", "mdr1", "k13", "csp", "msp1"}

fails = []


def check(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        fails.append(msg)


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def one(d, pat):
    hits = sorted(pathlib.Path(d).glob(pat))
    if not hits:
        sys.exit(f"MISSING: {d}/{pat}")
    return hits[0]


def base_id(s):
    return s.split("__", 1)[0]


def expected_cohort(samplesheet, runs_csv, sheet):
    """Cohort size = sample-sheet rows whose ont_multiplex_group is in --input.

    Derived, never hardcoded: the whole point is that this scales with the run.
    """
    import openpyxl
    runs = {r["run_name"] for r in read_csv(runs_csv)}
    wb = openpyxl.load_workbook(samplesheet, read_only=True, data_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    hdr = list(rows[0])
    gi, si = hdr.index("ont_multiplex_group"), hdr.index("sample_id")
    members = [(r[gi], r[si]) for r in rows[1:] if r[gi] in runs]
    return runs, members


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--sheet", default="samplesheet")
    ap.add_argument("--baseline", help="optional: an older per-run 04_downstream dir to compare scale against")
    a = ap.parse_args()

    res = pathlib.Path(a.outdir) / "04_resistance"
    runs, members = expected_cohort(a.samplesheet, a.input, a.sheet)
    n_expected = len(members)
    n_controls = sum(1 for _, s in members if s in CONTROL_IDS)
    n_field = n_expected - n_controls

    print("=" * 78)
    print(f"COHORT: {len(runs)} run(s), {n_expected} sequencing instances "
          f"({n_field} field + {n_controls} control)")
    print("=" * 78)

    # --- 1. coverage table -------------------------------------------------
    print("\n1. COVERAGE TABLE")
    cov = read_csv(one(res / "01_coverage_qc", "*coverage_by_run_sample*.csv"))
    ids = [r["sample_id"] for r in cov]
    print(f"  rows: {len(cov)}   e.g. {ids[:4]}")

    check(not (set(ids) & GENES),
          "sample_id holds sample names, not gene names")

    # Scales with the cohort: one row per sample, never more. Fewer is legitimate
    # (a barcode below --min_barcode_dir_size, e.g. an empty negative control).
    check(len(cov) <= n_expected,
          f"coverage rows <= cohort size, i.e. no merge fan-out "
          f"({len(cov)} <= {n_expected})")

    dupes = [k for k, v in collections.Counter(ids).items() if v > 1]
    check(not dupes, f"sample_id is unique in the coverage table {dupes[:3] or ''}")

    dropped = n_expected - len(cov)
    if dropped:
        print(f"  note: {dropped} instance(s) produced no coverage "
              f"(expected for an empty negative control)")

    # --- 2. frequencies ----------------------------------------------------
    print("\n2. FREQUENCIES")
    freq = read_csv(res / "04_summary/tables/summary_drug_resistance_frequencies_.csv")
    nc = [int(r["n_callable"]) for r in freq if r["n_callable"] not in ("", "NA")]
    check(bool(nc) and max(nc) > 0,
          f"at least one marker is callable (max n_callable={max(nc) if nc else 0})")

    names = {r["run_name"] for r in freq}
    check("run" not in names, f"run_name is the cohort label, not the literal 'run' -> {names}")

    totals = {int(r["n_total"]) for r in freq if r["n_total"]}
    check(totals == {n_field},
          f"n_total == field-sample count, i.e. controls excluded ({totals} vs {n_field})")

    by_gene = {r["gene"]: int(r["n_callable"])
               for r in freq if r["n_callable"] not in ("", "NA")}
    print(f"  n_callable by gene: {by_gene}")
    check(all(v <= n_field for v in by_gene.values()),
          "no gene reports more callable samples than exist")

    # --- 3. controls -------------------------------------------------------
    print("\n3. CONTROLS")
    k13 = read_csv(one(res / "02_genotype_calls", "*k13_variant_info.csv"))
    val = [r for r in k13 if r.get("k13_marker_classification") == "validated"]
    ctrl_hits = {r["sample_id"] for r in val if base_id(r["sample_id"]) in CONTROL_IDS}
    field_hits = {r["sample_id"] for r in val if base_id(r["sample_id"]) not in CONTROL_IDS}
    print(f"  k13 rows: {len(k13)};  WHO-validated calls: {len(val)}")
    for s in sorted(ctrl_hits)[:6]:
        print(f"    control {s}")

    # KH2 is the ART-R positive control: if it is in the cohort it MUST call C580Y.
    kh2_present = any(base_id(s) == "KH2" for _, s in members)
    if kh2_present:
        check(any(r["aa_mut"] == "C580Y" and base_id(r["sample_id"]) == "KH2" for r in val),
              "positive control KH2 calls k13 C580Y")
    else:
        print("  (no KH2 control in this cohort - skipping the positive-control check)")

    if field_hits:
        print(f"  {len(field_hits)} FIELD sample(s) carry a validated ART-R marker: "
              f"{sorted(field_hits)[:5]}")

    # --- 4. optional scale comparison against a known-good per-run baseline --
    if a.baseline:
        print("\n4. SCALE vs BASELINE")
        old = read_csv(pathlib.Path(a.baseline) /
                       "04_summary/tables/summary_drug_resistance_frequencies_.csv")
        old_gene = {r["gene"]: int(r["n_callable"])
                    for r in old if r["n_callable"] not in ("", "NA")}
        print(f"  baseline n_callable: {old_gene}")
        print("  NOTE: a clean integer multiple across runs is NOT proof of correctness -")
        print("        it is exactly what duplicated input data looks like.")
        for g in sorted(set(old_gene) & set(by_gene)):
            if old_gene[g]:
                print(f"    {g:5} baseline={old_gene[g]:>4}  cohort={by_gene[g]:>4}  "
                      f"x{by_gene[g] / old_gene[g]:.1f}")

    print()
    print("=" * 78)
    print("RESULT:", "ALL CHECKS PASSED" if not fails else f"{len(fails)} CHECK(S) FAILED")
    for f in fails:
        print("   FAILED:", f)
    print("=" * 78)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
