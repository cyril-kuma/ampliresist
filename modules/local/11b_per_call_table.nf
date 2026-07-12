// Flatten every variant call into one tidy table: per_call.tsv
//
// The R analysis stages collapse each call to a binary mutant/wild-type and then
// discard GT, AF and DP. Everything needed to QC the calls therefore vanishes
// before the analysis can see it -- which is how a systematic ~0.5-AF artifact
// (k13 C580Y, mdr1 Y184F, dhps K540N/A581G) survived as far as the results.
//
// The artifact signature is only visible in this table: a real variant is
// GT=1 (haploid) at AF ~0.9 and varies between specimens; a systematic error is
// GT=0/1 at AF ~0.5 in nearly EVERY specimen with a very tight cross-sample sd.
// Emit it as a first-class output so the QC is auditable, not accidental.

process PER_CALL_TABLE {
    tag "${run_name}"
    label 'process_low'
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(run_name), path(vcfs)

    output:
    tuple val(run_name), path("per_call.tsv"), emit: table

    script:
    """
    #!/usr/bin/env python3
    import glob, os, csv

    rows = []
    for path in sorted(glob.glob("*.vcf")):
        stem = os.path.basename(path)[:-4]              # <run>_<barcode>_<gene>
        gene = stem.rsplit("_", 1)[1]
        barcode = stem.rsplit("_", 2)[1]
        run = stem.rsplit("_barcode", 1)[0]

        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\\n").split("\\t")
                if len(f) < 10:
                    continue
                keys = f[8].split(":")
                vals = f[9].split(":")
                d = dict(zip(keys, vals))
                ad = d.get("AD", "")
                ref_dp, alt_dp = (ad.split(",") + ["", ""])[:2] if ad else ("", "")
                rows.append({
                    "run_name": run,
                    "barcode": barcode,
                    "gene": gene,
                    "pos": f[1],
                    "ref": f[3],
                    "alt": f[4],
                    "qual": f[5],
                    "filter": f[6],
                    "GT": d.get("GT", ""),
                    "DP": d.get("DP", ""),
                    "AF": d.get("AF", ""),
                    "ref_depth": ref_dp,
                    "alt_depth": alt_dp,
                })

    cols = ["run_name","barcode","gene","pos","ref","alt","qual","filter",
            "GT","DP","AF","ref_depth","alt_depth"]
    with open("per_call.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\\t")
        w.writeheader()
        w.writerows(rows)

    print(f"per_call.tsv: {len(rows)} calls from {len(glob.glob('*.vcf'))} VCFs")
    """

    stub:
    """
    printf 'run_name\\tbarcode\\tgene\\tpos\\tref\\talt\\tqual\\tfilter\\tGT\\tDP\\tAF\\tref_depth\\talt_depth\\n' > per_call.tsv
    """
}
