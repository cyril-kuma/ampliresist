// Concatenate each barcode's per-chunk FASTQs into one file per barcode.
//
// Barcode directories smaller than --min_barcode_dir_size (MB) are skipped: on a
// 96-barcode flowcell the unused/failed barcodes still produce a directory with a
// handful of stray reads, and calling variants on those is meaningless.

process SORT_FASTQS {
    tag "${run_name}"
    label 'process_low'

    input:
    tuple val(run_name), path(sequencing_dir, stageAs: 'seqdir'), path(sequencing_summary)

    output:
    tuple val(run_name), path("*.fastq.gz"), emit: fastqs
    path "versions.yml",                     emit: versions

    script:
    """
    threshold=${params.min_barcode_dir_size}
    fastq_pass_dir="seqdir/fastq_pass"

    if [ ! -d "\$fastq_pass_dir" ]; then
        echo "ERROR: fastq_pass directory not found in the sequencing dir for run ${run_name}" >&2
        exit 1
    fi

    found=false
    for dir in "\$fastq_pass_dir"/barcode*; do
        [ -d "\$dir" ] || continue
        found=true
        barcode=\$(basename "\$dir")

        shopt -s nullglob
        fastq_files=( "\$dir"/*.fastq.gz )
        shopt -u nullglob

        if [ \${#fastq_files[@]} -eq 0 ]; then
            echo "WARN: skipping \$barcode - no fastq.gz files" >&2
            continue
        fi

        disk_usage=\$(du -sm --dereference "\$dir" | awk '{ print \$1 }')
        if [ "\$disk_usage" -lt "\$threshold" ]; then
            echo "WARN: skipping \$barcode - \${disk_usage}MB < \${threshold}MB threshold" >&2
            continue
        fi

        cat "\${fastq_files[@]}" > "${run_name}_\${barcode}.fastq.gz"
    done

    if [ "\$found" = false ]; then
        echo "ERROR: no barcode directories found in \$fastq_pass_dir" >&2
        exit 1
    fi

    # A run where every barcode was filtered out is a configuration error, not an
    # empty result. Without this the process would emit nothing, every downstream
    # channel would be empty, and the pipeline would report SUCCESS having
    # processed zero samples.
    shopt -s nullglob
    produced=( *.fastq.gz )
    shopt -u nullglob
    if [ \${#produced[@]} -eq 0 ]; then
        echo "ERROR: no barcode in run ${run_name} met --min_barcode_dir_size=\${threshold}MB." >&2
        echo "       Every barcode directory was smaller than the threshold. Lower it and re-run." >&2
        exit 1
    fi
    echo "Merged \${#produced[@]} barcodes for run ${run_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${run_name}_barcode01.fastq.gz
    touch versions.yml
    """
}
