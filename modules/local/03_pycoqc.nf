// Run-level sequencing QC from the ONT sequencing_summary file.

process PYCOQC {
    tag "${meta.run_name}"
    label 'process_medium'
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"

    input:
    tuple val(meta), path(sequencing_summary)

    output:
    tuple val(meta), path("*.html"), path("*.json"), emit: report
    path "versions.yml",                             emit: versions

    script:
    """
    pycoQC -f ${sequencing_summary} \\
        -o ${meta.run_name}_pycoqc.html \\
        -j ${meta.run_name}_pycoqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version | sed 's/^pycoQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.run_name}_pycoqc.html ${meta.run_name}_pycoqc.json versions.yml
    """
}
