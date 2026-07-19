// Run-level sequencing QC from the ONT sequencing_summary file.

process PYCOQC {
    tag "${meta.run_name}"
    label 'process_medium'
    container "quay.io/biocontainers/pycoqc@sha256:ea0a084751a0b48b5ffe90e9d3adfa8f57473709a1b0a95c9cb38d434ee3a9a2"  // pycoqc 2.5.2

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
