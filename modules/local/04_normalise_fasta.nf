// Rewrite each reference FASTA with consistent line wrapping.
//
// The references are hand-curated amplicon/CDS sequences whose line lengths vary;
// samtools faidx requires uniform line lengths within a record, so normalise
// before indexing.

process NORMALISE_FASTA {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/biocontainers/biopython:1.78"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("normalised/${fasta.name}"), emit: fasta
    path "versions.yml",                               emit: versions

    script:
    """
    mkdir -p normalised
    python3 - <<'PY'
from Bio import SeqIO
records = SeqIO.parse("${fasta}", "fasta")
SeqIO.write(records, "normalised/${fasta.name}", "fasta")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p normalised && touch normalised/${fasta.name}
    touch versions.yml
    """
}
