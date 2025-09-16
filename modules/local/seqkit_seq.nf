process SEQKIT_SEQ {
    tag "${meta.alias}"
    label 'process_low'

    conda "bioconda::seqkit=2.9.0"
    container "quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0"

    input:
    tuple val(meta), path(fastq), val(q_score)

    output:
    tuple val(meta), path("*_filtered_Q${q_score}.fastq.gz")    , emit: filtered_fastq
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    seqkit \\
        seq \\
        --threads ${task.cpus} \\
        -Q ${q_score} \\
        ${args} \\
        ${fastq} \\
        | gzip > ${prefix}_filtered_Q${q_score}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    touch ${prefix}_filtered_${q_score}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
