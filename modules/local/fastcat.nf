process FASTCAT {
    tag "$meta.alias"
    label 'process_low'

    conda "nanoporetech::fastcat=0.20.0"
    container "docker://ontresearch/wf-common:shad28e55140f75a68f59bbecc74e880aeab16ab158"

    input:
    tuple val(meta), path(fastq_input)

    output:
    tuple val(meta), path("*.fastq.gz")     , emit: concat_fastq
    tuple val(meta), path("*.tsv")          , emit: stats
    tuple val(meta), path("*histograms")    , emit: histogram
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    fastcat \\
        -s ${prefix} \\
        -f ${prefix}_per_file_stats.tsv \\
        --histograms ${prefix}-histograms \\
        $args \\
        $fastq_input | \\
        bgzip > ${prefix}.fastq.gz
    
    
    if [ "\$(awk 'NR==1{for (i=1; i<=NF; i++) {ix[\$i] = i}} NR>1 {c+=\$ix["n_seqs"]} END{print c}' ${prefix}_per_file_stats.tsv)" = "0" ]; then
        mv ${prefix}.fastq.gz ${prefix}_empty.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcat: \$(fastcat --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    touch ${prefix}.fastq.gz
    touch ${prefix}_per_file_stats.tsv
    mkdir ${prefix}-histograms

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastcat: \$(fastcat --version)
    END_VERSIONS
    """
}
