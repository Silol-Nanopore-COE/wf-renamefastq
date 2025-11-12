process DORADO_DEMULTIPLEX {
    tag "$meta.alias"
    label 'process_medium'

    container "piroonj/dorado:hla-v1.2.0"

    input:
    tuple val(meta), path(fastq_input), val(is_dir), val(kit_name)

    output:
    tuple val(meta), path("dorado_demux_fastq/*.fastq") , emit: demux_fastq
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    def dir_input = is_dir ? "-r": ""
    def VERSION = '1.1.1+e72f1492+dirty' // WARN: Version information not provided by tool on CLI. Please update this string when changing container versions.
    """
    dorado \\
        demux \\
        $args \\
        $dir_input \\
        -t $task.cpus \\
        --kit-name $kit_name \\
        --output-dir dorado_demux_fastq \\
        $fastq_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    mkdir dorado_demux_fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: $VERSION
    END_VERSIONS
    """
}
