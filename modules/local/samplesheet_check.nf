process SAMPLESHEET_CHECK {
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "docker://ontresearch/wf-common:shad28e55140f75a68f59bbecc74e880aeab16ab158"

    input:
    path(samplesheet)

    output:
    path(samplesheet)     , emit: checked_sheet
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samplesheet_check.py --sample_sheet ${samplesheet} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
