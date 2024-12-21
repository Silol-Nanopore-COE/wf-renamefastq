/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check'
include { FASTCAT                } from '../modules/local/fastcat'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_renamefastq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RENAMEFASTQ {

    take:
    fastq // channel: path to FASTQ files from --fastq -> [ InputType, [ path/to/fastq/file or path/to/fastq/directory ] ]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // 
    // Create a meta map for three input FASTQ cases
    // 
    fastq
        .transpose(by: [1])
        .branch { 
            inputType, fastq ->
            singleFile: inputType == "SingleFile"
                def meta = ["alias": params.sample?: fastq.simpleName ]
                return [ meta, fastq ]
            topLevelDir: inputType == "TopLevelDir"
                def meta = ["alias": params.sample?: fastq.baseName ]
                return [ meta, fastq ]
            dirWithSubdirs: inputType == "DirWithSubDirs"
                def meta = ["alias": fastq.baseName, "barcode": fastq.baseName ]
                return [ meta, fastq ]
        }
        .set { ch_fastq_input }
    
    // Create a new meta map for the directory with subdirectory(s) type
    // If the sample sheet is supplied with --sample_sheet
    // 
    ch_dir_with_subdirs_new_meta = Channel.empty()
    if (params.sample_sheet) {
        // 
        // Check if the input sample sheet is valid
        // 
        ch_samplesheet = Channel.fromPath(file(params.sample_sheet, checkIfExists: true))
        SAMPLESHEET_CHECK {
            ch_samplesheet
        }
        ch_valid_samplesheet = SAMPLESHEET_CHECK.out.checked_sheet
        ch_versions          = ch_versions.mix(SAMPLESHEET_CHECK.out.versions.first())
        
        // Create a channel from a valid sample sheet
        ch_valid_samplesheet
            .splitCsv(sep: ',', skip: 1)
            .set { ch_samplesheet_for_joining }
        
        ch_fastq_input.dirWithSubdirs
            .map { meta, fastq -> [ meta.barcode, meta, fastq ] }
            .join(ch_samplesheet_for_joining, by: [0])
            .map { barcode, meta, fastq, alias ->
                def new_meta = [ "alias": alias, "barcode": barcode ]
                return [ new_meta, fastq ]
            }
            .set { ch_dir_with_subdirs_new_meta }
    } else {
        ch_dir_with_subdirs_new_meta = ch_fastq_input.dirWithSubdirs
    }

    // Combine channels of three input FASTQ cases 
    ch_fastq_with_meta = Channel.empty()
    ch_fastq_with_meta = ch_fastq_with_meta.mix(ch_fastq_input.singleFile)
    ch_fastq_with_meta = ch_fastq_with_meta.mix(ch_fastq_input.topLevelDir)
    ch_fastq_with_meta = ch_fastq_with_meta.mix(ch_dir_with_subdirs_new_meta)
    
    // MODULE: FASTCAT
    // Concatenating multiple FASTQ files into a single FASTQ file (if required)
    // Computing FASTQ statistics
    FASTCAT {
        ch_fastq_with_meta 
    }
    ch_fastcat_stats = FASTCAT.out.concat_fastq
    ch_versions      = ch_versions.mix(FASTCAT.out.versions.first())
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
