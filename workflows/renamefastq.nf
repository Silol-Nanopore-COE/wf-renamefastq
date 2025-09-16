/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if (params.sample_sheet) { ch_samplesheet = Channel.fromPath(file(params.sample_sheet), checkIfExists: true) } else { ch_samplesheet = Channel.empty() }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check'
include { DORADO_DEMULTIPLEX     } from '../modules/local/dorado_demux'
include { FASTCAT                } from '../modules/local/fastcat'
include { SEQKIT_SEQ             } from '../modules/local/seqkit_seq'
include { SEQKIT_STATS           } from '../modules/local/seqkit_stats'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
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
    
    // 
    // Check if the input sample sheet is valid
    //
    SAMPLESHEET_CHECK (
        ch_samplesheet
    )
    ch_valid_samplesheet = SAMPLESHEET_CHECK.out.checked_sheet
    ch_alias_for_map     = ch_valid_samplesheet.splitCsv(sep: ',', skip: 1) // skip header row
    ch_versions          = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    // 
    // Prepare the input FASTQ channel before renaming
    // 

    ch_not_demultiplexed_fastq = Channel.empty()
    ch_demultiplexed_fastq     = Channel.empty()
    
    if (!params.demultiplex) {
        ch_not_demultiplexed_fastq = ch_not_demultiplexed_fastq.mix(ch_fastq_input.singleFile)
        ch_not_demultiplexed_fastq = ch_not_demultiplexed_fastq.mix(ch_fastq_input.topLevelDir)
        
        // Map alias if the sample sheet is specified
        if (params.sample_sheet) {
            ch_fastq_input.dirWithSubdirs
                .map { meta, fastq -> [ meta.barcode, meta, fastq ] }
                .join(ch_alias_for_map, by: [0])
                .map { barcode, meta, fastq, alias ->
                    def new_meta = ["alias": alias, "barcode": meta.barcode ]
                        return [ new_meta, fastq ]
                }
                .set { ch_fastq_with_new_meta } 
            
            ch_not_demultiplexed_fastq = ch_not_demultiplexed_fastq.mix(ch_fastq_with_new_meta)
        } else {
            ch_not_demultiplexed_fastq = ch_not_demultiplexed_fastq.mix(ch_fastq_input.dirWithSubdirs)
        }
    } else {
        ch_fastq_for_demux  = Channel.empty()

        ch_fastq_for_demux
            .mix(ch_fastq_input.singleFile.map { meta, fastq -> [ meta, fastq, false ] })
            .mix(ch_fastq_input.topLevelDir.map { meta, fastq -> [ meta, fastq, true ]})
            .map { meta, fastq, is_dir -> [ meta, fastq, is_dir, params.kit_name ] }
            .set { ch_input_for_demux }
            
        //
        // MODULE: DORADO, Demultiplex FASTQ file
        //  
        DORADO_DEMULTIPLEX (
            ch_input_for_demux 
        )
        ch_demux_fastq = DORADO_DEMULTIPLEX.out.demux_fastq
        ch_versions    = ch_versions.mix(DORADO_DEMULTIPLEX.out.versions)

        if (params.sample_sheet) {
            // Separate demultiplexed and unclassified reads to different channels
            ch_demux_fastq
                .map { meta, fastq -> 
                    if (fastq instanceof List) {
                        return [ meta, fastq ]
                    } else {
                        return [ meta, [ fastq ] ]
                    }
                }
                .transpose(by: [1])
                .branch { meta, fastq ->
                    def fq_name = fastq.simpleName
                    demultiplexed: fq_name =~ /barcode/
                        def barcode = fq_name.split('_')[-1]
                        def new_meta = meta + [ "barcode": barcode, "demux_name": fq_name ]
                        return [ new_meta, fastq ]
                    unclassified: fq_name =~ /unclassified/
                        def new_meta = [ "alias": "unclassified", "demux_name": fq_name ]
                        return [ new_meta, fastq ]
                }
                .set { ch_reads }
        
            // Map alias to individual FASTQ file based on barcode
            ch_reads.demultiplexed
                .map { meta, fastq -> [ meta.barcode, meta, fastq ] }
                .join(ch_alias_for_map, by: [0])
                .map { barcode, meta, fastq, alias ->
                    def new_meta = [ "alias": alias, "demux_name": meta.demux_name ]
                    return [ new_meta, fastq ]
                }
                .mix(ch_reads.unclassified)
                .set { ch_demultiplexed_fastq }
        } else {
            // Create a new meta map for demultiplexed reads
            ch_demux_fastq
                .map { meta, fastq -> 
                    if (fastq instanceof List) {
                        return [ meta, fastq ]
                    } else {
                        return [ meta, [ fastq ] ]
                    }
                }
                .transpose(by: [1])
                .branch { meta, fastq ->
                    def fq_name = fastq.simpleName
                    demultiplexed: fq_name =~ /barcode/
                        def alias    = fastq.simpleName.split('_')[-1]
                        def new_meta = [ "alias": alias, "demux_name": fastq.simpleName ]
                        return [ new_meta, fastq ]
                    unclassified: fq_name =~ /unclassified/
                        def new_meta = [ "alias": "unclassified", "demux_name": fastq.simpleName ]
                        return [ new_meta, fastq ]
                }
                .set { ch_reads }

            ch_demultiplexed_fastq = ch_reads.demultiplexed.mix(ch_reads.unclassified)
        }
    }

    // Create the input channel for Fastcat
    ch_fastq_for_fastcat = Channel.empty()
    ch_fastq_for_fastcat = ch_fastq_for_fastcat.mix(ch_not_demultiplexed_fastq)
    ch_fastq_for_fastcat = ch_fastq_for_fastcat.mix(ch_demultiplexed_fastq)

    // 
    // MODULE: FASTCAT, Concatenate/Rename FASTQ files
    // 
    FASTCAT (
        ch_fastq_for_fastcat
    )
    ch_renamed_fastq  = FASTCAT.out.concat_fastq
    ch_versions       = ch_versions.mix(FASTCAT.out.versions.first())

    // Create the input channel for Seqkit stats
    ch_fastq_for_seqkit = Channel.empty()

    if (params.quality_filter) {
        ch_renamed_fastq
            .map { meta, fastq -> [ meta, fastq, params.q_score ] }
            .set { ch_fastq_for_filter }
        
        // 
        // MODULE: SEQKIT, Filtering reads by Q-score
        // 
        SEQKIT_SEQ (
            ch_fastq_for_filter
        )
        ch_filtered_fastq  = SEQKIT_SEQ.out.filtered_fastq
        ch_versions        = ch_versions.mix(SEQKIT_SEQ.out.versions.first())

        ch_fastq_for_seqkit = ch_fastq_for_seqkit.mix(ch_filtered_fastq)
    } else {
        ch_fastq_for_seqkit = ch_fastq_for_seqkit.mix(ch_renamed_fastq)
    }

    ch_fastq_for_seqkit
        .map { meta, fastq -> fastq }
        .toList()
        .map { list -> 
            def new_meta = ["alias": "all_fastq"]
                return [ new_meta, list ]
        }
        .set { ch_list_fastq }

    //
    // MODULE: SEQKIT, Computing statistics of FASTQ files
    //  
    SEQKIT_STATS (
        ch_list_fastq
    )
    ch_fastq_stats  = SEQKIT_STATS.out.stats
    ch_versions     = ch_versions.mix(SEQKIT_STATS.out.versions)

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
