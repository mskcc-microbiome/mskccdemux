/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                         } from '../modules/nf-core/fastqc/main'
include { CAT_FASTQ as MERGE_RUNS        } from '../modules/nf-core/cat/fastq/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { SEQKIT_STATS                   } from '../modules/nf-core/seqkit/stats/main'
include { add_demultiplex_info as adi_f  } from '../modules/local/add_demultiplex_info/main'
include { add_demultiplex_info as adi_r  } from '../modules/local/add_demultiplex_info/main'
include { tag_empty as tag_empty_f       } from '../modules/local/tag_empty/main'
include { tag_empty as tag_empty_r       } from '../modules/local/tag_empty/main'
include { make_map                       } from '../modules/local/make_map/main'
include { demultiplex as demux_f         } from '../modules/local/demux/main'
include { demultiplex as demux_r         } from '../modules/local/demux/main'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../subworkflows/local/utils_nfcore_mskccdemux_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process remove_primers {
    tag 'remove_primers'
    input:
    tuple val(meta), path(readsF), path(readsR)

    output:
    path 'reads1.fastq', emit: reads1
    path 'reads2.fastq', emit: reads2
    path 'barcodes.fastq', emit: barcodes
    path 'primer_removal.log'
    container 'ghcr.io/vdblab/biopython:1.70a'
    cpus 1
    memory '12 GB'
    script:
    def primer_f = meta.primer_f
    def primer_r = meta.primer_r
    """
    strip_addons3_py3.py \\
      -fw_primer $primer_f -rev_primer $primer_r -remove_bar_primer \\
       ${readsF} ${readsR} \\
    >> primer_removal.log 2>&1
    """
}

process guess_encoding {
    tag 'guess_encoding'
    input:
    path reads_fq
    output:
      path 'encoding.txt' , emit: encoding
    container 'ghcr.io/vdblab/biopython:1.70a'
    cpus 1
    script:
    """
    guess-encoding.py ${reads_fq} encoding.txt 2>> guess_encoding.log
    """
}

process rename_for_multiqc{
    """ Multiqc requires the _mqc string in the file name for auto detection

    """
    input:
    path seqkit
    output:
    path "demultiplex_seqkit_stats_mqc.out" , emit: mqc

    script:
    """
    cp $seqkit demultiplex_seqkit_stats_mqc.out
    """
}

workflow MSKCCDEMUX {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    def outdir = file(params.outdir)
    //LOG(ch_fastqs)
    //
    // MODULE: Run FastQC
    // Get the unique fastqs from the sample sheet and execute pool/multiplex level fastqc
//    ch_samplesheet.view()
    fqc_inputs = ch_samplesheet.flatMap{
	meta, reads ->
	reads.withIndex().collect { read, idx ->
            def new_meta = ["id": params.poolid, "run_accession": meta["run_accession"], "single_end":  meta["single_end"], "read_idx": idx + 1]
            [new_meta, read]
	}}.unique()

    FASTQC (
        fqc_inputs
    )
    // we can't concatenate files if there is not a second run, so we branch
    // here to separate them out, and mix back in after for efficiency
    ch_reads_grouped = fqc_inputs
        .map { meta, reads -> [[id: "pool", single_end: meta["single_end"]], reads] }
        .groupTuple ()
	.branch { meta, reads ->
	    cat: (meta.single_end && reads.size() > 1) || (!meta.single_end && reads.size() > 2)
	    skip: true
	}
    MERGE_RUNS ( ch_reads_grouped.cat)
    ch_reads_runmerged = MERGE_RUNS.out.reads
	.mix(ch_reads_grouped.skip) // dont need extra flatten cause this pipeline will always run on either a merged set of libraries or single library, not mixing at the sample level between both
   // ch_versions = ch_versions.mix(MERGE_RUNS.out.versions)


    ch_samplesheet_unique = ch_samplesheet
	.unique { meta, reads ->
            meta.rawid  // if we have multiple libraries, we need to toss  duplicated  entries in the sample sheet
	}

    def merged_reads = ch_reads_runmerged.first()
    persample_inputs = ch_samplesheet_unique
	.combine(ch_reads_runmerged.first())
	.map{ meta, reads, newmeta, newreads ->
	meta = ["id": params.poolid, "primer_f": meta["primer_f"], "primer_r": meta["primer_r"]]
	[meta, newreads[0], newreads[1]]
	}.unique()

    persample_inputs

    remove_primers (
	persample_inputs
    )


    guess_encoding (
	remove_primers.out[0]
    )
    /////////////////////////////////////
    map1_path = outdir.resolve( params.poolid + ".map.1")
    map2_path = outdir.resolve( params.poolid + ".map.2")
    header = Channel.value("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription")
    header
    .concat(
	ch_samplesheet_unique
	.map{ meta, reads ->
            "${meta.id}\t${meta.barcode_f}\t${meta.primer_f}\t${meta.primer_r}\t${meta.rawid}"
	    }
    )
    .collectFile(name: map1_path, newLine: true, sort:false)
    header
    .concat(
	ch_samplesheet_unique
	.map{ meta, reads ->
            "${meta.id}\t${meta.barcode_r}\t${meta.primer_r}\t${meta.primer_f}\t${meta.rawid}"
	    }
    )
    .collectFile(name: map2_path, newLine: true, sort:false)

    sampledir = file("${params.outdir}/sampleids/")
    sampledir.mkdir()

    ch_samplesheet_unique
    .map{ meta, reads ->
            "${meta.id}"
	    }
    .concat(Channel.value("Unassigned"))
    .collectFile{ x ->
        [ sampledir.resolve("${x}.sample"), x ]
	}
    .set{ samplefiles }

    /////////////////////////////////////
    //make_map (
    // 	params.input
    //)

    adi_f (
	remove_primers.out.reads1,
	map1_path,
	remove_primers.out.barcodes,
	guess_encoding.out.encoding,
	1
    )
    adi_r  (
	remove_primers.out.reads2,
	map2_path,
	remove_primers.out.barcodes,
	guess_encoding.out.encoding,
	2
    )

    demux_f(
	adi_f.out.seqsfq.toList(),
	samplefiles.flatten(),
	1
    )
    demux_r(
	adi_r.out.seqsfq.toList(),
	samplefiles.flatten(),
	2
    )

    seqkit_input = demux_f.out.samplefq.mix(
	 demux_r.out.samplefq
    )
	.collect(sort:true)
	.map{ x ->
	    meta = ["id": "${params.poolid}"]
	    [meta, x]
	}
    SEQKIT_STATS (
	seqkit_input

    )
    tag_empty_f(
	demux_f.out.samplefq
    )
    tag_empty_r(
	demux_r.out.samplefq
    )

    //SEQKIT_STATS.out.stats.collectFile( name: 'demultiplex_seqkit_stats_mqc.out')
    // SEQKIT_STATS.out.stats
    // 	.map{x -> [x]}
    // 	.collectFile(
    // 	    name: 'demultiplex_seqkit_stats_mqc.out',
    // 	    storeDir: "${params.outdir}/")
    rename_for_multiqc(
	SEQKIT_STATS.out.stats.collect{it[1]}
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}).mix(rename_for_multiqc.out.mqc)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first()).mix(SEQKIT_STATS.out.versions)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'mskccdemux_software_'  + 'mqc_'  + 'versions.yml',
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
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
