/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                         } from '../modules/nf-core/fastqc/main'
include { CAT_FASTQ as MERGE_RUNS        } from '../modules/nf-core/cat/fastq/main'
include { CUTADAPT as CUTADAPT_READTHROUGH } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { SEQKIT_STATS                   } from '../modules/nf-core/seqkit/stats/main'
include { add_demultiplex_info as adi_f  } from '../modules/local/add_demultiplex_info/main'
include { add_demultiplex_info as adi_r  } from '../modules/local/add_demultiplex_info/main'
include { MAKE_MANIFEST                  } from '../modules/local/make_manifest/main'
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

// code from https://github.com/nf-core/ampliseq/
// Complement table taken from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
def makeComplement(seq) {
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]
    def comp = seq.toUpperCase().collect { base -> complements[ base ] ?: 'X' }.join()
    return comp
}

process remove_primers {
    tag 'remove_primers'
    container 'ghcr.io/vdblab/biopython:1.70a'
    cpus 1
    memory '12 GB'
    input:
    tuple val(meta), path(readsF), path(readsR)

    output:
    tuple val(meta), path('reads1.fastq'), emit: reads1
    tuple val(meta), path('reads2.fastq'), emit: reads2
    path 'barcodes.fastq', emit: barcodes
    path 'primer_removal.log'

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
    container 'ghcr.io/vdblab/biopython:1.70a'
    cpus 1
    input:
    tuple val(meta), path(reads_fq)

    output:
    path 'encoding.txt' , emit: encoding

    script:
    def dealwithgz = reads_fq[0].getName().endsWith("gz")
    def uncompress_str = dealwithgz ? "zcat ${reads_fq[0]} | head -n 400 > tmp.fq" : ""
    def guess_input = dealwithgz  ? "tmp.fq" : reads_fq[0]
    """
    set +o pipefail
    $uncompress_str
    set -o pipefail
    guess-encoding.py $guess_input encoding.txt 2>> guess_encoding.log
    """
}

process rename_for_multiqc{
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
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    outdir = file(outdir)
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
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
	.unique { meta, _reads ->
            meta.rawid  // if we have multiple libraries, we need to toss  duplicated  entries in the sample sheet
	}

    persample_inputs = ch_samplesheet_unique
	.combine(ch_reads_runmerged.first())
	.map{ meta, _reads, _newmeta, newreads ->
	    meta = ["id": params.poolid, "primer_f": meta["primer_f"], "primer_r": meta["primer_r"]] +
	    [fw_primer_revcomp: makeComplement(meta["primer_f"].reverse())] +
                [rv_primer_revcomp: makeComplement(meta["primer_r"].reverse())]
	[meta, newreads[0], newreads[1]]
	}.unique()

    remove_primers (
	persample_inputs
    )
    ch_noprimers = params.trim_readthrough ?
	CUTADAPT_READTHROUGH(remove_primers.out.reads1
			     .combine(remove_primers.out.reads2)
			     .map{ meta1, reads1, _meta2, reads2 ->
		[meta1, [reads1, reads2]] } ).reads :
        remove_primers.out.reads1
	.combine(remove_primers.out.reads2)
	.map{ meta1, reads1, _meta2, reads2 ->
	    [meta1, [reads1, reads2]] }
    guess_encoding (
	ch_noprimers
    )
    /////////////////////////////////////
    map1_path = outdir.resolve( params.poolid + ".map.1")
    map2_path = outdir.resolve( params.poolid + ".map.2")
    header = channel.value("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription")
    header
    .concat(
	ch_samplesheet_unique
	.map{ meta, _reads ->
            "${meta.id}\t${meta.barcode_f}\t${meta.primer_f}\t${meta.primer_r}\t${meta.rawid}"
	    }
    )
    .collectFile(name: map1_path, newLine: true, sort:false)
    header
    .concat(
	ch_samplesheet_unique
	.map{ meta, _reads ->
            "${meta.id}\t${meta.barcode_r}\t${meta.primer_r}\t${meta.primer_f}\t${meta.rawid}"
	    }
    )
    .collectFile(name: map2_path, newLine: true, sort:false)

    sampledir = file("${params.outdir}/sampleids/")
    sampledir.mkdir()

    ch_samplesheet_unique
    .map{ meta, _reads ->
            "${meta.id}"
	    }
    .concat(channel.value("Unassigned"))
    .collectFile{ x ->
        [ sampledir.resolve("${x}.sample"), x ]
	}
    .set{ samplefiles }


    adi_f (
	ch_noprimers.first().map{ _meta, reads  -> reads[0]},
	map1_path,
	remove_primers.out.barcodes,
	guess_encoding.out.encoding,
	1
    )
    adi_r  (
	ch_noprimers.first().map{ _meta, reads  -> reads[1]},
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
    def all_demux_files = demux_f.out.samplefq.mix(
	demux_r.out.samplefq
    ).collect(sort: true)

    def seqkit_input = all_demux_files.map { x ->
	def meta = ["id": "${params.poolid}"]
	[meta, x]
    }
    SEQKIT_STATS (
	seqkit_input

    )

    //SEQKIT_STATS.out.stats.collectFile( name: 'demultiplex_seqkit_stats_mqc.out')
    // SEQKIT_STATS.out.stats
    // 	.map{x -> [x]}
    // 	.collectFile(
    // 	    name: 'demultiplex_seqkit_stats_mqc.out',
    // 	    storeDir: "${params.outdir}/")
    rename_for_multiqc(
	SEQKIT_STATS.out.stats.collect{ it -> it[1] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it -> it[1]}).mix(rename_for_multiqc.out.mqc)
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first()).mix(SEQKIT_STATS.out.versions)
    //all_demux_files.filter( ~/_R1.fastq/ ).view()
 //   MAKE_MANIFEST(all_demux_files.filter( ~/_R1.fastq/ ), paired=true)
MAKE_MANIFEST(
    demux_f.out.samplefq.mix(demux_r.out.samplefq).collect(sort: true),
    true
)
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'mskccdemux_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    //
    // MODULE: MultiQC
    //

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))


    MULTIQC (
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'mskccdemux'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
	})


    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
