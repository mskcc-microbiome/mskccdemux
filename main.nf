#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    amplicon-demux: Amplicon demultiplexing pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two-round barcode demultiplexing using cutadapt.
    Round 1: match forward barcodes.
    Round 2: re-check orphans against reverse barcodes (catches flipped fragments).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSE_OLIGOS     } from './modules/local/parse_oligos'
include { CUTADAPT as CUTADAPT_REMOVEPRIMERS         } from './modules/nf-core/cutadapt/main'
include { CUTADAPT_DEMUX   } from './modules/local/cutadapt_demux'
include { CUTADAPT_DEMUX_RC} from './modules/local/cutadapt_demux_rc'
include { MERGE_DEMUX      } from './modules/local/merge_demux'
include { OUTPUT_MANIFEST  } from './modules/local/output_manifest'
include { FASTQC           } from './modules/nf-core/fastqc/main'
include { MULTIQC          } from './modules/nf-core/multiqc/main'
include { SEQKIT_STATS     } from './modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_MERGED } from './modules/nf-core/seqkit/stats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// taken from ampliseq
// Complement table taken from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
def makeComplement(seq) {
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]
    def comp = seq.toUpperCase().collect { base -> complements[ base ] ?: 'X' }.join()
    return comp
}


// process rename_for_multiqc {
//     input:
//     path seqkit
//     output:
//     path "demultiplex_seqkit_stats_mqc.out" , emit: mqc

//     script:
//     """
//     cp $seqkit demultiplex_seqkit_stats_mqc.out
//     """
// }


workflow {

    // -------------------------------------------------------------------------
    // 1. Parse samplesheet: each row is one library (one R1/R2 pair)
    //    Columns: pool,R1,R2,oligos
    // -------------------------------------------------------------------------
    ch_samplesheet = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [id: row.pool]
            def r1   = file(row.R1, checkIfExists: true)
            def r2   = file(row.R2, checkIfExists: true)
            [meta, r1, r2, file(row.oligos, checkIfExists: true)]
        }
    ch_pool_primers = ch_samplesheet
	.map { meta, _r1, _r2, oligos ->
	    def pool_prf_prr = oligos
		.splitCsv(sep: '\t', header: false)
		.first()
		// .map { row ->

		//     tuple(row[0], row[1], row[2])
	    // }
	    if( pool_prf_prr.size() != 3 )
		error "Expected 3 fields in first line of , got ${pool_prf_prr.size()}: ${pool_prf_prr}"

	    meta + [
	     "fw_primer": pool_prf_prr[1],
	     "rv_primer": pool_prf_prr[2],
	     "fw_primer_revcomp": makeComplement(pool_prf_prr[1].reverse()),
	     "rv_primer_revcomp": makeComplement(pool_prf_prr[2].reverse())
	    ]
	}

    // Group all library FASTQs by pool
    ch_by_pool = ch_samplesheet
        .map { meta, r1, r2, _oligos -> [meta, r1, r2] }
        .groupTuple()
        .map { meta, r1s, r2s -> [meta, r1s.flatten(), r2s.flatten()] }

    // Get oligos (one per pool, take first)
    ch_oligos = ch_samplesheet
        .map { meta, _r1, _r2, oligos -> [meta, oligos] }
        .unique()

    // -------------------------------------------------------------------------
    // 2. FastQC on raw input libraries
    // -------------------------------------------------------------------------
    ch_fastqc_input = ch_samplesheet
        .map { meta, r1, r2, _oligos ->
            def lib_meta = [id: "${meta.id}_${r1.baseName}", single_end: false]
            [lib_meta, [r1, r2]]
        }

    FASTQC(ch_fastqc_input)

    // -------------------------------------------------------------------------
    // 3. Parse oligos file → barcode FASTAs
    // -------------------------------------------------------------------------
    PARSE_OLIGOS(ch_oligos)
    // -------------------------------------------------------------------------
    // 4. Concatenate multi-library FASTQs and split for parallel demux
    //    Using Nextflow's native splitFastq for scatter/gather
    // -------------------------------------------------------------------------
    ch_reads_split = ch_by_pool
        .map { meta, r1s, r2s ->
            // Concat conceptually happens here; if single library, pass through
            [meta, r1s, r2s]
        }
        .flatMap { meta, r1s, r2s ->
            // Pair up R1 and R2 files for each library chunk
            [r1s, r2s].transpose().collect { r1, r2 ->
                [meta, [r1, r2]]
            }
        }

    // -------------------------------------------------------------------------
    // 5. Demux round 1: forward barcodes
    // -------------------------------------------------------------------------
    ch_barcodes_r1 = PARSE_OLIGOS.out.barcodes_r1.map { _meta, f -> f }
    ch_barcodes_r2 = PARSE_OLIGOS.out.barcodes_r2.map { _meta, f -> f }

    CUTADAPT_DEMUX(ch_reads_split, ch_barcodes_r1.first())

    // -------------------------------------------------------------------------
    // 6. Demux round 2: reverse barcodes on orphans (flipped fragments)
    // -------------------------------------------------------------------------
    CUTADAPT_DEMUX_RC(CUTADAPT_DEMUX.out.orphans, ch_barcodes_r2.first())

    // -------------------------------------------------------------------------
    // 7. Merge rounds: combine per-sample reads from round 1 + round 2
    //    Flatten the demuxed outputs, extract sample_id from filename,
    //    group by sample, then cat together.
    // -------------------------------------------------------------------------

    // Extract individual sample files from round 1
    ch_round1 = CUTADAPT_DEMUX.out.demuxed
        .flatMap { meta, r1_files, r2_files ->
            def r1_list = r1_files instanceof List ? r1_files : [r1_files]
            def r2_list = r2_files instanceof List ? r2_files : [r2_files]
            [r1_list, r2_list].transpose().collect { r1, r2 ->
                def sample_id = r1.name.replaceAll('_R1\\.fastq\\.gz$', '').replaceAll('f_', '')

                [["sample_id":sample_id] + meta, r1, r2]
            }
        }

    // Extract individual sample files from round 2
    ch_round2 = CUTADAPT_DEMUX_RC.out.demuxed
        .flatMap { meta, r1_files, r2_files ->
            def r1_list = r1_files instanceof List ? r1_files : [r1_files]
            def r2_list = r2_files instanceof List ? r2_files : [r2_files]
            [r1_list, r2_list].transpose().collect { r1, r2 ->
                def sample_id = r1.name.replaceAll('_R1\\.fastq\\.gz$', '').replaceAll('rc_', '')
                [["sample_id":sample_id] + meta, r1, r2]
            }
        }

    // Also collect orphans from round 2 as "orphans" sample
    ch_orphans = CUTADAPT_DEMUX_RC.out.orphans
        .map { meta, r1, r2 ->
	    [["sample_id": "orphans"] + meta, r1, r2]
	}

    // Combine all, group by sample_id, merge
    ch_all_demuxed = ch_round1
        .mix(ch_round2, ch_orphans)
        .groupTuple()
    MERGE_DEMUX(ch_all_demuxed)
    def seqkit_input_merged = MERGE_DEMUX.output.reads.map {meta, r1, r2 ->
	def meta_ = ["sample_id": "${meta.sample_id}"]
	[meta_, [r1, r2]]
    }

    SEQKIT_STATS_MERGED(seqkit_input_merged)
    ch_for_primer_removal = MERGE_DEMUX.output.reads
	.map { meta, r1, r2 ->
	    // change the order for merging
	    tuple(meta.id, [meta, r1, r2])
	}
	.combine(ch_pool_primers
	      .map{contents ->
		    tuple(contents.id, contents)
		}, by: 0)
	.map { joined ->
	    def (id,  meta1_raw, meta2) = joined
	    def (meta1, r1, r2) = meta1_raw
	    def newmeta =  ["sample_id": meta1.sample_id] + meta2
	    tuple(newmeta, [r1, r2])
	}


    CUTADAPT_REMOVEPRIMERS(ch_for_primer_removal)
    def seqkit_input_noprimer = CUTADAPT_REMOVEPRIMERS.output.reads.map {meta, reads ->
	def meta_ = ["sample_id": "${meta.sample_id}"]
	[meta_, reads]
    }
    SEQKIT_STATS(seqkit_input_noprimer)
    // -------------------------------------------------------------------------
    // 8. Output manifest
    // -------------------------------------------------------------------------
    ch_gathered_fastqs = CUTADAPT_REMOVEPRIMERS.out.reads
	.collect(flat:false)


    ch_sample_list = PARSE_OLIGOS.out.sample_list.map { _meta, f -> f }

    //OUTPUT_MANIFEST(ch_gathered_fastqs, ch_sample_list.first())
    checked = CUTADAPT_REMOVEPRIMERS.out.reads
	.branch { row ->
	with_reads:
        row[1] &&
            row[1].size() >= 2 &&
            row[1][0].exists() &&
            row[1][0].countFastq() > 0

	missing_or_empty:
        true
    }

    checked.with_reads
	.map { row ->
            def (meta, reads)  = row
            tuple(meta.id, "${meta.id},${meta.sample_id},${reads[0]},${reads[1]}\n")
	}
	.collectFile(
	                seed: 'pool_id,sample_id,read1,read2\n',
		storeDir: params.outdir,
		sort: true,
	){ pool, line ->
	    ["${pool}_manifest.tsv", "$line"]

	}

    checked.missing_or_empty
	.map { row ->
            def meta = row[0]
            tuple( meta.id, "${meta.id},${meta.sample_id}\n")
	}
	.collectFile(
            seed: 'pool_id,sample_id\n',
	    storeDir: params.outdir,
	    sort: false,
	){ pool, line ->
	["${pool}_missing_or_empty.tsv", "$line"]

    }
    // -------------------------------------------------------------------------
    // 9. MultiQC
    // -------------------------------------------------------------------------
    // rename_for_multiqc(
    // 	SEQKIT_STATS_MERGED.out.stats.collect{ it -> it[1] }
    // )
    ch_multiqc_input = FASTQC.out.zip.map { _meta, zip -> zip }
	.mix(SEQKIT_STATS_MERGED.out.stats.collect{ it -> it[1] })
	.mix(SEQKIT_STATS.out.stats.collect{ it -> it[1] })
        .collect()
        .map { files -> [[id: "multqc"], files, [], [], [], []] }

    MULTIQC(ch_multiqc_input)
}
