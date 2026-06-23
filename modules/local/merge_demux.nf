/*
 * Merge per-sample FASTQs from round 1 and round 2 demux.
 * Each sample's reads from both rounds are concatenated into a single pair.
 */
process MERGE_DEMUX {
    tag "${meta.sample_id}"
    label 'process_single'

    input:
    tuple val(meta), path(fastqs_r1), path(fastqs_r2)

    output:
    tuple val(meta), path("${meta.sample_id}_R1.fastq.gz"), path("${meta.sample_id}_R2.fastq.gz"), emit: reads

    script:
    """
    cat ${fastqs_r1} > ${meta.sample_id}_R1.fastq.gz
    cat ${fastqs_r2} > ${meta.sample_id}_R2.fastq.gz
    """
}
