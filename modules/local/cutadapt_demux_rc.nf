/*
 * Demultiplex round 2: check reverse-complement barcodes on orphan reads.
 * Input R1/R2 are swapped relative to standard orientation to detect
 * fragments flipped during library prep.
 */
process CUTADAPT_DEMUX_RC {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/vdblab/cutadapt:3.7'

    input:
    tuple val(meta), path(orphan_r1), path(orphan_r2)
    path barcodes_r2

    output:
    tuple val(meta), path("demuxed_rc/*_R1.fastq.gz"), path("demuxed_rc/*_R2.fastq.gz"), emit: demuxed
    tuple val(meta), path("rc_orphans_R1.fastq.gz"), path("rc_orphans_R2.fastq.gz"),     emit: orphans

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p demuxed_rc

    # Round 2: swap R1/R2 input, swap -o/-p output to preserve sequencing orientation
    cutadapt \\
        -Z \\
        --cores ${task.cpus} \\
        -e 2 \\
        --no-indels \\
        --minimum-length 150 \\
        ${args} \\
        -g file:${barcodes_r2} \\
        -o 'demuxed_rc/rc_{name}_R2.fastq.gz' \\
        -p 'demuxed_rc/rc_{name}_R1.fastq.gz' \\
        --untrimmed-output rc_orphans_R2.fastq.gz \\
        --untrimmed-paired-output rc_orphans_R1.fastq.gz \\
        ${orphan_r2} ${orphan_r1}

    """
}
