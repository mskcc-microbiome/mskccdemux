/*
 * Demultiplex paired-end reads using cutadapt with barcode FASTA.
 * Produces per-sample paired FASTQs plus orphan (unmatched) reads.
 */
process CUTADAPT_DEMUX {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/vdblab/cutadapt:3.7'

    input:
    tuple val(meta), path(reads)
    path barcodes

    output:
    tuple val(meta), path("demuxed/*_R1.fastq.gz"), path("demuxed/*_R2.fastq.gz"), emit: demuxed
    tuple val(meta), path("orphans_R1.fastq.gz"), path("orphans_R2.fastq.gz"),     emit: orphans

    script:
    def (r1, r2) = reads
    def args = task.ext.args ?: ''
    """
    mkdir -p demuxed

    cutadapt \\
        -Z \\
        --cores ${task.cpus} \\
        -e 2 \\
        --no-indels \\
        --minimum-length 150 \\
        ${args} \\
        -g file:${barcodes} \\
        -o 'demuxed/f_{name}_R1.fastq.gz' \\
        -p 'demuxed/f_{name}_R2.fastq.gz' \\
        --untrimmed-output orphans_R1.fastq.gz \\
        --untrimmed-paired-output orphans_R2.fastq.gz \\
        ${r1} ${r2}

    cat <<-END_VERSIONS
    versions:
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
