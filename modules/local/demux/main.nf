process demultiplex {
    publishDir "${params.outdir}/isolated_oligos", mode: 'copy'
    tag 'demultiplex'

    input:
    path reads
    path sample_file
    val readDir

    output:
    path "${sample_file.simpleName}_R${readDir}.fastq.gz" , emit: samplefq

    container 'ghcr.io/vdblab/qiime:1.9.1'
    cpus 1
    memory '16 GB'
    script:
    """
    filter_fasta.py -f $reads --sample_id_fp $sample_file -o ${sample_file.simpleName}_R${readDir}.fastq 2>> demultiplex_${sample_file.simpleName}.log
    gzip -f ${sample_file.simpleName}_R${readDir}.fastq
    """
}
