process demultiplex {
    tag 'demultiplex'

    input:
    path reads
    path sample_file
    val readDir

    output:
    path "${sample_file.getName().replaceFirst(/.sample/, '')}_R${readDir}.fastq.gz" , emit: samplefq

    container 'ghcr.io/vdblab/qiime:1.9.1'
    cpus 1
    memory '16 GB'
    script:
    def samplename = sample_file.getName().replaceFirst(/.sample/, '')
    """
    filter_fasta.py -f $reads --sample_id_fp $sample_file -o ${samplename}_R${readDir}.fastq 2>> demultiplex_${samplename}.log
    gzip -f ${samplename}_R${readDir}.fastq
    """
}
