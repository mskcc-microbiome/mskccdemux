process demultiplex {
    """ Demultiplex with qiime, and tag empty file


    Ampliseq complains (correctly) about empty fastqs; after seqkit stats
    , this adds _empty to the fastqs so they don't get picked up by ampliseq's
    glob input specification

    """
    tag 'demultiplex'

    input:
    path reads
    path sample_file
    val readDir

    output:
    path "*[.gz|_empty]" , emit: samplefq



    container 'ghcr.io/vdblab/qiime:1.9.1'
    cpus 1
    memory '16 GB'
    script:
    def samplename = sample_file.getName().replaceFirst(/\.sample/, '')
    def outfile    = "${samplename}_R${readDir}.fastq"

    """
    filter_fasta.py -f $reads --sample_id_fp $sample_file -o $outfile 2>> demultiplex_${samplename}_r${readDir}.log
    if [ -s "$outfile" ]; then
        gzip -f $outfile
    else
        mv  $outfile  ${outfile}_empty
    fi

    """
}
