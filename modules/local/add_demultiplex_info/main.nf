process add_demultiplex_info {
    tag 'add_demux_info'
    input:
    path reads
    path map
    path barcodes
    path encoding
    val dirNum
    output:
    path "demultiplex_r${dirNum}/seqs.fastq", emit: seqsfq
    path "demultiplex_r${dirNum}/seqs.fna",   emit: seqsfna
    container 'ghcr.io/vdblab/qiime:1.9.1'
    cpus 1
    memory '12 GB'
    script:
    """
    split_libraries_fastq.py \\
    -i ${reads} -b ${barcodes} \\
      -o demultiplex_r${dirNum} \\
    --store_demultiplexed_fastq --phred_offset \$(cut ${encoding} -f2) \\
    -m ${map} -q 0 -s 100000000 --barcode_type 12 \\
      --max_barcode_errors 4 -n 5000 -p 1e-5 -r 1000000 \\
      --retain_unassigned_reads \
      >> add_demux_r${dirNum}.log 2>&1
    """
}
