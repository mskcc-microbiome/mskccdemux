process CUTADAPT_DEMUX_RC {
    tag "$meta.id"
    label 'process_medium'
    cpus 4

    container 'ghcr.io/vdblab/cutadapt:3.7'

    input:
    tuple val(meta), path(orphan_r1), path(orphan_r2)

    output:
    tuple val(meta), path("demuxed_rc/*_R1.fastq.gz"), path("demuxed_rc/*_R2.fastq.gz"), emit: demuxed
    tuple val(meta), path("rc_orphans_R1.fastq.gz"), path("rc_orphans_R2.fastq.gz"),     emit: orphans


    script:
    def args = task.ext.args ?: ''
    def outstring = lambda wc, input, output: (
            " -o {} -p {} ".format(output.R1, output.R2)
            if len(input.reads) > 1
            else " -o {}".format(output.R1)
        ),
    def primerstring = lambda wc, input, output: (
            " -g correct={F} -g flipped={R} -G flipped={F} -G correct={R}".format(
                F=config["primer_F"], R=config["primer_R"]
            )
            if len(input.reads) > 1
            else " -g correct={F} -g flipped={R} ".format(
                F=config["primer_F"], R=config["primer_R"]
            )
        ),
    def preprocess_min_read_length = config["preprocess_min_read_length"],
    log:
        o=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}.o",
        e=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}.e",
    """
        cutadapt \
            -Z \
            --cores {threads} \
            {params.primerstring} \
            -e 1.5 \
            --minimum-length {params.preprocess_min_read_length} \
            {params.outstring} \
            --rename='{{header}}  {{adapter_name}}' \
            {input.reads} \
            > {log.o} 2>> {log.e}
        """
