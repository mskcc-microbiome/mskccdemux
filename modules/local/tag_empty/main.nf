process tag_empty{
    tag 'tag_empty'
    """ Ampliseq complains (correctly) about empty fastqs; after seqkit stats
    , this adds _empty to the fastqs so they don't get picked up by ampliseq's
    glob input specification

    """
    publishDir "${params.outdir}/isolated_oligos", mode: 'copy'
    input:
    path x

    output:
    path "*[.gz|*_empty]", includeInputs: true
    script:
    """
    if LC_ALL=C gzip -l $x | awk 'NR==2 {exit(\$2!=0)}';
    then
	mv $x ${x}_empty
    fi
    """
}
