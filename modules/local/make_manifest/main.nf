process MAKE_MANIFEST {
    tag  'make_manifest'
    label 'process_single'

    input:
    path demux_files
    val paired

    output:
    path "manifest.tsv", emit: manifest
    path "missing.tsv", emit: missing
    tuple val("${task.process}"), val('make_manifest'), val("1.0"), emit: versions_make_manifest, topic: versions
    script:"""

    tmp_manifest="tmp_manifest"
    tmp_missing="tmp_missing"


    # Create/truncate output files
    printf 'sample_id\tR1\tR2\n' > manifest.tsv
    : > "\$tmp_manifest"
    : > "\$tmp_missing"

    for file in ${demux_files}
    do
	name=\$(basename "\$file")

	case "\$name" in
	    *_R1.fastq*)
		# sample_id = everything before the last "_R1.fastq"
		sample_id=\${name%_R1.fastq*}

		# suffix = everything after the last "_R1.fastq"
		suffix=\${name##*_R1.fastq}

		if [ "\$suffix" = "_empty" ]; then
		    printf '%s\n' "\$sample_id" >> "\$tmp_missing"
		else
		    if [ "\$sample_id" != "Unassigned" ]; then
			r1_path="\$file"
			r2_path=\$(printf '%s\n' "\$r1_path" | sed 's/_R1/_R2/')

			printf '%s\t%s\t%s\n' "\$sample_id" "\$r1_path" "\$r2_path" >> "\$tmp_manifest"
		    fi
		fi
		;;
	esac
    done

    sort "\$tmp_manifest" >> manifest.tsv
    cat "\$tmp_missing" > missing.tsv
    rm  "\$tmp_manifest" "\$tmp_missing"

"""

}
