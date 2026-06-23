/*
 * Parse a mothur-style oligos file and emit per-direction barcode FASTA files
 */
process PARSE_OLIGOS {
    tag "$meta.id"
    label 'process_single'

    //container 'biocontainers/python:3.12'

    input:
    tuple val(meta), path(oligos)

    output:
    tuple val(meta), path("barcodes_R1.fasta"), emit: barcodes_r1
    tuple val(meta), path("barcodes_R2.fasta"), emit: barcodes_r2
    tuple val(meta), path("samples.txt"),       emit: sample_list

    script:
    """
    #!/usr/bin/env python3
    import re, sys

    sample_ids = []
    with open("${oligos}") as f, \
         open("barcodes_R1.fasta", "w") as r1, \
         open("barcodes_R2.fasta", "w") as r2, \
         open("samples.txt", "w") as sl:
        for line in f:
            if not line.upper().startswith("BARCODE"):
                continue
            pieces = line.split()
            sample_id = pieces[3].split("..")[0]
            sample_id = re.sub(r"[_%+;: -]+", ".", sample_id)
            r1.write(f">{sample_id}\\n{pieces[1]}\\n")
            r2.write(f">{sample_id}\\n{pieces[2]}\\n")
            sl.write(f"{sample_id}\\n")
            sample_ids.append(sample_id)

    assert len(sample_ids) > 0, "No barcodes found in oligos file"

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f"  python: {sys.version.split()[0]}\\n")
    """
}
