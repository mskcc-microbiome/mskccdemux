/*
 * Generate a manifest TSV listing all demuxed samples and their FASTQ paths,
 * plus a separate file listing missing or incomplete samples.
 */
process OUTPUT_MANIFEST {
    tag "$meta.id"
    label 'process_single'

//    container 'biocontainers/python:3.12'

    input:
    tuple val(meta), path(fastqs)
    path sample_list

    output:
    tuple val(meta), path("${meta.id}_manifest.tsv"),              emit: manifest
    tuple val(meta), path("${meta.id}_missing_or_incomplete.tsv"), emit: missing
    path "versions.yml",                                           emit: versions

    script:
    """
    #!/usr/bin/env python3
    import gzip, os, sys

    def gz_size(fname):
        try:
            with gzip.open(fname, 'rb') as f:
                f.seek(0, 2)
                return f.tell()
        except Exception:
            return 0

    # Read expected sample list
    with open("${sample_list}") as f:
        samples = [l.strip() for l in f if l.strip()]

    # Add orphans
    samples.append("orphans")

    manifest_rows = []
    missing_rows = []

    for s in samples:
        r1 = f"{s}_R1.fastq.gz"
        r2 = f"{s}_R2.fastq.gz"
        r1_ok = os.path.exists(r1) and gz_size(r1) > 0
        r2_ok = os.path.exists(r2) and gz_size(r2) > 0
        if r1_ok and r2_ok:
            manifest_rows.append(f"{s}\\t{os.path.abspath(r1)}\\t{os.path.abspath(r2)}")
        else:
            missing_rows.append(f"{s}\\t{os.path.abspath(r1) if r1_ok else ''}\\t{os.path.abspath(r2) if r2_ok else ''}")

    with open("${meta.id}_manifest.tsv", "w") as f:
        f.write("sample_id\\tR1\\tR2\\n")
        f.write("\\n".join(manifest_rows) + "\\n" if manifest_rows else "")

    with open("${meta.id}_missing_or_incomplete.tsv", "w") as f:
        f.write("sample_id\\tR1\\tR2\\n")
        f.write("\\n".join(missing_rows) + "\\n" if missing_rows else "")

    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f"  python: {sys.version.split()[0]}\\n")
    """
}
