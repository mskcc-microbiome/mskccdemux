process MAKE_MANIFEST {
    tag  'make_manifest'
    input:
    path demux_files
    val paired

    output:
    path "manifest.tsv", emit: manifest
    path "missing.tsv", emit: missing

    exec:
    def sampleIds = [] // Find all sample IDs from R1 files
    def completeSamples = []
    def incompleteSamples = []
    demux_files.each { file ->
        def matcher = file.name =~ /(.+)_R1\.fastq(.*)/
        if (matcher.find()) {
            def sampleId = matcher.group(1)
	    def isbad = matcher.group(2) == "_empty"
	    if (isbad){
		incompleteSamples.add(sampleId)
	    } else {
		if (sampleId != "Unassigned") {
		    def r1Path = file.resolve().toString()
		    def r2Path = paired ? r1Path.replace("_R1", "_R2") : ""
		    completeSamples.add([sampleId, r1Path, r2Path])
		}
	    }
        }
    }

    completeSamples = completeSamples.sort()

    task.workDir.resolve("manifest.tsv").withPrintWriter { w ->
        w.println("sample_id\tR1\tR2")
        completeSamples.each { row ->
            w.println(row.join('\t'))
        }
    }
    task.workDir.resolve("missing.tsv").withPrintWriter { w ->
        incompleteSamples.each { row ->
            w.println(row.join('\t'))
        }
    }
}
