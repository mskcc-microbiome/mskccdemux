process make_map {

    publishDir params.outdir, mode: 'copy'
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    input:
    path csv_file

    output:
    path "${csv_file.simpleName}.map.1", emit: map1
    path "${csv_file.simpleName}.map.2", emit: map2
    path 'sampleids/*.sample',  emit: samplefiles
    /*
     *   1. Skips the CSV header.
     *   2. Re‑orders columns -> SampleID, BarcodeSequence, LinkerPrimerSequence, ReversePrimer, Description
     *   3. Writes the QIIME header plus transformed rows to a .tsv file.
     */
    exec:
    csv_file = file(csv_file)
    outfile1 = task.workDir.resolve("${csv_file.simpleName}.map.1")
    outfile2 = task.workDir.resolve("${csv_file.simpleName}.map.2")
    outfile1.withPrintWriter { w ->
        w.println('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription')
        csv_file.eachLine { line, index ->
            if (index != 1) { // skip header row
		def cols = line.split(',')
		def (sample, primer_f, primer_r, barcode_f, barcode_r) = cols
		w.println([sample, barcode_f, primer_f, primer_r, sample].join('\t'))
	    }
        }
    }
    outfile2.withPrintWriter { w ->
        w.println('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription')
        csv_file.eachLine { line, index ->
            if (index != 1) { // skip header row
		def cols = line.split(',')
		def (sample, primer_f, primer_r, barcode_f, barcode_r) = cols
		w.println([sample, barcode_f,  primer_r, primer_f, sample].join('\t'))
	    }
        }
    }
    // write out the sample files for scattering filter_fasta.py
    sampledir = task.workDir.resolve("sampleids")
    sampledir.mkdir()
    csv_file.eachLine { line, index ->
        if (index != 1) { // skip header row
	    def cols = line.split(',')
	    def (sample, primer_f, primer_r, barcode_f, barcode_r) = cols
	    thisout = task.workDir.resolve("sampleids/${sample}.sample")
	    thisout.withPrintWriter { w ->
		w.println([sample])
	    }
	}
    }
    sample = "Unassigned.sample"
    thisout = task.workDir.resolve("sampleids/$sample")
    thisout.withPrintWriter { w ->
	w.println([sample])
    }
}
