#!/usr/bin/env python
#
#Dec/24/2019
#Script to identify phred scale decoding in fastq file.
#
#phred score input is necessary in split_libraries_fastq.py
#
#Usage:
#head -1000 your.fastq | awk 'NR % 4 == 0' | guess-encoding.py -n 1000 | | awk '$0~/phred33/ {a=33} $0~/phred64/ {a=a""64} {print a}'
#
#Script was modified from this:
#https://wiki.bits.vib.be/index.php/Identify_the_Phred_scale_of_quality_scores_used_in_fastQ
#
#Return:
#- 33 if compatible with Illumina-1.8
#- 64 if compatible with Illumina-1.3 or Illumina-1.5
#- 33/64 if compatiblewith both (This would make decoding unidentifiable.)
#

import sys
#from Bio import SeqIO
#import optparse

RANGES = {
     'phred33' : (33, 73),
     'phred64' : (64, 104)
#    'Sanger': (33, 73),
#    'Illumina-1.8': (33, 74),
#    'Solexa': (59, 104),
#    'Illumina-1.3': (64, 104),
#    'Illumina-1.5': (67, 104)
}

def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in list(ranges.items()):
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def main(infile, outfile):
    n=-1
    gmin, gmax  = 99, 0
    valid = []
    with open(infile, "r") as inf:
        #for rec in SeqIO.parse(inf, "fastq"):
        for i, line in enumerate(inf):
            if (i+1) % 4 != 0:
                continue
            print(line)
            #score=rec.letter_annotations["phred_quality"]
            #print(score)
            #outf.write(score + "\n")
            #sys.exit(0)
            lmin, lmax = get_qual_range(line.rstrip())
            if lmin < gmin or lmax > gmax:
                gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                valid = get_encodings_in_range(gmin, gmax)
                if len(valid) == 0:
                    print("no encodings for range: %s" % str((gmin, gmax)), file=sys.stderr)
                    sys.exit(1)
                if len(valid) == 1 and n == -1:
                    break

            if n > 0 and i > n:
                break
    print(valid)
    with open(outfile, "w") as outf:
        outf.write("\t".join([valid[0], str(RANGES[valid[0]][0]), str(gmin), str(gmax)]) + "\n")
        #outf.write("\t".join([valid[0], str(gmin), str(gmax)]) + "\n")

if __name__ == "__main__":
    main(infile = sys.argv[1], outfile = sys.argv[2])
else:
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        maina(infile = snakemake.input[0], outfile = snakemake.output[0])
