#! /usr/bin/env python3
import sys
import os
import regex as re
import argparse
import gzip
import itertools
from Bio import SeqIO

# system("strip_addons3.py ./pool218_complete_CAGATC_L001_R1_001.fastq.gz ./pool218_complete_CAGATC_L001_R2_001.fastq.gz
#        -remove_bar_primer -fw_primer AYTGGGYDTAAAGNG -rev_primer CCGTCAATTYHTTTRAGT")


def primer_pattern(primer,pdiffs=0):
	assert re.match('^[ACGTKMRYSWBVHDXN]+$',primer),'YTError: function primer_pattern encountered strange letters!'
	pattern = re.sub('K','[GT]',primer)
	pattern = re.sub('M','[AC]',pattern)
	pattern = re.sub('R','[AG]',pattern)
	pattern = re.sub('Y','[CT]',pattern)
	pattern = re.sub('S','[CG]',pattern)
	pattern = re.sub('W','[AT]',pattern)
	pattern = re.sub('B','[CGT]',pattern)
	pattern = re.sub('V','[ACG]',pattern)
	pattern = re.sub('H','[ACT]',pattern)
	pattern = re.sub('D','[AGT]',pattern)
	pattern = re.sub('X|N','[ACTG]',pattern)
	if pdiffs>0:
		pattern = pattern+'{e<='+str(pdiffs)+'}'
	return pattern



#strip barcode addons at beginning of fastq files
def strip_barcode_addons(seqfile1,seqfile2,pdiffs=0,remove_bar_primer=False,fw_primer=None,rev_primer=None,test=False):
	if fw_primer is None:
		fw_primer = 'AYTGGGYDTAAAGNG'
	if rev_primer is None:
		rev_primer = 'CCGTCAATTYHTTTRAGT'
	print('primer',fw_primer,rev_primer)
	#can be fasta, fastq, or fastq.gz
	seq_pattern = re.compile('.(fasta|fastq)(.gz)?$')
	assert seq_pattern.search(seqfile1) is not None, 'YTError:'+seqfile1+'needs to end in .fasta, .fastq, or .fastq.gz!'
	assert seq_pattern.search(seqfile2) is not None, 'YTError:'+seqfile2+'needs to end in .fasta, .fastq, or .fastq.gz!'
	seq_format = seq_pattern.search(seqfile1).group(1)
	seq_format2 = seq_pattern.search(seqfile2).group(1)
	gzipped = seq_pattern.search(seqfile1).group(2) is not None
	gzipped2 = seq_pattern.search(seqfile2).group(2) is not None
	assert seq_format==seq_format2 and gzipped==gzipped2, 'YTError: These 2 sequence files do not match!'
	print('Trimming off barcode addons...')
	fw_primer_regex = re.compile(primer_pattern(fw_primer,pdiffs=pdiffs))
	rev_primer_regex = re.compile(primer_pattern(rev_primer,pdiffs=pdiffs))
	len_barcode = 12 #golay barcodes
	max_len_addon = 8 #addon is 1-8 bp
	max_len_fw = len_barcode + max_len_addon + len(fw_primer)
	max_len_rev = len_barcode + max_len_addon + len(rev_primer)
	trim_seqfile1 = 'reads1.fastq'
	trim_seqfile2 = 'reads2.fastq'
	scrap_seqfile = 'scrap.fastq'
	bar_seqfile = 'barcodes.fastq'
	if gzipped: #gzipped
		seq1 = gzip.open(seqfile1,'rt')
		seq2 = gzip.open(seqfile2,'rt')
	else: #not gzipped
		seq1 = open(seqfile1,'rt')
		seq2 = open(seqfile2,'rt')
	newseq1 = open(trim_seqfile1,'wt')
	newseq2 = open(trim_seqfile2,'wt')
	badseq = open(scrap_seqfile,'wt')
	if remove_bar_primer:
		barcode = open(bar_seqfile,'wt')
	n_forward_reverse = 0 #counts forward>reverse seqs (p1-p2)
	n_reverse_forward = 0 #counts reverse>forward seqs (p2-p1)
	n_noP1 = 0 #counts forward non-hits
	n_noP2 = 0 #counts reverse non-hits
	i=0
	for (forward,reverse) in zip(SeqIO.parse(seq1,seq_format),SeqIO.parse(seq2,seq_format)):
		assert forward.id==reverse.id, 'YTError: sequence IDs do not match!'
		#these will erase the second half of the header (id=first part,description=whole line)
		#for compatibility with usearch
		forward.description=forward.id
		reverse.description=reverse.id
		#look for forward primer in seq1 and reverse primer in seq2
		r_front = fw_primer_regex.search(str(forward.seq[:max_len_fw]))
		if r_front is not None: #forward primer hit seq1, look for rev primer in seq2.
			r_end = rev_primer_regex.search(str(reverse.seq[:max_len_rev]))
			if r_end is not None: #seq1 is forward, seq2 is reverse
				primerstart1 = r_front.start()
				barstart1 = max(primerstart1-len_barcode,0)
				seqstart1 = primerstart1+len(fw_primer)
				primerstart2 = r_end.start()
				barstart2 = max(primerstart2-len_barcode,0)
				seqstart2 = primerstart2+len(rev_primer)
				if remove_bar_primer:
					SeqIO.write(forward[seqstart1:],newseq1,seq_format)
					SeqIO.write(reverse[seqstart2:],newseq2,seq_format)
					bar1 = forward[barstart1:primerstart1]
					bar2 = reverse[barstart2:primerstart2]
					bar = bar1 + bar2
					bar.description = bar1.description
					SeqIO.write(bar,barcode,seq_format)
				else:
					SeqIO.write(forward[barstart1:],newseq1,seq_format)
					SeqIO.write(reverse[barstart2:],newseq2,seq_format)
				n_forward_reverse += 1
			else: #seq1 is forward but seq2 nohit, ERROR
				forward.id = forward.id+'|noP2'
				reverse.id = reverse.id+'|noP2'
				SeqIO.write(forward,badseq,seq_format)
				SeqIO.write(reverse,badseq,seq_format)
				n_noP2 += 1
		else: #first primer nohit, look for rev primer in seq1
			r_front = rev_primer_regex.search(str(forward.seq[:max_len_rev]))
			if r_front is not None: #rev primer hit seq1, look for forward primer in seq2.
				r_end = fw_primer_regex.search(str(reverse.seq[:max_len_fw]))
				if r_end is not None: #seq1 is reverse primer, seq2 is forward primer.
					primerstart1 = r_front.start()
					barstart1 = max(primerstart1-len_barcode,0)
					seqstart1 = primerstart1+len(rev_primer)
					primerstart2 = r_end.start()
					barstart2 = max(primerstart2-len_barcode,0)
					seqstart2 = primerstart2+len(fw_primer)
					if remove_bar_primer:
						SeqIO.write(forward[seqstart1:],newseq2,seq_format)
						SeqIO.write(reverse[seqstart2:],newseq1,seq_format)
						bar1 = forward[barstart1:primerstart1]
						bar2 = reverse[barstart2:primerstart2]
						bar = bar2 + bar1
						bar.description = bar2.description
						SeqIO.write(bar,barcode,seq_format)
					else:
						SeqIO.write(forward[barstart1:],newseq2,seq_format)
						SeqIO.write(reverse[barstart2:],newseq1,seq_format)
					n_reverse_forward += 1
				else: #rev primer hit seq1, but seq2 is nohit for forward primer.
					forward.id = forward.id+'|noP2'
					reverse.id = reverse.id+'|noP2'
					SeqIO.write(forward,badseq,seq_format)
					SeqIO.write(reverse,badseq,seq_format)
					n_noP2 += 1
			else: #seq1 is nohit for both primers
				forward.id = forward.id+'|noP1'
				reverse.id = reverse.id+'|noP1'
				SeqIO.write(forward,badseq,seq_format)
				SeqIO.write(reverse,badseq,seq_format)
				n_noP1 += 1
		if test and n_reverse_forward+n_forward_reverse>=300:
			break
	seq1.close()
	seq2.close()
	newseq1.close()
	newseq2.close()
	badseq.close()
	if remove_bar_primer:
		barcode.close()
	totalseqs = n_forward_reverse+n_reverse_forward+n_noP1+n_noP2
	print('forward-reverse primer:'+str(n_forward_reverse)+'('+str(round(float(n_forward_reverse)/totalseqs*100,1))+'%)')
	print('reverse-forward primer: '+str(n_reverse_forward)+'('+str(round(float(n_reverse_forward)/totalseqs*100,1))+'%)')
	print('no P1 match: '+str(n_noP1)+'('+str(round(float(n_noP1)/totalseqs*100,1))+'%)')
	print('no P2 match: '+str(n_noP2)+'('+str(round(float(n_noP2)/totalseqs*100,1))+'%)')
	print('Total Seqs: '+str(totalseqs)+' (100%)')




def main():
	#def strip_barcode_addons(seqfile1,seqfile2,pdiffs=0,separate_barfile=False,test=False):
	parser = argparse.ArgumentParser(description='YT: Trims off barcode addons. Output will be named reads1.fastq and reads2.fastq. Sequences that do not match both forward and reverse are stored in scrap.fastq.')
	parser.add_argument('seqfile1',help='Forward fastq file.')
	parser.add_argument('seqfile2',help='Reverse fastq file.')
	parser.add_argument('-pdiffs',type=int,help='Number of primer bp mismatches allowed (default is 0).',metavar='<int>',default=0)
	parser.add_argument('-remove_bar_primer',action='store_true',help='In addition to stripping addons, remove barcode and primer, and store barcode separately in barcodes.fastq. Default=false, where barcode and primer stay on the sequences.')
	parser.add_argument('-fw_primer',type=str,help='Forward primer. Default is 16S V4: AYTGGGYDTAAAGNG')
	parser.add_argument('-rev_primer',type=str,help='Reverse primer. Default is 16S V4: CCGTCAATTYHTTTRAGT.')
	parser.add_argument('-test',action='store_true',help='Test mode. Only does first 300 seqs.')
	args = parser.parse_args()
	strip_barcode_addons(**vars(args))

if __name__ == '__main__':
	main()
