# import all the libraries needed
# rename fasta header to be in format like: TEM_1_length
from Bio import SeqIO
import sys


with open(sys.argv[1]) as f1, open(sys.argv[2],'w') as f2:
	sequences = SeqIO.parse(f1,'fasta')
	n = 1
	for record in sequences:
		length = len(record.seq)
		fname = 'TEM_' + str(n) + '_' + str(length)
		n += 1
		record.id = fname
		record.description = fname
		SeqIO.write(record,f2,'fasta')
