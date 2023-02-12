#!/usr/bin/env python
from Bio import SeqIO
genbank_file = open('./AY810830.gbk','r')
fasta_file = open('./AY810830.fa', 'w')
records = SeqIO.parse(genbank_file, 'genbank')
SeqIO.write(records,fasta_file,'fasta')
fasta_file.close()
