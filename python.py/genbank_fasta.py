#!/usr/bin/env python
from Bio import SeqIO
gbk_file = "./AY810830.gbk"
fna_file = "./AY810830.fasta"
input = open(gbk_file, "r")
output = open(fna_file, "w")

for seq in SeqIO.parse(input, "genbank") :
    print("Convert %s from GenBank to FASTA format" % seq.id)
    output.write(">%s %s\n%s\n" % (
        seq.id, seq.description, seq.seq))

input.close()
output.close()
