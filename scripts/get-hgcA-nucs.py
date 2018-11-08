#! /usr/bin/env python

from Bio import SeqIO
import sys 

with open("hgcA-bacteroidetes-locus-tags.txt") as locus_list:
	loclines=locus_list.read().splitlines()
gb_file = sys.argv[1]

for record in SeqIO.parse(gb_file, "genbank"):
	for f in record.features:
		if f.type=="CDS" and "protein_id" in f.qualifiers:
			locus = f.qualifiers["protein_id"][0]
			if locus in loclines:
				locus, product, protein = f.qualifiers["protein_id"][0], f.qualifiers["product"][0], f.qualifiers["translation"][0] 
				with open(locus+".fa", "a") as outfile:
					outfile.write(">"+locus+' '+product+"\n"+protein+"\n")
				span = f.location
				beg= span.start
				end = span.end
				nuc = record.seq[beg:end].reverse_complement()
				with open (locus+"-hgcA.fna", "a") as out:
					out.write(">"+locus+"\n"+str(nuc)+"\n")
