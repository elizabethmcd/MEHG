#! /Users/emcdaniel/anaconda3/bin/python3
import pandas as pd 

with open("2500-hgcA-fasta-headers.txt") as meth_fasta:
    for line in meth_fasta:
        scaffold = line.split(" ")[0]
        annotation = line.split("id")[0].split(" ")[1:]
        print(annotation)
