#! /usr/bin/env python

# Archaeal list
with open("arch-genomes.txt") as arch_list:
    for index, line in enumerate(arch_list):
        archid= str(index + 1).zfill(5)
        archname=line.replace(".fna", "").strip().splitlines()[0]
        with open("arch-genomes-new-names.txt", "a") as outfile:
            outfile.write(str(archname)+"\t"+str(archname)+".fna"+"\t"+"archaea"+archid+"\t"+"archaea"+archid+".fna"+"\n")

# Bacterial list
with open("bac-genomes.txt") as bac_list:
    for index, line in enumerate(bac_list):
        bacid= str(index + 1).zfill(5)
        bacname=line.replace(".fna", "").strip().splitlines()[0]
        with open("bac-genomes-new-names.txt", "a") as outfile:
            outfile.write(str(bacname)+"\t"+str(bacname)+".fna"+"\t"+"bacteria"+bacid+"\t"+"bacteria"+bacid+".fna"+"\n")

