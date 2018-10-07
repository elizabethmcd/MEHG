# Geobacter Methylating Pangenome

To look at the pangenome of methylators, I am going to hone in on Geobacter and the surrounding outgroups. This would be a good test case because there are a number of mehylating Geobacter, and not just sparse ones. Additionally, the methylation protein HgcA will me more similar among Geobacter and therefore hopefully end up in the same cluster. The problem with looking at the pangenome of all methylators is that the mercury methylation protein itself is so divergent that it ends up in several different clusters. I could use this as a way to pinpoint core Geobacter methylating clusters, and then search for those proteins in other groups. The archaea would probably be good to do this way too because the methylators seem to all be in the same genus of the Metholo something. 

## Download the specific genome sets

```
ncbi-genome-download --genus Geobacter --format fasta,assembly-report --assembly-level complete bacteria

ncbi-genome-download --genus Pelobacter --format fasta,assembly-report --assembly-level complete bacteria

ncbi-genome-download --genus Desulfuromonas --format fasta,assembly-report --assembly-level complete bacteria

ncbi-genome-download --genus Geoalkalibacter --format fasta,assembly-report --assembly-level complete bacteria
```

I've downloaded these specific sets based on when I did my Deltaproteobacteria analysis of all methylators. I can then also layer this by methylating vs non-methylating for the pangenomic analysis. 

## Cleanup and Prokka annotation 

Move the cleaned .fa files to a directory called `PROKKA-ANNOTATIONS` to keep this separate from the original genomic files. Then use these annotations for building the pangenome. 

```
# File cleanup
for filename in */*.fna; do mv $filename ${filename%/*}/${filename%/*}.fna; done

# Defline cleanup for Anvi'o
for file in */*.fna; do
    f=$(basename $file .fna);
    anvi-script-reformat-fasta $file -o $f.fa --min-len 0 --simplify-names; 
done

# Run prokka

for fasta in *.fa; do
    N=$(basename $fasta .fa);
    prokka --outdir $N --prefix $N --cpus 15 $fasta;
done
```

## Metadata Extraction 

Make a METADATA directory. Use the script meth-filenames-prep.py to parse the strain names out of the assembly metadata files, and create a masterfile of all the taxa in a particular directory. 

```
#!/usr/bin/env python

import re, glob, os
import pandas as pd

# Create master file with dictionary of dictionaries
all_reports = glob.glob("*_assembly_report.txt")
all_dicts = {}
for report in all_reports:
    test_assembly_report=open(report, "r")
    report_dict = {}
    for line in test_assembly_report:
        if line.startswith("#\n"):
            break
        try: 
            descrip = re.match('#\s(.*):', line).group(1).lstrip(' ').replace(' ', '_')
        except: 
            descrip = None
        item = line.split(":", 1)[1].lstrip(' ').replace(' ', '_').rstrip()
        report_dict[descrip] = item
        all_dicts[os.path.basename(report).split("_assembly")[0]] = report_dict
df = pd.DataFrame.from_dict(all_dicts, orient='index', dtype=None)
df.to_csv("Masterfile.csv")

# Filename with strain name for tree making
masterfile = pd.read_csv("Masterfile.csv", header=0, index_col=0)
key = masterfile[['Organism_name', 'RefSeq_assembly_accession']]
key.rename(columns={'Organism_name':'name','RefSeq_assembly_accession':'contigs_db_path'},inplace=True)
key['name'].replace(regex=True, inplace=True, to_replace=r'\([^)]*\)',  value=r'')
key['contigs_db_path']=key['contigs_db_path'].astype(str)
key['contigs_db_path']=key['contigs_db_path'].str.replace('.','_')
key['contigs_db_path']=key['contigs_db_path'].astype(str) + '_contigs.db'
key.to_csv('meth-genomes.txt', sep='\t', index=False)
names=key['name]
names.to_csv('meth-names.txt', sep='\t', index=False)
```

Then create the additional layers file that looks like this: 

```
name	methylator
Geobacter_sulfurreducens_PCA_	Y
Pelobacter_carbinolicus_DSM_2380_	N
Geobacter_metallireducens_GS-15_	Y
Pelobacter_propionicus_DSM_2379_	N
Geobacter_uraniireducens_Rf4_	Y
Geobacter_lovleyi_SZ_	N
Geobacter_bemidjiensis_Bem_	Y
Geobacter_daltonii_FRC-32_	Y
Geobacter_sp._M21_	Y
Geobacter_sp._M18_	Y
Geobacter_sulfurreducens_KN400_	Y
Geobacter_pickeringii_	Y
Geoalkalibacter_subterraneus_	N
Desulfuromonas_soudanensis_	Y
Desulfuromonas_sp._DDH964_	Y
Geobacter_anodireducens_	Y
Pelobacter_sp._SFB93_	N
Pelobacter_acetylenicus_	N
```

## Anvi'o Pangenomics Workflow 

First parse the genbank files for importing external gene calls and annotations into Anvi'o using my genbank parser: 

```
#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

# Arguments 

parser = argparse.ArgumentParser(description = "Parse Prokka annotations from a genbank file to add external gene calls and functions to Anvi'o")
parser.add_argument('gbk_file', metavar='GBK', help='Annotation file from Prokka in Genbank format')
parser.add_argument('--gene-calls', default='gene_calls.txt', help='Output: External gene calls (Default: gene_calls.txt)')
parser.add_argument('--annotation', default='gene_annot.txt', help="Output: Functional annotation for external gene calls (Default: gene_annot.txt)")

args = parser.parse_args()

# Input and output files
GBK = args.gbk_file
OUT_CDS = open(args.gene_calls, "w")
OUT_ANNO = open(args.annotation, "w")

# Headers for Anvi'o output
OUT_CDS.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n")
OUT_ANNO.write("gene_callers_id\tsource\taccession\tfunction\te_value\n")

# Gene ID and e-value
gene_id = 1
e_value = "0"

# Parse the GBK file
for record in SeqIO.parse(GBK, "genbank"):
    contig = record.name
    for f in record.features:
        if f.type == 'CDS':
            span, inference, function = (f.location, f.qualifiers["inference"][0], f.qualifiers["product"][0])
            source = inference.split(":")[1]
            version = inference.split(":")[2]
            beg = span.start # biopython already reads in the genbank file to start counting from 0, so don't have to subtract 1
            end = span.end 
            strand = str(span.strand)
            r = strand.replace("-1", "r")
            direction = r.replace("1", "f")
            if (float(beg - end)/float(3)).is_integer() == True:
                partial = str(0)
            else:
                partial = str(1)
            try:
                gene_acc = f.qualifiers["gene"][0]
            except KeyError:
                gene_acc = ""
            # Write out to gene calls and annotation files
            OUT_CDS.write('%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' %(gene_id, contig, beg, end, direction, partial, source, version))
            OUT_ANNO.write('%d\t%s:%s\t%s\t%s\t%s\n' % (gene_id, 'Prokka', source, gene_acc, function, e_value))
            gene_id = gene_id + 1
```

Cleanup the directory with `for file in *.fa; do base=$(basename $file .fa); mv $file $base; done`

Then parse the genbank file: 

```
for file in */*.gbk; do
    f=$(basename $file .gbk);
    python genbank-parser.py $file --gene-calls $f/$f-gene-calls.txt --annotation $f/$f-gene-annot.txt; 
done 
```

## Anvi'o Pangenomics Workflow

```
# Generate contigs databases 
for fasta in */*.fa; do
    f=$(basename $fasta .fa);
    anvi-gen-contigs-database -f $fasta -o $f-contigs.db --external-gene-calls $f/$f-gene-calls.txt -n "$f" 
done

# annotations clean up
for file in */*-gene-annot.txt; do
    f=$(basename $file -gene-annot.txt);
    newname="${f//./_}-gene-annot.txt"
    mv $file $newname;
done

# Import annotations
for contig in *.db; do
    f=$(basename $contig _contigs.db);
    anvi-import-functions -c $contig -i $f-gene-annot.txt;
done

# Run HMMs
for contig in *_contigs.db; do
    anvi-run-hmms -c $contig --num-threads 16;
done

# Cleanup the file names
for file in *.db; do
    f=$(basename $file -contigs.db);
    newname="${f//./_}_contigs.db"
    mv $file $newname;
done

# Generate genomes storage 
anvi-gen-genomes-storage --external-genomes meth-genomes.txt --gene-caller Prodigal -o GEOBACTER-METHS-GENOMES.db

# Run pangenomic analysis
anvi-pan-genome -g GEOBACTER-METHS-GENOMES.db --project-name "Geobacter_Methylators" --output-dir GEOMETHPANG --num-threads 16 --mcl-inflation 2 --sensitive --minbit 0.5 --min-occurrence 2
```

## Working with the additional layers 

Going to layer the data by methylation vs non-methylation activity in the Geobacter area of Deltaproteobacteria.

```
anvi-import-misc-data Geobacter-meth-layers.txt -p GEOMETHPANG/Geobacter_Methylators-PAN.db --target-data-table layers
```

Visualize the pangenome: 

```
anvi-display-pan -p GEOMETHPANG/Geobacter_Methylators-PAN.db -g GEOBACTER-METHS-GENOMES.db
```

Explore the data: 

```
anvi-get-enriched-functions-per-pan-group -p GEOMETHPANG/Geobacter_Methylators-PAN.db -g GEOBACTER-METHS-GENOMES.db --category methylator --annotation-source Prokka:Prodigal -o Geobacter-methylators-enrichments --functional-occurrence-table-output Geobacter-functions-occurrence.txt
```

Since the methyl-mercury gene is annotated as different things or as hypothetical, it's hard to find in the presence/absence matrix. So I'll manually find those locus tags and change the annotation to "mercury methylation" for now to try and make my analysis a bit simpler. I don't think I can do this with `find/sed` because I have to find by the locus tag and then change the annotation. So read in the list of locus tags, all the genbank files, and if it finds the specific locus tag, then change the annotation from hypothetical protein to mercury methylation. 

1. I ran the HMM against the Geobacter pangenome protein files. Got the locus tag
2. Created a locus tag list file
3. Read the locus tag list file and all the genbank files into my script `manual-annotate.py` to change "hypothetical protein" to "mercury methylation" 
4. Output the genbank file, re-run the genbank parser for Anvi'o annotations, and then re-run the genomes storage and pangenomic analysis steps since the contigs databases will be updated with the new annotations

### Sunday 2018-07-29 

I still think the branch of Pelobacter might be too divergent for looking at presence vs absence against the Geobacter methylators. So I will do the analysis with just the next branch up, which is the non-methylating Geobacter and close relative Pelobacter. Then I will look for things enriched in the methylators not in the other two. Not a lot of power with just 2 outgroups. I can then also try to tone down the Archaeal workflow as well. 



