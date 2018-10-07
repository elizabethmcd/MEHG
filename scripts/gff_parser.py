#!/usr/bin/env python

## Antti Karkman
## University of Gothenburg
## antti.karkman@gmail.com
## 2017

import gffutils
import argparse

#parse the arguments
parser = argparse.ArgumentParser(description="""Parse Prokka annotated genome/metagenome to add external gene calls and functions to anvi'o. 
                    Input annotation in GFF3 format, outputs are tab-delimited text files, one for gene calls and one for annotations""")
parser.add_argument('gff_file', metavar='GFF3', help='Annotation file from Prokka in GFF3 format')
parser.add_argument('--gene-calls', default='gene_calls.txt', help='Output: External gene calls (Default: gene_calls.txt)')
parser.add_argument('--annotation', default='gene_annot.txt', help='Output: Functional annotation for external gene calls (Default: gene_annot.txt)')
parser.add_argument('--process-all', default=False, action="store_true", help="Prodigal returns anything it finds, including tRNAs and\
                    other genetic structures that are not nessecarily translatable, and can cause downstream issues especially if you would like to\
                    use your annotations in contigs databases for pangenomic analyses. As a precation this script only recovers open reading frames\
                    identifie dby Prodigal. Using this flag, you can import everything.")

args = parser.parse_args()

#Input file and output files
GFF = args.gff_file
OUT_CDS = open(args.gene_calls, 'w')
OUT_ANNO = open(args.annotation, 'w')

#load in the GFF3 file
db = gffutils.create_db(GFF, ':memory:')

#Print headers for anvi'o
OUT_CDS.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n")
OUT_ANNO.write("gene_callers_id\tsource\taccession\tfunction\te_value\n")

#running gene ID and a trumped-up e-value for the gene calls.
gene_id = 1
e_value = "0"

#parse the GFF3 file and write results to output files
for feature in db.all_features():
    # determine source
    source, version = feature.source.split(':')

    # move on if not Prodigal, unless the user wants it badly
    if not args.process_all:
        if source != 'Prodigal':
            continue

    start = feature.start - 1
    stop = feature.stop

    if (float(start - stop)/float(3)).is_integer() == True:
        partial = str(0)
    else:
        partial = str(1)

    try:
        gene_acc = feature.attributes['gene'][0]
    except KeyError:
        gene_acc = ""

    try:
        product = feature.attributes['product'][0]
    except KeyError:
        product = feature.attributes['note'][0]

    # skip if hypotethical proiten:
    if product == 'hypothetical protein':
        product = ""
        gene_acc = ""

    # determine direction
    if feature.featuretype=='repeat_region':
        direction='f'
    else:
        if feature.strand==' + ':
            direction='f'
        else:
            direction='r'
    
    OUT_CDS.write('%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n' %(gene_id, feature.seqid, start, stop, direction, partial, source, version))
    OUT_ANNO.write('%d\t%s:%s\t%s\t%s\t%s\n' % (gene_id, 'Prokka', source, gene_acc, product, e_value))

    gene_id = gene_id + 1