#! /bin/bash

# Requirements: ncbi-genome-download, anvi'o, prokka, gffutils and gffparser.py script
# Example with bacterial genomes for downloading and searching

# Gff parser script
wget https://raw.githubusercontent.com/karkman/gff_parser/master/gff_parser.py -O gff_parser.py

# Retrieve all complete, publicly available bacterial and archaeal genomes from NCBI
nohup ncbi-genome-download --format "protein-fasta,fasta,assembly-report" --assembly-level complete bacteria,archaea --parallel 16 &

# Rename RefSeq files based on directory name
for filename in */*.faa; do mv $filename ${filename%/*}/${filename%/*}.faa; done

# Create master FASTA files with genome name included in header, example of bacteria
for filename in bacteria/*/*.faa; do GENNAME=`basename ${filename%.faa}`; sed "s|^>|>${GENNAME} |" $filename; done > all-bac-refseq-prots.faa

# Change multi-line fasta to single-line fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' all-bac-refseq-prots.faa > all-bac-prots-singleline.faa

# Run HMM to find methylators
nohup hmmsearch -E 1e-50 --tblout bac-hgcA-hits.out hgcA.hmm all-bac-refseq-prots.faa &

# From HMM output get genome name and protein tag
awk '{print $1,$19}' bac-hgcA-hits.out > bac-hgcA-hits-list.txt

# grep hits to create fasta of hit and protein
grep -A 1 -w -F -f bac-hgcA-hits-list.txt all-bac-prots-singleline.faa > hgcA-hits-prots.faa
# cleanup the file
sed -i 's/.*--.*//' hgcA-hits-prots.faa
sed -i '/^\s*$/d' hgcA-hits-prots.faa

# Get all directories of only the hits
for dir in $(cat bac-meth-genomes-list.txt); do cp -rf refseq/bacteria/"$dir" METH-NUC/; done

# Cleanup fasta deflines to cooperate wtih Anvi'o
for file in */*.fna; do
    f=$(basename $file .fna);
    anvi-script-reformat-fasta $file -o $f.fa --min-len 0 --simplify-names; 
done

# Run prokka
for fasta in *.fa; do
    N=$(basename $fasta .fa);
    prokka --outdir $N --prefix $N --cpus 15 $fasta;
done

# Prokka for archael genomes
for fasta in *.fa; do 
    N=$(basename $fasta .fa);
    prokka --kingdom Archaea --outdir $N --prefix $N --cpus 15 $fasta;
done

# Parse GFFs for Anvi'o formatting
for file in */*.gff; do python gff_parser.py $file --gene-calls $file-gene-calls.txt --annotation $file-gene-annot.txt --process-all; done 

for file in *.fa; do base=$(basename $file .fa); mv $file $base; done

# Generate contigs database for Anvi'o and import functions from Prokka
for fasta in */*.fa; do
    f=$(basename $fasta .fa);
    anvi-gen-contigs-database -f $fasta -o $f-contigs.db --external-gene-calls $f/$f.gbk-gene-calls.txt -n "$f" 
done
    # Might have to add the --process-all flag because of CDS prokka issues

# Remove periods from genome names
for file in *.db; do
    f=$(basename $file _contigs.db);
    newname="${f//./_}_contigs.db"
    mv $file $newname;
done

# Rename annotations file 
for file in */*.gbk-gene-annot.txt; do
    f=$(basename $file .gbk-gene-annot.txt);
    newname="${f//./_}.gbk-gene-annot.txt"
    mv $file $newname;
done

# Import annotations

for contig in *.db; do
    f=$(basename $contig _contigs.db);
    anvi-import-functions -c $contig -i $f.gbk-gene-annot.txt;
done

# Run HMMs
for contig in *_contigs.db; do
    anvi-run-hmms -c $contig --num-threads 16;
done

# Import the annotations
for contig in *.db; do
    f=$(basename $contig _contigs.db);
    anvi-import-functions -c $contig -i $f/$f.gbk-gene-annot.txt;
done

# Make a genomes storage and pangenome analysis
# Make file of list of genome names and paths to contigs database
echo name$'\t'contigs_db_path > meth-genomes.txt
for contig in *.db; do 
    f=$(basename $contig _contigs.db); 
    echo $f$'\t'$contig >> meth-genomes.txt;
done

# Make list of path to contigs database and names of strains
python meth-filenames-prep.py
 
# Generate a genomes storage
anvi-gen-genomes-storage --external-genomes meth-genomes.txt --gene-caller Prodigal -o METHS-GENOMES.db
# error of having 0 gene calls even though I imported my gene calls externally, possibly because I skipped the process all step for the gff parser? 
# Use the --gene-caller parameter

# Run pangenomic analysis
anvi-pan-genome -g METHS-GENOMES.db --project-name "MethylatorsB" --output-dir METHPANG --num-threads 16 --mcl-inflation 2 --sensitive --minbit 0.5 --min-occurrence 2
# Very distantly related genomes, so can't have high granularity of of clusters, MCL i of 10 creates over 100,000 clusters
# I of 5 creates 22,364 clusters and threw out singletons, try i of 2, seems to be one of the defaults, without letter was 10, A is i of 5, and B is i of 2. Good news is the DIAMOND blastp doesn't have to be rerun every time.
# i of 2 still has a little over 20,000 clusters, so use --enforce-hierarchical-clustering flag to bypass the soft limit of 20,000 gene clusters and hope for the best
anvi-pan-genome -g METHS-GENOMES.db --project-name "MethylatorsC" --output-dir METHPANG --num-threads 16 --mcl-inflation 2 --sensitive --minbit 0.5 --min-occurrence 2 --enforce-hierarchical-clustering 

# Summarizing functional enrichments and presence/absence data
anvi-get-enriched-functions-per-pan-group -p MethylatorsC-PAN.db -g ../METHS-GENOMES.db --category num_gene_clusters --annotation-source Prokka:Prodigal -o meths-pan-enriched-functions-light.txt --functional-occurrence-table-output meths-occurrence-table.txt

# To get a master file of the metadata on hits
python assembly-report-cleaning.py

# Check CDS gene calls in genbank files
for file in */*.gbk
do
cat $file | grep 'CDS' | wc -l
done

# Get genome names from fasta file
awk '{print $1}' cleaned-hgcA-hits.faa | grep '>' | sed 's/[>]//g'

# Genome ID parsed out from locus tag
cat cleaned-hgcA-hits-genome-list.txt | grep 'GCA'| awk -F_ '{print $1 FS $2}' > peat-ids.txt

# Get genome names with matching ID tags
grep -w -F -f peat-ids.txt Peat_taxonomy.txt > peat-methylators-taxonomy.txt






