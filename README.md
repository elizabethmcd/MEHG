# Expanded Phylogenetic and Metabolic Diversity of Microbial Mercury Methylation

This notebook contains scripts, workflows, and results for analyzing publicly available reference genomes and MAGs that contain the _hgcA_ marker, an indicator of putative mercury methylation. The overall goal of this project is to characterize novel methylators and their metabolic intracacies using environmental metagenomic datasets. 

## Dependencies: 

- [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download)
- [Prokka](https://github.com/tseemann/prokka) 
- [HMMer](http://hmmer.org/)
- Muscle
- [AliView](http://www.ormbunkar.se/aliview/)
- FastTree 
- RaxML 
- ANICalculator
- GTDBtk
- CheckM

## Datasets Analyzed 

- Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system, [Anantharaman et al. 2016](https://www.nature.com/articles/ncomms13219)
- Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life [Parks et al. 2018](https://www.nature.com/articles/s41564-017-0012-7)
- Genome-centric view of carbon processing in thawing permafrost [Woodcroft et al. 2018](https://www.nature.com/articles/s41586-018-0338-1)
- NCBI Genbank Publicly Available Genomes 
- Trout Bog MAGs binned by Sarah Stevens (2015)
- Lake Mendota Hypolimnia Bins collected/binned by Ben Peterson (2018)
- Lake Tanganyika Bins from multiple depths collected by Peter McIntyre and binned by Anantharaman/Patricia Tran (2018)

All publicly available datasets were accessed through genbank in November of 2018. 

### Downloading Publicly Available Reference Genomes and MAGs from Environmental Datasets

The easiest way to retrieve publicly available genomes is using Kai Blin's [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download). Easy to install and use. On a server with 16 cores, it took under an hour to download all the protein FASTA files of all publicly available bacterial and archaeal genomes and their assembly statistics. Then do some simple formatting for downstream analyses, and the combined protein FASTA file is ready to search with against the custom build _hgcA_ gene HMM profile. 

```
# Retrieve all complete, publicly available bacterial and archaeal genomes from NCBI
nohup ncbi-genome-download -s genbank --format "protein-fasta,fasta,assembly-report" --assembly-level complete bacteria,archaea --parallel 16 &

# Rename Genbank files based on directory name
for filename in */*.fna; do mv $filename ${filename%/*}/${filename%/*}.fna; done
```

This step also works for pulling down all genome bins from Woodcroft et al. 2018, Parks et al. 2017, and Anantharaman et al. 2016. All of their genomes have been deposited as Bioprojects, the numbers are:  

- UBA genomes: PRJNA348753
- Aquifer genomes: PRJNA288027
- Permafrost genomes: PRJNA386568

For each of these datasets, download the assembly details file, and concatenate them together. Get the assembly ID's for downloading them all at once. Get the assembly ID's using `tail -n +3 genbank-combined-assembly-details.txt | cut -f1 > combined-assembly-ids.txt`. Supply that list to ncbi-genome-download with `nohup ncbi-genome-download -s genbank -F fasta -A combined-assembly-ids.txt all --parallel 16 &`. Create taxonomy files for each of the datasets as well for creating a metadata file for all identified genomes. Repeat the processing steps. 

Lake Mendota, Trout Bog, and Lake Tanganyika bins are analyzed from unpublished datasets, and were binned by Ben Peterson, Sarah Stevens, and Patricia Tran, all members of the McMahon lab. Throughout the analysis, keep archaea and bacterial genomes separate for annotation purposes. 

### Metadata and Reorganization

For all available genomes, I will organize the metadata and rename all the genomes in a uniform way to make all downstream processes easier. The genomes are split into archaea and bacteria folders, and the naming scheme will be `archaea0001` and `bacteria0001` and so on and so forth. First, we will create a dataframe of the existing names and the corresponding taxonomy. For the reference genomes from genbank, all metadata are in `_assembly_details.txt` files. The script `genbank-reference-metadata-prep.py` will take a directory of the assembly detail files and create two files, one with all the metadata for each genome, and one with the genome name and corresponding taxonomy. The UBA, permafrost, and aquifer datasets' taxonomy can be accessed from the assembly information BioProject files that were downloaded to get the accession numbers. 

To keep all genomes named uniformly, create a list of all the FNA files of all genome bins in the archaea and bacteria folders. Use the scripts `create-new-genome-names.py` to create new, sequential names based on the order of the FNA files. Then to get the matching metadata of the original genome names to the new names, run `all-genome-metadata.py`, which will also give classifications for the new names, and all genomes in the database. To change the names in the two lists: 

```
# renaming files to organized, universal numbered names
awk '{print $4}' bac-genomes-new-names.txt > new-names.txt
for file in *.fna; do read line; mv -v "${file}" "${line}"; done < new-names.txt
```

Now all the genomes in the database (almost 25,000 complete reference genomes and publicly available MAGs plus our own lakes datasets), are named the same way and easier to find, link back to in the future with functional annotations.

### Gene Calling

To ensure genes are called the same way among the genome database and annotated the same way, we will use Prodigal as a first pass to get the proteins for all genomes. If you've already installed Prokka, Prodigal is also already installed. Functional annotation will only be performed for the identified set of methylators and related genomes for metabolic and gene presence analysis. Install Prokka with `conda install -c conda-forge -c bioconda prokka` if using Anaconda. To get predicted genes and their proteins, for all genomes run: 

```
for file in *.fna; do 
    N=$(basename $file .fna);
    prodigal -i $file -a $N.faa;
done
```

### Searching for the _hgcA_ protein

Using the _hgcA_ HMM profile in the `files/` directory, run the following. I generally run a normal HMMSearch with an E value cutoff of 1e-50, since that is what the Podar 2015 paper did for analyzing short read metagenomic sequences for the _hgcA_ gene. Curation of hits and pruning is described below. 

```
# Create master FASTA files with genome name included in header, for both archaea and bacteria
for filename in *.faa; do GENNAME=`basename ${filename%.faa}`; sed "s|^>|>${GENNAME}_|" $filename; done > all-arch-prots.faa

for filename in *.faa; do GENNAME=`basename ${filename%.faa}`; sed "s|^>|>${GENNAME}_|" $filename; done > all-bac-prots.faa

# Change multi-line fasta to single-line fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' all-bac-prots.faa > all-bac-prots-singleline.faa

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' all-arch-prots.faa > all-arch-prots-singleline.faa

# Run HMM to find methylators, run separately
nohup hmmsearch -E 1e-50 --tblout arch-hgcA-hits.out hgcA.hmm all-arch-prots-singleline.faa &

nohup hmmsearch -E 1e-50 --tblout bac-hgcA-hits.out hgcA.hmm all-bac-prots-singleline.faa &

# From HMM output get genome name and protein tag
awk '{print $1}' arch-hgcA-hits.out > arch-hgcA-hits-list.txt
awk '{print $1}' bac-hgcA-hits.out | awk -F "_" '{print $1}' > bac-hgcA-genomes-list.txt

# Combined list of genomes and locus tags for combining with taxonomy information
awk '{print $1}' arch-hgcA-hits.out | awk -F "_" '{print $1"\t"$1"_"$2"_"$3}' > arch-hgcA-list-genomes-hits.txt
awk '{print $1}' bac-hgcA-hits.out | awk -F "_" '{print $1"\t"$1"_"$2"_"$3}' > bac-hgcA-list-genomes-hits.txt

# grep hits to create fasta of hit and protein
for hit in $(cat arch-hgcA-hits-list.txt); do grep -A 1 $hit all-arch-prots-singleline.faa; done > arch-hgcA-hits-prots.faa
for hit in $(cat bac-hgcA-hits-list.txt); do grep -A 1 $hit all-bac-prots-singleline.faa; done > bac-hgcA-hits-prots.faa

# cleanup the file
sed -i 's/.*--.*//' hgcA-hits-prots.faa
sed -i '/^\s*$/d' hgcA-hits-prots.faa

# Get all directories of only the hits
for hit in $(cat arch-hgcA-genomes-list.txt); do cp -rf "$hit".fna destination/; done
for hit in $(cat bac-hgcA-genomes-list-qced.txt); do cp -rf "$hit" desination/; done
```

Then make a concatenated protein file of all archaeal and bacterial hits for making a full _hgcA_ tree. When making the list of hits, remove any that are not above the 1e-50 e-value and the 300 score. Some will hover around 1e-50 but not have a good score, so remove those to have a stringent cutoff. 

### Functional Annotation 

To get functional annotations associated with the proteins, and all other file formats, run Prokka. 
```
# Bacteria
for fasta in *.fna; do
    N=$(basename $fasta .fna);
    prokka --outdir $N --prefix $N --cpus 15 $fasta --centre X --compliant;
done

# Archaea
for fasta in *.fna; do
    N=$(basename $fasta .fna);
    prokka --kingdom Archaea --outdir $N --prefix $N --cpus 15 --centre X --compliant $fasta;
done
```

### Align Protein Hits 

Align the _hgcA_ protein hits with Mafft `mafft hits.faa > hits.aln`, or with Muscle `muscle -in hits.faa -out hits.aln`. Using AliView, import the alignment files and manually inspect alignments for anomolies or specific conserved domains of interest. 

### Make Phylogenetic Tree of _hgcA_ Protein Hits 

Using FastTree, make a tree of just the _hgcA_ protein within all the genomes that had the hit. Install FastTree: 

```
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
mv FastTree /usr/local/bin
```
Usage is `FastTree aln-file.aln > tree.tre`. Pretty-fy it in iTOL based upon manual classifications. 

### Dereplicate Bins from Same Datasets and Find Related Non-Methylators

To schedule mass ANI comparisons on CHTC, we will use Sarah Stevens' [DAG pipeline](https://github.com/sstevens2/ani_compare_dag). This can be used for comparing the methylating bins wihtin datasets to see nucleotide similarity and dereplicate where the original authors did not, and also find genetically similar organisms within the dataset or others of non-methylating organisms. Clone the repo, and change the `group.sub` script to where the ANIcalculator is installed. If you want comparisons of all bins within a directory, write that directory to `groupslist.txt`. If wanting to compare between two directories, in each directory list the genomes and end in `_genome_list.txt`. Then put both the lists and sets of genomes into one directory, such as `all-genomes/`. Then in `groupslist.txt`, list the `all-genomes/` directory, and this will make comparisons against genomes in the directories, and not all-v-al within a directory. 

### Full-Genome Phylogenies of Methylators

Using the `metabolisHMM` package, run HMMs against the ribosomal protein set, and create a phylogenetic tree with RaxML. Alternatively, run all genomes through the GTDB-tk and create a phylogenetic tree of the single copy markers. 

### Wood-Ljungdhal Pathway Characterization 

The HgcA protein is most similar to the acetyl CoA synthase subunits of the Wood-Ljungdhal pathway. To characterize presence of this protein and the subunits (gamma, beta, alpha), I've pulled down these protein subunits from the model acetogen _Acetobacterium woodii_. To perform the BLAST run, make a BLAST database of the three protein subunits: 

```
makeblastdb -dbtype prot -in acetobacterium-wlj-acetyl-coA-synthase-proteins.fa -input_type faste -parse_sequid s-out WLJ-acSynth.db
```

Concatenate all bacterial proteins to include the genome name in the output with:

```
for filename in *.faa; do GENNAME=`basename ${filename%.faa}`; sed "s|^>|>${GENNAME}_|" $filename; done > all-bac-proteins.faa
```

Then run BLAST with: 

```
#! /bin/bash 

blastp -query all-bac-proteins.faa -db WLJ-acSynth.db -out bac-WLJ-blast.out -outfmt 11 -max_target_seqs 5

blast_formatter -archive bac-WLJ-blast.out -outfmt "7 qacc sacc evalue qstart qend sstart send pident qcovs qcovhsp" -out WLJ-blast-formatted.out

awk ' ($3 <=0.00001) && ($10 >=85) ' WLJ-blast-formatted.out > WLJ-blast-top-hits.txt

```

Additionally there are curated HMM profiles for key steps in the Wood-Ljungdahl pathway. I downloaded the TIGRFAMs version 15 on 2018-12-19. The existing TIGRFAMs for the acetyl-CoA pathway in metabolisHMM are those for formaldehyde dehydrogenase, in which CO2 is converted to formaldehyde. Those codes are TIGR01591, TIGR01582, TIGR01583.  The running hypothesis is that organisms that have this still may not perform the full WLJ pathway, as formate can be used for other things. The real rate limiting step is the incorporation of THF in bacteria to later form acetyl CoA from the CO dehydrogenase/acetyl coA synthase. 

For incorporation of THF, a ligase is needed to added the cofactor to formate to form THF-CHO. This is the 	5-formyltetrahydrofolate cyclo-ligase with code TIGR02727.

### Metabolic Markers

Create a high level summary of metabolic characteristics with the `metabolisHMM` package using custom and curated HMM markers. The `metabolisHMM` package is still under active [development](https://github.com/elizabethmcd/metabolisHMM). 
