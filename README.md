# Expanded Phylogenetic and Metabolic Diversity of Microbial Mercury Methylation

This notebook contains scripts, workflows, and results for analyzing publicly available reference genomes and MAGs that contain the _hgcA_ gene, an indicator of putative mercury methylation. The overall goal of this project is to characterize novel methylators and their metabolic intracacies using environmental metagenomic datasets. 

## Dependencies: 

- [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download)
- HMMer
- [Prokka](https://github.com/tseemann/prokka) 
- Muscle
- [AliView](http://www.ormbunkar.se/aliview/)
- ANICalculator
- FastTree 
- RaxML 
- GTDBtk
- CheckM

## Datasets Analyzed 

- Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system, [Anantharaman et al. 2016](https://www.nature.com/articles/ncomms13219) (Genomes from sulfur group 2018 overlap)
- Expanded diversity of microbial groups that shape the dissimilatory sulfur cycle [Anantharaman et al. 2018](https://www.nature.com/articles/s41396-018-0078-0)
- Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life [Parks et al. 2018](https://www.nature.com/articles/s41564-017-0012-7)
- Genome-centric view of carbon processing in thawing permafrost [Woodcroft et al. 2018](https://www.nature.com/articles/s41586-018-0338-1)
- NCBI Refseq Publicly Available Genomes Accessed June 2018 
- Trout Bog MAGs binned by Sarah Stevens (2015)
- Lake Mendota Hypolimnia Bins collected/binned by Ben Peterson (2018)
- Lake Tanganyika Bins from multiple depths collected by **insert collaborator name I forgot** and binned by Anantharaman/Patricia Tran (2018)

### Downloading Genome and MAG Datasets

The most high quality, experimentally verified, genomes of methylating organisms have been deposited on NCBI. The easiest (and most fun) way to retrieve these genomes is using Kai Blin's [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download). Easy to install and use. On a server with 16 cores, it took under an hour to download all the protein FASTA files of all publicly available bacterial and archaeal genomes and their assembly statistics. Then do some simple formatting for downstream analyses, and the combined protein FASTA file is ready to search with against the custom build _hgcA_ gene HMM profile. 

```
# Retrieve all complete, publicly available bacterial and archaeal genomes from NCBI
nohup ncbi-genome-download --format "protein-fasta,fasta,assembly-report" --assembly-level complete bacteria,archaea --parallel 16 &

# Rename RefSeq files based on directory name
for filename in */*.faa; do mv $filename ${filename%/*}/${filename%/*}.faa; done

# Create master FASTA files with genome name included in header, example of bacteria
for filename in */*.faa; do GENNAME=`basename ${filename%.faa}`; sed "s|^>|>${GENNAME}_|" $filename; done > all-bac-refseq-prots.faa

# Change multi-line fasta to single-line fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' all-bac-refseq-prots.faa > all-bac-prots-singleline.faa
```

This step also works for pulling down all genome bins from Woodcroft et al. 2018, that studied MAGs from a permafrost gradient. The project accession number is PRJNA386568, and on that NCBI Accession page, download the Assembly Details file. From that file, get the assembly id's using `tail -n +3 PRJNA386568_AssemblyDetails.txt | cut -f1 > assembly_ids.txt`, and then supply that list to the ncbi-genome-download tool with `ncbi-genome-download --section genbank -F fasta,protein-fasta --assembly-accessions assembly_ids.txt bacteria,archaea`. The Parks dataset (as far as I know yet), has not been deposited on NCBI, and has to be manually downloaded wtih `wget https://data.ace.uq.edu.au/public/misc_downloads/uba_genomes/uba_bac_prokka.tar.gz `. The size of the compressed folder is 90 GB, and takes about an hour and a half to download, and another hour and a half to decompress. Anantharaman datasets are manually downloaded from ggKbase, for which I downloaded the proteins file for the entire dataset. From what I can tell, selected bins of contigs then have to manually downloaded one by one. All other datasets can pull down the protein FASTA files by genome. Lake Mendota, Trout Bog, and Lake Tanganyika bins are analyzed from unpublished datasets, and were binned by Ben Peterson, Sarah Stevens, and Patricia Tran, all members of the McMahon lab. 

### Searching for the _hgcA_ protein

Using the _hgcA_ HMM profile in the `files/` directory, run the following. I generally run a normal HMMSearch with an E value cutoff of 1e-50, since that is what the Podar 2015 paper did for analyzing short read metagenomic sequences for the _hgcA_ gene. Curation of hits and pruning is described below. 

```
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
```

The above example is done on the bacterial NCBI reference genomes, but generally works the same for all datasets *except* the Anantharaman datasets since all of the genome bins are not downloaded with the dataset, only the protein FASTA file of the entire dataset. Therefore, hits have to manually pulled down if you want the entire genome files. Once this has been done for all dataset, a concatentated FASTA file of the _hgcA_ protein hits can be made for alignment and tree making purposes. Based off of e-value, I usually throw out any hits that hover around the cutoff, because usually a decent hit has an e-value hovering around 300-400. The hits around the cutoff have really long, spurious branches, and might likely actually be the acetyl CoA pathway protein. 

### Functional Annotation 

For functional annotation, and getting files other than the genome nucleotide files (GBK, protein-fasta, etc.), the quickest way to get other file formats and functional annotation is with Prokka. Install Prokka with `conda install -c conda-forge -c bioconda prokka` if using Anaconda. For a given set of bacterial/archaeal genomes, run the following: 

```
# Bacteria
for fasta in *.fa; do
    N=$(basename $fasta .contigs.fa);
    prokka --outdir $N --prefix $N --cpus 15 $fasta --centre X --compliant;
done

# Archaea
for fasta in *.fa; do
    N=$(basename $fasta contigs.fa);
    prokka --kingdom Archaea --outdir $N --prefix $N --cpus 15 $fasta;
done
```

Therefore your bacterial and archaeal bins per dataset will need to be split up accordingly. 

### Align Protein Hits 

Align the _hgcA_ protein hits with Mafft `mafft hits.faa > hits.aln`, or with Muscle `muscle -in hits.faa -out hits.aln`. Using AliView, import the alignment files and manually inspect alignments for anomolies or specific conserved domains of interest. 

### Dereplicate Bins from Same Datasets and Find Related Non-Methylators

To schedule mass ANI comparisons on CHTC, we will use Sarah Stevens' [DAG pipeline](https://github.com/sstevens2/ani_compare_dag). This can be used for comparing the methylating bins wihtin datasets to see nucleotide similarity and dereplicate where the original authors did not, and also find genetically similar organisms within the dataset or others of non-methylating organisms. Clone the repo, and change the `group.sub` script to where the ANIcalculator is installed. If you want comparisons of all bins within a directory, write that directory to `groupslist.txt`. If wanting to compare between two directories, in each directory list the genomes and end in `_genome_list.txt`. Then put both the lists and sets of genomes into one directory, such as `all-genomes/`. Then in `groupslist.txt`, list the `all-genomes/` directory, and this will make comparisons against genomes in the directories, and not all-v-al within a directory. 

### Make Phylogenetic Tree of _hgcA_ Protein Hits 

Using FastTree, make a tree of just the _hgcA_ protein within all the genomes that had the hit. Install FastTree: 

```
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
mv FastTree /usr/local/bin
```
Usage is `FastTree prot-file.faa > tree.tre`. Pretty-fy it in iTOL based upon manual classifications. 

### Classification 

### Quality Statistics 

### Metabolic Pathway Prediction 
- Looking @ organisms that haven't been characterized as methylators before, really make sure the HMM hit is good, and characterize their metabolism. Identify closely related organisms in the dataset with ANI comparison, can do genome presence/absence studies to figure out if something weird in the genome beside hgcA

