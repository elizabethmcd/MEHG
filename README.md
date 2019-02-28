# Expanded Phylogenetic and Metabolic Diversity of Microbial Mercury Methylation

This notebook contains scripts, workflows, and results for analyzing publicly available reference genomes and MAGs that contain the _hgcA_ marker, an indicator of putative mercury methylation. The overall goal of this project is to characterize novel methylators and their metabolic intracacies using environmental metagenomic datasets. 

## Dependencies: 

- [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download)
- [Prokka](https://github.com/tseemann/prokka) 
- [HMMer](http://hmmer.org/)
- [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/)
- [AliView](http://www.ormbunkar.se/aliview/)
- [FastTree](http://microbesonline.org/fasttree/)
- RaxML 
- [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM)
- GTDBtk
- CheckM
- Kallisto
- DESeq2

## Datasets Analyzed 

- Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system, [Anantharaman et al. 2016](https://www.nature.com/articles/ncomms13219)
- Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life [Parks et al. 2018](https://www.nature.com/articles/s41564-017-0012-7)
- Genome-centric view of carbon processing in thawing permafrost [Woodcroft et al. 2018]
- Molecular evidence for novel mercury methylating microorganisms in sulfate-impacted lakes [Jones et al. 2019](https://www.nature.com/articles/s41396-019-0376-1?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+ismej%2Frss%2Fcurrent+%28The+ISME+Journal+-+Current%29)(https://www.nature.com/articles/s41586-018-0338-1)
- NCBI Genbank Publicly Available (Complete) Genomes 
- Trout Bog MAGs binned by Sarah Stevens (2015)
- Lake Mendota Hypolimnia Bins collected/binned by Ben Peterson (2018)
- Lake Tanganyika Bins from multiple depths collected by Peter McIntyre and binned by Anantharaman/Patricia Tran (2018)

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

Jones et al. 2019 genomes were downloaded from JGI using the GOLD ID Gs0130353. Lake Mendota, Trout Bog, and Lake Tanganyika bins are analyzed from unpublished datasets, and were binned by Ben Peterson, Sarah Stevens, and Patricia Tran, all members of the McMahon lab. Throughout the analysis, keep archaea and bacterial genomes separate for annotation purposes. 

### Metadata and Reorganization

For all available genomes, I will organize the metadata and rename all the genomes in a uniform way to make all downstream processes easier. The genomes are split into archaea and bacteria folders, and the naming scheme will be `archaea0001` and `bacteria0001` and so on and so forth. First, we will create a dataframe of the existing names and the corresponding taxonomy. For the reference genomes from genbank, all metadata are in `_assembly_details.txt` files. The script `genbank-reference-metadata-prep.py` will take a directory of the assembly detail files and create two files, one with all the metadata for each genome, and one with the genome name and corresponding taxonomy. The UBA, permafrost, and aquifer datasets' taxonomy can be accessed from the assembly information BioProject files that were downloaded to get the accession numbers. 

To keep all genomes named uniformly, create a list of all the FNA files of all genome bins in the archaea and bacteria folders. Use the scripts `create-new-genome-names.py` to create new, sequential names based on the order of the FNA files. Then to get the matching metadata of the original genome names to the new names, run `all-genome-metadata.py`, which will also give classifications for the new names, and all genomes in the database. To change the names in the two lists: 

```
# renaming files to organized, universal numbered names
awk '{print $4}' bac-genomes-new-names.txt > new-names.txt
for file in *.fna; do read line; mv -v "${file}" "${line}"; done < new-names.txt
```

Now all the genomes in the database (almost 25,000 complete reference genomes and publicly available MAGs plus our own lakes datasets), are named the same way and easier to find, link back to in the future with functional annotations. However, the additional 19 genomes from Jones et al. 2019 have their JGI numbers for now in all names because I added them at the last minute. At the end of downloading genomes and doing all of the metadata reorganization, I realized that there is a flag in `ncbi-genome-download` to get the metadata of any downloaded set of genomes, which is `-m METADATA-TABLE.txt` added onto the download command. 

### Gene Calling

To ensure genes are called the same way among the genome database and annotated the same way, we will use Prodigal as a first pass to get the proteins for all genomes. If you've already installed Prokka, Prodigal is also already installed. Functional annotation will only be performed for the identified set of methylators and related genomes for metabolic and gene presence analysis. Install Prokka with `conda install -c conda-forge -c bioconda prokka` if using Anaconda. To get predicted genes and their proteins, for all genomes run: 

```
for file in *.fna; do 
    N=$(basename $file .fna);
    prodigal -i $file -a $N.faa;
done
```

### Searching for the _hgcA_ protein

Using the constructed _hgcA_ HMM profile, use the [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM) package for screening genomes for the _hgcA_ marker. The script `search-single-marker.py` will search for the marker, create an alignment, and build a phylogeny of the single marker all in one step, and therefore I don't have to use the manual steps I used to do that took forever. 

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

### Full-Genome Phylogenies of Methylators in Full Archaeal/Bacterial Trees

Pull down Refseq representative genomes to populate a full genome phylogeny for archaeal and bacterial methylators. 

```
nohup ncbi-genome-download --format protein-fasta -m archaea-metadata.txt -R representative archaea &
nohup ncbi-genome-download --format protein-fasta -m bacteria-metadata.txt -R representative bacteria --parallel 10 &
``` 

Then use the [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM) package for making these trees with the specific set of archaeal/bacterial genomes. 
 
### Pathway Characterizations 

Create a high level summary of metabolic characteristics with the `metabolisHMM` package using custom and curated HMM markers. The `metabolisHMM` package is still under active [development](https://github.com/elizabethmcd/metabolisHMM). With this workflow and downloading KEGG groups, I no longer use BLAST for pathway/genes characterization. 

### Transcription of Putative Methylators in a Permafrost System 
