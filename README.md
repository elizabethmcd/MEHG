# Expanded Phylogenetic and Metabolic Diversity of Microbial Mercury Methylation

This respository contains scripts, workflows, and results for analyzing publicly available reference genomes and MAGs that contain the _hgcA_ marker, an indicator of putative mercury methylation. The overall goal of this project is to characterize novel methylators and their metabolic intracacies using environmental metagenomic datasets. 
Analyses were performed for the mSystems publication "Expanded Phylogenetic Diversity and Metabolic Flexibility of Mercury-Methylating Microorganisms":

**McDaniel E.A.**, Peterson B., Stevens S.L.R., Tran P.Q., Anantharaman K., McMahon K.D. **Expanded Phylogenetic Diversity and Metabolic Flexibility of Microbial Mercury Methylation.** mSystems. Aug 2020, 5 (4) e00299-20; DOI: 10.1128/mSystems.00299-20 [[preprint]](https://www.biorxiv.org/content/10.1101/2020.01.16.909358v1) [[publication]](https://msystems.asm.org/content/5/4/e00299-20)

--------------------------------------------------------------

## Dependencies: 

- [genomes-MAGs-database](https://github.com/elizabethmcd/genomes-MAGs-database)
- [NCBI Genome Download Tool](https://github.com/kblin/ncbi-genome-download)
- [Prokka](https://github.com/tseemann/prokka) 
- [AliView](http://www.ormbunkar.se/aliview/)
- [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM)
- [GTDBtk](http://gtdb.ecogenomic.org/)
- [Kallisto](https://pachterlab.github.io/kallisto/)
- [KofamKOALA](https://www.genome.jp/tools/kofamkoala/)


### Downloading Publicly Available Reference Genomes and MAGs from Environmental Datasets

To search all publicly available Genbank genomes as of August 2019, I created a high-throughput pipeline for use on UW-Madison's CHTC HTCondor system found in the [genomes-MAGs-database](https://github.com/elizabethmcd/genomes-MAGs-database) repository. This pipeline will search all ~200,000 genomes with a given HMM and bring back any corresponding hits. Genomes with the hgcA marker were then downloaded with the `ncbi-genome-download` toolkit, and reformatted and functionally annotated with Prokka using the `reformat-annotate.sh` script found in the [genomes-MAGs-database](https://github.com/elizabethmcd/genomes-MAGs-database) repository. Freshwater lake MAGs from Trout Bog, Lake Mendota, Lake Tanganyika, and the Jones 2019 ISMEJ Minnesota lakes genomes were also added. All genome information and sample/environmental medadata can be found in the `files/` folder.  

### Searching for the _hgcA_ protein

Using the constructed _hgcA_ HMM profile, use the [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM) package for screening genomes for the _hgcA_ marker. The script `search-single-marker.py` will search for the marker, create an alignment, and build a phylogeny of the single marker all in one step. 

### Full-Genome Phylogenies of Archaeal/Bacterial Methylators

Use the [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM) package for making concatenated riobosomal protein trees of archael and bacterial methylators. These phylogenies will need to be made separately, as the ribosomal markers for archaea/bacteria are different. 

### Pathway Characterizations 

Create a high level summary of metabolic characteristics with the `metabolisHMM` package using custom and curated HMM markers. Additionally, perform specific metabolic reconstructions with the KofamKOALA package to get annotations through the KEGG database with curated HMMs to have more confidence in the annotations. 

```
./exec_annotation ../Actinobacteria-metabolism/bacteria23279.faa -p profiles/ -k ko_list -o ../Actinobacteria-metabolism/bacteria23279-kofam-annots.txt --cpu 10 &
```

### Investigating Specific Methylators

For this paper, I specifically honed in on the putative Actinobacterial methylators, because they are somewhat weird and they haven't been identified as mercury methylators before. They are preliminarily classified as Coriobacteriia, Thermoleophilia, and some weird UBA classes. To download all genomes of a specific taxonomic type from NCBI: 

1. Go to the [NCBI Taxonomy Database](https://www.ncbi.nlm.nih.gov/taxonomy)
2. Type in your taxonomy identfier of interest, such as Coriobacteriia
3. Click on the taxonomy names twice
4. In the taxonomy brower, look at the table of Entrez records, and click the number in the Assembly field.
5. At the top of the page, click the "Send to" dropdown button
6. Click ID Table (text) and order however you want

Once you have saved the assembly-ids file, modify it to only get the list of genbank accession numbers by: `awk '{print $1}' coriobacteriia-class-accessions.txt | tail -n +2 > coriobacteriia-ids.txt`. Then feed this into `ncbi-genome-download` with `nohup ncbi-genome-download --section genbank --assembly-accessions coriobacteriia-ids.txt --format "fasta" -m metadata.txt bacteria &`. 

Then you have the genomes and can make trees with your MAGs, metabolic comparsions etc.

### Transcription of Putative Methylators in a Permafrost System 

Woodcroft et al. 2018 analyzed metatranscriptomes along a perfmafrost thawing gradient. I identified ~100 methylators in their system and want to understand the expression profiles of those organisms at different depths. 

1. Create a concatenated FASTA file of all predicted ORFs for methylating genomes from this system
2. Build a genome index with kallisto: `kallisto index -i ebpr-orfs fastafile`
3. Pseudoalign and quantify with kallisto: `kallisto quant -i index -o outdir fastqfiles`

The last step can be queued with the list of each paired end file with the sample name to create a directory with mapping results for each sample depth. 
