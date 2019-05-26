#! /bin/bash
# old bash commands for manually checking hgcA

# Create master FASTA files with genome name included in header, for both archaea and bacteria
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