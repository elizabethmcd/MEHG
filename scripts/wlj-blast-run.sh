#! /bin/bash 

for file in *.faa; do
    N=$(basename $file .faa);
    blastp -query $file -db WLJ-acSynth.db -out $N-blast.out -outfmt 11 -max_target_seqs 5; 
done

for file in *-blast.out; do
    N=$(basename $file -blast.out);
    blast_formatter -archive $file -outfmt "7 qacc sacc evalue qstart qend sstart send pident qcovs qcovhsp" -out $N-formatted.out;
done

for file in *-formatted.out; do
    N=$(basename $file -formatted.out);
    awk ' ($3 <=0.00001) && ($10 >=75) ' $file > $N-top-hits.txt;
done

