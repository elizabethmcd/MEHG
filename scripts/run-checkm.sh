#! /bin/bash

##########################
# Run checkM on bins on WEI/GLBRC server
##########################
export PATH=/opt/bifxapps/bin:$PATH
export PYTHONPATH=/opt/bifxapps/python/lib/python2.7/site-packages/
export PATH=/opt/bifxapps/hmmer/bin:$PATH
export PATH=/home/GLBRCORG/pcamejo/bin/pplacer:$PATH
export PATH=/opt/bifxapps/prodigal-2.6.3/:$PATH
cd /home/GLBRCORG/emcdaniel/mehg/all-meths/
if [ -d checkM ]; then
 echo “Removing old checkM folder”
 rm -r checkM
fi
mkdir checkM
checkm lineage_wf \
       -x .fna \
       -t 16 \
       nucs \
       checkM
checkm qa checkM/lineage.ms \
         checkM \
         -o 2 \
         -f checkM/checkm.out \
         --tab_table
awk -F ‘\t’ \
   -v OFS=‘,’ \
   ‘{ print $1,$6,$7,$8,$9,$11,$13,$15,$17,$19,$23 }’ \
   checkM/checkm.out \
   > checkM/checkM_stats.tsv