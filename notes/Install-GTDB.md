# Install the GTDB Tool for Classifying MAGs 

These instructions are for manually installing the [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) tool from Ecogenomics. I've created a Docker image of the package and dependencies. However, the way the supplementary database files are downloaded (25 GB total) and has to be pointed to manually in the configuration file, a Docker image installation method for this tool doesn't work quite yet. These steps are for installing the tool on an Ubuntu Linux system. More specifically, this applies to the VM server systems through UW, so sometimes a little more handholding is necessary. 

## Ubuntu Dependencies

```
sudo apt get install make build-essential autoconf zlib1g-dev libgsl2 libgsl-dev libmoose-perl libipc-run-perl bioperl
```

## Dependencies 

```
# Prodigal
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux &&
mv prodigal.linux /usr/local/bin/prodigal && chmod 755 /usr/local/bin/prodigal

# HMMER
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz &&
tar -xvzf hmmer-3.1b2-linux-intel-x86_64.tar.gz &&
cd hmmer-3.1b2-linux-intel-x86_64 && ./configure && make

# pplacer
wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && 
unzip pplacer-linux-v1.1.alpha19.zip &&
cp pplacer-Linux-v1.1.alpha19/pplacer /usr/local/bin/ && 
cp pplacer-Linux-v1.1.alpha19/guppy /usr/local/bin/ && 
cp pplacer-Linux-v1.1.alpha19/rppr /usr/local/bin/

# FastANI
git clone https://github.com/ParBLiSS/FastANI.git &&
cd FastANI &&
./bootstrap.sh &&
./configure --prefix=/usr/local/bin/

# FastTree 
wget http://www.microbesonline.org/fasttree/FastTree.c &&
gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm &&
cp FastTree /usr/local/bin/
```

### GTDB-TK 

Then install the GTDB-tk with `pip`: 
```
pip install gtdbtk
```

Install the necessary datafile by creating a directory to store the database. This will take some time to download and decompress. 
```
wget https://data.ace.uq.edu.au/public/gtdbtk/release_86/gtdbtk_r86_data.tar.gz
tar -xvzf gtdbtk_r86_archived_data.tar.gz
```

The configuration file of GTDB-tk has to be manually edited. The configuration file is in the directory like so depending on your machine/version of python `/home/emcdaniel/anaconda2/lib/python2.7/site-packages/gtdbtk/config`. Copy the configuration file to `cp config_template.py config.py`. 

### Example Run 

```
gtdbtk classify_wf --cpus 14 --genome_dir ./LAKE-TANG-GENOMES --out_dir output --extension fasta
```