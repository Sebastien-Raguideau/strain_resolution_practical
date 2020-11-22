
# Practical : Strain resolution
##  STRONG - Strain Resolution ON Graphs
STRONG resolves strains on assembly graphs by resolving variants on core COGs using co-occurrence across multiple samples.

###  STRONG : quickstart
#### Installation
STRONG is already installed on your VM, you will find the git repos at : 
`/usr/local/STRONG`
It is a clone of the one available online and installation steps are detailed [there].(https://github.com/chrisquince/STRONG)

**Packages needed**:
Current installation depend entirely on [conda](https://docs.conda.io/en/latest/), by creating an environment where all packages and their dependencies are installed. 

Look at : `/home/student/repos/STRONG/conda_env.yaml`

Bowtie2 is a dependencies, it has been installed, try :

    bowtie2 -h
<details><summary>not working?</summary>
<p>
Try activating the relevant conda environment :

    conda env list
    conda activate STRONG
    bowtie2 -h

</p>
</details>


**Databases**

 - [COG database](ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian) , you will find it installed at
    `/home/student/databases/rpsblast_cog_db`

 - (optional) [GTDB](https://pubmed.ncbi.nlm.nih.gov/30148503/) , used
   with gtdb-tk, (77Gb) takes way more than 20Gb ram on execution.

#### Dataset
[Sharon&al 2012](https://pubmed.ncbi.nlm.nih.gov/22936250/): Infant gut metagenomes, first metagenomic time series
Dataset path : 

    /home/student/datasets/infant_gut

Check the number of reads :

       zcat sample1_R1.fastq.gz |wc -l


STRONG will look for samples in the directory specified by the **data** parameter, inside this directory subdirectories should be present named sample1, ..., sampleN these correspond to different samples. The program will expect sequencing reads present in each subdirectory 'sampleX' with file names 'sampleX_R1.fq.gz' and 'sampleX_R2.fq.gz' for the forward and reverse reads. These are assumed paired, other file formats e.g. not gzipped should also work. In the sample test data used above eight samples are used.

This is a bit unwieldy and we are changing this in the next version of STRONG.


#### Launching STRONG :  config file
To use STRONG, you need to create a config file which is  used to store the parameters of a run. It is divided into  sections with parts corresponding to the different steps of the pipeline. Here is a template config file. 

```
# ------ Samples ------
samples: '*' # specify a list samples to use or '*' to use all samples
data:  ?  # path to data folder

# ------ Resources ------
threads : 8 # single task nb threads

# ----- Annotation database -----
cog_database: ? # COG database

# ----- Binning parameters ------
concoct:
    contig_size: 1000

# ------ Assembly parameters ------ 
read_length: ?
assembly: 
    assembler: spades
    k: [?]
    mem: ?
    threads: 24

# ----- BayesPaths parameters ------
bayespaths:
    nb_strains: 5
    nmf_runs: 1
    max_giter: 1
    min_orf_number_to_merge_bins: 10
    min_orf_number_to_run_a_bin: 10
    percent_unitigs_shared: 0.1

# ----- DESMAN parameters ------
desman:
    execution: 1
    nb_haplotypes: 10
    nb_repeat: 5
    min_cov: 1
```
Let's create a working directory for STRONG:

    mkdir -p "/home/student/Strain_resolution/STRONG_run"
    cd /home/training/Strain_resolution/
And an empty config file

    nano config.yaml

Copy previous template and fill the "?"

What are the consequence of a small or high K value for the assembly? [answer](https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size)

<details><summary>command line help</summary>
<p>

nb of cpu : `lscpu`

ram availlable : `free -h`

length of reads : `zgrep -v "@" -m 1 sample1_R1.fastq.gz |wc -c`

</p>
</details>

##### STRONG command lines
We are now ready to run STRONG, the minimal command needed is : 

    STRONG /home/student/Strain_resolution/STRONG_run --threads 4 --config config.yaml --verbose --dryrun

To launch it for real, try without --dryrun
