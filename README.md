
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
To use STRONG, you need to create a config file which is  used to store the parameters of a run. It is divided into  sections with parts corresponding to the different steps of the pipeline.

Here is a template config file. 

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

<details><summary> All config file options</summary>
<p>

### ------ Samples ------
**samples**:  specify a list samples to use or \"\*\" to selects all samples
**data:**   path to data folder

### ------ Resources ------
**threads :** 8  single task nb threads

### ----- Annotation database -----
**cog_database:**  # COG database

### ----- Binning parameters ------
**concoct**:
&nbsp;&nbsp;&nbsp;&nbsp;**contig_size**: mininum contig length for the CONCOCT binning defaults to 1000 
&nbsp;&nbsp;&nbsp;&nbsp;**fragment_size**: size at which contigs are fragmented for CONCOCT inputs defaults to 10000
&nbsp;&nbsp;&nbsp;&nbsp;**bin_multiplier**: the program calculates the median SCG number and then multiplies 
by this value to get the initial bin number for CONCOCT, defaults to 3 
&nbsp;&nbsp;&nbsp;&nbsp;**bin_max**: maximum initial bin number for CONCOCT defaults to 2000 reduce this to speed up the CONCOCT binning

***mag_quality_threshold***: fraction of SCGs in single-copy for a bin to be considered a MAG, should be set between 0 and 1, defaults to 0.75, a higher value will give higher quality MAGs

### ------ Assembly parameters ------ 
**read_length**:  read length used for sequencing 
**assembly:** 
&nbsp;&nbsp;&nbsp;&nbsp;    **assembler**: program for coassembly currently only metaSPAdes is supported specify as ***spades*** which is also the default
&nbsp;&nbsp;&nbsp;&nbsp;    **k**:  kmer length for assembly 77 is a good choice for 150 bp reads. It is possible to use a list of 
kmers here but they should all be odd so for instance '[33,55,77]' defaults to [21, 33, 55]
&nbsp;&nbsp;&nbsp;&nbsp;    ***mem***: This is the maximum memory allocated to metaSPAdes in Gb defaults to 120 
for complex data sets:
&nbsp;&nbsp;&nbsp;&nbsp;    ***threads***: The number of threads used by metaSPAdes defaults to 16



### ----- BayesPaths parameters ------
**bayespaths:**
&nbsp;&nbsp;&nbsp;&nbsp; ***nb_strains****: initial strain number in a MAG prior to automatic relevance determination, this is the maximum number that can be resolved per MAG, defaults to 16
&nbsp;&nbsp;&nbsp;&nbsp; ***nmf_runs***: number of initialisation using NMF 
&nbsp;&nbsp;&nbsp;&nbsp; ***min_orf_number_to_merge_bins***: The number of overlapping ORFs that triggers bin merging defaults to 10
&nbsp;&nbsp;&nbsp;&nbsp; **percent_unitigs_shared***: Fraction of unitigs shared that cause ORFs to be flagged as overlapping defaults to 0.1
&nbsp;&nbsp;&nbsp;&nbsp; ***min_orf_number_to_run_a_bin***:  After removing shared orfs, what minimum number of orfs is needed to run BayesPath 
&nbsp;&nbsp;&nbsp;&nbsp;  ***max_giter***: number of iterations of SCG filtering defaults to 4
   

### ----- DESMAN parameters ------
**desman:**
&nbsp;&nbsp;&nbsp;&nbsp; ***execution***: determines whether the DESMAN steps will be run 0 or 1, defaults to 0
&nbsp;&nbsp;&nbsp;&nbsp; ***nb_haplotypes***: maximum number of DESMAN haplotypes defaults to 10
&nbsp;&nbsp;&nbsp;&nbsp; ***nb_repeat***: repeats for the Gibbs sampler per haplotype defaults to 5
&nbsp;&nbsp;&nbsp;&nbsp; ***min_cov***: minimum coverage for a sample to be used defaults to 1

</p>
</details>

Let's create a working directory for STRONG:

    mkdir -p "/home/student/Strain_resolution/STRONG_run"
    cd /home/training/Strain_resolution/
And an empty config file

    nano config.yaml

Copy previous template and fill the "?"

What are the consequence of a small or high K value for the assembly? How should we chose K? [answer](https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size)

<details><summary>command line help</summary>
<p>

nb of cpu : `lscpu`

ram availlable : `free -h`

length of reads : `zgrep -v "@" -m 1 sample1_R1.fastq.gz |wc -c`

</p>
</details>

##### STRONG command lines
STRONG should be run from within the STRONG repository, however you can symlink it into your path.

    ll $(which STRONG)

We are now ready to run STRONG, let's look at the help:

    STRONG -h

The required argument are an output directory, and a config file. You can also specify a step to run. Try the following command:

    STRONG /home/student/Strain_resolution/STRONG_run assembly --threads 4 --config config.yaml --verbose --dryrun

What is happening?

To launch it for real, remove "assembly" and "--dryrun"
<details><summary>In depth STRONG commands</summary>
<p>
Basic STRONG commands : 
```
STRONG outputdir --config config.yaml
```
This will run all steps generate output in ***outputdir*** and search for the config yaml in ***outputdir*** if the config file is not specified : you can shorten the command line if you create the config file inside ***outputdir*** .

It is also possible to run just part of the pipeline:

1. ***assembly***: Runs just the [assembly steps](#assembly)  
2. ***graphextraction***: Runs just the [graph extraction](#graphextraction)
3. ***bayespaths***: Runs just the [bayespaths](#bayespaths)
4. ***evaluation***:  Runs just the [evaluation steps](#evaluation)
5. ***results***:  Runs just the [results steps](#results)
6. ***desman***:  Runs just the [desman steps](#desman)

By specifying desired step e.g.:
```
STRONG outputdir --config config.yaml bayespaths
```

This is useful if for example you wish to rerun strain resolution with alternate parameters then similarly remove the bayespaths directory and rerun as above.

The program also takes the following optional parameters:

```
  --threads THREADS, -t THREADS
```

Specify the maximum thread number to be used.

```
  --verbose, -v         Increase verbosity level
```

Useful to obtain more info from SnakeMake

```
  --dryrun, -n          Show tasks, do not execute them
```

Again snakemake command to list commands that would be run not to actually execute them.


```
  --unlock, -u          Unlock the directory
```

Will unlock directories if snakemake fails.

```
  --dag DAG, -d DAG     file where you want the dag to be stored
```

Generate diagram of jobs.


```
  -s ...                Pass additional argument directly to snakemake
```

Enables any additional commands to be passed to snakemake.

</p>
</details>

### Let's restart
As the assembly on this dataset may take up to 45 min, we are going to skip it.
Interrupt your current STRONG run, (ctrl+c)
Symlink a pre-generated assembly folder to your STRONG_run folder : 

    cd /home/student/Strain_resolution/STRONG_run
    rm -r assembly
    ln -s /home/student/repos/strain_resolution_practical/assembly .

Try and relaunch STRONG.


## Pipeline : detailed description

<a name="assembly"/>

### Assembly, COG annotation and binning 

The first step of the pipeline is a coassembly of all samples followed by binning. This involves multiple steps, including mapping with bowtie2 to get coverages and annotations 
to COGs with RPS-Blast, this is summarised in the figure:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/Dag_rules1.png)



This part of the pipeline produces a number of intermediate output files. We detail the key ones here:
1. ***assembly/spades/***: This directory contains the standard metaSPAdes run including ***assembly.fasta*** the contigs used in MAG construction 
2. ***assembly/high_res/***: This directory contains the high resolution assembly graph pre- ***graph_pack.gfa*** and post-simplication ***simplified.gfa*** 
and also ***simplified.mult_prof*** the unitig kmer coverages of the simplified graph across samples
3. ***annotation***: This directory contains contains the contig ORF predictions and COG annotations with RPS-BLAST
4. ***binning***: Contains the CONCOCT bins post refinement and merging these are given in ***clustering_gt1000_merged.csv*** as a csv file of contig names 
with bin assignments together with a list of MAGs satisfying 75% single-copy core genes in single copy ***list_mags.tsv*** 

The list of single-copy core genes are given as COGs in the data file ***SnakeNest/scg_data/scg_cogs_to_run.txt*** as default but this file can be changed.


Now let's have a look at the assembly graph, low resolution and high resolution.

    mkdir AssemblyGraphs
    cd AssemblyGraphs 
    
    wget https://desmantutorial.s3.climb.ac.uk/contigs_colorM.gfa
    wget https://desmantutorial.s3.climb.ac.uk/high_res_colorM.gfa


Start up Bandage:

    Bandage

Open up the first of these files contigs_colorM.gfa in Bandage. You should see something like this:

![Bandage contigs](https://github.com/chrisquince/Ebame5/blob/master/Figures/contigs_bandage.png)

The colors correspond to four MAGs we obtained from clustering the Spades contigs: 
* Bin3 magenta, E. faecalis
* Bin7 blue 
* Bin12 red
* Bin19 green, Staph. epidermidis

Why are some of the bins fragmented?

Can you find any contigs that are misassigned

Locate these two contigs with the search feature:

NODE_55_length_32977_cov_19.323249, NODE_327_length_8496_cov_5.646014

These correspond to contigs annotated to the single-copy core gene COG0060 in Bin19. Why are there two of them? 
Try blasting the sequences against the NCBI.

Now open up the file high_res_colorM.gfa

And find these nodes. Corresponding to COG0016 in Bin19:
* start 2816027
* end 2524601

You should be able to determine that at least two strains are present from the variant bubbles in the graph.

This is what the STRONG pipeline resolves into strains.


<a name="graphextraction"/>
<a name="bayespaths"/>

### Graph extraction and BayesPaths

This part of the pipeline extracts the single-copy core genes for each MAG from the simplified HRAG together with the unitig coverage profiles. These then undergo another round of simplification using the MAG coverages as well as potential merging of bins that share unitigs on the SCG subgraphs. These SCGs for each MAG are then run in the BayesPaths strain resolution algorithm. In detail:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/Dag_rules2.png)

This section produces a number of outputs:

1. ***subgraphs/bin_merged/bin_name/simplif***: directory contains the simplified SCG graphs for each MAG
2. ***bayespaths/selected_bins.txt***: text file listing MAGs that are run by BayesPaths
3. ***bayespaths/bin_name***: directory contains the BayesPaths output for each MAG with id bin_name. This contains:
   1. bin_nameF_Haplo_X.fa: SCG haplotypes for the X strains predicted for this MAG
   2. bin_nameF_Intensity.csv: the intensities for each strain in each sample (coverage depth/read length) 
   3. bin_nameF_varIntensity.csv: the variance of the intensities for each strain in each sample (coverage depth/read length) 
   4. bin_nameF_Divergence.csv: the divergences for each strain, these are proportional to 
path uncertainties
   5. bin_nameF_maxPath.tsv: most likely unitig paths for each strain
   6. bin_nameF_geneError.csv: errors associated with individual SCGs
   7. bin_nameF_Bias.csv: unitig biases 
   8. bin_nameF_Precision.csv: unitig precisions

<a name="evaluation"/>

### Evaluation

This section of the pipeline should only be run if known reference genomes are available because this is a benchnarking run with synthetic reads or in vitro mock communities. 

<a name="results"/>

### Results

This part of the pipeline generates results summaries of the MAG strain inference. The steps are detailed below:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/Dag_rules7.png)

A sub-directory is generated for each bin. These contain:

1. A phylogentic trees of strains created on the single-copy genes with a combined heat map of percent sequence identity for the bin, for example:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/TreeExample24.png)

This will include the Bin consensus contig sequence (Bin_Name) (and alternatives if multiple COGs are present in bin - Bin_Name_nb) and evaluation strains when available. 

2. In the ***graph*** sub-directory gfa files coloured by haplotype. These are viewable with [Bandage](https://rrwick.github.io/Bandage/) the file ***joined_SCG_graph.gfa*** contains all scgs in a single graph and the individual graphs are in the  ***cogs*** subdirectory

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/GraphExample.png)

3. A barchart of the strain intensities in each sample:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/IntensityExample.png)

There are also combined pdfs in the top level of results ***haplotypes_coverage.pdf***
and ***haplotypes_tree.pdf***. Finally ***summary.tsv*** contains some info on the assembly and number of strains resolved.

<a name="desman"/>

### DESMAN

This section runs the [DESMAN](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1309-9) program on each MAG to infer strains using variant frequencies across samples obtained from read mapping. This can be viewed as a complement to the BayesPaths algorithm and a validation of its predictions. The detailed pipeline is:

![alt tag](https://github.com/chrisquince/STRONG/blob/master/Figures/Dag_rules4.png)

The DESMAN steps are parametrised by the following options in the ***desman*** subsection of the config yaml:



This produces outputs in the ***desman*** directory for each MAG there is a sub-directory 
***bin_name*** which contains:

1. ***best_run.txt***: the predicted strain number for the bin format 
(G,H,R,Err,TauFile)
where G is the number of haplotypes, H the number that are reliable, Err their mean error and the TauFile points to the best haplotype encodings
2. ***Run_m_n***: DESMAN run repeat n for m haplotypes
3. ***Deviance.pdf***: Deviance plot of fit with haplotype number

<a name="Synthetic"/>










