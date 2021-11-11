
# ------- assembly  -------

mkdir -p  ~/Projects/AD_snakemake
cd ~/Projects/AD_snakemake

# create file for megahit
ls ~/Data/AD_small/*/*R1.fastq | tr "\n" "," | sed 's/,$//' > R1.csv
ls ~/Data/AD_small/*/*R2.fastq | tr "\n" "," | sed 's/,$//' > R2.csv

# launch megahit
megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 4 -o Assembly


# ------- map reads  -------

# index assembly
bwa index Assembly/final.contigs.fa

# loop over all fastq.files
mkdir -p Map
for file in ~/Data/AD_small/*/*R1.fastq
do 
   
   stub=${file%_R1.fastq}
   name=${stub##*/}
   
   echo $name

   file2=${stub}_R2.fastq

   # do a succession of task, by feeding output of command into input of the next one
   # map reads to assembly, filter to only maping reads, sort reads in bam files

   bwa mem -t 4 Assembly/final.contigs.fa $file $file2 | samtools view -b -F 4 - | samtools sort - > Map/$name.mapped.sorted.bam
done

# ------- Metabat2  -------
# create the depth file
cd ~/Projects/AD_snakemake/Map
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

# launch metabat2
cd ~/Projects/AD_snakemake
mkdir Binning
mv Map/depth.txt Binning/depth.txt
metabat2 -i Assembly/final.contigs.fa -a Binning/depth.txt -t 4 -o Binning/Bins/Bin





