
# Bioinformatics analysis
## Dependencies

* sratoolkit
* fastqc
* multiqc
* trimmomatic
* hisat2
* samtools
* htseq
## Create a conda environment & load modules
```
conda create --env Anopheles
conda activate Anopheles
module load <name of the module>
```

## Fetching raw data from NCBI

**fastq-dump** is a tool in SRAtoolkit used for downloading sequenced reads from NCBI Sequence Read Archive(SRA).The data is dowloaded in fastq format.Here we are using the two options *--gzip* for compressing the sequence reads and *--split-files* to give both forward and reverse reads since the reads are from Illumina platform.

Getting the data one data-set at a time


```
fastq-dump --gzip --split-files <Accession number>
fastq-dump --gzip --split-files SRR9987839
```
Getting all data at once
```
mkdir raw-data
cd raw-data
for i in $(cat …/SraAccLis.txt);    #SraAcclis.txt contains a list of the accesion numbers
do
    echo $i
    fastq-dump --gzip --split-files $i  #fastq-dump gets data in fastq format
done
```
## Quality analysis
**fastqc** checks the quality of the raw reads from high throughput sequencing platforms like Illumina.It provides a html report summary of the quality of reads in graphs and tables.This makes you aware of the quality of reads for downstream analysis.
**multiqc** makes fastqc output more manageable by compiling them and generating one report.

```
mkdir rawfastqc
cd rawfastqc
for i in `../*.gz`;
do
  	fastqc $i
done
multiqc ./
```

## Trimming
low quality reads together with adapters are removed after quality assessment.**Trimmomatic** is a commandline tool used to trim and crop reads. 

```
#basename is a command that strips trailing suffix in a file name
for i in *_1.fastq.gz;
do
    name=$(basename ${i} _1.fastq.gz)
    trimmomatic PE ${i} ${name}_2.fastq.gz\
                   ${name}_1.trim.fastq.gz ${name}_1un.trim.fastq.gz \
                    ${name}_2.trim.fastq.gz ${name}_2un.trim.fastq.gz \
                    HEADCROP:11
done
```
The trimmed reads are then checked for quality.
## Mapping
RNA Seq reads are mapped to the reference genome.**hisat2** is a fast and sensitive splice-aware aligner that compresses the genome using an
indexing scheme to reduce the amount of space needed to store the genome. This also makes the genome quick to search, using a whole-genome index.We use samtools to convert the output file from mapping to bam format and to index the bam files.Indexing creates a searchable index of sorted bam files required in some programs.

```
#Building a reference genome index
#hisat2-build -p25 ./VectorBase-53_AgambiaePEST_AnnotatedCDSs.fasta  ../hisat-index/VectorBase-53_AgambiaePEST.idx

#Run hisat2 using indexed reference
mkdir Alignment_hisat
for i in $(cat ../SraAccList.txt);
do
   echo ${i}
   hisat2 -p25 -x ../hisat-index/VectorBase-53_AgambiaePEST.idx\
               -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz\
               -S Alignment_hisat/${i}_hisat.sam
   samtools view -Sb Alignment_hisat/${i}_hisat.sam  | samtools sort  > Alignment_hisat/${i}_hisat_sorted.bam
   samtools index Alignment_hisat/${i}_hisat_sorted.bam
   rm Alignment_hisat/${i}_hisat.sam 
done
```
## Abundance estimation

Once you have your aligned reads and a gff file ,**htseq** is used to give counts of reads mapped to each feature.A feature is an interval on a chromosome.
```
htseq-count -t CDS -i ID -f bam SRR9987840_hisat_sorted.bam VectorBase-53_AgambiaePEST.gff
```

## Differential analysis

Differential analysis involves using read counts to perform statistical analysis to discover quantitative changes in gene expression levels between experimental groups;exposed and non-exposed.**DESeq2** is used for differential analysis. 

 

