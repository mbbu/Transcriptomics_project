
# Codes


Fetching rawdata from NCBI

Getting the data one data-set at a time

```
fastq-dump --gzip --split-files <Accession number>
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
Quality analysis
looping through the rawdata to do quality check and combining the reports using multiqc
```
mkdir rafastqc
cd rafastqc
for i in `../*.gz`;
do
  	fastqc $i
done
multiqc ./
```
Trimming
```
for i in *_1.fastq.gz
do
    name=$(basename ${i} _1.fastq.gz)
    trimmomatic PE ${i} ${name}_2.fastq.gz\
                   ${name}_1.trim.fastq.gz ${name}_1un.trim.fastq.gz \
                    ${name}_2.trim.fastq.gz ${name}_2un.trim.fastq.gz \
                    HEADCROP:11
done
```
Mapping
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
