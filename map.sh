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
