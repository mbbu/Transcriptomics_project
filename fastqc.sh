mkdir rawfastqc
cd rawfastqc
for i in `../*.gz`;
do
  	fastqc $i
done
multiqc ./
