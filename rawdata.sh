mkdir raw-data
cd raw-data
for i in $(cat …/SraAccLis.txt);    #SraAcclis.txt contains a list of the accesion numbers
do
    echo $i
    fastq-dump --gzip --split-files $i  #fastq-dump gets data in fastq format
done
