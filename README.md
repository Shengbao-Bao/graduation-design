1.the method to get refernce:
for i in $(seq 1 22) X Y M;do echo $i;wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;done
gunzip *.gz
for i in $(seq 1 22) X Y M;do cat chr${i}.fa >> hg19.fasta;done

1.python is used to remove peaks apear in less than 10% cells
2.R is used to yield a file with chr and its bases
