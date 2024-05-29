1.the method to get refernce:
for i in $(seq 1 22) X Y M;do echo $i;wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;done
gunzip *.gz
for i in $(seq 1 22) X Y M;do cat chr${i}.fa >> hg19.fasta;done
2.the homo file can be got from https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme



1.python is used to remove peaks apear in less than 10% cells.
2.R is used to yield a file beginning with >, flowing chromosome place and bases.
3.python is used to scanning TF motif with fimo.
