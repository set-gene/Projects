#!/bin/bash

####SRR download
cat /HDD8T2/jeeh99/PRJNA607321_mm_WAT/SRR.txt| while read id

do
fastq-dump ${id} --gzip
done

####fastq mapping
cat /HDD8T2/jeeh99/PRJNA607321_mm_WAT/SRR.txt| while read id

do 
hisat2 -p 8 -x /HDD8T/jeeh99/GRCm39/GRCm39 -U ${id}.fastq.gz -S ./${id}.sam 
done  

####sam to sorting bam
cat /HDD8T2/jeeh99/PRJNA607321_mm_WAT/SRR.txt| while read id

do 
samtools view -Sb ./${id}.sam  > ./${id}.bam 

rm ./${id}.sam 

samtools sort -@ 10 ./${id}.bam > ./${id}_sort.bam
done  

####counting bam
cat /HDD8T2/jeeh99/PRJNA607321_mm_WAT/SRR.txt | while read id 

do 
featureCounts -T 6 -p -t exon -g gene_name -a /HDD8T/jeeh99/GRCm39/gencode.vM32.annotation.gtf -o ${id}.txt ${id}_sort.bam
done
