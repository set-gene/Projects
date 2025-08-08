#!/bin/bash

####SRR download
cat /HDD8T2/jeeh99/GSE123310_Muscle_mm/SRR.txt | while read id

do
fastq-dump ${id} --split-3 --gzip
done

####fastq mapping
cat /HDD8T2/jeeh99/GSE123310_Muscle_mm/SRR.txt | while read id
do 
hisat2 -p 8 -x /HDD8T/jeeh99/GRCm39/GRCm39 -1 ${id}_1.fastq.gz -2 ${id}_2.fastq.gz -S ${id}.sam 
done  

####sam to sorting bam
cat /HDD8T2/jeeh99/GSE123310_Muscle_mm/SRR.txt | while read id
do 
samtools view -Sb ${id}.sam > ${id}.bam 

rm ${id}.sam

samtools sort -@ 10 ${id}.bam > ${id}_sort.bam

done  

####counting bam
cat /HDD8T2/jeeh99/GSE123310_Muscle_mm/SRR.txt | while read id
do 
featureCounts -T 6 -p -t exon -g gene_name -a /HDD8T/jeeh99/GRCm39/gencode.vM32.annotation.gtf -o ${id}.txt ${id}_sort.bam
done

