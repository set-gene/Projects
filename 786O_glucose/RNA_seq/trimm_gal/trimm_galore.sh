#!/bin/bash

cat sample_trimm.txt | while read id 
do
	fq1=./raw_fq/${id}_1.fastq.gz
	fq2=./raw_fq/${id}_2.fastq.gz
	echo "start trim_galore for ${id}" 'date'
	trim_galore --paired --phred33 --length 30 --gzip --cores 8 -o ./trimm_gal $fq1 $fq2 >> ./trim_gal/${id}_trimm.log 2>&1
	echo "end trim_galore for ${id}" 'date'
done
