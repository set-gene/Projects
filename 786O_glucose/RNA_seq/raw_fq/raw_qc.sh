#!/bin/bash

cat /HDD8T2/jeeh99/240529_manvi/raw_fq/sample_raw.txt | while read id
do
	fastqc --outdir ./raw_qc/ --threads 16 ./${id}*.fastq.gz >> ./raw_qc/${id}_fastqc.log 2>&1 
done 
