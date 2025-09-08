cat sample_trimm.txt | while read id
do
	fq1=./trimm_gal/${id}_1_val_1.fq.gz
	fq2=./trimm_gal/${id}_2_val_2.fq.gz
	echo "start hisat2 for ${id}" 'date'
	hisat2 -p 8 -x /HDD8T/jeeh99/GRCV_hisat2/gencode44 -1 ${fq1} -2 ${fq2} -S ./hisat2/${id}.sam >> ./hisat2/${id}.log 2>&1
	echo "end hisat2 for ${id}" 'date'
  echo "start samtools for ${id}" 'date'
  samtools view -Sb ./hisat2/${id}.sam > ./bam/${id}.bam
  samtools sort -@ 10 ./bam/${id}.bam  > ./sort/${id}_sort.bam
	rm ./hisat2/${id}.sam
  echo "end samtools for ${id}" 'date'
done

