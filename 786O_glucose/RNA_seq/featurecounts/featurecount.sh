cat sample_trimm.txt | while read id 
do 
	echo "start featurecounts for ${id}" 'date'
featureCounts -M -T 6 -p -t exon -g gene_name -a /HDD8T/jeeh99/GRCV_hisat2/gencode.v44.annotation.gtf -o ./featurecounts/${id}.txt ./sort/${id}_sort.bam >> ./featurecounts/${id}.log 2>&1
	echo "end featurecounts for ${id}" 'date'
done
