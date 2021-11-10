# Initial processing steps #

Most of these processes are incorporated into loops over all samples.

Trim raw fastqs

```
java -jar trimmomatic-0.33.jar
PE -threads 7 -phred33
sample_read1.fastq.gz
sample_read2.fastq.gz
sample.trimmo.R1.fq.gz sample.trimmo.unpairF.fq.gz sample.trimmo.R2.fq.gz
sample.trimmo.unpairR.fq.gz
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Demultiplex using process_radtags from Stacks

```
process_radtags -1 sample.trimmo.R1.fq.gz -2
sample.trimmo.R2.fq.gz -o ./SAMPLES/ -b ./BARCODES/list.of.barcodes
-c -q -r --renz_1 ASEI --renz_2 BSTBI -i gzfastq --disable_rad_check
```

Align to reference using BWA-MEM:
```
bwa mem -t 7 ../GENOME/pdae_genome.v1.fa ./BST_PREPROC/$ind.1.fq
./BST_PREPROC/$ind.2.fq | samtools view -bS - | samtools sort -
./BST_ALIGN/$ind.sort.bam
```

Add read groups to bams and filter for mapping quality >= 30
```
samtools view -q 30 -@ 7 -bh ./BST_ALIGN/$ind.sort.bam.bam |
samtools addreplacerg -r 'ID:id'$ind -r 'LB:lb'$ind -r 'SM:'$ind - |
samtools sort -@ 7 - > ./BST_ALIGN/Q30RG/$ind.RG.q30.sort.bam
```

Use freebayes to call genotypes
```
freebayes command
```

Filter freebayes vcf based on RAD focused filters according to dDocent pipeline. See http://www.ddocent.com/filtering/ for filtering approach and discussion of individual filters.


