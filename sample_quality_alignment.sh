#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

#! /bin/bash

SMP=$1
ACC=$2

## Accesing to working directory 
cd /home/pedro/arabidopsis/rnaseq_studies/RNAseq_mpipz_2021/method2/samples/$SMP

## Sample quality control and read mapping to reference genome
#fastqc ${ACC}_1.fq.gz
#fastqc ${ACC}_2.fq.gz
STAR --runThreadN 4 --genomeDir ../../genome/ --readFilesIn ${ACC}_1.fq.gz ${ACC}_2.fq.gz --outFileNamePrefix $SMP --readFilesCommand gunzip -c

## Generating sorted bam file
samtools sort -o $SMP.bam ${SMP}Aligned.out.sam
rm ${SMP}Aligned.out.sam
samtools index $SMP.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam $SMP.bam -o $SMP.bw

## Transcript assembly
stringtie -G ../../annotation/arabidopsis_thaliana.gtf -o $SMP.gtf -l $SMP $SMP.bam

## Preparing merge list file for transcriptome merging
echo /home/pedro/arabidopsis/rnaseq_studies/RNAseq_mpipz_2021/method2/samples/$SMP/$SMP.gtf >> ../../results/merge_list.txt

## Moving bam file to another folder in order to perform featureCounts
mv $SMP.bam ../bams


