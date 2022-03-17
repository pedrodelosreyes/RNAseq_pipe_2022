#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

#! /bin/bash

#SMP=$1

## Accesing to working directory 
cd /home/pedro/arabidopsis/rnaseq_studies/co_rnaseq_2022/co_rnaseq_14dag_v2/samples/bams

## Gene Expression Quantification

#single end
#featureCounts -O -a ../../annotation/annotation.gtf -o counts.txt $(ls *.bam) -s 2 --verbose

#Paired end
featureCounts -O -p --countReadPairs -a ../../annotation/annotation.gtf -o counts.txt $(ls *.bam) -s 2 --verbose

#You have to modify the -s parameter depending on the strandness of your library


