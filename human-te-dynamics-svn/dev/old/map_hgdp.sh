#!/bin/bash --norc

wd=/home/a1222423/scratch/TEs
ref=${wd}/ref/humrep25.ref
study='HGDP-Pilot-TEs'
fastq=${wd}/fastq
mapped=${wd}/bam

samples=(HGDP00450_Mbuti HGDP00463_Mbuti HGDP00458_Biaka HGDP00464_Biaka HGDP00525_French HGDP00534_French HGDP00550_PapuanHighlands HGDP00556_PapuanHighlands HGDP00779_Han HGDP00822_Han HGDP00992_San HGDP01029_San HGDP01408_BantuKenya HGDP01416_BantuKenya)


i=1
for sample in ${samples[@]}
do
    echo =========================================
    echo Sample $i/${#samples[@]} $sample
    echo =========================================

    bwa mem $ref $fastq/*${sample}*R1*gz $fastq/*${sample}*R2*gz -t 54 -T 0 -R "@RG\tID:$sample\tSM:$study\tPL:Illumina" | samtools view -F4 -OBAM > $mapped/${sample}.bam
    ((i=i+1))
done