#!/bin/bash

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz 
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 30 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@  30 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 30 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 30 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 30 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 30 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 30 ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 20, frecuencia de nucleotido predominante: 55%#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix} -q 20 -t 0.55 -m 15 ;
samtools mpileup -aa -A -d 600000 -B -Q 0 ${prefix}.bam | ivar variants -p ${prefix} -q 20 -t 0.55 -r reference.fasta ;
done ;

#8# estimar la profunidad de cobertura y el porcentaje de Ns#
for r1 in *bam
do
prefix=$(basename $r1 .bam)
DEP=( `samtools depth $r1 | awk '{sum+=$3}END{print sum/10735}' `)
NPE=( `seqtk comp ${prefix}.fa | awk '{x=$3+$4+$5+$6;y=10735;print 1-(y-x)/y}'`)
echo "${prefix} ${DEP}x $NPE" >> profundidad_ns.txt ;
done ;

rm *.fastq.gz.bam *.fastq.gz.bam.bai *.fastq.gz.fa *.fastq.gz.qual.txt *.fastq.gz.tsv ; 
mkdir assembly_15 ;   
cat *.fa > secuencias.fasta ; 
cat secuencias.fasta reference.fasta > assembly_15/in.fasta ; 
mafft --thread 30 --auto --inputorder assembly_15/in.fasta > assembly_15/out.fasta ; 
mv *.fa *.tsv *.txt assembly_15/ ; 
ls -lh ; 

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 25, frecuencia de nucleotido predominante: 60%, profundidad: 3#
for r1 in *bam
do
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix} -q 25 -t 0.6 -m 3 ;
samtools mpileup -aa -A -d 600000 -B -Q 0 ${prefix}.bam | ivar variants -p ${prefix} -q 25 -t 0.6 -r reference.fasta ;
done ;

#############################
########### BLAST ###########
#############################

## BLAST ##
makeblastdb -in /home/vjimenez/Documentos/DENV/DENV_2023/BLAST/NCBI.fasta -dbtype nucl ;
blastn -db /home/vjimenez/Documentos/DENV/DENV_2023/BLAST/NCBI.fasta -query out_edition.fasta -perc_identity 97 -max_target_seqs 10 -outfmt 6 -num_threads 15 > blast.csv ;
ls -lh ;

R 
dir()
getwd()
library(tidyr)
blast <- read.csv("blast.csv", header=FALSE, sep="\t") ; 
names(blast) <- c("query.acc.ver","subject.acc.ver","perc.identity","alignment.length","mismatches","gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score") ; 
head(blast) ; 
write.table(blast,"blast_renamed2.tsv", sep="\t", row.names=F)

m <- read.csv("metadata.csv")
names(blast) <- c(names(blast)[1],"Accession",names(blast)[3:12])
bm <- merge(blast,m,by="Accession", all.x=F)
lt <- unique(bm$Accession)

mafft --addfragments sequences.fasta --adjustdirection --6merpair --thread 15 out_edition.fasta > ins2.fasta ;
aliview ins2.fasta ;
ls -lh ; 

#############################

sed 's/Consensus_//g' ins3.fasta > a.fasta ;
sed 's/_L001_.*//g' a.fasta > b.fasta ;
aliview b.fasta ; 

#############################

a <- read.csv("denv3.txt", header=F) ; 
f <- rep(c("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"),((length(a$V1))/2)) ; 
a1 <- paste0(a$V1,f) ; 
write.table(a1,"denv3_new.txt",col.names=F, quot=F, row.names=F) ; 

## alternative ##
a <- read.csv("denv2.txt", header=F) ; 
a1 <- c(paste0(a$V1,"_L001_R1_001.fastq.gz"),paste0(a$V1,"_L001_R2_001.fastq.gz")) ; 
write.table(sort(a1),"denv2_new.txt",col.names=F, quot=F, row.names=F) ; 

mkdir denv1 ; 
for file in $(cat denv1_new.txt)
do 
mv ${file} denv1/
done

mkdir denv2 ;
for file in $(cat denv2_new.txt)
do 
mv ${file} denv2/
done

mkdir denv3 ;
for file in $(cat denv3_new.txt)
do 
mv ${file} denv3/
done









































 




