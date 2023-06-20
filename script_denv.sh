# 1 # fastqc #
fastqc -t 25 *
mkdir fastqc ; 
mv *.html *.zip fastqc/ ; 
ls -lh ; 

# 2 # kraken viral #
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
kraken2 --paired --use-names --gzip-compressed --db /home/administrador/Documentos/KRAKENVIRDB/ --threads 28 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
rm *.fastq.gz_report.txt ; 
mkdir kraken_out ;
mv *.out  kraken_out/ ;
mkdir kraken_txt ; 
mv *.txt kraken_txt/ ;  
cd kraken_txt/ ; 
ls -lh ; 

# 3 # PAVIAN #
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

pavian::runApp(port=5000)

# 4 # trimming #
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
java -jar trimmomatic-0.39.jar PE -threads 28 $r1 $r2 ${prefix}_f_paired.fq.gz ${prefix}_f_unpaired.fq.gz ${prefix}_r_paired.fq.gz ${prefix}_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:45 ;
done ; 

mkdir trimm ; 
mv *_paired.fq.gz trimm/ ; 
rm *_unpaired.fq.gz ;
cd trimm/ ;
fastqc *.gz -t 25 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ;

# 5 # de-novo assembly SPADES #
for r1 in *.fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
spades --pe1-1 $r1 --pe1-2 $r2 --careful -t 4 --phred-offset 33 -m 15 -o ${prefix}_spades ;
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
cp ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rm -r *.fq.gz_assembly ;
ls -lh ;

#5 # de-novo assembly MEGAHIT #
for r1 in *.fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
megahit -1 $r1 -2 $r2 -t 4 -o ${prefix}_assembly ;
mv ${prefix}_assembly/final.contigs.fa ${prefix}_assembly/${prefix}_final.contigs.fa ;
cp ${prefix}_assembly/${prefix}_final.contigs.fa . ;
done ;

#6# compress#
rar a Sample_01_X_2023.fa.rar Sample_01_X_2023.fa ;
rar a Sample_02_X_2023.fa.rar Sample_02_X_2023.fa ;
rar a Sample_03_X_2023.fa.rar Sample_03_X_2023.fa ;
rar a Sample_04_X_2022.fa.rar Sample_04_X_2022.fa ;
rar a Sample_05_X_2022.fa.rar Sample_05_X_2022.fa ;
rar a Sample_07_Y_2023.fa.rar Sample_07_Y_2023.fa ;
rar a Sample_08_Y_2023.fa.rar Sample_08_Y_2023.fa ;
ls -lh ;

# install.packages("seqinr") #
library(seqinr) ;
r <- dir() ;
head <- gsub(".fa","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fas")) ;
} ;
q("no") ;

cat *.fas > contigs.fasta

cd .. ;
mkdir scaffolds ; 
mv *scaffolds.fasta script_2.R scaffolds/ ;
cd scaffolds/ ;

## BLAST ##
makeblastdb -in DB.fasta -dbtype nucl ;
blastn -db /home/vjimenez/Documentos/polio/DB_all/DB.fasta -query contigs.fasta -perc_identity 65 -max_target_seqs 20 -outfmt 6 -num_threads 15 > blast.csv ;
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
ins <- sort(unique(bm$query.acc.ver))
bn <- separate(bm,query.acc.ver,c("sample","remover"),sep="_NODE_")
bo <- data.frame(bn[1:2],bm[2:37])
bed <- blast[ , c(1,7,8)]

write.table(bed,"seq.bed",row.names=F, col.names=F,sep="\t",quot=F)
write.table(ins, "extraerdecontigs.fasta.txt", sep="\t", row.names=F, col.names=F, quot=F)
write.table(lt, "extraerdeDB.fasta.txt", sep="\t", row.names=F, col.names=F, quot=F)
write.table(bo, "blast_total2.tsv", sep="\t", row.names=F)

quit("no") ; 
samtools faidx /home/vjimenez/Documentos/polio/DB_all/DB.fasta ;
samtools faidx contigs.fasta ; 
seqtk subseq contigs.fasta extraerdecontigs.fasta.txt > selectos_ins.fasta ; 
seqtk subseq /home/vjimenez/Documentos/polio/DB_all/DB.fasta extraerdeDB.fasta.txt > selectos.fasta ; 
bedtools getfasta -fi contigs.fasta -bed seq.bed > ins.fasta ; 
bedtools getfasta -fi contigs.fasta -bed seq.bed > ins.fasta ; 
mafft --addfragments ins.fasta --adjustdirection --6merpair --thread 15 sabin.fasta > ins2.fasta ;
mafft --addfragments selectos.fasta --adjustdirection --6merpair --thread 15 ins2.fasta > ins3.fasta ;
aliview ins3.fasta ;
ls -lh ; 
 





#6 # looping mapping and estimate abundance #
for r1 in *fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
r3=${prefix}.fas

bwa index $r3 ;
bwa mem -t 4 $r3 $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T $r3 ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_count.bam ;
samtools index -@ 4 ${prefix}_count.bam ;

rm ${prefix}_uno.sam ${prefix}_unoa.bam ;
done ;

for p1 in *_count.bam
do
prefix=$(basename $p1 _count.bam)
samtools view -@ 4 $p1 | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c > ${prefix}_readcounts.txt ;
done ;
rm *.amb *.ann *.bwt *.fai *.pac *.sa ;
cat *_readcounts.txt > counts.txt ;
grep "Sample" counts.txt > counts2.txt ;
ls -lh ; 

#10# CAT/BAT#
CAT contigs -c all.fasta -d /home/administrador/Documentos/CAT_BAT/CAT_database.2023-04-06 -t /home/administrador/Documentos/CAT_BAT/CAT_taxonomy.2023-04-06 -o all2 ;
CAT add_names -i all2.contig2classification.txt -o all2.contig2classification.official_names.txt -t /home/administrador/Documentos/CAT_BAT/CAT_taxonomy.2023-04-06 --only_official ;
CAT summarise -c all.fasta -i all2.contig2classification.official_names.txt -o all2.summary.txt ;
ls -lh ;

#16# merge odentification and abundances#
Rscript script_3.R ;
ls -lh ;





###############################################################################################

# 7 # next #
## download DB : https://ccb.jhu.edu/software/kraken/ ##
kraken --db /home/administrador/Documentos/kraken/minikraken_20171019_8GB/ --fasta-input all.fasta --threads 25 --unclassified-out unclassified --classified-out classified --output iih_files ; 
kraken-translate --db /home/administrador/Documentos/kraken/minikraken_20171019_8GB/ iih_files > iih_files.labels.csv ; 

# 8 # CAT/BAT #
CAT prepare --fresh -d /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_CAT_database2/ -t /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_taxonomy2/ ; 
CAT contigs -c all.fasta -d /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_CAT_database -t /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_taxonomy -o all2 ;

## mapping and extracting ##
bowtie2-build --threads 25 DB.fasta DB

for r1 in *_paired.fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
bowtie2 --no-unal -x DB --threads 25 -1 $r1 -2 $r2 -S ${prefix}.sam ; 
samtools view -bS -T DB.fasta ${prefix}.sam > ${prefix}.bam ;
samtools sort -n ${prefix}.bam -o ${prefix}_sorted.bam ;
samtools index ${prefix}_sorted.bam ;

#4# obtencion de: 1)archivos bam condormado por solo reads No-mapeados y 2) fastq files "f" y "r" de estos reads mapeados#
samtools view -@ 15 -b -f 3 ${prefix}_sorted.bam > ${prefix}.mapped.bam ;
samtools index -@ 15 ${prefix}.mapped.bam ;
samtools fastq -1 ${prefix}_f.fastq -2 ${prefix}_r.fastq -0 /dev/null -s /dev/null -n ${prefix}.mapped.bam ;
done ;
mkdir aligned ; 
mv *.sam *.mapped.bam *.fastq aligned/ ; 
cd aligned/ ; 
ls -lh ; 

## iva assembly#
for r1 in *.fastq.gz
do
prefix=$(basename $r1 _f.fastq.gz)
r2=${prefix}_r.fastq.gz
iva -t 25 -f $r1 -r $r2 ${prefix}_iva ;
mv ${prefix}_iva/contigs.fasta ${prefix}_iva/${prefix}_iva_contigs.fasta ;
mv ${prefix}_iva/${prefix}_iva_contigs.fasta . ;
done ;




