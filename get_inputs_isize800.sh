#!/bin/bash

if [ "$1" == "-h" ] ; then
    echo -e "Usage: "$(basename $0)" [fasta file genome] [PE-read 1] [PE-read 2] [threads]"
    exit 1
fi


f=$1
R1=$2
R2=$3
threads=$4


mkdir bowtie 
cd bowtie

cp $f ./fasta.fa

bowtie2-build-s fasta.fa fasta 
bowtie2-align-s -p $threads --very-fast -q -x fasta -1 $R1 -2 $R2 --no-unal -S fasta.sam 
samtools view --threads $threads -h -b -S fasta.sam | samtools sort > fasta0.bam

## get only NM0 reads 
samtools view fasta0.bam | grep "NM:i:0" > fasta0_filter.tmp   
samtools view -H fasta0.bam > fasta0_header.tmp    
cat fasta0_header.tmp fasta0_filter.tmp > fasta0_filter.sam
samtools depth -a fasta0_filter.sam > depth.txt     
samtools view --threads $threads -h -b -S fasta0_filter.sam | samtools sort > fasta0_filter.bam     

rm *.sam *.tmp *.bt2

# prepare fasta.bam

samtools view fasta0.bam | perl -pe 's/\t\w\w:i:/\t/g' | perl -pe 's/\t\w\w:\w:/\t/g' | cut -f1-18 > fasta0.txt

cd ..
 
# get repeats by blastn (save as blastn.txt)

makeblastdb -in $f -dbtype nucl -parse_seqids
blastn -word_size 7 -perc_identity 80 -evalue 0.001 -db $f -out blastn.tmp -outfmt "6 qseqid sseqid qstart qend sstart send length bitscore evalue pident qcovs qcovhsp btop sstrand" -query $f -parse_deflines
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tlength\tbitscore\tevalue\tpident\tqcovs\tqcovhsp\tbtop\tsstrand" > head.tmp
cat head.tmp blastn.tmp > blastn.txt
rm "${f}".* *.tmp

# get identical regions

mkvtree -db $f -dna -pl -allout -v
vmatch -d -p -l 20 -showdesc 10 $f | perl -pe s'/^\s+//g' |  perl -pe 's/\s+/\t/g' | perl -pe 's/100.00\t/100.00\n/g' | awk -v OFS='\t' '{print $2,$3,$3+$1-1,$6,$7,$7+$1-1,$1}' | tail  -n +2 > identical_regions.txt
rm "${f}".*

# get chromosomal chr lengths


awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $f | perl -pe 's/>(.*)\n(.*)/$1\t$2/g' > contigs_lengths.txt

