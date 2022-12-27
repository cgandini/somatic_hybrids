
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(seqinr)
library(Biostrings)
library(SamSeq)
library(ape)
library(bedtorch)
library(valr)

p1_fasta_file = "./inputs/nt_mtDNA.fa"
p2_fasta_file = "./inputs/po_mtDNA.fa"
NM_max = 0
max_gap_repeats = 550 ## inner size for PE-reads with inner size > 300, if not 300 bp
cut_size = c(0) ## 0 for PE-reads with large isize (>500 bp) isize, c(0,50,75) for PE reads with small isize (<500 bp)
isize = 800
read_length = 125
innersize = isize - 2*read_length
max_isize = isize + isize/2
min_isize = isize - isize/2
min_depth = 3
min_size_depth = 5
min_flanking = 5
p1_complete_name = "N. tabacum"
p2_complete_name = "P. orientalis"
hybrid_name = "ntpoA01"

revcomp = function(strand,read_seq) {
  if (strand == "minus") {
    read_seq2 = toupper(c2s(rev(comp(s2c(as.character(read_seq))))))
  } else {
    read_seq2 = read_seq
  }
    return(read_seq2)
}

### Add fasta files

p1_fasta=read.fasta(p1_fasta_file)
p2_fasta=read.fasta(p2_fasta_file)

p1_name = getName(p1_fasta)
p2_name = getName(p2_fasta)
p1_length = getLength(p1_fasta)
p2_length = getLength(p2_fasta)

### Add genes table

p1_genes = read.delim(file="./inputs/p1_gene_regions.txt", header=F,  col.names = c("p", "start", "end", "gene", "strand")) %>% 
  mutate(length = end-start+1)
p2_genes = read.delim(file="./inputs/p2_gene_regions.txt", header=F, col.names = c("p", "start", "end", "gene", "strand")) %>% 
  mutate(length = end-start+1)

### Get parental putative regions by coverage (>= min_depth, join if they are separated by < min_size_depth, min_length 1000 bp)
 
depth_file = read.delim("./inputs/depth.txt", col.names= c("p","pos","depth"), sep = "\t")

depth_p1 = depth_file %>% 
  filter(p == p1_name, depth >= min_depth)
depth_p1 = tapply(depth_p1$pos, cumsum(c(TRUE, diff(depth_p1$pos) > min_size_depth)), range) %>% 
  data.frame(matrix(unlist(.), ncol=2, byrow=T)) %>% 
  mutate(p=p1_name, length=X2-X1+1) %>% select(4,start=2,end=3,5)

depth_p2 = depth_file %>% 
  filter(p == p2_name, depth >= min_depth)
depth_p2 = tapply(depth_p2$pos, cumsum(c(TRUE, diff(depth_p2$pos) > min_size_depth)), range) %>% 
  data.frame(matrix(unlist(.), ncol=2, byrow=T)) %>% 
  mutate(p=p2_name, length=X2-X1+1) %>% select(4,start=2,end=3,5)

depth_all = rbind(depth_p1, depth_p2) 

parental_regions_by_cov = rbind(depth_p1, depth_p2) %>% 
  filter(length >= read_length) %>% 
  mutate(filter=ifelse(length >=1000, 1,0)) %>% 
  mutate(filter2= ifelse(lag(filter)==1 & start-lag(end) <= 500 | lead(filter)==1 & lead(start)-end <= 500,1,0)) %>%
  filter(filter == 1 | filter2 == 1)

rm(depth_p1, depth_p2) 

write.table(parental_regions_by_cov, "./outputs/parental_regions_by_cov.txt", sep="\t", row.names = F, quote=F)
write.table(parental_regions_by_cov %>% select(p, start, end), "./outputs/parentalregions4circos.txt", sep="\t", row.names = F, quote=F, col.names = F)
write.table(parental_regions_by_cov %>% mutate(start=start-1) %>% mutate(name=paste(start,"_",end, sep = "")), "./outputs/parental_regions_by_cov.bed", sep="\t", row.names = F, quote=F, col.names = F)

## plot coverage

depth_file = depth_file %>% mutate(cov=0)

for(i in 1:dim(parental_regions_by_cov)[1]) {
  start=parental_regions_by_cov$start[i]
  end=parental_regions_by_cov$end[i]
  name=parental_regions_by_cov$p[i]
  
  depth_file = depth_file %>% mutate(cov=ifelse(between(pos,start,end) & p==name,1,cov))
}

p1_plot_cov = ggplot(depth_file %>% filter(p == p1_name), aes(x=pos, y= depth)) + 
    ggtitle(paste(p1_complete_name, "genome coverage by", hybrid_name, sep=" ")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, p1_length), labels = scales::comma) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + 
    geom_area(data = depth_file %>% filter(p == p1_name) %>% filter(cov == 1), mapping = aes(depth = ifelse(cov == 1, depth, 0)), fill = "#E88056", alpha=0.3) +
    geom_line(color = "#E88056", size=0.03) +
    xlab("genome (bp)") +
    theme_classic()  
  
ggsave(file="./outputs/ntabacum_cov.pdf", plot = last_plot(), device = NULL, path = NULL,
  scale = 1, width = 10, height = 3.5, units = "in",
  dpi = 1200, limitsize = F) 

p2_plot_cov = ggplot(depth_file %>% filter(p == p2_name), aes(x=pos, y= depth)) + 
    ggtitle(paste(p2_complete_name,"genome coverage by", hybrid_name, sep=" ")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, p2_length),labels = scales::comma) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) + 
    geom_area(data = depth_file %>% filter(p == p2_name) %>% filter(cov == 1), mapping = aes(depth = ifelse(cov == 1, depth, 0)), fill = "#5BA864", alpha=0.3) +
    geom_line(color = "#5BA864", size=0.03) +
    xlab("genome (bp)") +
    theme_classic()  
  
ggsave(file="./outputs/porientalis_cov.pdf", plot = last_plot(), device = NULL, path = NULL,
  scale = 1, width = 10, height = 3.5, units = "in",
  dpi = 1200, limitsize = F) 
    

## get repeat file, modify subject start and end in minus pairs, eliminate duplicate repeats

repeat_file = read.delim(file="./inputs/blastn.txt", header=T) 

repeat_regions = repeat_file %>% 
  filter(qseqid == p1_name & sseqid == p2_name) %>%
  mutate(sstart2 = if_else(sstart < send, sstart, send)) %>% 
  mutate(send2 = if_else(send > sstart, send, sstart )) %>% 
  select(qseqid, sseqid, qstart, qend, sstart=sstart2, send=send2, length, pident, sstrand)  %>%
  arrange(qseqid, qstart) %>%
  tibble::rowid_to_column("i")

p1_repeats = repeat_regions %>% mutate(qstart=qstart-1000, qend=qend+1000) %>% select(chrom=qseqid, start=qstart, end=qend, i)
p1_genes_i = p1_genes %>% select(chrom=p, start, end, gene)

p1_genes_i = bed_intersect(p1_repeats, p1_genes_i) %>% group_by(i.x) %>% select(i=i.x, gene = gene.y) %>% distinct() %>% summarize(genes = str_c(gene, collapse = ";")) 

p2_repeats = repeat_regions %>% mutate(sstart=sstart-1000, send=send+1000) %>% select(chrom=sseqid, start=sstart, end=send, i)
p2_genes_i = p2_genes %>% select(chrom=p, start, end, gene)

p2_genes_i = bed_intersect(p2_repeats, p2_genes_i) %>% group_by(i.x) %>% select(i=i.x, gene = gene.y) %>% distinct() %>% summarize(genes = str_c(gene, collapse = ";")) 

repeat_regions = left_join(repeat_regions, p1_genes_i, by="i") %>% left_join(., p2_genes_i, by="i") %>% dplyr::rename(genes_p1 = genes.x, genes_p2 = genes.y)

## repeats2bedfile

repeat_regions_first = repeat_regions %>% unite(id, c(i,pident), remove = T, sep = "_") %>% mutate(cero = 0) %>% select(name = qseqid, start=qstart, end=qend, id, cero, sstrand)
repeat_regions_second = repeat_regions %>%  unite(id, c(i,pident), remove = T, sep = "_") %>% mutate(cero = 0) %>% select(name = sseqid, start=sstart, end=send, id, cero, sstrand)
repeat_regions_bed = bind_rows(repeat_regions_first, repeat_regions_second) %>% arrange(id) %>% mutate(sstrand = ifelse(sstrand == "forward", "+", "-"))
write.table(repeat_regions_bed, "./outputs/repeat_regions.bed", sep="\t", row.names = F, quote=F, col.names = F)
write.table(repeat_regions, "./outputs/repeat_regions.txt", sep="\t", row.names = F, quote=F, col.names = F)

rm(repeat_regions_first, repeat_regions_second)
 
## join repeats separated by less than "max_gap_repeat", and SNPs tables if so.

repeat_regions_joined = repeat_regions[FALSE,]
repeat_regions_all = repeat_regions
repeat_regions_all$qseqid = as.character(repeat_regions_all$qseqid)
repeat_regions_joined$sseqid = as.character(repeat_regions_joined$qseqid)
repeat_regions_all$sseqid = as.character(repeat_regions_all$sseqid)
repeat_regions_joined$qseqid = as.character(repeat_regions_joined$sseqid)
repeat_regions_all$i = as.character(repeat_regions_all$i)
repeat_regions_joined$i = as.character(repeat_regions_joined$i)
repeat_regions_all$length = as.character(repeat_regions_all$length)
repeat_regions_joined$length = as.character(repeat_regions_joined$length)
repeat_regions_all$pident = as.character(repeat_regions_all$pident)
repeat_regions_joined$pident = as.character(repeat_regions_joined$pident)


while (dim(repeat_regions_all)[1]>0) {
   qseqid=repeat_regions_all$qseqid[1]
   qstart=as.numeric(repeat_regions_all$qstart[1])
   qend=as.numeric(repeat_regions_all$qend[1])
   sseqid=repeat_regions_all$sseqid[1]
   sstart=as.numeric(repeat_regions_all$sstart[1])
   send=as.numeric(repeat_regions_all$send[1])
   strand=as.character(repeat_regions_all$sstrand[1])

   row = repeat_regions_all[1,]
   repeat_regions_all = repeat_regions_all[-c(1), ]
   df = repeat_regions_all 
   
  if(dim(df)[1]>0){
      for (j in 1:(dim(df)[1])) {
      if(qseqid == df$qseqid[j] & sseqid == df$sseqid[j] & strand == df$sstrand[j] & (between(df$qstart[j],qstart-max_gap_repeats,qend+max_gap_repeats) |
                                                            between(df$qend[j],qstart-max_gap_repeats,qend+max_gap_repeats)) & 
                                                            (between(df$sstart[j],sstart-max_gap_repeats,send+max_gap_repeats) |
                                                            between(df$send[j],sstart-max_gap_repeats,send+max_gap_repeats)))  {
                                                          
                                                            row = bind_rows(row, df[j,])
                                                            id = df$i[j]
                                                            repeat_regions_all = repeat_regions_all %>% filter(i != id)
                                                            qstart=min(as.numeric(row$qstart))
                                                            qend=max(as.numeric(row$qend))
                                                            sstart=min(as.numeric(row$sstart))
                                                            send=max(as.numeric(row$send))
      }
    }
  }
  
   if(dim(row)[1] > 1) {
     row_join = row %>% group_by() %>% mutate(length=paste0(length, collapse = ";")) %>% mutate(pident=paste0(pident, collapse = ";")) %>% mutate(i=paste0(i, collapse = ";")) %>% mutate(genes_p1 = apply(rbind(unique(genes_p1)), 1, function(x) paste(x[!is.na(x)], collapse = ";"))) %>%
         mutate(genes_p2 = apply(rbind(unique(genes_p2)), 1, function(x) paste(x[!is.na(x)], collapse = ";"))) %>%
         summarise(i=last(i), 
         qseqid = last(qseqid), 
         sseqid =  last(sseqid), 
         qstart = min(qstart), 
         qend = max(qend),
         sstart = min(sstart), 
         send = max(send), 
         length = last(length), 
         pident = last(pident), 
         sstrand = last(sstrand),
         genes_p1 = last(genes_p1),
         genes_p2 = last(genes_p2))
     
     repeat_regions_joined = bind_rows(repeat_regions_joined, row_join)
     row_join_i = row_join[1] %>% as.vector()
     

   } else {
     repeat_regions_joined = bind_rows(repeat_regions_joined, row)
   }
}

rm(row, row_join, row_join_i, i, repeat_regions_all, df,qseqid, qstart, qend, sseqid, sstart, send, strand)

## repeats_joined2bedfile

repeat_regions_joined_first = repeat_regions_joined %>% select(name = qseqid, start=qstart, end=qend, id=i)
repeat_regions_joined_second = repeat_regions_joined %>% select(name = sseqid, start=sstart, end=send, id=i)
repeat_regions_joined_bed = bind_rows(repeat_regions_joined_first, repeat_regions_joined_second) %>% arrange(id)
write.table(repeat_regions_joined_bed, "./outputs/repeat_regions_joined.bed", sep="\t", row.names = F, quote=F, col.names = F)
write.table(repeat_regions_joined, "./outputs/repeat_regions_joined.txt", sep="\t", row.names = F, quote=F, col.names = T)

rm(repeat_regions_joined_first, repeat_regions_joined_second)

## load identical regions (100% identical)
  
identical_regions = read.delim("./inputs/identical_regions.txt", header=F, sep = "\t")  %>% select(qseqid = 1,qstart = 2,qend = 3, sseqid = 4,sstart = 5, send =6, length = 7) %>%  
  filter(qseqid != sseqid) %>%
   tibble::rowid_to_column("i") 

identical_regions_first = identical_regions %>% select(name = qseqid, start=qstart, end=qend, i) 
identical_regions_second = identical_regions %>%  select(name = sseqid, start=sstart, end=send, i)
identical_regions_bed = bind_rows(identical_regions_first, identical_regions_second) %>% arrange(i) %>% mutate(cero = 0) %>% select(name, start, end, i, cero)
write.table(identical_regions_bed, "./outputs/identical_regions.bed", sep="\t", row.names = F, quote=F, col.names = F)

rm(identical_regions_first, identical_regions_second)

## Parse Bowtie SAM files
 
paired_parental = data_frame()

for(f in 1:length(cut_size)) {
  
 cut = cut_size[f]
 sam_file_discordant = read.delim(file=paste("./inputs/fasta", cut, ".txt", sep = ""), col.names= c("READ","FLAG","REFERENCE", "POS", "MAPQ", "CIGAR", "MATE", "MPOS", "LENGTH", "SEQ", "QUAL", "AS", "XN", "XM", "XO", "XG", "NM", "MD"), sep = "\t") %>% filter(MATE != "=")
 
 sam_file_concordant = read.delim(file=paste("./inputs/fasta", cut, ".txt", sep = ""), col.names= c("READ","FLAG","REFERENCE", "POS", "MAPQ", "CIGAR", "MATE", "MPOS", "LENGTH", "SEQ", "QUAL", "AS", "XN", "XM", "XO", "XG", "NM", "MD"), sep = "\t") %>% filter(MATE == "=") %>% mutate(span = MPOS+read_length-POS) %>% filter(span > min_isize & span < max_isize)
 
 assign(paste("sam_file_discordant",cut,sep = "_"), sam_file_discordant)  
 assign(paste("sam_file_concordant",cut,sep = "_"), sam_file_concordant)  
 
 sam_flags = sam_file_discordant %>% select(FLAG) %>% distinct()
 
 sam_flags_explained = data.frame(matrix(ncol = 13, nrow = 0))
 row=as.data.frame(t(as.matrix(samFlags(sam_flags$FLAG[1])))) %>% mutate(FLAG=sam_flags$FLAG[1]) %>% select(13, 1:12)
 colnames(sam_flags_explained) = colnames(row) %>% as.vector()
 rm(row)
 
 for (i in 1:dim(sam_flags)[1]) {
   FLAG = sam_flags$FLAG[i]
   flag_explained = samFlags(FLAG) %>% as.data.frame() 
   row = as.data.frame(t(as.matrix(flag_explained))) %>% mutate(FLAG=FLAG) %>% select(13, 1:12)
   sam_flags_explained = bind_rows(sam_flags_explained, row)
 }
 
 rm(row, flag_explained, FLAG)
 
 sam_flags2keep_R1 = sam_flags_explained %>% 
   filter(FIRST_OF_PAIR == T, MATE_UNMAPPED == F, READ_UNMAPPED == F, READ_PAIRED == T) %>% select(FLAG) 
 
 sam_flags2keep_R2 = sam_flags_explained %>% 
   filter(SECOND_OF_PAIR == T,MATE_UNMAPPED == F, READ_UNMAPPED == F, READ_PAIRED == T) %>% select(FLAG)
 
 sam_file_discordant_R1 = left_join(sam_flags2keep_R1, sam_file_discordant, by="FLAG") %>% 
   select(FLAG1 = FLAG, READ, REFERENCE, POS,  MATE, MPOS, NM1=NM, MD1=MD)
 
 sam_file_discordant_R2 = left_join(sam_flags2keep_R2, sam_file_discordant, by="FLAG")  %>% 
   select(FLAG2 = FLAG, READ, REFERENCE, POS,  MATE, MPOS, NM2=NM, MD2=MD)
 
  paired_parental_tmp = full_join(sam_file_discordant_R1, sam_file_discordant_R2, by = c("READ" = "READ", "REFERENCE" = "MATE", "MATE"= "REFERENCE", "POS" = "MPOS", "MPOS" = "POS")) %>%
   filter(NM1 <= NM_max, NM2 <= NM_max) %>% 
   mutate(POSEND = POS + read_length-1, MPOSEND = MPOS + read_length-1, READ_LENGTH = read_length - cut) 
 
  paired_parental = rbind( paired_parental,  paired_parental_tmp)
 
}

 paired_parental_p1 =  paired_parental %>% filter(REFERENCE == p1_name)
 paired_parental_p2 =  paired_parental %>% filter(REFERENCE == p2_name) %>% select(FLAG1=FLAG2, READ, REFERENCE=MATE, POS=MPOS, MATE=REFERENCE, MPOS=POS, NM1=NM2, MD1=MD2, FLAG2=FLAG1, NM2=NM1, MD2=MD1, POSEND=MPOSEND, MPOSEND=POSEND, READ_LENGTH)

paired_parentals = bind_rows( paired_parental_p1,  paired_parental_p2)

 
## eliminate all reads that falls entirely inside identical regions

 df = paired_parentals
 identical_regions_read_length = identical_regions %>% filter(length >= read_length)
 
 for(j in 1:dim(identical_regions_read_length)[1]) {
       Qseqid=identical_regions_read_length$qseqid[j]
       Qstart=identical_regions_read_length$qstart[j]
       Qend=identical_regions_read_length$qend[j]
       Sseqid=identical_regions_read_length$sseqid[j]
       Sstart=identical_regions_read_length$sstart[j]
       Send=identical_regions_read_length$send[j]
   
       df = df %>% mutate(out = ifelse(Qseqid == REFERENCE & Sseqid == MATE & between(POS, Qstart-1, Qend+1) & between(POSEND, Qstart-1, Qend+1) |
                                       Qseqid == REFERENCE & Sseqid == MATE & between(MPOS, Sstart-1, Send+1) & between(MPOSEND, Sstart-1, Send+1) |
                                       Qseqid == MATE & Sseqid == REFERENCE & between(MPOS, Qstart-1, Qend+1) & between(MPOSEND, Qstart-1, Qend+1) |
                                       Qseqid == MATE & Sseqid == REFERENCE & between(POS, Sstart-1, Send+1) & between(POSEND, Sstart-1, Send+1), "out", "in")) %>%
         filter(out == "in")
       
 }
 
 paired_parental_filter = df %>% distinct(POS, MPOS, .keep_all = T) ## remove PCR duplicates
 
 rm(df, Qseqid, Qstart, Qend, Sstart, Send, Sseqid, j)

 ## get recombining regions by overlapping reads
 
 not_properly_paired = paired_parental_filter %>%
  arrange(REFERENCE, MATE, POS) %>%
  tibble::rowid_to_column("i") 

recombining_regions = data.frame(matrix(ncol = 8, nrow = 0))
colnames(recombining_regions) = c('i', 'REFERENCE', 'POS', 'POSEND', 'MATE', 'MPOS', 'MPOSEND', 'READSJOINED')
recombining_regions$REFERENCE = as.character(recombining_regions$REFERENCE)
recombining_regions$MATE = as.character(recombining_regions$MATE)
x=0
  
while (dim(not_properly_paired)[1]>0) {
   read = as.character(not_properly_paired$READ[1])
   ref = as.character(not_properly_paired$REFERENCE[1])
   rstart = as.numeric(not_properly_paired$POS[1])
   rend = as.numeric(not_properly_paired$POSEND[1])
   mate = as.character(not_properly_paired$MATE[1])
   mstart = as.numeric(not_properly_paired$MPOS[1])
   mend = as.numeric(not_properly_paired$MPOSEND[1])

   row = not_properly_paired[1,]
   row_prev = not_properly_paired[0,]
   not_properly_paired = not_properly_paired[-c(1), ]
   df = not_properly_paired
   x = x+1
   
   if(dim(df)[1]>0){
     while(dim(row_prev)[1] < dim(row)[1] ) {
              
      row_prev = row
      df2 = df
      
    if(dim(df2)[1]>0){ 
      for (j in 1:(dim(df2)[1])) {
        
        if(ref == df2$REFERENCE[j] & mate == df2$MATE[j] & (between(df2$POS[j], rstart - isize, rend + isize) & between(df2$MPOS[j], mstart - isize, mend + isize) | 
                                                    between(df2$POS[j], rstart - isize, rend + isize) & between(df2$MPOSEND[j], mstart - isize, mend + isize) |     
                                                    between(df2$POSEND[j], rstart - isize, rend + isize) & between(df2$MPOS[j], mstart - isize, mend + isize) |
                                                    between(df2$POSEND[j], rstart - isize, rend + isize) & between(df2$MPOSEND[j], mstart - isize, mend + isize))) {
            
        row = bind_rows(row, df2[j,])  
        id = df2$i[j]
        not_properly_paired = not_properly_paired  %>% filter(i != id)
        df = df  %>% filter(i != id)
        rstart=min(row$POS)
        rend=max(row$POSEND)
        mstart=min(row$MPOS)
        mend=max(row$MPOSEND)
        
        }
       }
      }
     }
    }
  
   if (dim(row)[1] > 1) {
     n = row %>% distinct(READ) %>% summarise(n=n())
     assign(paste("reads_joined",x,sep = "_"), row)
       if(dim(row %>% filter(grepl(read_length,READ_LENGTH)))[1] > 0) {
           row = row %>% group_by() %>% summarise(i= x, REFERENCE = last(REFERENCE), POS = min(POS), POSEND = max(POSEND), MATE = last(MATE), MPOS = min(MPOS), MPOSEND = max(MPOSEND)) %>% mutate(READSJOINED = n[1] %>% as.numeric())
           recombining_regions = bind_rows(recombining_regions, row) 
      
     }
   }
}

rm(x, not_properly_paired, id, j, row, row_prev, mate, mstart, mend, ref, rstart, rend, df, df2, read)

## add information of overlapping repeats

recombining_regions = recombining_regions %>% mutate(LENGTH = POSEND - POS) %>% mutate(repeat_id = 0, repeat_length = 0, repeat_pident = 0, repeat_strand = 0, overlap_repeat_left = 0, overlap_repeat_right =0) 

for(j in 1:dim(recombining_regions)[1]) {
  i=recombining_regions$i[j]
  REFERENCE=recombining_regions$REFERENCE[j]
  MATE=recombining_regions$MATE[j]
  POS=recombining_regions$POS[j]
  POSEND=recombining_regions$POSEND[j]
  MPOS=recombining_regions$MPOS[j]
  MPOSEND=recombining_regions$MPOSEND[j]
  LENGTH=recombining_regions$LENGTH[j]
  
  df1 = repeat_regions_joined %>%  mutate(recombined1 = ifelse(REFERENCE == qseqid & MATE == sseqid & POS < qstart & POSEND > qend, 
                                         "R1R2", 
                                         ifelse(REFERENCE == qseqid & MATE == sseqid & between(POSEND, qstart - (isize-read_length*2), qend + (isize+read_length*2)),
                                          "R1",
                                         ifelse(REFERENCE == qseqid & MATE == sseqid & between(POS, qstart - (isize-read_length*2), qend + (isize+read_length*2)),
                                          "R2", 0)))) %>%
                                  mutate(recombined2 = ifelse(REFERENCE == qseqid & MATE == sseqid & MPOS < sstart - (isize-read_length*2) & MPOSEND > send + (isize+read_length*2), 
                                         "R1R2", 
                                         ifelse(REFERENCE == qseqid & MATE == sseqid & between(MPOSEND, sstart - (isize-read_length*2), send + (isize+read_length*2)),
                                          "R2",
                                         ifelse(REFERENCE == qseqid & MATE == sseqid & between(MPOS, sstart - (isize-read_length*2), send + (isize+read_length*2)),
                                          "R1", 0)))) %>%
                                  filter(recombined1 != 0 & recombined2 != 0)
  
  if(dim(df1)[1]>0){
  recombining_regions$repeat_id[j] = df1 %>% summarize(text = str_c(i, collapse = ";")) %>% as.character()
  recombining_regions$repeat_length[j] = df1 %>% summarize(text = str_c(length, collapse = ";")) %>% as.character()
  recombining_regions$repeat_pident[j] = df1 %>% summarize(text = str_c(pident, collapse = ";")) %>% as.character()
  recombining_regions$repeat_strand[j] = df1 %>% summarize(text = str_c(sstrand, collapse = ";")) %>% as.character()
  recombining_regions$overlap_repeat_left[j] = df1 %>% summarize(text = str_c(recombined1, collapse = ";")) %>% as.character()
  recombining_regions$overlap_repeat_right[j] = df1 %>% summarize(text = str_c(recombined2, collapse = ";")) %>% as.character()
  } 
}

write.table(recombining_regions %>% as_tibble(), "./outputs/recombining_regions.txt", sep="\t", row.names = F, quote=F)

## evaluate evidence for recombination across repeats

blast_parental = repeat_regions_joined

blast_parental$AB2 = 0
blast_parental$A4 = 0
blast_parental$B4 = 0
blast_parental$p1top2 = 0
blast_parental$p2top1 = 0

for (i in 1:(dim(blast_parental)[1])) {
  k=as.character(blast_parental$i[i])
  p1=as.character(blast_parental$qseqid[i])
  p1start=as.numeric(blast_parental$qstart[i])
  p1end=as.numeric(blast_parental$qend[i])
  p2=as.character(blast_parental$sseqid[i])
  p2start=as.numeric(blast_parental$sstart[i])
  p2end=as.numeric(blast_parental$send[i])
  strand=as.character(blast_parental$sstrand[i])
  length = as.numeric(blast_parental$qend[i] - blast_parental$qstart[i] + 1)
  
          p1_seq = toupper(c2s(p1_fasta[[1]][p1start:p1end])) ## get sequences by coordinates
          
          if(strand == "plus") {
              p2_seq = toupper(c2s(p2_fasta[[1]][p2start:p2end]))
          } else {
              p2_seq = toupper(c2s(rev(comp(p2_fasta[[1]][p2start:p2end])))) ## get reverse complement
          }
  
          msA = DNAStringSet(c(p1_seq, p2_seq)) %>% as.DNAbin()
          dfA = muscle(msA) %>% as.character() %>% t() %>% as.data.frame()
          colnames(dfA) = c("p1", "p2")
          
          rm(msA)
          
          ## add real coordinates
          
          dfN = data.frame(matrix("N",   
                          nrow = 1000,
                          ncol = 2))
          colnames(dfN) = c("p1", "p2")
          
          df = bind_rows(dfN,dfA, dfN) %>% mutate(p1_coord=0, p2_coord=0)
          
          counT = p1start-1000
          for(l in 1:nrow(df)) {
              if(df[l,1] != "-") {
                df$p1_coord[l] = counT
                counT = counT + 1
              } else {
                df$p1_coord[l] = counT
            }
          }
          
          if(strand == "plus"){
            counT = p2start-1000
            for(l in 1:nrow(df)) {
                if(df[l,2] != "-") {
                  df$p2_coord[l] = counT
                  counT = counT + 1
                } else {
                  df$p2_coord[l] = counT
            }}} else {
            counT = p2end+1000
            for(l in 1:nrow(df)) {
                if(df[l,2] != "-") {
                  df$p2_coord[l] = counT
                  counT = counT - 1
                } else {
                  df$p2_coord[l] = counT
                }}}
          
        
          
          ## add coverage to alignment
          
          df = left_join(df, depth_file %>% filter(p==p1_name), by=c("p1_coord"="pos")) %>% select(p1, p2, p1_coord, p2_coord, cov_p1 = depth) %>%
            left_join(., depth_file %>% filter(p==p2_name), by=c("p2_coord"="pos")) %>% select(p1, p2, p1_coord, p2_coord, cov_p1, cov_p2 = depth) %>% 
            replace(is.na(.), 0)
          
          ## add gene regions to alignment
          
          df$p1_genes = 0
          df$p2_genes = 0
          
          for(l in 1:nrow(df)) {
            p1_coord=df$p1_coord[l]
            p2_coord=df$p2_coord[l]
            
              gene_p1_align = p1_genes %>% mutate(cov = ifelse(between(p1_coord, start, end), 1, 0))
              gene_p2_align = p2_genes %>% mutate(cov = ifelse(between(p2_coord, start, end), 1, 0))
            
              if(sum(gene_p1_align$cov) > 0){
                gene_name = gene_p1_align %>% filter(cov > 0) %>% select(gene) %>% group_by() %>% aggregate(gene ~ gene, data = ., FUN = paste, collapse = "; ") %>% as.character()
                df$p1_genes[l] = gene_name
              }           
              if(sum(gene_p2_align$cov) > 0){
                gene_name = gene_p2_align %>% filter(cov > 0) %>% select(gene) %>% group_by() %>% aggregate(gene ~ gene, data = ., FUN = paste, collapse = "; ") %>% as.character()
                df$p2_genes[l] = gene_name
              } 
          }
          
          SNP = df %>% 
                tibble::rowid_to_column("alg_pos") %>%
                mutate(SNP = ifelse(p1 != p2 | p1 == "N" & p2 == "N", "X", "-")) %>% filter(SNP=="X")
          
          assign(paste("align_",k,sep = ""), tibble::rowid_to_column(df, "alg_pos")) 
          assign(paste("SNP_",k,sep = ""), SNP) 
          rm(l, counT, p1_seq,  p2_seq, df, SNP, p1_coord)
          
          ## analyze intergenic connections
        
          df = paired_parental_filter %>% filter(between(POS, p1start - 1000, p1end + 1000) & between(MPOS, p2start - (1000 + max_isize), p2end + 1000 + max_isize)) %>% mutate(r1end=POS+read_length-1,r2end=MPOS+read_length-1) %>% select(rname=READ,r1=REFERENCE, r1start=POS, r1end, r2=MATE, r2start=MPOS, r2end) %>% mutate(A2_SNP=0, B2_SNP=0, p1_SNP=0, p2_SNP=0)

          if (dim(df)[1]>0) {
          
          ## A is from p1 to p2, B from p2 to p1
          ## 1, first read before R start
          ## 2, both reads within repeat
          ## 3, second read after R end
          ## 4, both reads outside the repeat (max repeat size ~450 bp)
  
    
          for (j in 1:(dim(df)[1])) {
            rname=df[j,1]
            r1 = df[j,2]
            r1start = df[j,3]
            r1end = df[j,4]
            r2 = df[j,5]
            r2start = df[j,6]
            r2end = df[j,7]
            
            A2 = 0
            B2 = 0
            A4 = 0
            B4 = 0
            lA2 = c()
            lB2 = c()
            
              if(length > innersize) { ## A2, B2 idem plus and minus
                    
                lA2_tmp = c()
                lB2_tmp = c()
                
                    for (l in 1:(dim(get(paste("SNP_",k,sep = "")))[1])) {
                    
                    p1SNP = get(paste("SNP_",k,sep = ""))[l,4]
                    p2SNP = get(paste("SNP_",k,sep = ""))[l,5]
                    alg_pos = get(paste("SNP_",k,sep = ""))[l,1]
                    
                            if (between(p1SNP,r1start,r1end)) {
                              lA2_tmp = c(lA2_tmp, alg_pos) ## A2 p1 SNPs, B2 p2 SNPs
                            } else if (between(p2SNP,r2start,r2end)) {
                              lB2_tmp = c(lB2_tmp, alg_pos)
                            }
                        }
                    
                    if(length(lA2_tmp) >= 1 & length(lB2_tmp) >= 1) {
                          lA2 = lA2_tmp
                          lB2 = lB2_tmp
                          A2 = 1
                          B2 = 1
                    }

              }
            
                
                rm(lA2_tmp, lB2_tmp)
                
            
            if(length <= innersize) {
              
                if (strand == "plus" & 
                    p1 == r1 & 
                    between(r1start, p1end - max_isize, p1start - min_flanking) & 
                    between(r2start, p2end - (read_length-min_flanking), p2start + max_isize - read_length)) { ## A4 plus
                      
                  A4 = 1
                  
              } else if (strand == "minus" & 
                         p1 == r1 & 
                         between(r1start, p1end - max_isize, p1start - min_flanking) & 
                         between(r2start, p2end - max_isize, p2start - min_flanking)) { ## A4 minus
                      
                  A4 = 1
                  
              } else if (strand == "plus" & 
                         p1 == r1 & 
                         between(r2start, p2end - max_isize, p2start - min_flanking) & 
                         between(r1start, p1end - (read_length-min_flanking), p1start + max_isize - read_length)) { ## B4 plus
                      
                  B4 = 1

              } else if (strand == "minus" & 
                         p1 == r1 &
                         between(r2start, p2end - (read_length-min_flanking), p2start + max_isize - read_length) & 
                         between(r1start, p1end - (read_length-min_flanking), p1start + max_isize - read_length)) { ## B4 minus
                      
                  B4 = 1

              }
            }
          
          df$A2_SNP[j] = toString(unique(lA2))
          df$B2_SNP[j] = toString(unique(lB2))
          df$AB2[j] = A2
          df$A4[j] = A4
          df$B4[j] = B4
          
             if(A2>0) {
  
                  if(lA2[1] > lB2[1]) {
                  
                  df$p1_SNP[j] = lA2[1]
                  df$p2_SNP[j] = rev(lB2)[1]
                  
                  } else {
                  
                  df$p1_SNP[j] = rev(lA2)[1]
                  df$p2_SNP[j] = lB2[1]
                  
                  }
              
            } else if (A4>0) {
  
              df$p1_SNP[j] = ifelse(r1end < p1start, 1000-(p1start-r1end), 1000)
              df$p2_SNP[j] = ifelse(r2start > p2end, 1000 + length + (r2start-p2end), 1000 + length)
              
            } else if (B4>0) {
  
              df$p1_SNP[j] = ifelse(r1start > p1end, 1000 + length + (r1start-p1end), 1000 + length)
              df$p2_SNP[j] = ifelse(r2end < p2start, 1000-(p2start-r2end), 1000)
              
            }
          }
          
          df = df %>% 
            filter(A2_SNP > 0 | B2_SNP > 0 | A4 > 0 | B4 > 0 ) 
          
          if(dim(df)[1]>0) {
          
          df = df %>%
            mutate(direction = ifelse(p1_SNP < p2_SNP, "p1 > p2", "p1 < p2")) %>% 
            mutate(p1top2 = ifelse(p1_SNP < p2_SNP, 1, 0)) %>% 
            mutate(p2top1 = ifelse(p1_SNP > p2_SNP, 1, 0)) %>%
            mutate(insert_size = abs(p1_SNP - p2_SNP) + read_length) %>% 
            filter(between(insert_size, -max_isize, max_isize)) %>% 
            select(-insert_size)
          
          }

          if(dim(df)[1]>0) {
            
            assign(paste("reads_",k,sep = ""), df)
            
            df2 = df %>% select("AB2", "A4", "B4", "p1top2", "p2top1") %>% apply( 2, function(x) sum(x>0))
            
            blast_parental$AB2[i] = df2[1]
            blast_parental$A4[i] = df2[2]
            blast_parental$B4[i] = df2[3]
            blast_parental$p1top2[i] = df2[4]
            blast_parental$p2top1[i] = df2[5]
              
            }
          }
}

          

rm(i, j, lA1, lA2, lA3, lB1, lB2, lB3, A1, A2, A3, A4, B1, B2, B3, B4, length, strand, p1, p1end, p1start, p2, p2end, p2start , p1SNP, p2SNP, r1, r1end, r1start, r2, r2end, r2start, df, df2)

blast_parental[is.na(blast_parental)] = 0
cols_names = c("p1top2","p2top1")
blast_parental[cols_names] = sapply(blast_parental[cols_names],as.numeric)


blast_parental_evidence = blast_parental %>% mutate(total = AB2+A4+B4) %>% 
  mutate(total_p1top2 = p1top2 + A4) %>%
  mutate(total_p2top1 = p2top1 + B4) %>%
  filter(total_p1top2>0 | total_p2top1>0) %>%
  select(i, qseqid, sseqid, qstart, qend, sstart, send, length, pident, sstrand)
  

write.table(blast_parental_evidence, "./outputs/blast_parental_evidence.txt", sep="\t", row.names = F, quote=F)

## join recombinations over repeats

## by p1

blast_parental_evidence_joined = blast_parental_evidence[FALSE,]
blast_parental_evidence_tmp = blast_parental_evidence
blast_parental_evidence_joined$i = as.character(blast_parental_evidence_joined$i)
blast_parental_evidence_tmp$i = as.character(blast_parental_evidence_tmp$i)
blast_parental_evidence_joined$sstart = as.character(blast_parental_evidence_joined$sstart)
blast_parental_evidence_tmp$sstart = as.character(blast_parental_evidence_tmp$sstart)
blast_parental_evidence_joined$send = as.character(blast_parental_evidence_joined$send)
blast_parental_evidence_tmp$send = as.character(blast_parental_evidence_tmp$send)

while (dim(blast_parental_evidence_tmp)[1]>0) {
   qseqid=blast_parental_evidence_tmp$qseqid[1]
   qstart=as.numeric(blast_parental_evidence_tmp$qstart[1])
   qend=as.numeric(blast_parental_evidence_tmp$qend[1])
   id1=blast_parental_evidence_tmp$i[1]

   row = blast_parental_evidence_tmp[1,]
   blast_parental_evidence_tmp = blast_parental_evidence_tmp[-c(1), ]
   df = blast_parental_evidence_tmp
   
  if(dim(df)[1]>0){
      for (j in 1:(dim(df)[1])) {
      if(qseqid == df$qseqid[j] & qstart == df$qstart[j] & qend == df$qend[j])  {
                                                          
                                                            row = bind_rows(row, df[j,])
                                                            id2 = df$i[j]
                                                            blast_parental_evidence_tmp = blast_parental_evidence_tmp %>% filter(i != id2)
                                                            
       align = full_join(get(paste0("align_", id1, sep = "")), get(paste0("align_", id2, sep = "")), by = c("alg_pos" , "p1", "p2", "p1_coord", "cov_p1", "p1_genes"),  keep = F) %>% 
       unite(p2_coord,c("p2_coord.x", "p2_coord.y"), sep = "/") %>%
       mutate(cov_p2 = cov_p2.x + cov_p2.y) %>%
       mutate(p2_genes = ifelse(p2_genes.x != 0, p2_genes.x, p2_genes.y)) %>%
       select(alg_pos, p1, p2, p1_coord, p2_coord, cov_p1, cov_p2, p1_genes, p2_genes)
       
       assign(paste("align_",id1, "/", id2,sep = ""), align)
       
       reads = bind_rows(get(paste0("reads_", id1, sep = "")), get(paste0("reads_", id2, sep = "")))
       
       assign(paste("reads_",id1, "/", id2,sep = ""), reads)
         
       SNP = bind_rows(get(paste0("SNP_", id1, sep = "")), get(paste0("SNP_", id2, sep = ""))) %>% distinct(alg_pos, .keep_all = T)
       
       assign(paste("SNP_",id1, "/", id2,sep = ""), SNP)
       
       id1 = paste(id1, id2, sep="/")
       
      }
    }
  }
  
   if(dim(row)[1] > 1) {
     row_join = row %>% group_by() %>% mutate(i=paste0(i, collapse = "/")) %>% mutate(sstart=paste0(sstart, collapse = "/")) %>% mutate(send=paste0(send, collapse = "/")) %>% summarise(i=last(i), 
                   qseqid = last(qseqid), 
                   sseqid =  last(sseqid), 
                   qstart = last(qstart), 
                   qend = last(qend),
                   sstart = last(sstart), 
                   send = last(send), 
                   length = last(length), 
                   pident = last(pident), 
                   sstrand = last(sstrand))
     
     blast_parental_evidence_joined = bind_rows(blast_parental_evidence_joined, row_join)
     row_join_i = row_join[1] %>% as.vector()
     

   } else {
     blast_parental_evidence_joined = bind_rows(blast_parental_evidence_joined, row)
   }
}

rm(row, row_join, row_join_i, id1, id2, j, blast_parental_evidence_tmp, df,qseqid, qstart, qend, sstart, send)

## by p2

blast_parental_evidence_joined2 = blast_parental_evidence_joined[FALSE,]
blast_parental_evidence_tmp = blast_parental_evidence_joined
blast_parental_evidence_joined2$i = as.character(blast_parental_evidence_joined2$i)
blast_parental_evidence_tmp$i = as.character(blast_parental_evidence_tmp$i)
blast_parental_evidence_joined2$sstart = as.character(blast_parental_evidence_joined2$sstart)
blast_parental_evidence_tmp$sstart = as.character(blast_parental_evidence_tmp$sstart)
blast_parental_evidence_joined2$send = as.character(blast_parental_evidence_joined2$send)
blast_parental_evidence_tmp$send = as.character(blast_parental_evidence_tmp$send)
blast_parental_evidence_joined2$qstart = as.character(blast_parental_evidence_joined2$qstart)
blast_parental_evidence_tmp$qstart = as.character(blast_parental_evidence_tmp$qstart)
blast_parental_evidence_joined2$qend = as.character(blast_parental_evidence_joined2$qend)
blast_parental_evidence_tmp$qend = as.character(blast_parental_evidence_tmp$qend)

while (dim(blast_parental_evidence_tmp)[1]>0) {
   sseqid=blast_parental_evidence_tmp$sseqid[1]
   sstart=blast_parental_evidence_tmp$sstart[1]
   send=blast_parental_evidence_tmp$send[1]
   id1=blast_parental_evidence_tmp$i[1]

   row = blast_parental_evidence_tmp[1,]
   blast_parental_evidence_tmp = blast_parental_evidence_tmp[-c(1), ]
   df = blast_parental_evidence_tmp
   
  if(dim(df)[1]>0){
      for (j in 1:(dim(df)[1])) {
      if(sseqid == df$sseqid[j] & sstart == df$sstart[j] & send == df$send[j])  {
                                                          
                                                            row = bind_rows(row, df[j,])
                                                            id2 = df$i[j]
                                                            blast_parental_evidence_tmp = blast_parental_evidence_tmp %>% filter(i != id2)
                                                            
                                                            
       align = full_join(get(paste0("align_", id1, sep = "")), get(paste0("align_", id2, sep = "")), by = c("alg_pos" , "p1", "p2", "p2_coord", "cov_p2", "p2_genes"),  keep = F) %>% 
       unite(p1_coord,c("p1_coord.x", "p1_coord.y"), sep = "/") %>%
       mutate(cov_p1 = cov_p1.x + cov_p1.y) %>%
       mutate(p1_genes = ifelse(p1_genes.x != 0, p1_genes.x, p1_genes.y)) %>%
       select(alg_pos, p1, p2, p1_coord, p2_coord, cov_p1, cov_p2, p1_genes, p2_genes)
       
       assign(paste("align_",id1, "/", id2,sep = ""), align)  
       
       reads = bind_rows(get(paste0("reads_", id1, sep = "")), get(paste0("reads_", id2, sep = "")))
       
       assign(paste("reads_",id1, "/", id2,sep = ""), reads)
       
       SNP = bind_rows(get(paste0("SNP_", id1, sep = "")), get(paste0("SNP_", id2, sep = ""))) %>% distinct(alg_pos, .keep_all = T)
       
       assign(paste("SNP_",id1, "/", id2,sep = ""), SNP)
       
       id1 = paste(id1, id2, sep="/")
       
      }
    }
  }
  
   if(dim(row)[1] > 1) {
     row_join = row %>% group_by() %>% mutate(i=paste0(i, collapse = "/")) %>% mutate(qstart=paste0(qstart, collapse = "/")) %>% mutate(qend=paste0(qend, collapse = "/")) %>% summarise(i=last(i), 
                   qseqid = last(qseqid), 
                   sseqid =  last(sseqid), 
                   qstart = last(qstart), 
                   qend = last(qend),
                   sstart = last(sstart), 
                   send = last(send), 
                   length = last(length), 
                   pident = last(pident), 
                   sstrand = last(sstrand))
     
     blast_parental_evidence_joined2 = bind_rows(blast_parental_evidence_joined2, row_join)
     row_join_i = row_join[1] %>% as.vector()
     

   } else {
     blast_parental_evidence_joined2 = bind_rows(blast_parental_evidence_joined2, row)
   }
}


blast_parental_evidence_joined = blast_parental_evidence_joined2

rm(row, row_join, row_join_i, id1, id2, j, blast_parental_evidence_tmp, blast_parental_evidence_joined2, df,qseqid, qstart, qend, sstart, send)


## get breakpoints ranges and filter blast parental evidence

blast_parental_evidence_filter = blast_parental_evidence_joined 
blast_parental_evidence_filter$total = 0
blast_parental_evidence_filter$genes = 0

all_breakpoints = as_tibble()

for(l in 1:dim(blast_parental_evidence_filter)[1]) {
        
    i=blast_parental_evidence_filter$i[l]
    length = as.numeric(dim(get(paste("align_",i,sep = "")))[1])
    genes_p1 = get(paste("align_",i,sep = "")) %>% select(genes = p1_genes) 
    genes_p1$genes = as.character(genes_p1$genes)
    genes_p2 = get(paste("align_",i,sep = "")) %>% select(genes = p2_genes) 
    genes_p2$genes = as.character(genes_p2$genes)
    genes = bind_rows(genes_p1, genes_p2) %>%
      filter(genes != 0) %>%
      distinct() %>%
      arrange(genes) %>%
      summarise(str_c(genes, collapse = ";")) %>% as.character()
    
    genes = ifelse(genes == "", "no", genes) %>% as.character()
    align = get(paste("align_",i,sep = ""))
    
    df3 = as_tibble()
    
    df_p1 = get(paste("reads_",i,sep = ""))  %>% select(p1_SNP, p2_SNP, direction) %>% arrange(direction, p1_SNP) %>% ungroup() %>% filter(direction == "p1 > p2")
    df_p2 = get(paste("reads_",i,sep = ""))  %>% select(p1_SNP, p2_SNP, direction) %>% arrange(direction, p2_SNP) %>% ungroup() %>% filter(direction == "p1 < p2")

    if(dim(df_p1)[1]>0) {
      
        p1_max = df_p1$p1_SNP[1]
        p2_min = df_p1$p2_SNP[1]
        df_p1$group = 0
        count = 0
        
        for (j in 1:dim(df_p1)[1]) {
          p1_SNP = df_p1$p1_SNP[j]
          p2_SNP = df_p1$p2_SNP[j]
          
          if(p1_SNP >= p2_min) {
            df_p1$group[j] = count + 1
            count = count + 1
            p2_min = p2_SNP
            p1_max = p1_SNP
            
          } else if(p1_SNP < p2_min) {
              
              if(p2_SNP < p2_min){
                p2_min = p2_SNP
              }
              
              if(p1_SNP > p1_max) {
                p1_max = p1_SNP
              }
              
              df_p1$group[j] = count
            }
        } 
        
        df1 = df_p1 %>% 
          group_by(group) %>% 
          summarise(start = max(p1_SNP), end = min(p2_SNP), reads = n(), direction = last(direction)) 
        
        df3 = bind_rows(df3,df1)
    }
    
    if(dim(df_p2)[1]>0) {
    
          p1_min = df_p2$p1_SNP[1]
          p2_max = df_p2$p2_SNP[1]
          df_p2$group = 0
          count = 0
          
          for (j in 1:dim(df_p2)[1]) {
            p1_SNP = df_p2$p1_SNP[j]
            p2_SNP = df_p2$p2_SNP[j]
            
            if(p2_SNP >= p1_min) {
              df_p2$group[j] = count + 1
              count = count + 1
              p1_min = p1_SNP
              p2_max = p2_SNP
              
            } else if(p2_SNP < p1_min) {
                
                if(p1_SNP < p1_min){
                  p1_min = p1_SNP
                }
                
                if(p2_SNP > p2_max) {
                  p2_max = p2_SNP
                }
                
                df_p2$group[j] = count
              }
          }
          
          df2 = df_p2 %>% 
            group_by(group) %>% 
            summarise(start = min(p1_SNP), end = max(p2_SNP), reads = n(), direction = last(direction)) 
          
          df3 = bind_rows(df3,df2)
    }

  df3 = df3 %>% mutate(i=i) 
  
  for(f in 1:dim(df3)[1])  {
    bk_start = df3$start[f]
    bk_end = df3$end[f]
    
    df3$p1_bk_start[f] = align %>% filter(alg_pos == bk_start) %>% select(p1_coord) %>% as.character()
    df3$p1_bk_end[f] = align %>% filter(alg_pos == bk_end) %>% select(p1_coord) %>% as.character()
    df3$p2_bk_start[f] = align %>% filter(alg_pos == bk_start) %>% select(p2_coord) %>% as.character()
    df3$p2_bk_end[f] = align %>% filter(alg_pos == bk_end) %>% select(p2_coord) %>% as.character()
  }

  df = bind_rows(df_p1,df_p2)  
  df4 = df %>% inner_join(., df3, by=c("group", "direction")) %>% select(-group) 

  df5 = df3 %>%  mutate(start2 = ifelse(direction == "p1 > p2" & start-1000+1 > 1, start-1000+1,
                                                              ifelse(direction == "p1 < p2" & end-1000+1 > 1, end-1000+1, 1)),
                                       end2=ifelse(direction == "p1 > p2" & end-1000+1 < length, end-1000+1, 
                                                   ifelse(direction == "p1 < p2" & start-1000+1 < length, start-1000+1, length))) %>% 
                                select(-group) %>% arrange(start2)
  
  all_breakpoints = bind_rows(all_breakpoints, df5)
  
  assign(paste("R_breakpoints_",i,sep = ""), df3 %>% select(-group))
  
  assign(paste("R_sites_",i,sep = ""), df4)
  
  total = df3 %>% summarize(sum(reads)) %>% as.numeric()
  
  blast_parental_evidence_filter$total[l] = total
  blast_parental_evidence_filter$genes[l] = genes
    
  rm(df, df1, df2, df3, df4, df5, df_p1, df_p2, genes_p1, genes_p2, genes)
  
}

all_breakpoints =  all_breakpoints %>%
  select(i,bk_start=start2, bk_end=end2, reads, direction, p1_bk_start, p1_bk_end, p2_bk_start, p2_bk_end)
write.table(all_breakpoints, "./outputs/all_breakpoints.txt", sep="\t", row.names = F, quote=F)

blast_parental_evidence_filter=full_join(blast_parental_evidence_filter, all_breakpoints, by="i") %>% filter(reads > 2)
write.table(blast_parental_evidence_filter, "./outputs/blast_parental_evidence_filter.txt", sep="\t", row.names = F, quote=F)


## basic stats

## calculate parental length without LR (>1000 bp)

LR_regions = read.delim(file="./inputs/LR.txt", header=T) 

LR_2_substract = LR_regions %>% filter(status == "repeated") %>% select(chrom=p, start, end)

genome_lengths_woR = LR_regions %>% 
  filter(status == "repeated") %>% 
  group_by(p) %>% 
  summarise(repeat_length = sum(length)) %>% 
  mutate(length=ifelse(p == p1_name, p1_length, p2_length)) %>% 
  mutate(length_woR = length-repeat_length)

## get homologous regions not repeated (and exclusive by substracting)

hom_regions_p1 = repeat_regions %>% 
  select(chrom=qseqid, start=qstart, end=qend, pident) %>% 
  arrange(start) %>% 
  distinct()

hom_regions_p2 = repeat_regions %>% 
  select(chrom=sseqid, start=sstart, end=send, pident) %>% 
  arrange(start) %>% 
  distinct()

hom_regions_merged = bind_rows(hom_regions_p1, hom_regions_p2) %>% 
  valr::bed_merge(., avg_pident = mean(pident)) %>% 
  valr::bed_subtract(., LR_2_substract) 

hom_regions = bind_rows(hom_regions_p1, hom_regions_p2) 

total_hom_regions_count = hom_regions %>% 
  mutate(length = end - start +1) %>%
  group_by(chrom) %>% 
  tally() %>%
  select(total_hom_regions_count =n)

total_hom_regions_count_100 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length<100) %>%
  group_by(chrom) %>% 
  tally() %>%
  select(total_hom_regions_count_100 =n)

hom_regions_pident_100 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length<100) %>%
  group_by(chrom) %>% 
  summarise(hom_regions_100_pident = mean(pident)) %>% select(2)

total_hom_regions_count_100_1000 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length>=100 & length <1000) %>%
  group_by(chrom) %>% 
  tally() %>%
  select(total_hom_regions_count_100_1000 =n)

hom_regions_pident_100_1000 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length>=100 & length<1000) %>%
  group_by(chrom) %>% 
  summarise(hom_regions_100_1000_pident = mean(pident))  %>% select(2)

total_hom_regions_count_1000 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length >= 1000) %>%
  group_by(chrom) %>% 
  tally() %>%
  select(total_hom_regions_count_1000 =n)

hom_regions_pident_1000 = hom_regions %>% 
  mutate(length = end - start +1) %>%
  filter(length>=1000) %>%
  group_by(chrom) %>% 
  summarise(hom_regions_1000_pident = mean(pident))  %>% select(2)

total_length_hom_regions = hom_regions_merged %>%
  mutate(length = end - start +1) %>%
  group_by(chrom) %>% 
  summarise(hom_regions = sum(length))  %>% select(2)

genome_lengths_woR = bind_cols(genome_lengths_woR, total_length_hom_regions , total_hom_regions_count, total_hom_regions_count_100, hom_regions_pident_100, total_hom_regions_count_100_1000, hom_regions_pident_100_1000,  total_hom_regions_count_1000, hom_regions_pident_1000) %>%
  mutate(excl_regions =  length_woR - hom_regions)


## get identical regions not repeated 

id_regions_p1 = identical_regions %>% 
  select(chrom=qseqid, start=qstart, end=qend) %>% 
  arrange(start) %>% 
  distinct()

id_regions_p2 = identical_regions %>% 
  select(chrom=sseqid, start=sstart, end=send) %>% 
  arrange(start) %>% 
  distinct()

id_regions = bind_rows(id_regions_p1, id_regions_p2) %>% 
  valr::bed_merge(.) %>% 
  valr::bed_subtract(., LR_2_substract) 

total_id_regions = id_regions %>% mutate(length = end - start +1) %>% 
  group_by(chrom) %>% 
  summarise(id_regions = sum(length))%>%
  select(id_regions)

genome_lengths_woR = bind_cols(genome_lengths_woR, total_id_regions) %>%
  mutate(id_perc = id_regions/hom_regions*100)

## get parental regions (greater than 1000 bp) retained in the hybrid (subtracting LR)

parental_regions_by_cov_1000 = parental_regions_by_cov %>% 
  filter(length >= 1000)  %>% 
  select(chrom=p, start, end)

retained_regions = valr::bed_subtract(parental_regions_by_cov_1000, LR_2_substract) 

total_retained_regions = retained_regions %>% 
  mutate(length = end - start +1) %>% 
  group_by(chrom) %>% 
  summarise(retained = sum(length)) %>% 
  select(retained)

genome_lengths_woR = bind_cols(genome_lengths_woR, total_retained_regions) %>% 
  mutate(retained_perc = retained/length_woR*100) %>%
  mutate(lost =  length_woR - retained) %>%
  mutate(lost_perc = lost/length_woR*100) 

retained_hom_regions = valr::bed_intersect(retained_regions, hom_regions_merged) %>% 
  group_by(chrom) %>% 
  summarise(retained_hom = sum(.overlap)) %>% 
  select(retained_hom) 

genome_lengths_woR = bind_cols(genome_lengths_woR, retained_hom_regions) %>% 
  mutate(retained_hom_perc = retained_hom/hom_regions*100) %>%
  mutate(retained_excl = retained-retained_hom) %>%
  mutate(retained_excl_per = retained_excl /excl_regions*100)

write.table(genome_lengths_woR , "./outputs/stats_parental_contribution.txt", sep="\t", row.names = F, quote=F)

## get mapped reads

mapped_reads = data.frame(concordant_p1 = dim(sam_file_concordant %>% filter(REFERENCE==p1_name))[1], 
                          concordant_p2 = dim(sam_file_concordant %>% filter(REFERENCE==p2_name))[1]) %>%
  mutate(total_concordant = concordant_p1 +concordant_p2) %>%
  mutate(discordant = dim(sam_file_discordant)[1], filter = dim(paired_parental_filter)[1]) %>%
  mutate(total = total_concordant + discordant) %>% 
  mutate(filter_perc = filter/discordant*100)

write.table(mapped_reads, "./outputs/stats_mapped_reads.txt", sep="\t", row.names = F, quote=F)

## get genes coverage

p1_genes_cov = p1_genes %>%
  select(chrom=p, start, end, gene)

p2_genes_cov = p2_genes %>%
  select(chrom=p, start, end, gene)

genes = bind_rows(p1_genes_cov, p2_genes_cov)

genes_cov = valr::bed_intersect(genes, parental_regions_by_cov_1000) %>% 
  mutate(gene_length = end.x-start.x) %>% 
  select(chrom, start=start.x, end=end.x, gene=gene.x, gene_length, o_start = start.y, o_end=end.y, overlap=.overlap) %>%
  mutate(overlap_perc = overlap/gene_length*100)

write.table(genes_cov, "./outputs/gene_coverage.txt", sep="\t", row.names = F, quote=F)

## plot R_sites for recombination spots with >= min_depth reads evidence

library(jcolors)
library(ggtext)

blast_parental_evidence4plot = blast_parental_evidence_filter %>% 
  filter(total >= min_depth) %>%
  distinct(i,qseqid,sseqid,qstart, qend, sstart, send, length,pident,sstrand, total)

for(l in 1:dim(blast_parental_evidence4plot)[1]) {
   
  i = blast_parental_evidence4plot$i[l]
  qstart = blast_parental_evidence4plot$qstart[l]
  qend = blast_parental_evidence4plot$qend[l]
  sstart = blast_parental_evidence4plot$sstart[l]
  send = blast_parental_evidence4plot$send[l]
  length =  blast_parental_evidence4plot$length[l]
  length2 =  dim(get(paste("align_",i,sep = "")))[1]
  pident =  blast_parental_evidence4plot$pident[l]
  sstrand = blast_parental_evidence4plot$sstrand[l]
  total = blast_parental_evidence4plot$total[l]
  gene_cov_p1 = get(paste("align_",i,sep = "")) %>% select(alg_pos,p1_genes) %>% filter(p1_genes!="0") %>% group_by(p1_genes) %>% summarise(min = min(alg_pos), max = max(alg_pos)) %>% mutate(parent="p1")
  gene_cov_p2 = get(paste("align_",i,sep = "")) %>% select(alg_pos,p2_genes) %>% filter(p2_genes!="0") %>% group_by(p2_genes) %>% summarise(min = min(alg_pos), max = max(alg_pos)) %>% mutate(parent="p2")
  
  
   align_cov = get(paste("align_",i,sep = "")) %>% 
     mutate(cov_p1=ifelse(cov_p1==0, 0.25, ifelse(cov_p1>0 & cov_p1<100, (cov_p1*0.25)/100+0.25, 0.50))) %>% 
     mutate(cov_p2=ifelse(cov_p2==0, -0.25, ifelse(cov_p2>0 & cov_p2<100, -((cov_p2*0.25)/100+0.25), -0.50)))
   
  ## start and end of repeats

   if(grepl(";", i) == T & grepl("/", i) == T) {
     repeats_p1 = i %>% str_split("/", simplify = T) %>% str_split(";", simplify = T) %>% t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% select(i) %>% mutate(start=0, end=0, parent="p1") 
   } else if (grepl(";", i) == T) {
     repeats_p1 = i %>% str_split(";", simplify = T) %>%  t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% mutate(start=0, end=0, parent="p1")
   } else if (grepl("/", i) ==T) {
     repeats_p1 = i %>% str_split("/", simplify = T) %>%  t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% mutate(start=0, end=0, parent="p1")
   } else {
     repeats_p1 = i %>% as.tibble() %>% dplyr::rename(i=value) %>% mutate(start=0, end=0, parent="p1")
   }
   
  
  for(k in 1:dim(repeats_p1)[1]){
    i2=repeats_p1$i[k]
    df =  repeat_regions %>% filter(i==i2) 
    repeats_p1$start[k] =  align_cov %>% 
      filter(stringr::str_detect(as.character(p1_coord), as.character(df$qstart[1]))) %>%
      select(alg_pos) %>% 
      unlist()
    repeats_p1$end[k] =  align_cov %>% 
      filter(stringr::str_detect(as.character(p1_coord), as.character(df$qend[1]))) %>% 
      select(alg_pos) %>% 
      unlist()
  }
 
 
   if(grepl(";", i) == T & grepl("/", i) == T) {
     repeats_p2 = i %>% str_split("/", simplify = T) %>% str_split(";", simplify = T) %>% t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% select(i) %>% mutate(start=0, end=0, parent="p2") 
   } else if (grepl(";", i) == T) {
     repeats_p2 = i %>% str_split(";", simplify = T) %>%  t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% mutate(start=0, end=0, parent="p2")
   } else if (grepl("/", i) ==T) {
     repeats_p2 = i %>% str_split("/", simplify = T) %>%  t() %>% as.tibble() %>% dplyr::rename(i=V1) %>% mutate(start=0, end=0, parent="p2")
   } else {
     repeats_p2 = i %>% as.tibble() %>% dplyr::rename(i=value) %>% mutate(start=0, end=0, parent="p2")
   }
  
  for(k in 1:dim(repeats_p2)[1]){
    i2=repeats_p2$i[k]
    df =  repeat_regions %>% filter(i==i2) 
    repeats_p2$start[k] = align_cov %>% 
      dplyr::filter(stringr::str_detect(as.character(p2_coord), as.character(df$sstart[1]))) %>%
      select(alg_pos) %>% 
      unlist() 
    repeats_p2$end[k] = align_cov %>% 
      dplyr::filter(stringr::str_detect(as.character(p2_coord), as.character(df$send[1]))) %>%
      select(alg_pos) %>% 
      unlist() 
  }
 
 repeats_p2 = repeats_p2 %>% 
      mutate(start2=ifelse(start>end, end, start), end2=ifelse(start>end, start, end)) %>%
      select(i, start=start2, end=end2, parent) 
 
 repeats=bind_rows(repeats_p1, repeats_p2) %>% distinct(start, end, parent, .keep_all = T)

  breakpoints = get(paste("R_breakpoints_", i, sep="")) %>% 
  filter(reads >= min_depth) 
  
  SNP = get(paste("SNP_", i, sep="")) %>% filter(p1 != "N")

  R_sites = get(paste("R_sites_",i,sep = "")) %>% as_tibble() %>% 
  filter(reads >= min_depth) 
  
  df1 = R_sites %>%
    ungroup() %>%
    select(p1_SNP,direction,reads) %>% 
    mutate(p=p1_name) %>%
    mutate(y=0.25) %>%
    tibble::rowid_to_column("pair") %>%
    dplyr::rename(SNP = "p1_SNP")
  
  df2 = R_sites %>%
    ungroup() %>%
    select(p2_SNP,direction,reads)  %>% 
    mutate(p=p2_name) %>%
    mutate(y=-0.25) %>%
    tibble::rowid_to_column("pair") %>%
    dplyr::rename(SNP = "p2_SNP")
  
  df = rbind(df1,df2) 
  df$SNP = df$SNP %>% unlist() %>% as.integer()
  df$reads = df$reads %>% unlist() %>% as.integer()
  df$direction = df$direction %>% unlist() %>% as.character()
  
if(dim(breakpoints)[1]>0) {
  
    p = ggplot(df, aes(x=SNP, y=y)) + 
    ggtitle(paste("HR #",i, ", ", length, " bp, ", pident, "<br/>", "*N. tabacum*", ", ", qstart, ":", qend, " bp <br/>", "*P. orientalis*", ", ", sstart, ":", send, " bp <br/>", sep="")) +
    theme_void() +
    theme(legend.position = "none", plot.title = element_markdown()) + 
    scale_x_continuous(expand = c(0, 0), limits = c(1, length2)) 
    
    p = p + geom_area(data = align_cov, aes(x=alg_pos, y=cov_p1), fill= "#E88056")
    p = p + geom_area(data = align_cov, aes(x=alg_pos, y=cov_p2), fill = "#5BA864")

    p= p + geom_rect(xmin=0,xmax=length2,ymin=-0.25,ymax=0.25, fill = "white") +
        geom_hline(yintercept = 0.25, color = "darkgray") +
        geom_hline(yintercept = -0.25, color = "darkgray") 

    p= p + geom_segment(data = repeats %>% filter(parent=="p1"), aes(x=start,xend=end,y=0.05,yend=0.05), color = "gray", size=2.5, lineend = "butt")
  
    p= p + geom_segment(data = repeats %>% filter(parent=="p2"), aes(x=start,xend=end,y=-0.05,yend=-0.05), color = "gray", size=2.5, lineend = "butt")
    
    p = p + geom_line(aes(group=pair, color=direction)) + 
    scale_color_manual(values=c("#E88056", "#5BA864")) +
    geom_point(size = 2, colour = "gray40")
    
    if(dim(SNP)[1] > 0) {
    p = p + geom_segment(data = SNP, aes(x=alg_pos,xend=alg_pos, y=0.028, yend=-0.028), color = "firebrick1", size=0.1)
    }
  
    p = p + geom_segment(data = breakpoints, aes(x=start,xend=end, y=0, yend=0), color = "gray25", size=0.5, lineend="butt") +
    geom_segment(data = breakpoints, aes(x=start,xend=start, y=0.02, yend=-0.02), color = "gray25", size=0.5) +
    geom_segment(data = breakpoints, aes(x=end,xend=end, y=0.02, yend=-0.02), color = "gray25", size=0.5)
    
    if(dim(breakpoints %>% filter(direction == "p1 > p2" ))[1] > 0) {
    p = p + geom_text(breakpoints %>% filter(direction == "p1 > p2"), mapping = aes(x=end+50, y=0.08, label=reads, vjust=1), show.legend = FALSE)
    }
    
    if(dim(breakpoints %>% filter(direction == "p1 < p2" ))[1] > 0) {
    p = p + geom_text(breakpoints %>% filter(direction == "p1 < p2"), mapping = aes(x=end+50, y=-0.05, label=reads, vjust=1), show.legend = FALSE)
    }
    
    if(dim(gene_cov_p1)[1] > 0){
        p = p + geom_segment(data = gene_cov_p1, aes(x=min,xend=max, y=0.25, yend=0.25), color = "magenta4", size=2.5, lineend = "butt") + 
          geom_text(gene_cov_p1, mapping = aes(x=max-(max-min)/2, y=0.20, label=p1_genes, vjust=1), show.legend = FALSE,  fontface = "italic")
    }
    
    if(dim(gene_cov_p2)[1] > 0){
    p = p + geom_segment(data = gene_cov_p2, aes(x=min,xend=max, y=-0.25, yend=-0.25), color = "magenta4", size=2.5, lineend = "butt") + 
          geom_text(gene_cov_p2, mapping =aes(x=max-(max-min)/2, y=-0.15, label=p2_genes, vjust=1), show.legend = FALSE,  fontface = "italic")
    }


  ggsave(file=paste("./outputs/R_sites_", str_replace_all(i, "/", "."), ".pdf",sep=""), plot = last_plot(), device = NULL, path = NULL,
  scale = 1, width = 10, height = 3.5, units = "in",
  dpi = 1200, limitsize = F) 
  
  }
}