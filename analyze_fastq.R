library(parallel)
library(ShortRead)
library(data.table)
library(Biostrings)
library(seqinr)
library(foreach)
library(doParallel)


OVERLAP_CUTOFF <- 30 # minimal alignment length of V region

write.FASTA <- function(file.name, out.dir, head, seq) {
  sequences <- seq
  names(sequences) <- head
  writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, file.name, '.fasta'), width=1000)
}

remove.low.quality <- function(dt.fastq){

  no.clusters <- detectCores()-1
  cl <- makeCluster(no.clusters)
  MIN.PERC <- 0.25 # % of low quality nucletoides 
  MIN.VAL <- 40   
  good <- parSapply(cl, dt.fastq[,QUALITY], function(x) length(which(as.numeric(charToRaw(x)) < MIN.VAL)) <  MIN.PERC*nchar(x))
  
  stopCluster(cl)
  return(dt.fastq[good,])
}


read.clean.fastq <- function(dir, file.name, trim.header=F){
  
  tmp <- readFastq(dir, pattern=file.name)
  dt <- data.table(ID=as.character(id(tmp)), SEQ=as.character(sread(tmp)), QUALITY=as.character(quality(quality(tmp))), stringsAsFactors=F)
  # remove low quality reads
  dt <- remove.low.quality(dt)
  
  if(trim.header)
    dt <- dt[, ID:=sapply(ID, function(x) substr(x, 1, nchar(x)-9))]
  
  return(dt)
  
}

merge.seq <- function(seq1, seq2, Q1, Q2){

  # !!! always reverse sequences from R1 !!!
  seq1<- toString(reverseComplement(DNAString(seq1)))
  # reverse quality sequence
  Q1 <- paste(rev(strsplit(Q1, NULL)[[1]]), collapse='')
  # local alignment
  nw.res <- pairwiseAlignment(seq1, seq2, type='local', gapOpening = 10 , gapExtension=10)
  # get alignment indices
  s1.start <- start(pattern(nw.res))
  s2.start <- start(subject(nw.res))
  s1.end <- end(pattern(nw.res))
  s2.end <- end(subject(nw.res))
  # check if overlap is long enough
  if( s1.end-s1.start < 30)
    return('')
  # get overlapping region indices
  overlap.region1 <- s1.start:s1.end
  overlap.region2 <- s2.start:s2.end
  # which position are of better quality in 1st sequence
  pos <- which(as.numeric(charToRaw(Q1))[overlap.region1]  >= as.numeric(charToRaw(Q2))[overlap.region2])
  char1 <- unlist(strsplit(seq1, split=''))[overlap.region1]
  char2 <- unlist(strsplit(seq2, split=''))[overlap.region2]
  
  # replace good nts from seq1 in seq2 
  char2[pos] <- char1[pos] 
  
  # concatenate sequences
  if(s1.start < s2.start) # S1 is more 3'
    merged.seq <-  paste0(substr(seq2, 1, overlap.region2[1]-1), paste0(char2, collapse=''), substr(seq1, overlap.region1[length(overlap.region1)]+1, nchar(seq1)), collapse='')
  else  # S1 is more 5'
    merged.seq <- paste0(substr(seq1, 1, overlap.region1[1]-1), paste0(char2, collapse=''), substr(seq2, overlap.region2[length(overlap.region2)]+1, nchar(seq2)), collapse='')
  
  return(merged.seq)
  
}

check.strand <- function(dt){
  # check first 100 sequences
  ind <- sample(1:nrow(dt), 100)
  res <- rep(F, 100)
  for(i in 1:length(ind)){
    seq1 <- dt[ind[i], SEQ]
    seq2 <- dt[ind[i], SEQ2]
    # seq1 inverted 
    seq11<- toString(reverseComplement(DNAString(seq1)))
    nw.res1 <- pairwiseAlignment(seq11, seq2, type='local', gapOpening = 10 , gapExtension=10)
    # seq2 inverted 
    seq22<- toString(reverseComplement(DNAString(seq1)))
    nw.res2 <- pairwiseAlignment(seq1, seq22, type='local', gapOpening = 10 , gapExtension=10)
    if(score(nw.res1) > score(nw.res2)){
      res[i] <- T
    }
  }
  if(sum(res) > 50) # R1 is inverted 
    return(1)
  else
    return(2)
}


main <- function(in.dir, dir.name, organism, chain, merge=F){
  
  
  # working directory
  w.dir <- getwd()

  # directory of FASTQfiles
  fastq.dir <- paste0(in.dir, 'raw_fastq/')
  
  # output directory for FASTA files
  fasta.dir <- paste0(in.dir, 'raw_fasta/')
  dir.create(fasta.dir, showWarnings = F, recursive = T)
  
  # get file(s) in directory
  fastq.files <- list.files(fastq.dir, pattern='fastq$')
  
  if( merge == F ){
    
    for (i in 1:length(fastq.files)){
      
      # read fastQ file and remove low quality reads
      df.fastq <- read.clean.fastq(fastq.dir, fastq.files[i], trim.header = F)
      
      # save sequences in FASTA format
      write.FASTA(substr(fastq.files[i], 1, nchar(fastq.files[i])-1), fasta.dir, df.fastq$ID, df.fastq$SEQ) 
       
    }
  }else{ # merge reads
    
    uni.files <- unique(sapply(fastq.files, function(x) substr(x, 1, nchar(x)-13)))
   
    for(i in 1:length(uni.files)){
      print(paste0(i, '/', length(uni.files)))
      ind <- which(attributes(regexpr(uni.files[i], fastq.files))$match.length!=-1)
      
      dt1 <- read.clean.fastq(fastq.dir, fastq.files[ind[1]], trim.header = T)
      dt2 <- read.clean.fastq(fastq.dir, fastq.files[ind[2]], trim.header = T)
      
      setkey(dt1,ID)
      setkey(dt2,ID)
      
      # merge data tables
      int <-intersect(dt1[,ID], dt2[,ID])
      dt1 <- dt1[match(int, dt1[,ID]),]
      dt2 <- dt2[match(int, dt2[,ID]),]
      dt <- dt1
      dt <- dt[, SEQ2:=dt2[,SEQ]]
      dt <- dt[, QUAL2:=dt2[,QUALITY]]
      
      # merge sequences
      no.clusters <- detectCores()-1
      cl <- makeCluster(no.clusters)
       clusterExport(cl, "merge.seq")
      clusterCall(cl, function() library(Biostrings))
      res <- check.strand(dt)
      if(res==1) # R1 is inverted
        system.time(dt[, MERGED := parRapply(cl, dt, function(x) merge.seq(x[2], x[4], x[3], x[5]))])
      else # R2 is inverted
        system.time(dt[, MERGED := parRapply(cl, dt, function(x) merge.seq(x[4], x[2], x[5], x[3]))])
      stopCluster(cl)
       # remove short sequences
      len <- sapply(dt[, MERGED], nchar)
      dt <- dt[len!=0,]
      # save as FASTA file
      write.FASTA(uni.files[i], fasta.dir, dt[,ID], dt[, MERGED])
   
    }
  }
}
