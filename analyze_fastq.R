
library(ShortRead)

library(Biostrings)
library(seqinr)

IGBLAST.PATH <- '/IgBLAST/ncbi-igblast-1.4.0/'
GERMLINE.PATH <- '/germline/'

MIN.PERC <- 0.25 # % of low quality nucletoides 
MIN.VAL <- 40    # minimal quality score
ALIGN_CUTOFF <- 30 # minimal alignment length of V region


parse.igblast <- function( df.fastq, igblastRes.dir, Vgerm, Jgerm, chain ){
  
  # load igblast output file
  igblast.out <- scan(paste0(igblastRes.dir, 'igblast.txt'), what='character', sep='\n', 
                      blank.lines.skip=F, strip.white=F)
  
  # get start of each query
  start.ind <- grep('# IGBLASTN\\s', igblast.out, perl=T, fixed=F) 
  n.queries <- length(start.ind)
  
  
  # parse each query results separatly
  for(i in 1:n.queries){
    
    # get current lines in igblast output file
    if (i==n.queries)
      q.data <- igblast.out[c(start.ind[i]:length(igblast.out))]
    else
      q.data <- igblast.out[c(start.ind[i]:(start.ind[i+1]-1))]
    
    # get sequence name
    q.name <- sub('# Query: ', '', q.data[2]) 
    # get index in original data frame
    q.ind <- which(df.fastq$ID==q.name)
    
    # check if correct chain was found
    if (length(grep(paste0('\\s',chain,'\\s'), q.data, perl=T, fixed=F))==0) { next }
    
    # check if alignment is long enough
    V.lines <- grep('^V\\s', q.data, perl=T, fixed=F)
    Valign <- strsplit(q.data[V.lines[1]], '\t')[[1]]
    if (as.integer(Valign[[5]]) < ALIGN_CUTOFF) { next }
    
    # check strand
    rearrangment.info <- strsplit(q.data[grep('\\srearrangement summary\\s', q.data, perl=T, fixed=F)+1], '\t')[[1]] 
    if (rearrangment.info[length(rearrangment.info)] == '-') {
      # reverse orginal sequence
      df.fastq[q.ind, 'SEQ'] <- toString(reverseComplement(DNAString(df.fastq[q.ind, 'SEQ'])))
      # reverse quality sequence
      df.fastq[q.ind, 'QUALITY'] <- paste(rev(strsplit(df.fastq[q.ind, 'QUALITY'], NULL)[[1]]), collapse='')
    }
    
    # get FR/CDR regions (according to IMGT)
    s.ind <- grep('# Alignment summary\\s', q.data, perl=T, fixed=F)
    e.ind <- grep('# Hit table\\s', q.data, perl=T, fixed=F)
    if( length(s.ind)==0 | length(s.ind)==0 ) { next() }
    IMGT.regions <- q.data[c((s.ind+1):(e.ind-3))]
    
    # search for FR2-IMGT
    reg.ind <- grep('FR2-IMGT\\s', IMGT.regions, perl=T, fixed=F)
    if (length(reg.ind)!=0){
      if (reg.ind > 1){ # FR2 is not the first and is thus complete
        split.res <- strsplit(IMGT.regions[[reg.ind]],'\t')[[1]]
        df.fastq[q.ind, 'FR2'] <- as.integer(split.res[[2]])
      }
    }
    
    # search for CDR2-IMGT
    reg.ind <- grep('CDR2-IMGT\\s', IMGT.regions, perl=T, fixed=F)
    if (length(reg.ind)!=0){
      split.res <- strsplit(IMGT.regions[[reg.ind]],'\t')[[1]]
      df.fastq[q.ind, 'CDR2'] <- as.integer(split.res[[2]])
    }  
  }
  return(df.fastq)
}

run.igblast <- function(igblast.path, igblastRes.dir, organism){
  
  curr.dir <- getwd()
  # IGBLAST must be run from directory in which the database directory is located
  setwd(igblast.path)
  print("Run IgBLAST")
  pt <- proc.time()
  system(paste0('./bin/igblastn -num_threads 16 -germline_db_V database/', organism, '_db/', organism, '_gl_V -germline_db_J database/', organism, '_db/', organism, '_gl_J ',
                '-germline_db_D database/', organism, '_db/', organism, '_gl_D -organism ', organism, ' -domain_system imgt -query ',
                curr.dir, '/tmp.fasta  -out ', igblastRes.dir, '/igblast.txt -auxiliary_data optional_file/', organism, '_gl.aux -show_translation -outfmt 7' ))
  pt2 <- proc.time()-pt
  print(pt2)
  setwd(curr.dir)
}

check.V.position <- function(df.fastq, organism, chain, germline.path, igblast.path){
  
  # load germline sequences
  Vgerm <- read.fasta(paste0(germline.path, organism, '/', chain, '/Vgermline_', organism, '_', chain, '.fasta'), 
                      seqtype = "DNA", as.string = T, forceDNAtolower = F, set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)  
  Jgerm <- read.fasta(paste0(germline.path, organism, '/', chain, '/Jgermline_', organism, '_', chain, '.fasta'), 
                      seqtype = "DNA", as.string = T, forceDNAtolower = F, set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)
  
  igblastRes.dir = paste0(getwd(),'/')
  
  
  df.fastq$CDR2 <- -1
  df.fastq$FR2 <- -1
  
  # divide file and run only 1000 sequences at a time
 # num.seq <- nrow(df.fastq)
#  num.iter <- ceiling(num.seq/1000)
#  ind <- 1
#  for( i in 1:num.iter){
#    ptm <- proc.time()
#    print(paste0(i, '/', num.iter))
#    if( i < num.iter )
#      tmp.fasta <- df.fastq[c(ind:(ind+1000-1)),c('ID', 'SEQ')] ### ind+1000-1
#    else
#      tmp.fasta <- df.fastq[c(ind:num.seq), c('ID', 'SEQ')]
    
   # write.fasta(as.list(tmp.fasta$SEQ),tmp.fasta$ID, 'tmp.fasta', open = "w", nbchar = 60)
    write.fasta(as.list(df.fastq$SEQ),df.fastq$ID, 'tmp.fasta', open = "w", nbchar = 60)
    
    #run IgBlast
    run.igblast(igblast.path, igblastRes.dir, organism)
    
    # parse Igblast
    df.fastq <- parse.igblast(df.fastq, igblastRes.dir, Vgerm, Jgerm, chain )
  #  if( i < num.iter )
  #    df.fastq[c(ind:(ind+1000-1)),] <- tmp[c(ind:(ind+1000-1)),] ### ind+1000-1
  #  else
  #    df.fastq[c(ind:num.seq), ] <- tmp[c(ind:num.seq), ]
  #  ind <- ind + 1000###
  
   # proc.time()-ptm
  #}
  
  return(df.fastq)
}

remove.low.quality <- function(df.fastq){
  df.fastq$GOOD <- F
  for (i in 1:nrow(df.fastq)){
    qual <- df.fastq[i,'QUALITY']
    bad <- length(which(as.numeric(charToRaw(qual)) < MIN.VAL))
    # check if there are more than MIN.PERC low quality nts
    if (bad < MIN.PERC*nchar(qual))
      df.fastq[i, 'GOOD'] <- T
  }
  return(df.fastq[df.fastq$GOOD==T,])
}

merge.sequences <- function(good1, good2){
  
  seq.df <- data.frame(HEAD=character(min(nrow(good1), nrow(good2))), SEQ=character(min(nrow(good1), nrow(good2))), stringsAsFactors=F)
  df.ind <- 1
  
  # merge sequences
  for (j in 1:nrow(good1)){
    ind <- which(good1[j, 'ID']==good2$ID)
    
    if (length(ind)!=0){
      reg <- which(good1[j,c('CDR2', 'FR2')]!=-1 & good2[ind,c('CDR2', 'FR2')]!=-1)
      if ( length(reg) ==0 )
        next
      reg <- reg[1] # take first if both are good
      seq1 <- good1[j, 'SEQ']
      seq2 <- good2[ind, 'SEQ']
      
      if (reg==1){ # CDR2
        overlap.length <- min(nchar(seq1)-good1[j,'CDR2'], nchar(seq2)-good2[ind,'CDR2'])
        overlap.region1 <- good1[j,'CDR2']:(good1[j,'CDR2']+overlap.length)
        overlap.region2 <- good2[ind,'CDR2']:(good2[ind,'CDR2']+overlap.length)
      }else{
        overlap.length <- min(nchar(seq1)-good1[j,'FR2'], nchar(seq2)-good2[ind,'FR2'])
        overlap.region1 <- good1[j,'FR2']:(good1[j,'FR2']+overlap.length)
        overlap.region2 <- good2[ind,'FR2']:(good2[ind,'FR2']+overlap.length)
      }
      
      # which position are of better quality in 1st sequence
      pos <- which(as.numeric(charToRaw(good1[j,'QUALITY']))[overlap.region1]  >= as.numeric(charToRaw(good2[ind,'QUALITY']))[overlap.region2])
      char1 <- unlist(strsplit(seq1, split=''))[overlap.region1]
      char2 <- unlist(strsplit(seq2, split=''))[overlap.region2]
      
      # replace good nts from seq1 in seq2 
      char2[pos] <- char1[pos] 
      
      # concatenate sequences
      seq.df[df.ind, 'HEAD'] <- good1[j, 'ID']
      if ((nchar(seq1)-good1[j,'CDR2']) > (nchar(seq2)-good2[ind,'CDR2'])) # R1 1 is longer
        seq.df[df.ind, 'SEQ'] <- paste0(substr(seq2, 1, overlap.region2[1]-1), paste0(char2, collapse=''), substr(seq1, overlap.region1[length(overlap.region1)]+1, nchar(seq1)), collapse='')
      else
        seq.df[df.ind, 'SEQ'] <- paste0(substr(seq1, 1, overlap.region1[1]-1), paste0(char2, collapse=''), substr(seq2, overlap.region2[length(overlap.region2)]+1, nchar(seq2)), collapse='')
      df.ind <- df.ind + 1
    }
  }
  
  # remove empty rows
  seq.df <- seq.df[-(df.ind:nrow(seq.df)),]
  
  return(seq.df)
}

merge.files <- function(csv.files, fasta.dir){
  
  uni.files <- unique(sapply(csv.files, function(x) substr(x, 1, nchar(x)-4)))
  
  summ <- data.frame(FILE.NAME=character(length(uni.files)), 
                     R1.tot=numeric(length(uni.files)), R1.perc=numeric(length(uni.files)), 
                     R2.tot=numeric(length(uni.files)), R2.perc=numeric(length(uni.files)), 
                     OVERLAP = numeric(length(uni.files)),
                     stringsAsFactors=F)
  
  for (i in 1:length(uni.files)){
    print(paste0(i, '/', length(uni.files)))
    
    summ[i, 'FILE.NAME'] <- uni.files[i]
    
    ind <- which(attributes(regexpr(uni.files[i], csv.files))$match.length!=-1)
    
    df1 <- read.csv(paste0(out.dir, csv.files[ind[1]]), header = T, stringsAsFactors=F)
    df2 <- read.csv(paste0(out.dir, csv.files[ind[2]]), header = T, stringsAsFactors=F)
    
    good1 <- df1[(df1$CDR2!=-1 | df1$FR2!=-1),]
    good2 <- df2[(df2$CDR2!=-1 | df2$FR2!=-1),]
    
    good1$ID <- sapply(good1$ID, function(x) substr(x, 1, nchar(x)-9))
    good2$ID <- sapply(good2$ID, function(x) substr(x, 1, nchar(x)-9))
    
    summ[i, 'R1.tot'] <- nrow(df1)
    summ[i, 'R2.tot'] <- nrow(df2)
    
    summ[i, 'R1.perc'] <- nrow(good1)/nrow(df1)
    summ[i, 'R2.perc'] <- nrow(good2)/nrow(df2)
    
    summ[i, 'OVERLAP'] <- length(which(good1$ID%in%good2$ID))/max(nrow(good1), nrow(good2))
    
    seq.df <- merge.sequences(good1, good2)
    
    write.fasta(as.list(seq.df$SEQ),seq.df$HEAD, paste0(merge.dir, uni.files[i], '.fasta'), open = "w", nbchar = 60)
  }
  
  return(summ)
}

main <- function(in.dir, dir.name, organism, chain, merge=F){
  
  # working directory
  w.dir <- getwd()
  germline.path <- paste0(w.dir, GERMLINE.PATH)
  igblast.path <- paste0(w.dir, IGBLAST.PATH)
  
   # directory of FASTQfiles
  fastq.dir <- paste0(in.dir, 'raw_fastq/')
  
  # output directory for FASTA files
  fasta.dir <- paste0(in.dir, 'raw_fasta/')
  dir.create(fasta.dir, showWarnings = F, recursive = T)
  
  # get file(s) in directory
  fastq.files <- list.files(fastq.dir, pattern='fastq$')
  
  if( merge == F ){
    
    for (i in 1:length(fastq.files)){
      
      # read fastQ file
      tmp <- readFastq(fastq.dir, pattern=fastq.files[i])
      df.fastq <- data.frame(ID=as.character(id(tmp)), SEQ=as.character(sread(tmp)), QUALITY=as.character(quality(quality(tmp))), stringsAsFactors=F)
      
      # remove low quality reads
      df.fastq <- remove.low.quality(df.fastq)
      
      # save sequences in FASTA format
      write.fasta(as.list(df.fastq$SEQ),df.fastq$ID, paste0(fasta.dir, substr(fastq.files[i], 1, nchar(fastq.files[i])-1), 'a'), open = "w", nbchar = 60)
      
    }
  }else{ # merge reads
    
    # output directory for igblast output files
    igblast.dir  <- paste0(in.dir, 'raw_csv/')
    dir.create(igblast.dir, showWarnings = F, recursive = T)
    
    for (i in 1:length(fastq.files)){
      print(i)
      
      # read fastQ file
      tmp <- readFastq(fastq.dir, pattern=fastq.files[i])
      df.fastq <- data.frame(ID=as.character(id(tmp)), SEQ=as.character(sread(tmp)), QUALITY=as.character(quality(quality(tmp))), stringsAsFactors=F)
      
      # remove low quality reads
      df.fastq <- remove.low.quality(df.fastq)
      #ptm <- proc.time()
      # run Igblast to detect V position
      df.fastq <- check.V.position(df.fastq,organism, chain, germline.path, igblast.path)
      #print(proc.time()-ptm)
      # save output
      out.name <- paste0(igblast.dir, substr(fastq.files[i], 1, nchar(fastq.files[i])-6), '.csv')
      write.table(df.fastq, file = out.name, col.names = T, sep=",", row.names = F)
    }
    
    # merge files
    csv.files <- list.files(igblast.dir)
    summ <- merge.files(csv.files)
    
    write.table(summ, file = paste0(in.dir, 'fastq_summary_', dir.name, '.csv'), col.names = T, sep=",", row.names = F)
  }
  
}