
require(seqinr)
require(Biostrings)
require(data.table)


# FLAGS:
REGIONS <- T 
START_END <- T

ALIGN_CUTOFF <- 30 # minimal alignment length of V region
MATCH_CUTOFF <- 60 # mininmal percentage of matching nucleotides

# Read FASTA file and conert to data frame
#
# Params:   FASTA file name
#           file directory
#
# Returns:  data frame of sequences
read.FASTA <- function(file.name, in.dir){
  
  fastaFile <- readDNAStringSet(paste0(in.dir, file.name))
  HEAD <- names(fastaFile)
  SEQ <- paste(fastaFile)
  return(data.table(HEAD, SEQ))
}

# Get rearrangment information: germline assignments and sequence functionality 
#
# Params: res.df = a data frame for parsed IgBLAST output
#         q.ind = index of current sequence query in res.df
#         q.data = IgBLAST output for current sequence query
#
# Returns: updated res.df
check.strand <- function( q.data ){
  
  # get rearrangment information
  rearrangment.info <- strsplit(q.data[grep('\\srearrangement summary\\s', q.data, perl=T, fixed=F)+1], '\t')[[1]] 
  
  ifelse(rearrangment.info[8] == '+', return(T), return(F)) 
  
}


# Fix additions and deletions in V region (according to germline)
#
# Params: q.seq - input query sequence
#         s.seq - germline sequence
#
# Returns: list of - 
#               q.aligned = fixed query sequence
#               ins.ind = insertion indices
#               del.ind = deletion indices
fix.sequence <- function( q.aligned, s.aligned ){
  
  # fix gaps in gemline = insertions in query sequence - remove nt in query sequence
  ins.ind <- unlist(gregexpr("(?=-)", s.aligned, perl=TRUE))
  ins.ind <- ins.ind[ins.ind < nchar(q.aligned)-20 & ins.ind > 0]
  
  # fix gaps in query = deletions in germline sequence - replace gap with matching ny in germline
  del.ind <- unlist(gregexpr("(?=-)", q.aligned, perl=TRUE))
  del.ind <- del.ind[del.ind < nchar(s.aligned)-20 & del.ind > 0]
  
  ins.ind <- rev(ins.ind) # fix insertions from the last to the first
  for( i in ins.ind){
    # remove gap in germline
    s.aligned <- paste0(substr(s.aligned, 1, i-1), substr(s.aligned,i+1, nchar(s.aligned)))
    # remove nt in query
    q.aligned <- paste0(substr(q.aligned, 1, i-1), substr(q.aligned,i+1, nchar(q.aligned)))
  }
  
  del.ind2 <- unlist(gregexpr("(?=-)", q.aligned, perl=TRUE)) # get indices again (in case nts were removed above)
  del.ind2 <- del.ind2[del.ind2 < nchar(q.aligned)-20 & del.ind2 > 0]
  for( i in del.ind2)
    q.aligned <- paste0(substr(q.aligned, 1, i-1 ), substr(s.aligned, i, i), 
                        substr(q.aligned, i+1, nchar(q.aligned) ))
  
  # fix Ns in V region (up to 20 nucleotides before end of germline)
  N.ind <- unlist(gregexpr("(?=N)", q.aligned, perl=TRUE))
  N.ind <- N.ind[N.ind < nchar(q.aligned)-20 & N.ind >0]
  for( i in N.ind)
    q.aligned <- paste0(substr(q.aligned, 1, i-1 ), substr(s.aligned, i, i), 
                        substr(q.aligned, i+1, nchar(q.aligned) ))
  
  # returns fixed sequence, ins.ind and del.ind to fix regions
  return(list(q.aligned, ins.ind, del.ind))
}


# Get alignment information and fix indels 
#
# Params: res.df = a data frame for parsed IgBLAST output
#         q.ind = index of current sequence query in res.df
#         q.data = IgBLAST output for current sequence query
#         chain = chain type (VH,VK etc.)
#         Vgerm = germline V sequences
#         Jgerm = germline J sequences
#
# Returns: updated res.df 
get.alignment <- function( res.df, q.ind, q.data, chain, Vgerm, Jgerm ){
  
  # get alignment information
  # V gene
  V.line <- grep('^V\\s', q.data, perl=T, fixed=F)
  Valign <- strsplit(q.data[V.line], '\t')[[1]]
  Vgene <- Valign[[3]]
  s.Vseq <- Vgerm[[Vgene]]
  if(is.null(s.Vseq)) # if best V is not in database 
    return()
  # Check if alignment is long enough
  Valign.len <- as.numeric(Valign[[5]]) 
  Vmatch <- as.numeric(Valign[[4]])
  if( Valign.len < ALIGN_CUTOFF | Vmatch < MATCH_CUTOFF )
    return() 
  # svae V gene
  res.df <- res.df[q.ind, V_CALL := Vgene]
  # get alignment and indices
  q.Vstart <- as.numeric(Valign[[6]]) 
  q.Vend <- as.numeric(Valign[[7]])
  s.Vend <- as.numeric(Valign[[9]])
  q.aligned <- Valign[[10]]
  s.aligned<- Valign[[11]]
  
  # J gene
  J.line <- grep('^J\\s', q.data, perl=T, fixed=F)
  if( length(J.line)==0)
    return() # no J was found - sequence too short
  
  Jalign <- strsplit(q.data[J.line], '\t')[[1]]
  Jgene <- Jalign[[3]]
  s.Jseq <- Jgerm[[Jgene]]
  if(is.null(s.Jseq)) # if best V is not in database 
    return()
  # save J gene
  res.df <- res.df[q.ind, J_CALL := Jgene]
  # get indices
  q.Jstart <- as.numeric(Jalign[[6]])
  s.Jstart <- as.numeric(Jalign[[8]])
  
  # V start index
  # q.seq <- res.df[q.ind, SEQUENCE] 
  fixed.indel <- fix.sequence(q.aligned, s.aligned) 
  
  # return fixed sequence, insertion and deletion indices
  q.seq <- res.df[q.ind, SEQUENCE]
  fix.seq <- paste0(fixed.indel[[1]],substr(q.seq, q.Vend+1, nchar(q.seq)))
  fix.seq <- gsub('-', '', fix.seq)
  
  res.df <- res.df[q.ind, SEQUENCE := fix.seq]
  
  # fix region indices if there are indels 
  diff <- length(fixed.indel[[3]])-length(fixed.indel[[2]])
  #diff <- ifelse(fixed.indel[[2]]<0, 0, x) + ifelse(fixed.indel[[3]]<0, 0, x)
  
  # fix and save indices - only insertions in V region are removed - substract diff from all indices
  q.Vend <- q.Vend+diff
  q.Jstart <- q.Jstart+diff
  
  # check germline ends
  if(nchar(s.Vseq)-s.Vend > 0)
    q.Vend <- q.Vend + nchar(s.Vseq)-s.Vend 
  if(s.Jstart > 1)
    q.Jstart <- q.Jstart - (s.Jstart-1)
  res.df <- res.df[q.ind, V_END := q.Vend]
  res.df <- res.df[q.ind, J_START := q.Jstart]
  res.df <- res.df[q.ind, VJ_DIST := q.Jstart-q.Vend-1]
  
  return(res.df) 
}

create.res.table <- function(N){
  # create data.frame for results
  res.df <- data.table(SEQUENCE_ID=character(N), # Unique sequence identifier
                       SEQUENCE='', # Input sequence
                       CHAIN='',
                       STRAND='', # Input sequence
                       V_CALL='', # V gene assignment(s)
                       J_CALL='', # J gene assignment(s)
                       V_END=0,
                       J_START=0,
                       VJ_DIST=0,
                       stringsAsFactors=F)
  return(res.df)
}

get.name <- function(res.df, q.ind, q.data){
  q.name <- sub('# Query: ', '', q.data[2]) 
  res.df <- res.df[q.ind, SEQUENCE_ID := q.name ]
  return(res.df)
}

check.chain <- function(q.data, chain){
  if (length(grep(paste0('\\s',chain,'\\s'), q.data, perl=T, fixed=F))==0)
    return(F)
  return(T)
}

parse.igblast <- function( input.sequences, igblastRes.dir, in.name, Vgerm, Jgerm, chain ){
  
  # load input FASTA file
  # input.sequences <- read.fasta(paste0(inSeq.dir,in.name), seqtype = "DNA", as.string = T, forceDNAtolower = F,
  #                               set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)
  setkey(input.sequences, HEAD)
  # load igblast output file
  igblast.out <- scan(paste0(igblastRes.dir, 'igblast.txt'), what='character', sep='\n', 
                      blank.lines.skip=F, strip.white=F)
  
  # get start of each query
  start.ind <- grep('# IGBLASTN\\s', igblast.out, perl=T, fixed=F) 
  n.queries <- length(start.ind)
  
  # create logical vector (to remove bad sequences at the end)
  good.seq <- logical(n.queries)
  
  # create data.frame for results
  res.df <- create.res.table(n.queries)
  
  # parse each query results separatly
  for(i in 1:n.queries){
    
    # get current lines in igblast output file
    if (i==n.queries)
      q.data <- igblast.out[c(start.ind[i]:length(igblast.out))]
    else
      q.data <- igblast.out[c(start.ind[i]:(start.ind[i+1]-1))]
    
    # get sequence name
    res.df <- get.name(res.df, i, q.data)
    
    # check if correct chain was found
    if(!check.chain(q.data, chain))
      next 
    
    # get sequence 
    q.seq <- input.sequences[grep(res.df[i,SEQUENCE_ID], input.sequences[,HEAD]), SEQ]
    
    
    # check stand - if '-' - reverse it
    res.df <- res.df[i, SEQUENCE := ifelse(check.strand(q.data), q.seq, toString(reverseComplement(DNAString(q.seq))))]
    
    
    # get alignment and fix sequencing errors, fix RF, 
    tmp <- get.alignment( res.df, i, q.data, chain, Vgerm, Jgerm )
    if( is.null(tmp) ) { next }# alignment is not good enough
    res.df <- tmp 
    good.seq[i] <- T
  }
  res.df <- res.df[good.seq,]
  return(res.df)  
}

run.igblast <- function(file.path, igblast.path, igblastRes.dir, organism){
  
  curr.dir <- getwd()
  # IGBLAST must be run from directory in which the database directory is located
  setwd(igblast.path)
  print("Run IgBLAST")
  n.cores <- detectCores()
  system(paste0('./bin/igblastn -germline_db_V database/', organism, '_db/', organism, '_gl_V -germline_db_J database/', organism, '_db/', organism, '_gl_J ',
                '-germline_db_D database/', organism, '_db/', organism, '_gl_D -organism ', organism, ' -domain_system imgt -query ',
                file.path, ' -num_threads ',  n.cores, ' -out ', igblastRes.dir, '/igblast.txt -auxiliary_data optional_file/', organism, 
                '_gl.aux -show_translation -outfmt "7 qseqid sseqid pident length qstart qend sstart send qseq sseq" -num_alignments_V 1 -num_alignments_D 0 -num_alignments_J=1' ))
  setwd(curr.dir)
}

# Get alignment information and fix indels in V gene
#
# Params: in.name = input FASTA file name
#         inSeq.dir = input directory
#         out.dir = IGBLAST output directory
#         igblast.path = IgBLAST program directory
#         Vgerm = germline V sequences
#         Jgerm = germline J sequences
#         chain = chain type (VH,VK etc.)
#         organism 
#
# Returns: res.df = a data frame for parsed IgBLAST output
igblast <- function(in.name, inSeq.dir, out.dir, igblast.path, Vgerm, Jgerm, chain, organism){
  
  # check if file is empty
  if( file.info(paste0(inSeq.dir,in.name))$size == 0 ){
    print(paste0(in.name, " is empty!"))
    return(F)
  }
  
  igblastRes.dir = paste0(getwd(),'/')
  
  # load input FASTA file
  input.sequences <- read.FASTA(in.name, inSeq.dir)
  
  # run IgBLAST
  run.igblast(paste0(inSeq.dir,in.name), igblast.path, igblastRes.dir, organism)
  
  # parse Igblast
  res.df <- parse.igblast(input.sequences, igblastRes.dir, in.name, Vgerm, Jgerm, chain )
  
  # remove bad sequences
  res.df <- res.df[VJ_DIST>=-40,]
  # write output 
  out.name <- paste0(out.dir, substr(in.name, 1, nchar(in.name)-6), '.csv')
  # new
  fwrite(res.df, file = out.name, col.names = T, sep=",", row.names = F)
  
  # remove temporary files
  file.remove(paste0(igblastRes.dir, 'igblast.txt'))
  return(T)
}

# Check maximun length of V gene and J gene in sequences
#
# Params: igblast.out = IGBLAST output directory
#         Vgerm = germline V sequences
#         Jgerm = germline J sequences

#
# Returns: Vlen = maximum length of V (numeric)
#          Jlen = maximum length of J (numeric)
check.VJ.len <- function(igblast.out, Vgerm, Jgerm ){
  
  # max Vgerm and Jegerm length
  max.V.len <- floor(quantile(sapply(Vgerm, function(x) nchar(x)), .01))
  max.J.len <- floor(quantile(sapply(Jgerm, function(x) nchar(x)), .01))
  
  V.dist <- numeric(max.V.len)
  J.dist <- numeric(max.J.len)
  
  files <- list.files(igblast.out, pattern = 'csv$')
  for(i in 1:length(files)){
    df <- fread(paste0(igblast.out, files[i]), stringsAsFactors = F)
    Vs <- df[,V_END]
    Vs[Vs>max.V.len] <- max.V.len
    V.dist <- V.dist + hist(Vs, breaks=0:(max.V.len), plot=F)$counts
    Js <- nchar(df[,SEQUENCE]) - df[,J_START] + 1
    Js[Js>max.J.len] <- max.J.len
    Js[Js<0] <- 0
    J.dist <- J.dist + hist(Js, breaks=0:(max.J.len), plot=F)$counts
  }
  # get 5 percentile for each
  V.dist <- cumsum(V.dist/sum(V.dist))
  Vlen <- length(V.dist[V.dist <=0.05])
  J.dist <- cumsum(J.dist/sum(J.dist))
  Jlen <- length(J.dist[J.dist <=0.05])
  
  return(c(Vlen, Jlen))
}
