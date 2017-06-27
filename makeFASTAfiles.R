
# flag - 
MIN.DIS <- -50 # minimum V-J distance allowed 

# Create fasta file from clone data.frame
#
# Params:  clone.df = data.frame of clone with [taxa, seq] columns
#          out.file = the file name to write to
#
# Returns: NULL
write.FASTA <- function(file.name, out.dir, head, seq, append=F) {
  sequences <- seq
  names(sequences) <- head
  writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, file.name, '.fasta'), append=append, width=1000)
}

# Remove non-unique sequences 
#
# Params: head = a character vector of sequence IDs 
#         seq = a character vector of sequences of current V-J-distance combination
#         VJL.str = string of V, J and distance (as numbers)
#         cp.num.del = delimiter for copy number (if exists)
#
# Returns: list - [head, seq, del ]
deleteEqual <- function(head, seq, file.name, seq.ID, VJL.str, cp.num.del){
  
  # get unique sequences
  uni.seq <- unique(seq)
  for(i in 1:length(uni.seq)){

    ind <- which(seq==uni.seq[i])
    num <- length(ind)
    
    # if raw sequences already contain copy number
    if (!is.null(cp.num.del))
      num <- sum(as.numeric(sapply(strsplit(head[ind], cp.num.del), function(x) { x[[length(x)]]})))
    
    # fix header of first occurence
    #head[ind[1]] <- paste0('Seq_', seq.ID, '_',file.name, '_', VJL.str, '_', as.character(num))
    head[ind[1]] <- paste0('Seq_', seq.ID, '_',file.name, '_', as.character(num))
    seq.ID <- seq.ID + 1
    if(length(ind) > 1){
      ind <- ind[-1]
      head <- head[-ind]
      seq <- seq[-ind]
    }
  }
  
  return(list(head, seq, seq.ID))
}

# Add V-J outgroup
#
# Params: Vlen = maximum length of V (numeric)
#         Jlen = maximum length of J (numeric)
#         Vseq = germline sequence of current V gene
#         Jseq = germline sequence of current J gene
#         dis = distance between end of V and start of J
#
# Returns: VJ outgroup sequence (char)
addVJ <- function(Vlen, Jlen, Vseq, Jseq, dis){
  
  # create V-J outgroup sequence
  # V segment
  if(nchar(Vseq) < Vlen) # if Vseq is too short - add gaps at beginning
    Vseq <- paste0(paste(rep('-',Vlen-nchar(Vseq)), collapse = ""), Vseq)
  else
    Vseq <- substr(Vseq, nchar(Vseq)-Vlen+1, nchar(Vseq))
  
  # J segment
  Jseq <- substr(Jseq, 1, Jlen)
  
  # distance region
  if(dis < 0){ # if V-J distance is negative, removes symetrically NTs from V and J segments
    if(ceiling(abs(dis/2))>Jlen | ceiling(abs(dis/2))+1>Jlen) # remove more nts from V
      # remove Jlen nts from J -> J is empty
      # remove abs(dis)-Jlen nts from V
      VJseq <- substr(Vseq,1, nchar(Vseq)-(abs(dis)-Jlen))
    else if(abs(dis) %% 2 == 0) # distance is even
      VJseq <- paste0(substr(Vseq,1, nchar(Vseq)-ceiling(abs(dis/2))), substr(Jseq, ceiling(abs(dis/2))+1, nchar(Jseq)))
    else # distance is odd
      VJseq <- paste0(substr(Vseq,1, nchar(Vseq)-ceiling(abs(dis/2))), substr(Jseq, ceiling(abs(dis/2)), nchar(Jseq)))
  }else
    VJseq <- paste0(Vseq, paste(rep('-',dis), collapse = ""), Jseq)
  
  return(VJseq)
}
# Cut all sequences to have same V and J lengths and add V-J outgroup
#
# Params: head = a character vector of sequence IDs 
#         seq = a character vector of sequences of current V-J-distance combination
#         VJ.ind = a data.frame with V end and J start indices
#         dis = an integer of the distance in current V-J-distance combination
#         Vlen = maximum length of V (numeric)
#         Jlen = maximum length of J (numeric)
#         Vseq = germline sequence of current V gene
#         Jseq = germline sequence of current J gene
#
# Returns: list - [del, head, seq]
cutSeq <- function(head, seq, VJ.ind, Vlen, Jlen){
  
  # get sequence lengths
  seq.len <- unname(sapply(seq, nchar))
  
  # lengths of V and J regions
  Vlen.seq <- VJ.ind[,1]
  Jlen.seq <- seq.len-VJ.ind[,2]
  
  # remove short sequences
  bad <- which(Vlen.seq<Vlen | Jlen.seq<Jlen)
  del <-  length(bad) # keep number of removed sequences
  if(del>0){
    head <- head[-bad] 
    seq <- seq[-bad]
    VJ.ind <- VJ.ind[-bad,]
  }
  
  # if no sequences are left for this V-J-distance combination
  if(length(head)==0) 
    return(del) 
  
  # cut start and end of sequence to have Vlen NTs in V segment and Jlen NTs in J segment
  if(length(seq)==1){
    seq <- substr(seq, VJ.ind[1]-Vlen+1, VJ.ind[2]+Jlen-1 )
  }else{
    for(i in 1:length(seq))
      seq[i] <- substr(seq[i], VJ.ind[i,1]-Vlen+1, VJ.ind[i,2]+Jlen-1 )
  }
  

  return(list(del,head,seq))
}

get.VJdis <- function (df, Vgerm, Jgerm){
  
  # convert V and J germline names to numbers
  VJdis <- matrix(0,nrow(df),3)
  VJdis[,3] <- df[,VJ_DIST]
  Vs <- names(Vgerm)
  VJdis[,1] <- sapply(df[,V_CALL], function(x) which(Vs==x))
  Js <- names(Jgerm)
  VJdis[,2] <- sapply(df[,J_CALL], function(x) which(Js==x))
  
  
  # get unique V-J-distance combinations
  uni.VJdis <- unique(VJdis)
  
  # remove extremely short sequences, resulting in distance < -50
  uni.VJdis <- uni.VJdis[uni.VJdis[,3]>=MIN.DIS,,drop=FALSE]
  
  return(list(uni.VJdis, VJdis))
}
# MAIN FUNCTION - Create a FASTA file with each V-J-distance combination and add the V-J outgroup sequence
#
# Params: in.name = input csv file name 
#         in.dir = directory of input csv file
#         out.dir = output directory for FASTA files
#         Vgerm = a list of germline V sequences
#         Jgerm = a list of germline J sequences
#         Vlen = maximum length of V (numeric)
#         Jlen = maximum length of J (numeric)
#
# Returns: NULL
makeFASTAfiles <- function(dir.name, file.name, in.dir, out.dir, Vgerm, Jgerm, Vlen, Jlen, cp.num.del = NULL){
  
  seq.ID <- 1
  # remove fasta extensions
  file.name <- substr(file.name, 1, nchar(file.name)-6)
  # load IgBLAST output file
  if(!file.exists(paste0(in.dir, file.name, '.csv')))
    return(0)
  df <- fread(paste0(in.dir, file.name, '.csv'), header=T, stringsAsFactors=F)
  
  # convert V and J germline names to numbers and get unique combinations
  tmp <- get.VJdis(df, Vgerm, Jgerm)
  uni.VJdis <- tmp[[1]]
  VJdis <- tmp[[2]]
  
  # for each combination - cut sequences, add V-J outgroup and remove identical sequences
  tot.del <- 0
  for(i in 1:nrow(uni.VJdis)){
    
    # get sequences with current combintaion
    ind <- which(VJdis[,1]==uni.VJdis[i,1] & VJdis[,2]==uni.VJdis[i,2] & VJdis[,3]==uni.VJdis[i,3])
    tmp.df <- df[ind,]
    
    # save V-J-distance as string of IF numbers
    VJL.str <- sprintf('%03d.%03d.%03d', uni.VJdis[i,1], uni.VJdis[i,2], uni.VJdis[i,3]);
    
    # get sequences and IDs
    head <- tmp.df[,SEQUENCE_ID]
    seq <- tmp.df[,SEQUENCE]
    VJ.ind <- as.matrix(tmp.df[,c('V_END','J_START')])
    
    # cut sequences to match Vlen and Jlen
    res <- cutSeq(head, seq, VJ.ind, Vlen, Jlen)
    del <- res[[1]]
    tot.del = tot.del + del
    
    # no sequences are left from this V-J-distance combination
    if(length(res)==1) { next }
    head <- res[[2]]
    seq <- res[[3]]
    
    # collapse duplicate sequences and add VJ outgroup sequence
    res <- deleteEqual(head, seq, file.name, seq.ID, VJL.str, cp.num.del)
    head <- res[[1]]
    seq <- res[[2]]
    seq.ID <- res[[3]]
      
    # write FASTA file
    if(file.exists(paste0(out.dir, dir.name, '_', VJL.str, '.fasta' ))){
      # append
      write.FASTA(paste0(dir.name, '_', VJL.str ), out.dir, head, seq, append = T)
    }else{
     # head <- c('VJ_1', head)
    #  seq <- c(addVJ(Vlen, Jlen, Vgerm[[uni.VJdis[i,1]]], Jgerm[[uni.VJdis[i,2]]], uni.VJdis[i,3]), seq)
      write.FASTA(paste0(dir.name, '_', VJL.str ), out.dir, head, seq, append = F)
    } 
  }
  return(del/nrow(df))
}
