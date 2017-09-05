library(seqinr)
library(data.table)
library(Biostrings)
library(stringr)
# Read FASTA file and convert to data frame
#
# Params:   FASTA file name
#           file directory
#
# Returns:  data frame of sequences
read.FASTA <- function(file.name, in.dir){
  
  fastaFile <- readDNAStringSet(paste0(in.dir, file.name))
  head <- names(fastaFile)
  seq <- paste(fastaFile)
  return(data.table(head, seq))
}

#cutting the tree and removing vertex with only one sun
cutting_the_tree<-function(edge.dt, numOfChange, n.leaves, smallerthantresh, validEdge){
  #sotring the data frame according to the from vertex
  ord <- order(edge.dt[,from])
  edge.dt<- edge.dt[ord, ]
  smallerthantresh <- smallerthantresh[ord]
  validEdge <- validEdge[ord]
  #creating a new data frame that consist data about edges that have one or zero suns and if there have only one suns the 
  #number of the line of this edge.
  vertex<- c(edge.dt[,from], edge.dt[,to])
  max.num <- max(as.numeric(vertex[!grepl("Seq", vertex)]))
  lineNum<- rep.int(-2, max.num-n.leaves)#max(as.numeric(vertex), na.rm = T))
  numOfChild<- rep.int(-2, max.num-n.leaves)#max(as.numeric(vertex), na.rm = T))
  # validvertex<- data.frame(lineNum, numOfChild)
  validvertex <- matrix(-2, nrow=max.num-n.leaves, ncol=2)
  colnames(validvertex) <- c('lineNum', 'numOfChild')
  #counting the suns for the first vertex
  countSuns<-0
  #saving the index edge of the only sun
  index<-0
  #the first from vertx
  oneVertex <- edge.dt[1,from]
  
  if(smallerthantresh[1]){
    countSuns<- 1
    index<-1
  }
  j<-edge.dt[1,from]
  for (i in 2:nrow(edge.dt)){
    secondVertex <- edge.dt[i,from]
    #if im running and im still on the same vertex then i keep on counting the suns  
    if (oneVertex == secondVertex){
      countSuns <- countSuns + smallerthantresh[i]
      #updating the index - the line if i only have one sumn
      if((index==0) && (countSuns==1)){
        index<-i
      } 
    } else{
      # i reached to a new vertex
      #if i have problomatic vertex
      if(countSuns<=1){
        #updating the line number and the number of suns
        validvertex[j-n.leaves,'lineNum'] <- index
        validvertex[j-n.leaves, 'numOfChild'] <- countSuns
      } else {
        validvertex[j-n.leaves,'lineNum'] <- -1
        validvertex[j-n.leaves, 'numOfChild'] <- -1
      }
      #a new vertex
      j<- edge.dt[i, from]
      #starsing to count the suns from the begining
      countSuns<-0
      index<-0
      oneVertex <- edge.dt[i,from]
      if(smallerthantresh[i]){
        countSuns<-1
        index<-i
      }
    }
  }
  if (countSuns<=1){
    validvertex[j-n.leaves, 'lineNum'] <- index
    validvertex[j-n.leaves, 'numOfChild']<-countSuns
  } else {
    validvertex[j-n.leaves,'lineNum'] <- -1
    validvertex[j-n.leaves, 'numOfChild'] <- -1
  }
  #running over the column of the to vertex
  for (m in 1:nrow(edge.dt)){
    vertexTo <- edge.dt[m, to]
    if(!(grepl("_",vertexTo))){
      index<- as.numeric(vertexTo)
      if(validvertex[index-n.leaves, 'numOfChild']==0){
        validEdge[m] <- F
        numOfChange <- numOfChange + 1
      }
      else if(validvertex[index-n.leaves, 'numOfChild']==1){
        index2<- as.numeric(validvertex[index-n.leaves, 'lineNum'])
        edge.dt <- edge.dt[m, to := edge.dt[index2, to]]
        edge.dt <- edge.dt[m, distance := edge.dt[index2, distance] + edge.dt[m, distance]]
        validEdge[index2] <- F
        numOfChange <- numOfChange + 1
      } else if(validvertex[index-n.leaves, 'numOfChild']==-2){
        validEdge[m] <- F
        numOfChange <- numOfChange + 1
      }
    }
  }
  
  edge.dt <- edge.dt[validEdge,]
  
  return(list("numOfChange" = numOfChange, "edge.dt" = edge.dt))
  
}
#this is the main function - you need to call only it
cutting <- function(edge.dt, tree.threshold){
  n.leaves <-min(edge.dt[,from])-1
  smallerthantresh <- edge.dt[,distance]<=tree.threshold
  validEdge <- edge.dt[,distance]<=tree.threshold
  numOfChange<-100
  while (numOfChange!=0){
    numOfChange<-0
    list_objects<-cutting_the_tree(edge.dt, numOfChange, n.leaves, smallerthantresh, validEdge)
    numOfChange<- list_objects$numOfChange
    edge.dt<- list_objects$edge.dt
    smallerthantresh <- edge.dt[,distance]<=tree.threshold
    validEdge <- edge.dt[,distance]<=tree.threshold
    if(nrow(edge.dt)==0){
      numOfChange<-0
    } else {
      for (t in 1:nrow(edge.dt)){
        if(edge.dt[t, distance] > tree.threshold){
          validEdge[t] <- F
          numOfChange<- numOfChange + 1
        }
      }
    }
    if (nrow(edge.dt) == 1){
      numOfChange <- 0
    }
  }
  return(edge.dt)
}

cut.clones <- function(nameSeq.dt, edge.dt, tree.threshold, file.name, curr.tree){
  
  nameSeq.dt <- nameSeq.dt[,head2 := sapply(strsplit(nameSeq.dt[,head], paste0(file.name, '_', curr.tree)), function(x) substr(x[2],2, nchar(x[2])))]
  setkey(nameSeq.dt, head2)
  # add column for clone ID
  nameSeq.dt <- nameSeq.dt[, cloneID := -1 ]
  
  # Cut edges longer than threshold
  trim.edge.dt <- cutting(edge.dt, tree.threshold)
  #trim.edge.dt <- edge.dt[edge.dt[,distance]<=tree.threshold,]
  clone.id <- 1
  if(nrow(trim.edge.dt)>0){
    
    # Get roots of all subtrees
    sub.roots <- setdiff(trim.edge.dt[,from], trim.edge.dt[,to])
    
    # For each root, get subtree edges and nodes
    for(root in sub.roots) {
      at.bottom <- F
      nodes <- nameSeq.dt[head2==root]
      # Trace down tree until at all leaves
      while(at.bottom != T) {
        # Add children of all nodes in tree thus far
        new.nodes <- unique(rbind(nodes, nameSeq.dt[nameSeq.dt[,head2] %in% trim.edge.dt[trim.edge.dt[,from] %in% nodes[,head2], to],]))
        
        # If no children are to be added (tree is complete)
        if(nrow(nodes) == nrow(new.nodes)) {
          nodes <- new.nodes
          # delete internal sequences
          nodes <- nodes[grep('Seq_', nodes[,head2])] 
          if(nrow(nodes)!=0){
            nameSeq.dt[nodes[,head2], cloneID := as.numeric(clone.id)]
            clone.id <- clone.id + 1
          }
          at.bottom <- T
        }else{
          nodes <- new.nodes
        }
      }
    }
    
    # save sub-roots, sequences and clone IDs
    if( length(sub.roots)!=0){
      roots <- data.table(root=sub.roots, seq='', clone.ID=1:length(sub.roots))
      roots <- roots[,seq := sapply(sub.roots, function(x) nameSeq.dt[nameSeq.dt[head2==x],'seq'])]
    }else
      roots <- -1
  }else
    roots <- -1
  # get singletons
  sing <- setdiff(edge.dt[,to], trim.edge.dt[,to])
  
  # remove VJ and internal nodes
  sing <- sing[grep('Seq_', sing)]
  
  # remove internal sequences from nameSeq.dt
  nameSeq.dt <- nameSeq.dt[grep('Seq_', nameSeq.dt[,head2]),] 
  
  # add singletons as clones
  if(length(sing)>0){
    nameSeq.dt[sing, cloneID := as.numeric(clone.id:(clone.id+length(sing)-1))]
  }
  
  # remove internal sequences from nameSeq.dt
  nameSeq.dt <- nameSeq.dt[grep('Seq_', nameSeq.dt[,head2]),] 
  
  # order by cloneID
  nameSeq.dt <- nameSeq.dt[order(nameSeq.dt[,cloneID]),]
  
  return(list(nameSeq.dt, roots))
}

create.out.table <- function(N){
  
  # N = 10000
  out.dt <- data.table(FILE_ID=character(N), # original file name( before merging)
                       SEQUENCE_ID='', # Unique sequence identifier
                       CLONE_ID='', # Unique clone identifier
                       SEQUENCE='', # Input sequence
                       AA_SEQUENCE='', # Input sequence
                       FUNCTIONAL=F, # T: productive, F: unproductive 
                       STRAND='', # Input sequence
                       INFRAME=F, # T: junction is in-frame, F: junction is out-of-frame
                       STOP=F, # T: stop codon and unproductive sequence, F: no stop codon
                       CP_NUM=-1,# 
                       V_CALL='', # V gene assignment(s)
                       D_CALL='', # D gene assignment(s)
                       J_CALL='', # J gene assignment(s)
                       FR1=-1,# 
                       CDR1=-1,# 
                       FR2=-1,# 
                       CDR2=-1,# 
                       FR3=-1,# 
                       CDR3=-1,# 
                       V_END=-1,# 
                       D_START=-1,# 
                       D_END=-1,# 
                       J_START=-1,# 
                       VJ_DIST=-1,
                       stringsAsFactors=F)# 
  return(out.dt)
}


get.VJdis <- function(out.dt, tree.name,Vgerm, Jgerm){
  #name.split <- unlist(strsplit(tree.name, "_", fixed = TRUE))
  #class <- name.split[length(name.split)]
  VJdis.num <- as.numeric(unlist(strsplit(tree.name, ".", fixed = TRUE)))
  out.dt <- out.dt[,V_CALL := names(Vgerm)[VJdis.num[1]]]
  out.dt <- out.dt[,J_CALL := names(Jgerm)[VJdis.num[2]]]
  out.dt <- out.dt[,VJ_DIST := VJdis.num[3]]
  
  return(out.dt)
}


get.RF <- function(Vlen, Vgerm){
  
  # germline V sequence lengths
  germ.len <- nchar(Vgerm)
  
  # trim beginning of germline V sequences to match Vlen
  trim.germ <- lapply(Vgerm, function(x) substr(x, nchar(x)-Vlen+1, nchar(x)))
  
  # modulu of full germline V end
  rem <- germ.len%%3
  
  # modulu of trimmed germline V end (without extra NTs at the end of the full sequences)
  rem2 <- nchar(mapply(function(x,y) substr(x, 1, nchar(x)-y), x=trim.germ, y=rem))%%3
  
  # reading frame of trimmed germline V sequences
  V.RF <- rem2+1;
  
  # if there are germline sequences shorter than Vlen - set RF to be NA
  V.RF[germ.len<Vlen] <- NA
  
  return(V.RF)
}


get.regions <- function(regions, Vlen, V.RF){
  
  # for too short sequences - regions will be automatically NAs
  # calculate variable region length according to Vlen
  diff <- regions[,V_END]-Vlen+V.RF-1
  
  # fix region indices
  reg.names <- colnames(regions)[-1]
  regions[, (reg.names) := regions[, reg.names, with=F]-diff ]
  
  # for Vlen which does not include all regions, set first indices to -1
  # set first index to 1
  for(i in 1:nrow(regions)){
    neg.ind <- which(regions[i,reg.names, with=F]<0)
    if(length(neg.ind)==1)
      regions <- regions[i, reg.names[neg.ind] := 1]
    else if(length(neg.ind)>1){
      regions <- regions[i, reg.names[1:(length(neg.ind)-1)] := -1]
      regions <- regions[i, reg.names[length(neg.ind)] := 1]
    }    
  }
  
  return(regions)
}

fix.seq.RF <- function(out.dt, nameSeq.dt, V.RF){
  
  RF <- V.RF[out.dt[1,V_CALL]]
  out.dt[,SEQUENCE := sapply(nameSeq.dt[,seq], function(x) substr(x, RF, nchar(x)))]
  out.dt[,SEQUENCE_ID := nameSeq.dt[,head] ]
  return(out.dt)
}

set.regions <-  function(out.dt, regions, Jlen){
  
  # current V gene
  V.name <- out.dt[1,'V_CALL']
  
  # set CDR/FR and V end indices
  out.dt[,c('FR1','CDR1','FR2','CDR2','FR3','CDR3','V_END')] <- regions[regions$Vgene==V.name,-1]
  
  # set J start index
  out.dt$J_START <- nchar(out.dt[1,'SEQUENCE'])-Jlen+1
  
  return(out.dt)
}


translate.seq <- function(out.dt){
  
  out.dt <- out.dt[, AA_SEQUENCE := sapply(out.dt[, SEQUENCE], function(x) paste(seqinr::translate(s2c(substr(x, 1, nchar(x)-(nchar(x)%%3))),
                                                                                                   numcode = 1, NAstring = "X", ambiguous = F), collapse=''))]
  return(out.dt)  
}

check.functionality <- function( out.dt, Jlen, Jgerm ){
  
  # check presence of stop codons in each sequence
  SC <- which(sapply(out.dt[,AA_SEQUENCE], function(x) length(grep("\\*", x)))==1)
  out.dt[SC,'STOP'] <- T
  
  # all sequence in current tree are aligned, so only one sequence needs to be tested
  # check if J is in same RF as V
  Jseq <- Jgerm[[out.dt[1,J_CALL]]]
  Jstart <- out.dt[1,J_START]
  if( nchar(paste0(substr(out.dt[1,SEQUENCE],1,Jstart-1), substr(Jseq,1,nchar(Jseq)-1)))%%3 == 0 )
    out.dt <- out.dt[,INFRAME := T]
  
  # sequence is functional if is in frame and has no stop codons
  func <- out.dt[, INFRAME] == T & out.dt[,STOP] == F
  out.dt <- out.dt[func,FUNCTIONAL := T ]
  
  return(out.dt)
}

set.regions <-  function(out.dt, regions, Jlen){
  
  # current V gene
  V.name <- out.dt[1,V_CALL]
  reg.names <- colnames(regions)[-1]
  # set CDR/FR and V end indices
  out.dt <- out.dt[, (reg.names) := regions[Vgene==V.name, reg.names, with=F] ]
  
  # set J start index
  out.dt <- out.dt[,J_START := nchar(out.dt[1,SEQUENCE])-Jlen+1]
  
  return(out.dt)
}

get.cp.num <- function(out.dt){
  
  # get copy number form sequence headers
  cp.num <- sapply(strsplit(out.dt[,SEQUENCE_ID], '_'), function(x) x[[length(x)]])
  
  out.dt <- out.dt[,CP_NUM := as.numeric(cp.num)]
  
  return(out.dt)
}


change.names <- function(out.dt, names.df){
  
  # match between original and temporary names 
  out.dt$SEQUENCE_ID <- as.character(sapply(out.dt$SEQUENCE_ID, function(x) names.df[x==names.df$head, 'head2']))
  
  return(out.dt)
}

set.clone.IDs <- function(out.dt, nameSeq.dt, tree.name){
  
  # paste tree name and clone ID
  out.dt <- out.dt[,CLONE_ID := paste0(tree.name, '_clone', nameSeq.dt[, cloneID])]
  
  return(out.dt)
}

fix.sing.format <- function(nameSeq.dt){
  
  # add clone number
  nameSeq.dt <- nameSeq.dt[,cloneID := 1]
  
  return(nameSeq.dt)
}
fix.pair.format <- function(nameSeq.dt, tree.threshold){
  # check distance between the sequences
  seq.1.char <- unlist(strsplit(nameSeq.dt[1,seq], ''))
  seq.2.char <- unlist(strsplit(nameSeq.dt[2,seq], ''))
  # get mutation position
  distance <- length(which(seq.1.char!= seq.2.char))
  # add clone number
  if(distance <= tree.threshold){ # same clone
    nameSeq.dt <- nameSeq.dt[,cloneID := 1]
  }else{ # different clones
    nameSeq.dt <- nameSeq.dt[,cloneID := as.list(1:2)]
  }
  
  return(nameSeq.dt)
}
get.consensus <- function(seq){
  con.seq <- consensusString(DNAStringSet(seq), ambiguityMap='N')
  amb.pos <- gregexpr(con.seq, pattern='N')
  if(sum(amb.pos!=-1)>1){
    split.com <- strsplit(con.seq, '')
    split.seq <- strsplit(seq, '')
    for(i in 1:length(amb.pos)){
      nts <- sapply(split.seq, function(x) x[amb.pos[i]])
      split.con[[1]][amb.pos[i]] <- nts[[samples(1:length(nts),1)]]
    }
    con.seq <- c2s(split.con[[1]])
  }
  return(con.seq)
}

merge.clones <- function( out.dt, nameSeq.dt, roots){
  setkey(nameSeq.dt, cloneID)
  setkey(out.dt, SEQUENCE_ID)
  n.clones <- max(nameSeq.dt[,cloneID])
  out.dt.clones <- create.out.table( n.clones )
  for(i in 1:n.clones){
    # get all current clone sequences
    tmp.dt <- out.dt[ nameSeq.dt[.(i), head],]
    out.dt.clones <- out.dt.clones[i, names(tmp.dt) := tmp.dt[1,]]
    
    if(nrow(tmp.dt)>1){ # not a singleton - get the ancestral sequence
      if(is.numeric(roots)) # in case of a tree with 2 sequences
        out.dt.clones[i, SEQUENCE := get.consensus(nameSeq.dt[,seq])]
      else
        out.dt.clones[i, SEQUENCE := roots[i, seq]]
    }
    # update copy number (clone size)
    out.dt.clones <- out.dt.clones[i, CP_NUM := sum(tmp.dt[,CP_NUM])] 
    
  }
  return(out.dt.clones)
}

split.clones <- function(out.dt.clones, out.dt, out.seq.dir, out.clone.dir, dir.name, V){
  # for each clone - split sequences into files, and update clone size
  setkey(out.dt, CLONE_ID)
  for(i in 1:nrow(out.dt.clones)){
    tmp.dt <- out.dt[out.dt.clones[i, CLONE_ID],]
    setkey(tmp.dt, FILE_ID)
    files <- unique(tmp.dt[,FILE_ID])
    for( j in 1:length(files)){
      tmp.dt.file <- tmp.dt[files[j],]
      
      tmp.dt.clones <- out.dt.clones[i,]
      
      # update copy number (clone size)
      tmp.dt.clones <- tmp.dt.clones[,CP_NUM := sum(tmp.dt.file[,CP_NUM])]
      
      # write files
      if(!file.exists(paste0(out.seq.dir, dir.name, '_', files[j], '_', V, '_out.csv'))){
        fwrite(tmp.dt.file, file = paste0(out.seq.dir, dir.name, '_', files[j], '_', V, '_out.csv'), sep=",", append = F, row.names = F, col.names = T)
        fwrite(tmp.dt.clones, file = paste0(out.clone.dir, dir.name, '_', files[j], '_', V, '_out.csv'), sep=",", append = F, row.names = F, col.names = T)   
      }else{
        fwrite(tmp.dt.file, file = paste0(out.seq.dir, dir.name, '_', files[j], '_', V, '_out.csv'), sep=",", append = T, row.names = F, col.names = F)
        fwrite(tmp.dt.clones, file = paste0(out.clone.dir, dir.name, '_', files[j], '_', V, '_out.csv'), sep=",", append = T, row.names = F, col.names = F)   
      } 
    }
  }
}

get.file.ID <- function(out.dt, dir.name, curr.tree){
  
  tmp <- strsplit(sapply(strsplit(out.dt[,SEQUENCE_ID], 'Seq_'), function(x) x[2]),'_')
  out.dt <- out.dt[,FILE_ID := sapply(tmp, function(x) paste0(x[2:(length(x)-1)], collapse = '_'))]
  if(length(which(is.na(out.dt[,FILE_ID])))){
    print('ss')
  }
  return(out.dt)
}

get.clone.info <- function(dir.name, out.seq.dir, out.clones.dir){
  
  all.files <- list.files(out.seq.dir)
  n.sub.files <- length(all.files)-1
  # get index of merged file
  m.ind <- which(all.files==paste0(dir.name, '_out.csv'))
  clone.info <- data.frame(file=c(all.files[m.ind], all.files[-m.ind]), tot.seq=0, unique.seq=0, unique.clone=0, stringsAsFactors = F)
  
  # add info for merged file
  dt.clone <-  fread(paste0(out.clones.dir, all.files[m.ind]))
  dt.seq <-  fread(paste0(out.seq.dir, all.files[m.ind]))
  clone.info[1,'tot.seq'] <- sum(dt.clone[,CP_NUM])
  clone.info[1,'unique.seq'] <- nrow(dt.seq)
  clone.info[1,'unique.clone'] <- nrow(dt.clone)
  
  # add info for each file
  ind <- 2
  for(file in all.files[-m.ind]){
    df.clone <-  read.csv(paste0(out.clones.dir, file))
    df.seq <-  read.csv(paste0(out.seq.dir, file))
    clone.info[ind,'tot.seq'] <- sum(df.clone$CP_NUM)
    clone.info[ind,'unique.seq'] <- nrow(df.seq)
    clone.info[ind,'unique.clone'] <- nrow(df.clone)
    ind <- ind + 1
  }
  return(clone.info)
}

compute.mutations <- function(dir.name, curr.tree, nameSeq.dt, edge.dt, regions){
  setkey(nameSeq.dt, head)
  # get V index
  V.ind <- as.numeric(substr(curr.tree,1,3))
  V.regions <- as.numeric(regions[V.ind, 2:8])
  V.regions.names <- colnames(regions)[2:8][V.regions!=-1] # remove 'V_END'
  V.regions.names <- V.regions.names[-length(V.regions.names)]
  # remove regions not in range of sequences
  V.regions <- V.regions[V.regions!=-1]
  V.regions.aa <- ceiling(V.regions/3)
  # mutations.df <- data.frame(NS.FR1=nrow(edge.dt), NS.CDR1=0, NS.FR2=0, NS.CDR2=0, NS.FR3=0, NS.CDR3=0, S.FR1=0, S.CDR1=0, S.FR2=0, S.CDR2=0, S.FR3=0, S.CDR3=0)
  for(i in 1:nrow(edge.dt)){
    print(i)
    # get sequences
    seq.father <- nameSeq.dt[paste0(dir.name, '_', curr.tree, '_', edge.dt[i, from]),seq]
    seq.son <- nameSeq.dt[paste0(dir.name, '_', curr.tree, '_', edge.dt[i, to]),seq]
    # convert strings to chars
    seq.father.char <- unlist(strsplit(seq.father, ''))
    seq.son.char <- unlist(strsplit(seq.son, ''))
    # get mutation position
    mut.pos.nt <- which(seq.father.char!= seq.son.char)
    edge.dt[i,distance := length(mut.pos.nt)]
    # translate sequences
    seq.father.aa <- seqinr::translate(tolower(seq.father.char))
    seq.son.aa <- seqinr::translate(tolower(seq.son.char))
    # find mutation positions at amino acid level
    mut.pos.aa <- ceiling(mut.pos.nt/3)
    # if mutations are found in the end of the sequence which does not have a full codon
    if(sum(mut.pos.aa>length(seq.son.aa))>0)
      mut.pos.aa <- mut.pos.aa[mut.pos.aa<=length(seq.son.aa)]
    NS.mut <- mut.pos.aa[seq.father.aa[mut.pos.aa]!=seq.son.aa[mut.pos.aa]]
    S.mut <- mut.pos.aa[seq.father.aa[mut.pos.aa]==seq.son.aa[mut.pos.aa]]
    
    # divided by regions
    if(length(NS.mut)>0){
      a <- hist(NS.mut,  breaks=c(1, V.regions.aa[2:(length(V.regions.aa)-1)], length(seq.father.aa)), plot=F)
      edge.dt[i,  (paste0('NS.', V.regions.names)) := as.list(a$counts)]
    }
    if(length(S.mut)>0){
      a <- hist(S.mut,  breaks=c(1, V.regions.aa[2:(length(V.regions.aa)-1)], length(seq.father.aa)), plot=F)
      edge.dt[i,  (paste0('S.', V.regions.names)) := as.list(a$counts)]
    }
  }
  return(edge.dt)
}


compute.edge <- function(nameSeq.dt, edge.dt){
  
  n.edge <- nrow(edge.dt)
  edge.dt[,distance := 0]
  setkey(nameSeq.dt, head)
  for(i in 1:n.edge){
    from.seq <- unlist(strsplit(nameSeq.dt[edge.dt[i,from], seq],''))
    to.seq <- unlist(strsplit(nameSeq.dt[edge.dt[i,to], seq],''))
    edge.dt[i,distance := length(which(from.seq!=to.seq))]
  }
  
  return(edge.dt)
}

convert.format <- function(V, dir.name, in.dir, out.seq.dir, out.clones.dir, Vlen, Jlen, Vgerm, Jgerm, regions, tree.threshold){
  
  out.seq.name <- paste0(out.seq.dir, dir.name, '_', V, '_out.csv')
  out.clones.name <- paste0(out.clones.dir, dir.name, '_', V, '_out.csv')
  
  V.RF <- get.RF(Vlen, Vgerm)
  regions <- get.regions(regions, Vlen, V.RF)
  
  if(!file.exists(paste0(in.dir, dir.name, '_sequences_', V, '.fasta'))){
    print('no trees!')
    return()
  }
  
  TREES <- F
  all.nameSeq.dt <- read.FASTA(paste0(dir.name, '_sequences_', V, '.fasta'), in.dir)
  if(file.exists(paste0(in.dir, dir.name, '_edges_', V, '.tab'))){
    TREES <- T
    all.edge.dt <- fread( paste0(in.dir, dir.name, '_edges_', V, '.tab'), sep="\t", header=T)
    all.edge.dt[,c('NS.FR1', 'NS.CDR1', 'NS.FR2', 'NS.CDR2', 'NS.FR3', 'NS.CDR3', 'S.FR1', 'S.CDR1', 'S.FR2', 'S.CDR2', 'S.FR3', 'S.CDR3') := 0] 
  }
  FIRST <- T
  start.ind <- 1
  end.ind <- 2
  edge.start <- 1
  while(start.ind < nrow(all.nameSeq.dt)){
    STOP <- F
    tmp<-gregexpr(pattern='_', all.nameSeq.dt[start.ind, head])[[1]]
    curr.tree.long <- substr(all.nameSeq.dt[start.ind, head], 1, tmp[length(tmp)]-1)
    while(STOP == F){
      if(length(grep(curr.tree.long, all.nameSeq.dt[end.ind, head]))==0 | end.ind > nrow(all.nameSeq.dt)) # new tree
        STOP <- T
      end.ind <- end.ind + 1
    }
    nameSeq.dt <- all.nameSeq.dt[start.ind:(end.ind-2),]
    curr.tree <- strsplit(strsplit(curr.tree.long, dir.name)[[1]][2], '_')[[1]][2]
    print(curr.tree)
    
    start.ind <- end.ind-1
    if(nrow(nameSeq.dt) == 1){ # singleton
      nameSeq.dt <- fix.sing.format(nameSeq.dt)
    }else if (nrow(nameSeq.dt) == 2){
      nameSeq.dt <- fix.pair.format(nameSeq.dt, tree.threshold)
      roots <- -1
    }else{
      edge.dt <- all.edge.dt[edge.start:(edge.start+nrow(nameSeq.dt)-2), ]
      
      edge.dt <- compute.mutations(dir.name, curr.tree, nameSeq.dt, edge.dt, regions)
      all.edge.dt <- all.edge.dt[edge.start:(edge.start+nrow(nameSeq.dt)-2), names(all.edge.dt):= edge.dt]
      edge.start <- edge.start+nrow(nameSeq.dt)-1
      
      tmp <- cut.clones(nameSeq.dt,edge.dt,tree.threshold, dir.name, curr.tree) 
      nameSeq.dt <- tmp[[1]]
      roots <- tmp[[2]]
    }
    
    # create temporary output table for current tree
    out.dt <- create.out.table(nrow(nameSeq.dt))   
    # get V/J/distance
    out.dt <- get.VJdis(out.dt, curr.tree, Vgerm, Jgerm)
    # cut start of sequences to fix reading frame
    out.dt <- fix.seq.RF(out.dt, nameSeq.dt, V.RF)
    # get copy number
    out.dt <- get.cp.num(out.dt)
    # set clone ID's in output table
    out.dt <- set.clone.IDs(out.dt, nameSeq.dt, curr.tree)
    # get FR/CDR regions
    out.dt <- set.regions(out.dt, regions, Jlen)
    # translate NT sequences
    out.dt <- translate.seq(out.dt)
    # check functionality
    out.dt <- check.functionality( out.dt, Jlen, Jgerm )
    # get file name for each sequence
    out.dt <- get.file.ID(out.dt, dir.name, curr.tree)  
    # merge clone sequences
    if(nrow(out.dt)>1){
      out.clones.dt <- merge.clones( out.dt, nameSeq.dt, roots)
    }else{
      out.clones.dt <- out.dt
    }
    split.clones(out.clones.dt, out.dt, out.seq.dir, out.clones.dir, dir.name, V)
    
    # save clones in output table
    if(!file.exists(out.seq.name)){
      fwrite(out.dt, file = out.seq.name, sep=",", append = F, row.names = F, col.names = T)
      fwrite(out.clones.dt, file = out.clones.name, sep=",", append = F, row.names = F, col.names = T)
    }else{
      fwrite(out.dt, file = out.seq.name, sep=",", append = T, row.names = F, col.names = F)
      fwrite(out.clones.dt, file = out.clones.name, sep=",", append = T, row.names = F, col.names = F)
    }
  }  
  if(TREES)
    fwrite(all.edge.dt, file = paste0(in.dir, dir.name, '_edge_mutations_', V, '.tab'), sep="\t", row.names = F, col.names = T)
  
  return()
}

split.tree.files <- function(dir.name, tree.out){
  
  all.nameSeq.dt <- read.FASTA(paste0(dir.name, '_sequences.fasta'), tree.out)
  all.edge.dt <- fread( paste0(tree.out, dir.name, '_edges.tab'), sep="\t", header=T)
  
  # split by Vs
  tmp.seq <- str_extract(all.nameSeq.dt[,head], pattern='[0-9]{3}\\.[0-9]{3}')
  tmp.edge <- str_extract(all.edge.dt[,Tree], pattern='[0-9]{3}\\.[0-9]{3}')
  Vs.seq <- str_sub(tmp.seq, 1, 3)
  Vs.edge <- str_sub(tmp.edge, 1, 3)
  uni.Vs <- unique(Vs.seq)
  # save files
  for(i in 1:length(uni.Vs)){
    tmp.nameSeq.df <- all.nameSeq.dt[Vs.seq==uni.Vs[i],]
    write.FASTA(paste0(dir.name, '_sequences_', uni.Vs[i]), tree.out, tmp.nameSeq.df, append=F)
    tmp.edge.dt <- all.edge.dt[Vs.edge==uni.Vs[i],]
    fwrite(tmp.edge.dt, file=paste0(tree.out, dir.name, '_edges_', uni.Vs[i], '.tab'), quote=F, sep='\t', col.names=T, row.names=F) 
  }
  return(uni.Vs)
}

merge.out.files <- function(dir.name, seq.out, clones.out, seq.out.merged, clones.out.merged){
  
  all.files <- list.files(clones.out)
  big.files.ind <-grep(all.files, pattern=paste0(dir.name, '_[0-9]{3}_out\\.csv'))
  big.files <- all.files[big.files.ind]
  all.files <- all.files[-big.files.ind]
  # merge big files
  out.seq.name <- paste0(seq.out.merged, dir.name, '_out.csv')
  out.clones.name <- paste0(clones.out.merged, dir.name, '_out.csv')
  for(i in 1:length(big.files)){
    out.dt <- fread(paste0(seq.out, big.files[i]))
    out.clones.dt <- fread(paste0(clones.out, big.files[i]))
    if(!file.exists(out.seq.name)){
      fwrite(out.dt, file = out.seq.name, sep=",", append = F, row.names = F, col.names = T)
      fwrite(out.clones.dt, file = out.clones.name, sep=",", append = F, row.names = F, col.names = T)
    }else{
      fwrite(out.dt, file = out.seq.name, sep=",", append = T, row.names = F, col.names = F)
      fwrite(out.clones.dt, file = out.clones.name, sep=",", append = T, row.names = F, col.names = F)
    }
  }
  # merge all files
  files <- sapply(str_split(all.files, pattern='_[0-9]{3}_out\\.csv'), function(x) x[[1]])
  uni.files <- unique(files)
  for(i in 1:length(uni.files)){
    curr.files <- all.files[grep(files, pattern=uni.files[i])]
    out.seq.name <- paste0(seq.out.merged, uni.files[i], '_out.csv')
    out.clones.name <- paste0(clones.out.merged, uni.files[i], '_out.csv')
    for(j in 1:length(curr.files)){
      out.dt <- fread(paste0(seq.out, curr.files[j]))
      out.clones.dt <- fread(paste0(clones.out, curr.files[j]))
      if(!file.exists(out.seq.name)){
        fwrite(out.dt, file = out.seq.name, sep=",", append = F, row.names = F, col.names = T)
        fwrite(out.clones.dt, file = out.clones.name, sep=",", append = F, row.names = F, col.names = T)
      }else{
        fwrite(out.dt, file = out.seq.name, sep=",", append = T, row.names = F, col.names = F)
        fwrite(out.clones.dt, file = out.clones.name, sep=",", append = T, row.names = F, col.names = F)
      }
    }
  }
}

merge.tree.mutation.files <- function(dir.name, tree.out, tree.out.merged){
  
  mut.files <- list.files(tree.out, pattern='mutations')
  tree.files <- list.files(tree.out, pattern='edges')
  seq.files <- list.files(tree.out, pattern='sequences')
  
  mut.name <- paste0(tree.out.merged, dir.name, '_edge_mutations.tab')
  tree.name <- paste0(tree.out.merged, dir.name, '_edges.tab')
  seq.name <- paste0(tree.out.merged, dir.name, '_sequences')
  
  # merge tree and mutation files
  for(i in 1:length(mut.files)){
    mut.dt <- fread(paste0(tree.out, mut.files[i]))
    tree.dt <- fread(paste0(tree.out, tree.files[i]))
    if(!file.exists(mut.name)){
      fwrite(mut.dt, file = mut.name, sep="\t", append = F, row.names = F, col.names = T)
      fwrite(tree.dt, file = tree.name, sep="\t", append = F, row.names = F, col.names = T)
    }else{
      fwrite(mut.dt, file = mut.name, sep="\t", append = T, row.names = F, col.names = F)
      fwrite(tree.dt, file = tree.name, sep="\t", append = T, row.names = F, col.names = F)
    }
  }
  # merge sequence files
  for(i in 1:length(seq.files)){
    seq.dt <- read.FASTA(seq.files[i], tree.out)
    if(!file.exists(paste0(seq.name, '.fasta'))){
      write.FASTA(paste0(dir.name, '_sequences'), tree.out.merged, seq.dt, append=F)
    }else{
      write.FASTA(paste0(dir.name, '_sequences'), tree.out.merged, seq.dt, append=T)
    }
  }
  
}