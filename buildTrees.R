library(ape)
library(igraph)
library(phangorn)
library(seqinr)
library(Biostrings)
library(parallel)
library(data.table)
library(stringr)


MIN.SEQ <- 2
MAX.SEQ <- 100

# Create fasta file from clone data.frame
#
# Params:  clone.df = data.frame of clone with [taxa, seq] columns
#          out.file = the file name to write to
#
# Returns: NULL
write.FASTA <- function(file.name, out.dir, nameSeq.dt, append = F) {
  sequences <- nameSeq.dt[,seq]
  names(sequences) <- nameSeq.dt[,head]
  if(append)
    writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, file.name, '.fasta'), append=T, width=1000)
  else
    writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, file.name, '.fasta'), width=1000)
}

# run neighbor joining
run.neighbor <- function(dna.seq){
  
  # make distance matrix
  dist.mat <- dist.dna(dna.seq, model='K80')
  
  # run NJ 
  tree <- NJ(dist.mat)
  
  return(tree)
}

# convert matrix to DNA sequence
mat.to.seq <- function(mat){
  
  nt <- c('A', 'C', 'G', 'T')
  
  return(paste0(unlist(apply(mat, 1, function(x) nt[x==1])), collapse=''))
}

rand.nt <- function(r){
  ind <- which(max(r)==r)
  s.nt <- sample(1:length(ind),1)
  r[ind[s.nt]] <- 1
  r[ind[-s.nt]] <- 0
  return(r)
}
choose.nt <- function(r.node, r.father){
  # check if father's nt is in son's possibilities
  ind <- which(max(r.node)==r.node)
  f.ind <- which(r.father==1)
  if(f.ind%in%ind){ # choose father's nt
    r.node[f.ind] <- 1
    r.node[-f.ind] <- 0
  }else{
    r.node <- rand.nt(r.node)
  }
  return(r.node)
}

fix.amb <- function(anc.seq, edges, node, father){
  # check if node is not a leaf
  if(node%in%edges[,1]){
    # check if root
    if(father==0){
      mat <- anc.seq[[node]]
      ind <- which(rowSums(mat==1)==0) # which rows have ambiguities
      if(length(ind)==1){
        mat[ind,] <- rand.nt(mat[ind,])
      }else if(length(ind)>1){
        mat[ind,] <- t(apply(mat[ind,], 1, rand.nt))
      }
    }else{
      mat <- anc.seq[[node]]
      mat.father <- anc.seq[[father]]
      ind <- which(rowSums(mat==1)==0) # which rows have ambiguities
      if(length(ind)==1)
        mat[ind,] <- choose.nt(mat[ind,], mat.father[ind,])
      else if(length(ind)>1)
        mat[ind,] <- t(sapply(ind, function(x) choose.nt(mat[x,], mat.father[x,])))
    }
    anc.seq[[node]] <- mat
    sons <- edges[edges[,1]%in%node,2]
    
    # recursive call for each son
    for(i in 1:length(sons)){
      anc.seq <- fix.amb(anc.seq, edges, sons[i], node)
    }
    
  }
  return(anc.seq)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# get internal sequences
get.internal.nodes <- function(tree, dna.seq){
  
  # convert sequences to phyDat format
  phy.seq <- as.phyDat(dna.seq)
  
  # get internal states
  anc.seq <- ancestral.pars(tree, phy.seq, type = 'ACCTRAN')
  seq.names <- names(anc.seq)
  seq.ind <- setdiff(1:length(seq.names),c(grep('VJ', seq.names),grep('Seq', seq.names)))
  # find root --> node with 3 sons
  edges <- tree$edge
  root <- Mode(edges[,1])
  anc.seq <- fix.amb(anc.seq, edges, root, 0)
  
  # take only inner nodes
  inn.nodes <- subset(anc.seq, seq.ind)
  # get pattern indices
  ind <- attr(anc.seq, 'index')
  # convert matrix of patterns to characters
  patt.chars <- strsplit(unlist(lapply(inn.nodes, mat.to.seq)), '')
  # get inner node full sequences
  inn.seq <- unlist(lapply(patt.chars, function(x) paste0(x[ind], collapse='')))
  # get leave full sequences
  leave.seq <- toupper(sapply(dna.seq, paste, collapse=""))
  
  return(c(leave.seq, inn.seq))
}


create.fromTo <- function(tree, seq.char, file.name){
  
  nam <- names(seq.char)
  dt <- data.table(Tree=rep(file.name, nrow(tree$edge)), from='', to='')
  
  dt <- dt[,  c('from', 'to') := as.data.table(t(apply(tree$edge, 1, function(x) nam[x])))]
  
  
  return(dt)
}

# compute edge lengths
#
# Params:  nameSeq.dt - data frame with headers and sequences
#          edge.dt - data frame with columns - parent(from), child(to), edge weight (weight)
#
# Returns: edge.dt - data frame with columns - parent(from), child(to), edge weight (weight) , edge length(distance)
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

# Main function - Build phylogenetic tree with maximum parsomony or neigbor joining 
#
# Params:  file.name - original FASTA file name
#          input directory for 
#
# Returns: -

build.trees <- function( dir.name, in.dir, out.dir, V, VJdis.files ){
  
  #VJdis.files <- list.files(in.dir)
  Vs <- str_sub(str_extract(VJdis.files, pattern='[0-9]{3}\\.[0-9]{3}'),1, 3)
  VJdis.files <- VJdis.files[Vs==V]
  if(length(VJdis.files)==0)
    return()
  for(i in 1:length(VJdis.files)){
    FASTA.file <- VJdis.files[i]
    # read FASTA file as DNAbin format
    fasta.dna_total <- read.dna(paste0(in.dir, FASTA.file), format='fasta')
    
    # skip small or big files
    n.seq <- nrow(fasta.dna_total)
    
    # make distance matrix
    if (n.seq>1){
      dist.mat <- as.matrix(dist.hamming(as.phyDat(fasta.dna_total)))#dist.dna(fasta.dna, model='K80', as.matrix=T)
      bad <- which(rowSums(is.nan(dist.mat))>(nrow(dist.mat)/4))
      if(length(bad)>0){
        fasta.dna_total <- fasta.dna_total[-bad,]
        dist.mat <- as.matrix(dist.hamming(as.phyDat(fasta.dna_total)))#dist.dna(fasta.dna, model='K80', as.matrix=T)
        bad2 <- unique(which(is.nan(dist.mat), arr.ind=T)[,1])
        if(length(bad2)>0){
          fasta.dna_total <- fasta.dna_total[-bad2,]
          dist.mat <- dist.dna(fasta.dna_total, model='K80', as.matrix=T)
        }
      }
    }
    is.Prim<-TRUE
    numOfComponents<- n.seq
    if( n.seq > MIN.SEQ) { 
      print ("Clustering: RUN PRIM")
      tresh_old_clustering<- 8/ncol(fasta.dna_total)
      dist.mat[dist.mat>tresh_old_clustering]<- 0
      if(sum(dist.mat)==0){
        is.Prim<-FALSE
      } else {
      g<- graph.adjacency(dist.mat, weighted = TRUE)
      g_mst<- mst(g, algorithm = 'prim')
      c<- components(g_mst)
      numOfComponents<- count_components(g_mst)
      }
    } else{
      is.Prim<-FALSE
    }
    
    
    for (j in 1: numOfComponents){
      
      if(is.Prim){
        indexes<- as.vector(which(c$membership==j))
        fasta.dna<- fasta.dna_total[indexes,]
      } else{
        fasta.dna<- fasta.dna_total[j,]
      }
      
      # skip small or big files
      n.seq <- nrow(fasta.dna)
      
      file.name <- paste0(substr(FASTA.file, 1, nchar(FASTA.file)-6),'.',sprintf('%03d',j))
      
      if( n.seq <= MIN.SEQ) { 
        print(paste0(FASTA.file,' - Not enough sequences to make tree with Phylip'))
        # only 1 sequence - fix header
        head <- paste0(file.name, '_', dimnames(fasta.dna)[[1]])
        seq <- toupper(sapply(fasta.dna, paste, collapse=""))
        fasta.dt <- data.table(head=head, seq=seq)
        # save 
        if(file.exists(paste0(out.dir, dir.name, '_sequences_', V, '.fasta')))
          write.FASTA(paste0(dir.name, '_sequences_', V), out.dir, fasta.dt, append = T)
        else
          write.FASTA(paste0(dir.name, '_sequences_', V), out.dir, fasta.dt, append = F)
        
        next()
      }
      
      # if number of sequence is between MIN.SEQ and MID.SEQ - build tree with Maximum Parsimony
      if( n.seq > (MIN.SEQ + 1) & n.seq < MAX.SEQ ){
        print(paste0(FASTA.file, ' - maximum parsimony'))
        
        # run maximum parsimony
        tree <- pratchet(as.phyDat(fasta.dna))
        
        # if number of sequence is between MID.SEQ and MAX.SEQ - build tree with Neigbor joining
      }else if( n.seq >= MAX.SEQ | n.seq == (MIN.SEQ+1)){
        print(paste0(FASTA.file, ' - neighbor joining'))
        
        # make distance matrix
        dist.mat <- as.matrix(dist.hamming(as.phyDat(fasta.dna)))#dist.dna(fasta.dna, model='K80', as.matrix=T)
        bad <- which(rowSums(is.nan(dist.mat))>(nrow(dist.mat)/4))
        if(length(bad)>0){
          fasta.dna <- fasta.dna[-bad,]
          dist.mat <- as.matrix(dist.hamming(as.phyDat(fasta.dna)))#dist.dna(fasta.dna, model='K80', as.matrix=T)
          bad2 <- unique(which(is.nan(dist.mat), arr.ind=T)[,1])
          if(length(bad2)>0){
            fasta.dna <- fasta.dna[-bad2,]
            dist.mat <- dist.dna(fasta.dna, model='K80', as.matrix=T)
          }
        } 
        
        # run NJ 
        tree <- NJ(dist.mat)
        
      }
      print('Get internal nodes:')
      seq.char <- get.internal.nodes(tree, fasta.dna)
      
      edge.dt <- create.fromTo(tree, seq.char, file.name)
      nameSeq.dt <- data.table(head= sapply(names(seq.char), function(x) paste0(file.name, '_', x)), 
                               seq=seq.char)
      
      # compute edge lengths
      edge.dt <- compute.edge(nameSeq.dt, edge.dt)
      
      # save sequences as FASTA file
      if(file.exists(paste0(out.dir, dir.name, '_sequences_', V, '.fasta')))
        write.FASTA(paste0(dir.name, '_sequences_', V), out.dir, nameSeq.dt, append = T)
      else
        write.FASTA(paste0(dir.name, '_sequences_', V), out.dir, nameSeq.dt, append = F)
      # save Fome/To/distances table (edges)
      if(file.exists(paste0(out.dir, dir.name, '_edges_', V, '.tab')))
        fwrite(edge.dt, file=paste0(out.dir, dir.name, '_edges_', V, '.tab'), quote=F, sep='\t', col.names=F, row.names=F, append = T)
      else
        fwrite(edge.dt, file=paste0(out.dir, dir.name, '_edges_', V, '.tab'), quote=F, sep='\t', col.names=T, row.names=F) 
    }
  }
  return()
}