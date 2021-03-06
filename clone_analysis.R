
# imports
#-------------------
require(seqinr)
require(Biostrings)
require(stringr)

# directories
#-------------------
IGBLAST.PATH <- '/IgBLAST/ncbi-igblast-1.7.0/' # contains IgBLAST program
GERMLINE.PATH <- '/germline/' # contains germline sequences and IMGT CDR/FWR regions

# flags
#---------------
TREE.THRESHOLD <- 4

# Summarize project run and output to file
#
# Params:   in.dir = FASTA file directory
#           dir.name = project name
#           Vlen = maximum length of V (numeric)
#           Jlen = maximum length of J (numeric)
#           del = number of removed sequences after igblast
#           clone.info = clone and sequence counts
#
# Returns:  -
print.info <- function(in.dir, dir.name, Vlen, Jlen, del, clone.info){
  out.file <- paste0(in.dir, 'RUN_INFO_', dir.name, '.tab')
  
  # print project name
  cat(paste0('Project: ', dir.name, '\n'), file=out.file)
  
  # print V and J cuts
  cat(paste0('Max V: ', Vlen, ', Max J: ', Jlen, '\n'), file=out.file, append=T)
  
  # print number of deleted sequences
  cat(paste0('Number of sequences removed:', del, '\n'), file=out.file, append=T)
  
  # print clone info table
  cat(colnames(clone.info), file=out.file, append=T, sep='\t')
  cat('\n', file=out.file, append=T)
  for(i in 1:nrow(clone.info)){
    cat(unlist(clone.info[i,]), file=out.file, append=T, sep='\t')
    cat('\n', file=out.file, append=T)
  }
  
}

# main function - run entire pipeline
#
# Params:   in.dir = FASTA file directory
#           dir.name = project name
#           organism = human or mouse
#           chain = VH/VK/VL
#           cp.num.del = if sequence header contains copy number - 
#                        then it must appear at the end of the header with a specified delimiter
#
# Returns:  -
main <- function(in.dir, dir.name, organism, chain, IGBLAST = T, VJ.FASTA = T, TREES = T, FORMAT = T, cp.num.del = NULL){
  
  # set directories
  w.dir <- getwd()
  germline.path <- paste0(w.dir, GERMLINE.PATH)
  igblast.path <- paste0(w.dir, IGBLAST.PATH)
  
  # load germline sequences
  Vgerm <- read.fasta(paste0(germline.path, organism, '/', chain, '/Vgermline_', organism, '_', chain, '.fasta'), 
                      seqtype = "DNA", as.string = T, forceDNAtolower = F, set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)  
  Jgerm <- read.fasta(paste0(germline.path, organism, '/', chain, '/Jgermline_', organism, '_', chain, '.fasta'), 
                      seqtype = "DNA", as.string = T, forceDNAtolower = F, set.attributes = F, legacy.mode = T, seqonly = F, strip.desc = F)
  
  # directory of FASTA files
  fasta.dir <- paste0(in.dir, 'raw_fasta/')
  
  # get file(s) in directory
  fasta.files <- list.files(fasta.dir, pattern='fasta$')
  
  #---------------------------------------------------------------------------------------------------
  # PART I - Runs IgBLAST, parse output file and create csv file
  #---------------------------------------------------------------------------------------------------
  
  # load IgBLAST functions
  source(paste0(w.dir, '/runIgblast.R'))
  
  # create IgBLAST output directory
  igblast.out <- paste0(in.dir,'IGBLAST/') 
  dir.create(igblast.out, recursive = T, showWarnings = F)
  
  if(IGBLAST){
    print('Run and parse IgBLAST')
    print (Sys.time())

    #run and parse IgBLAST for each FASTA file
    for(i in 1:length(fasta.files)){
      res <- igblast(fasta.files[i], fasta.dir, igblast.out , igblast.path, Vgerm, Jgerm, chain, organism)
    }
    #run the fix allele 
    igblast.files<- list.files(igblast.out, pattern = 'csv$')
    source(paste0(w.dir,'/fixing_alleles.R'))
    runFixAlleles(igblast.files, igblast.out)
  }
  
  # determine Vlen and Jlen from data
  print('Determine Vlen and Jlen')
  print (Sys.time())
  len <- check.VJ.len(igblast.out, Vgerm, Jgerm )
  Vlen <- len[[1]]
  Jlen <- len[[2]]
  
  #---------------------------------------------------------------------------------------------------
  # PART II - Read csv files and create FASTA files for each V-J-distance group
  #---------------------------------------------------------------------------------------------------
  source(paste0(w.dir, '/makeFASTAfiles.R'))
  
  # create VJdis output directory
  VJdis.out <- paste0(in.dir, Vlen, '.', Jlen, '/VJdis/')
  dir.create(VJdis.out, recursive = T, showWarnings = F)
  
  del <- 0
  if(VJ.FASTA){
    print('Create V-J-distance files')
    print (Sys.time())
    for(i in 1:length(fasta.files)){
      print(paste0(i, '/', length(fasta.files), ': ', fasta.files[i]))
      tmp.del <- makeFASTAfiles(dir.name, fasta.files[i], igblast.out, VJdis.out, Vgerm, Jgerm, Vlen, Jlen, cp.num.del)
      del <- del + tmp.del
    }
  }
  
  #---------------------------------------------------------------------------------------------------
  # PART III - Build phylogenetic tree for each V-J-distance group and parse output tree files
  #---------------------------------------------------------------------------------------------------
  source(paste0(w.dir, '/buildTrees.R'))
  
  # create output directory
  tree.out <- paste0(in.dir, Vlen, '.', Jlen, '/Trees_TMP/')
  dir.create(tree.out, recursive=T, showWarnings = FALSE)
  
  # get all V-J-distance FASTA files
  VJdis.files <- list.files(path = VJdis.out, pattern='*.fasta', full.names = F)
  # get Vs
  Vs <- str_sub(str_extract(VJdis.files, pattern='[0-9]{3}\\.[0-9]{3}'),1, 3)
  uni.Vs <- unique(Vs)
  # separate big files
  VJdis.files.size <- file.size(paste0(VJdis.out, VJdis.files))
  VJdis.files.big<- VJdis.files[VJdis.files.size>=1000000 & VJdis.files.size < 5000000]# throw huge trees
  VJdis.files <- VJdis.files[VJdis.files.size<1000000]
  
  # run "small" files in parrallel
  if(TREES){
    print('Build phylogenetic trees')
    print (Sys.time())
    # run buildTrees in parrallel
    print('start cluster')
    cl <- makeCluster(4)
    clusterCall(cl,function() {source('buildTrees.R')})
    clusterExport(cl, list('dir.name', 'VJdis.out', 'tree.out', 'VJdis.files'), envir=environment())
    print('parSapply starts')
    parSapply(cl, uni.Vs, function(x) build.trees(dir.name, VJdis.out, tree.out, x, VJdis.files))
    stopCluster(cl)
    print (Sys.time())
    
    # for(i in 1: length(uni.Vs)){
    #   if (i==9){
    #     print (i)
    #   }
    #   build.trees(dir.name, VJdis.out, tree.out,uni.Vs[i], VJdis.files)
    # }
    # print(Sys.time())

    # run big trees sequencially
    print('run big trees')
    for(i in 1:length(uni.Vs))
      build.trees(dir.name, VJdis.out, tree.out, uni.Vs[i], VJdis.files.big)
    print (Sys.time())
    
  }
  
  #---------------------------------------------------------------------------------------------------
  # PART IV - Cut trees into clones and convert all sequences into new format 
  #---------------------------------------------------------------------------------------------------
  source(paste0(w.dir, '/convertFormat.R'))
  
  # create temporary output directory
  seq.out <- paste0(in.dir, Vlen, '.', Jlen, '/Sequences_TMP/')
  clones.out <- paste0(in.dir, Vlen, '.', Jlen, '/Clones_TMP/')
  
  dir.create(seq.out, recursive=T, showWarnings = FALSE)
  dir.create(clones.out, recursive=T, showWarnings = FALSE)
  
  # load germline CDR and FWR regions info
  regions <- fread(paste0(germline.path, organism, '/', chain, '/CDR_FR_regions_IgBLAST_', organism, '_', chain, '.csv'),header=T)
  
  if(FORMAT){
    print('Format clone files')
    #if(!file.exists(paste0(tree.out, dir.name, '_sequences.fasta'))){
    #  print('no trees!')
    #  return()
    #}
    
    ## split files by V gene
    #uni.Vs <- split.tree.files(dir.name, tree.out)
    
    # run convert.format in parrallel
    print('start cluster')
    cl <- makeCluster(3)
    clusterCall(cl,function() {source('convertFormat.R')})
    clusterExport(cl, list('dir.name', 'tree.out', 'seq.out', 'clones.out', 'Vlen', 'Jlen', 'Vgerm', 'Jgerm', 'regions', 'TREE.THRESHOLD'), envir=environment())
    print('parSapply starts')
    parSapply(cl, uni.Vs, function(x) convert.format(x, dir.name, tree.out, seq.out, clones.out, Vlen, Jlen, Vgerm, Jgerm, regions, TREE.THRESHOLD))
    stopCluster(cl)
    # for(i in 1: length(uni.Vs)){
    # convert.format(uni.Vs[i],dir.name, tree.out, seq.out, clones.out, Vlen, Jlen, Vgerm, Jgerm, regions, TREE.THRESHOLD)
    # }
    
    print (Sys.time())
    
    # merge output files
    seq.out.merged <- paste0(in.dir, Vlen, '.', Jlen, '/Sequences/')
    clones.out.merged <- paste0(in.dir, Vlen, '.', Jlen, '/Clones/')
    tree.out.merged <- paste0(in.dir, Vlen, '.', Jlen, '/Trees/')
    
    dir.create(seq.out.merged, recursive=T, showWarnings = F)
    dir.create(clones.out.merged, recursive=T, showWarnings = F)
    dir.create(tree.out.merged, recursive=T, showWarnings = F)
    
    merge.out.files(dir.name, seq.out, clones.out, seq.out.merged, clones.out.merged)
    merge.tree.mutation.files(dir.name, tree.out, tree.out.merged)
    
    # get total clone/sequence number for summary
    clone.info <- get.clone.info(dir.name, seq.out.merged, clones.out.merged)
    
    # save all info about current run
    print.info(in.dir, dir.name, Vlen, Jlen, del, clone.info)
    
    # remove temporary files
    #unlink(seq.out, recursive = T)
    #unlink(clones.out, recursive = T)
    #unlink(tree.out,  recursive = T)
  }
  
}
