BCellAnalysis<- function(in.dir,dir.name,organism,
                         chain,fastq, merge,IGBLAST, VJ.FASTA,TREES, FORMAT){
  if(fastq){
    source('analyze_fastq.R')
    main(in.dir, dir.name, organism, chain, merge)
  }
  # if input are fasta files - put raw fasta files in paste0(in.dir, 'raw_fasta')
  source('clone_analysis.R')
  main(in.dir, dir.name, organism, chain, IGBLAST, VJ.FASTA, TREES, FORMAT)
}

# organism <- 'human'
# chain <- 'VH'
# fastq <- T # if input files are fastq
# merge <- T # if input files are fastq but in paired end - need to merge
# 
# IGBLAST = T
# VJ.FASTA = T
# TREES = T
# FORMAT = T

#dir.name <- 'M29_BM_PBL' # project name (provided by user)

#DATA.PATH <- '/media/raid10/jennifer/Analyzed_Bcells/' # directory of analyzed datasets

#in.dir <- paste0(DATA.PATH, organism, '/', dir.name, '/') # home directory of project
#dir.create(in.dir, recursive = T, showWarnings = F)

#w.dir <- '/home/adar/clonal_analysis/' # working directory (contains all scripts)


# if input are fastq files - put raw fastq files in paste0(in.dir, 'raw_fastq')

