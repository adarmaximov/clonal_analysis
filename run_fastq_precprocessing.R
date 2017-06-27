organism <- 'mouse'
chain <- 'VK'
merge <- F

dir.name <- 'lamina_propria'

DATA.PATH <- '/media/raid10/jennifer/Analyzed_Bcells/'
in.dir <- paste0(DATA.PATH, organism, '/', dir.name, '/') # change if fastq...
# create this direcrotry and put files in paste0(in.dir, 'raw_fasta')
dir.create(in.dir, recursive = T, showWarnings = F)


w.dir <- '/media/raid10/jennifer/BuildTrees/merged_version3/'
setwd(w.dir)

# put raw fastq files in paste0(in.dir, 'raw_fastq')

source('analyze_fastq.R')
#main(in.dir, dir.name, organism, chain, merge)

source('clone_analysis.R')
main(in.dir, dir.name, organism, chain)
