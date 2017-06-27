setwd('/home/adar/Documents/App1')

in.dir <- '/media/raid10/adar/data/human/BM_PBL/'
dir.name='BM_PBL'
organism ='human'
chain='VH'

source('clone_analysis.R')
main(in.dir, dir.name, organism, chain)

