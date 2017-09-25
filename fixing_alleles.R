# runFixAlleles - find the dominant allele of each gene after running the igblast
#
# params - 
# igblastFile - the csv files that created after running the igblast.
# out.dir - the directory of the csv file
runFixAlleles<- function(igblastFile, out.dir){
  out.name<-NULL
  df.igblast<- NULL
  
  #merage all the igblast file to one df
  for (m in 1:length(igblastFile)){
    out.name[m] <- paste0(out.dir, igblastFile[m])
    df.temp<- read.csv(out.name[m], header = T, stringsAsFactors = F)
    df.temp$file<-igblastFile[m]
    df.igblast<-rbind(df.igblast,df.temp)
  }
  # extract all the different J genes
  allelesJ<- unique(substring(df.igblast$J_CALL,1, nchar(df.igblast$J_CALL)-3))
  
  #update the dominant allele in the J genes and the Jstart
  for (i in 1:length(allelesJ)){
    j<- df.igblast$J_CALL
    countAllelesJ<- table(j[which(substring(j,1,nchar(j)-3)==allelesJ[i])])
    names<- names(countAllelesJ)
    #get the dominant allele for a specific J gene
    dominant_allele<-names[which.max(as.vector(countAllelesJ))]
    dominant_allele_Jstart<- df.igblast$J_START[df.igblast$J_CALL==dominant_allele][1]
    indexes<- which(substring(df.igblast$J_CALL,1,nchar(df.igblast$J_CALL)-3)== substring(dominant_allele, 1, nchar(dominant_allele)-3))
    df.igblast$J_CALL[indexes]<- dominant_allele
    df.igblast$J_START[indexes]<- dominant_allele_Jstart 
  }
  # extract all the different v genes
  allelesV<- unique(substring(df.igblast$V_CALL,1, nchar(df.igblast$V_CALL)-3))
  
  #update the dominant allele in the V genes and the Vend
  for (i in 1:length(allelesV)){
    v<- df.igblast$V_CALL
    countAllelesV<- table(v[which(substring(v,1,nchar(v)-3)==allelesV[i])])
    names<- names(countAllelesV)
    #get the dominant allele for a specific V gene
    dominant_allele<-names[which.max(as.vector(countAllelesV))]
    dominant_allele_Vend<- df.igblast$V_END[df.igblast$V_CALL==dominant_allele][1]
    indexes<- which(substring(df.igblast$V_CALL,1,nchar(df.igblast$V_CALL)-3)== substring(dominant_allele, 1, nchar(dominant_allele)-3))
    df.igblast$V_CALL[indexes]<- dominant_allele
    df.igblast$V_END[indexes]<- dominant_allele_Vend  
  }
  #update the distance between V and J 
  df.igblast$VJ_DIST<- df.igblast$J_START-df.igblast$V_END-1
  
  #updating the csv files after fixing the allels 
  for (m in 1: length(igblastFile)){
    df.temp<- subset(df.igblast, df.igblast$file == igblastFile[m])
    df.temp$file<-NULL
    write.csv(df.temp, out.name[m], row.names = FALSE)
  }
  
  
  
}