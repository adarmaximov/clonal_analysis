library(shiny)
library(mailR)
options(shiny.maxRequestSize=300*1024^2)

shinyServer(function(input,output, session){
  
  observe({
    
    if (input$do[1] == 0 || is.null(input$file))  {
      return()
    }else{
      DATA.PATH <-'/media/raid10/run_user3/'
      w.dir = '/media/raid10/run_user3/BuildTrees/'
      setwd(w.dir)
       in.dir <- paste0(DATA.PATH, input$organism, '/', input$directory, '/') # change if fastq...
     # if (!(dir.exists(in.dir))){

        #create this direcrotry and put files in paste0(in.dir, 'raw_fasta')
        dir.create(in.dir, recursive = T, showWarnings = F)
        temp<-paste0(in.dir, 'raw_fasta', '/')
        dir.create(temp, recursive = T, showWarnings = F)
        #create paste0(in.dir, 'raw_fasta') and upload all fasta files that the user choose into it
        for (i in 1:(nrow(input$file))){
          file.copy(input$file[[i,'datapath']], paste0(temp, input$file[[i,'name']]))
        }
        source('run_pipeline.R')
        BCellAnalysis(in.dir, input$directory, input$organism, input$chain, input$fastq, input$merge, input$IgBLAST, input$VJ.FASTA, input$TREES, input$FORMAT)
        sender <- "adarmaximov@gmail.com"
        recipients <- c(input$email)
        msg <- paste0('The analysis finshed! the directory to your reasult is: /home/run_user3/Documents/', input$organism, '/', input$directory)
        send.mail(from = sender,
                  to = recipients,
                  subject="Reasult for B_cell analysis",
                  body = msg,
                  smtp = list(host.name = "smtp.gmail.com", port = 465, 
                              user.name="adarmaximov@gmail.com", passwd="9093668gal", ssl=TRUE),
                  authenticate = TRUE,
                  send = TRUE)
  #    }else {
   #     stop("This directory is already exists")
    #  }
    }
  })
})