library(shiny)
# ui.R

shinyUI(fluidPage(
  titlePanel(strong("B cells analysis:", style = "color:blue")),
  
  sidebarLayout(
    sidebarPanel(
      
      fileInput("file",strong("Upload the file", style = "color:navy"), multiple = TRUE), # fileinput() function is used to get the file upload contorl option
      
      textInput("directory", label = strong("Please insert directory name", style = "color:navy")),
      
      # Copy the line below to make a text input box
      textInput("email", label = strong("Please insert your e-mail address", style = "color:navy")),
      helpText(strong("please choose input files format", style = "color:navy")),
               checkboxInput("fastq", label = "FASTQ", FALSE),
               checkboxInput("merge", label = "merge", FALSE),
               selectInput("organism", label = strong("Select organism", style = "color:navy"), 
                           choices = list("human" = "human", "mouse" = "mouse")),
               
               radioButtons("chain", label = strong("Please choose a chain option", style = "color:navy"),
                            choices = list("VH" = "VH", "VK" = "VK", "VL" = "VL", "VK/VL" = "VK/VL"), 
                            selected = 1),
               helpText(strong("Please choose which session you would like to run:", style = "color:navy")),
                        checkboxInput("IgBLAST", label = "IgBLAST", FALSE),
                        checkboxInput("VJ.FASTA", label = "VJ.FASTA", FALSE),
                        checkboxInput("TREES", label = "TREES", FALSE),
                        checkboxInput("FORMAT", label = "FORMAT", FALSE),
                        
                        
                        actionButton("do", "Click Me")
               ), 
               mainPanel(
                 tableOutput('table1')
               )
      )
    ))