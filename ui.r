library(shinyjs)
library(shiny)

#filenames <- list.files(path = "data",pattern="\\.csv$")
#names(filenames) <- gsub(pattern = "\\.csv$", "", filenames)
shinyUI(fluidPage(theme = "bootstrap.css",
                  (navbarPage("MAMGED BETA Version",
                              position = c("fixed-top"),
                              fluid=TRUE,
                              navbarMenu("Help",
                                         tabPanel(
                                           a("Packages Required",
                                             target="_blank",href="Stand-alone.pdf")
                                           ),
                                         tabPanel(
                                           a("Reference Manual",
                                             target="_blank", href="MAMGEDManual.pdf")
                                           )
                                         ),
                              navbarMenu("Sample Data",
                                         tabPanel(
                                           downloadButton("AffymetrixData", "Affymetrix")
                                           ),
                                         tabPanel(
                                           downloadButton("CodelinkData", "Codelink")
                                           ),
                                         tabPanel(
                                           downloadButton("IlluminaData", "Illumina")
                                           )
                                         ),
                              navbarMenu("Stand-Alone Version",
                                         tabPanel("Code", "mamged"),
                                         tabPanel(
                                           a("Stand-alone Manual",
                                             target = "_blank", href= "Stand-alone.pdf")
                                         )
                              )
                  )
                  ),
                  
                  br(),
                  titlePanel(
                    headerPanel( 
                      h3("Meta-Analysis of Microarray Gene Expression Data",
                         align="center", style="bold")
                      )
                    ),
                  
                  br(),
                  br(),
                  useShinyjs(),  ## initialize shinyjs to reset input files.
                  sidebarLayout(
                    sidebarPanel(
                      h5("Upload Data Files",style="bold"),
                      fileInput("files", 
                                "Choose CSV/txt processed files or raw files",
                                multiple = "TRUE",
                                accept=c('text/csv',
                                         'text/comma-separated-values,
                                         text/plain', '.csv','.cel','.TXT','.txt')),
                      #selectInput('dataset',"Choose platform annotation file", c("Please select a file" ='',filenames)),
                      fluidRow(
                        column(5,
                               radioButtons("radio", label = h5("Data uploaded"),
                                            choices = list("Affymetrix" = 1, "Codelink" = 2,
                                                           "Illumina" = 3),selected = NULL)
                               )
                        ),
                      fluidRow(
                        column(10,
                               h5("Differential Expression Call", style = "bold"),
                               checkboxInput("checkbox",
                                             label = "Differential Expression", value = FALSE))),
                      br(),
                      uiOutput('AnnotationFile'),
                      uiOutput('annotation_range'),
                      uiOutput('hide'),
                      uiOutput('FDRAffymetrix'),
                      uiOutput('FDRCodlink'),
                      br(),
                      
                      tags$head(tags$script(src = "message-handler.js")),
                      fluidRow(
                        column(4,
                               actionButton("Submit", label = "Submit")),
                        column(4,
                               actionButton("Reset_Input", label = "Reset"))
                        ),
                      br(),
                      br(),
                      titlePanel(
                        fluidRow(
                          column(5,
                                 h5("Downloads", style="bold"),
                                 downloadButton('downloadData1', 'Summary')),
                          br(),
                          column(7,
                                 downloadButton("downloadData2",'Supplementary')
                                 )
                          )
                        ),
                      br(),
                      br(),
                      h4("Development:",style = "color:#0000"),
                      p("The application is developed in R.",
                        "This version deals with",
                        a("Affymetrix", 
                          href = "http://en.wikipedia.org/wiki/Affymetrix"),
                        ",",
                        a("Codelink",
                          href="http://www.appliedmicroarrays.com/back_up/pdf/CodeLink_Human_Whole_Genome_Bioarray.pdf"),
                        "and",
                        a("Illumina",
                          href= "http://www.illumina.com/applications/transcriptome-analysis
                          /gene-expression-analysis/gene-expression-arrays.html"),
                        "gene expression data ."),
                      p("Each of the three data platform type has three functionalities, which cover meta-analysis of
                         processed data, processing raw data to perform meta-analysis
                         and finally differential expression"),
                      br(),
                      br(),
                      
                      img(src = "bigorb.png", height = 100, width = 100),
                      br(),
                      h5("Powered by:", 
                         span("HLS and team", colour="Yellow"))),
                    
                    mainPanel(
                      tabsetPanel(id = "MamgedTabs",
                        #tabPanel("File-list",dataTableOutput("file")),
                        tabPanel("Source-data", dataTableOutput("sourced")),
                        tabPanel("Annotation-data",dataTableOutput("annotation")),
                        tabPanel("Summary",dataTableOutput("final")),
                        tabPanel("Supplementary file", dataTableOutput("full"))
                      )
                    )
                    )))