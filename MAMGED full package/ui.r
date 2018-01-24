library(shiny)
library(shinyjs)
library(simpleaffy)
library(affy)
library(oligo)
library(markdown)
library(plyr)
library(data.table)
library(dplyr)
library(codelink)
library(lumi)
library(h20kcod.db)
library(h10kcod.db)
library(hwgcod.db)
library(limma)
library(tibble)
library(dtplyr)
library(tools)
library(affyio)
library(gtools)

filenames <- list.files(path = "data",pattern="\\.txt$")
names(filenames) <- gsub(pattern = "\\.txt$", "", filenames)
shinyUI(fluidPage(theme = "bootstrap.css",
                  (navbarPage("MAMGED BETA Version",
                              position = c("fixed-top"),
                              fluid=TRUE,
                              navbarMenu("Help",
                                         tabPanel(
                                           a("Reference Manual",
                                             target="_blank", href="MAMGEDManual.pdf")
                                           ),
                                         tabPanel(
                                           a("GPLs Supported",
                                             target="_blank", href="gpl.pdf")
                                         ),
                                         tabPanel(
                                           a("Video Tutorials",
                                             downloadLink("AbsoluteExpression", " Absolute Expression", class=" fa fa-cloud-download"),
                                             downloadLink("DifferentialExpression", " Differential Expression", class=" fa fa-cloud-download")
                                             )
                                           )
                                         ),
                              navbarMenu("Sample Data",
                                         tabPanel(
                                           downloadLink("AffymetrixData", " Affymetrix", class=" fa fa-cloud-download")
                                           ),
                                         tabPanel(
                                           downloadLink("CodelinkData", " Codelink", class=" fa fa-cloud-download")
                                           ),
                                         tabPanel(
                                           downloadLink("IlluminaData", " Illumina", class=" fa fa-cloud-download")
                                           )
                                         ),
                              navbarMenu("Stand-Alone Version",
                                         tabPanel(
                                           downloadLink("CodeandData", " MAMGED", class=" fa fa-cloud-download")
                                           ),
                                         tabPanel(
                                           a("Stand-alone Manual",
                                             target = "_blank", href= "Stand-alone.pdf")
                                           )
                                         )
                                         
                              )
                   ),
                  
                  br(),
                  br(),
                  titlePanel(
                    headerPanel( #title=div(img(src="dna.jpg",  height = 120, width = 120),
                     h2( "Meta-Analysis of Microarray Gene Expression Data",
                         align="center", style="bold"
                      )
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
                      
                      uiOutput('Display_source_data'),
                    
                      #selectInput('dataset',"Choose platform annotation file", c("Please select a file" ='',filenames)),
                      fluidRow(
                        column(5,
                               radioButtons("radio", label = h5("Data uploaded"),
                                            choices = list("Affymetrix" = 1, "Codelink" = 2,
                                                           "Illumina" = 3),selected = 1)
                               )
                        ),
                      fluidRow(
                        column(10,
                               h5("Differential Expression Call", style = "bold"),
                               checkboxInput("checkbox",
                                             label = "Differential Expression", value = FALSE))),
                      br(),
                      
                      conditionalPanel(
                        condition = "input.checkbox != true | (input.checkbox == true && input.radio != 2)",
                        
                        # conditionalPanel(
                        #   condition = "input.checkbox != true",
                          wellPanel(
                            h5("Upload Annotation File"),
                            selectInput('dataset', "Choose platform annotation file",
                                        c("Please select a file" = '', filenames), multiple = TRUE))
                        ),
                        
                        
                        # conditionalPanel(
                        #   condition = "input.checkbox == true && input.radio != 2",
                        #   wellPanel(
                        #     h5("Upload Annotation File"),
                        #     selectInput('dataset1', "Choose platform annotation file",
                        #                 c("Please select a file" = '', filenames)))
                        # ),
                        
                       
                        
                        uiOutput('Display_annotation_file'),
                      # ),
                      
                      uiOutput('annotation_range'),
                      
                      ## Display Differential Expression Parameters *******************************
                      conditionalPanel(
                        condition = "input.checkbox == true && input.radio == 1",
                        wellPanel(
                          h5("Affymetrix Differential Expression"),
                          list(numericInput("num", label = h5("Fold Change"),
                                            value = NULL, min = 0.0, max = 10, step = 0.1),
                               numericInput("num1", label = h5("P-Value"),
                                            value = NULL, min = 0.001, max = 0.09, step = 0.01),
                               textInput("text", label = h5("Make Contrasts"),
                                         value = "Prostate-Normal")))
                      ),
                      
                      conditionalPanel(
                          condition = "input.checkbox == true && input.radio == 2",
                          wellPanel(
                            h5("Codelink Differential Expression"),
                            list(
                              textInput("text1", label = h5("Make Contrasts"), 
                                        value = "spermatogenesis2-spermatogenesis5"),
                              selectizeInput("select", label = h5("Choose Annotation"),
                                             choices = list("h10kcod.db" = "h10kcod",
                                                            "h20kcod.db" = "h20kcod",
                                                            "hwgcod.db" = "hwgcod"),
                                             options = list(placeholder = 'Choose DB',
                                                            onInitialize = I('function()
                                                         { this.setValue(""); }'))),
                              numericInput("num3", label = h5("Fold Change"),
                                           value = NULL, min = 0.0, max = 10, step = 0.1),
                              numericInput("num4", label = h5("P-Value"),
                                           value = NULL, min = 0.001, max = 0.09, step = 0.01)))
                        ),
                      conditionalPanel(
                        condition = "input.checkbox == true && input.radio == 3",
                        wellPanel(
                          h5("Illumina Differential Expression"),
                          list(
                            textInput("text3", label = h5("Make Contrasts"), value = "Enter conditions.."),
                            numericInput("num5", label = h5("Fold Change"),
                                         value = NULL, min = 0.0, max = 10, step = 0.1),
                            numericInput("num6", label = h5("P-Value"),
                                         value = NULL, min = 0.001, max = 0.09, step = 0.01)))
                      ),
                      
                      conditionalPanel(
                        condition = "input.checkbox == true && (input.radio == 1 | input.radio == 3)",
                        wellPanel(
                          h5("Additional information you want in supplementary file"),
                          list(
                            checkboxInput("checkbox1", label = "t-statistic", value = FALSE),    
                            checkboxInput("checkbox2", label = "adjusted p-value", value = FALSE),
                            checkboxInput("checkbox3", label = "B-statistic", value = FALSE),
                            actionButton("selectAll", label = "Select All"),
                            actionButton("deselectAll", label = "Deselect All"),
                            br(),
                            br()
                            #checkboxInput("checkbox4", label = "All", value = FALSE),
                            # radioButtons("radio1", label = h5("Choose False Discovery Control"),
                            #              choices = list("Benjamini & Hochberg (BH)" = "BH",
                            #                             "Benjamini & Yekutieli (BY)" = "BY",
                            #                             "Holm" = "holm",
                            #                             "Hochberg" = "hochberg",
                            #                             "Hommel" = "hommel",
                            #                             "Bonferroni" = "bonferroni",
                            #                             "FDR" = "fdr",
                            #                             "None" = "none"),
                            #              selected = NULL)
                            ))
                      ),
                        
                      conditionalPanel(
                        condition = "input.checkbox == true && input.radio == 2",
                        wellPanel(
                          h5("Additional information you want in supplementary file"),
                          list(
                            checkboxInput("checkbox4", label = "meanSNR", value= FALSE),
                            checkboxInput("checkbox5", label = "t-statistic", value = FALSE),
                            checkboxInput("checkbox6", label = "adjusted p-value", value = FALSE),
                            checkboxInput("checkbox7", label = "B-statistic", value = FALSE),
                            actionButton("selectAll1", label = "Select All"),
                            actionButton("deselectAll1", label = "Deselect All"),
                            br(),
                            br()
                            #checkboxInput("checkbox5", label = "All", value = FALSE),
                            # radioButtons("radio2", label = h5("Choose False Discovery Control"),
                            #              choices = list("Benjamini & Hochberg (BH)" = "BH",
                            #                             "Benjamini & Yekutieli (BY)" = "BY",
                            #                             "Holm" = "holm",
                            #                             "Hochberg" = "hochberg",
                            #                             "Hommel" = "hommel",
                            #                             "Bonferroni" = "bonferroni", 
                            #                             "FDR" = "fdr",
                            #                             "None" = "none"),
                            #              selected = NULL)
                            ))
                      ),
                      
                      conditionalPanel(
                        condition = "input.checkbox ==true",
                        wellPanel(
                          radioButtons("radio1", label = h5("Choose False Discovery Control"),
                                       choices = list("Benjamini & Hochberg (BH)" = "BH",
                                                      "Benjamini & Yekutieli (BY)" = "BY",
                                                      "Holm" = "holm",
                                                      "Hochberg" = "hochberg",
                                                      "Hommel" = "hommel",
                                                      "Bonferroni" = "bonferroni", 
                                                      "FDR" = "fdr",
                                                      "None" = "none"),
                                       selected = NULL))
                      ),
                        
                      # uiOutput('AnnotationFile'),
                      # uiOutput('annotation_range'),
                      # uiOutput('hide'),
                      # uiOutput('FDRAffymetrix'),
                      # uiOutput('FDRCodlink'),
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
                          href="https://www.fda.gov/downloads/scienceresearch/bioinformaticstools/microarrayqualitycontrolproject/ucm132970.pdf"),
                        "and",
                        a("Illumina",
                          href= "https://www.illumina.com/techniques/microarrays/gene-expression-arrays.html"),
                        "gene expression data ."),
                      p("Each of the three data platform type has three functionalities, which cover meta-analysis of
                         processed data, processing raw data to perform meta-analysis
                         and finally differential expression"),
                      br(),
                      br(),
                      
                      # img(src = "bigorb.png", height = 100, width = 100),
                      br(),
                      h5("Powered by:", 
                         span("HLS and team", colour="Yellow"))),
                    
                    mainPanel(
                      tabsetPanel(id = "MamgedTabs",
                        # tabPanel("File-list",dataTableOutput("file")),
                        tabPanel("Source-data", dataTableOutput("sourced")),
                        tabPanel("Annotation-data",dataTableOutput("annotation")),
                        tabPanel("Summary",dataTableOutput("final")),
                        tabPanel("Supplementary file", dataTableOutput("full"))
                      )
                      
                      # uiOutput("alltabs")
                    )
                    )))