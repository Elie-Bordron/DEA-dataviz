library(shiny)
library(shinybusy)
library(tidyverse)
library(DT)

working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
results_dir = file.path(working_dir_shiny, "gendex_results")
if(!dir.exists(results_dir))dir.create(results_dir)
setwd(working_dir_shiny)
options("max.print"=100)
rstudioapi::filesPaneNavigate(working_dir_shiny)

#### load server information
source(file.path(working_dir_shiny, "gendex_server.R"))

ui <- fluidPage(
    titlePanel("GenDex - Genomic Index visualization tool"),
    add_busy_spinner(spin = "fading-circle"),
    tabsetPanel(selected="Home",

        tabPanel(title = "Home",
            sidebarLayout(
                sidebarPanel(
                    fileInput("probeset_txt", "Choose probeset.txt file", accept = ".txt"),
                    checkboxInput("probeSetHeader", "Header", TRUE),
                    fileInput("OSCHP", "Choose .OSCHP file", accept = ".OSCHP"),
                    textInput("prefix", "insert prefix for sample")
                ),
                mainPanel(
                    # tableOutput("contents")
                    DT::dataTableOutput("probesetContent")
                )
            )
        ),
        
        
        tabPanel(title = "CGHcall",
            mainPanel(width = 12,align = "center",
                verticalLayout(
                    fillPage(
                        div(style = "float:left;width:90%;",plotOutput("testPlot")),
                        div(style = "float:right;width:10%;position:relative;top:70px;",
                            radioButtons("profile_plot_choice","Profile display",choices = c("Call probability","segmentation results"),
                            selected = "segmentation results")
                        ),
                    ),
                    sidebarLayout(
                        sidebarPanel(width = "2",
                            sliderInput("undoSD",
                                        "This represents the distance under which two segments are fused together:",
                                        min = 0.01,
                                        max = 10,
                                        value = 3
                            ),
                            selectInput("prior",
                                        "How call probabilities should be calculated:", 
                                        choices = c('Per chromosome arm'='not all', 'On whole genome'='all')
                            ),
                            checkboxInput ("correctCell", 
                                           "Correct data using proportion of tumoral cells", 
                                           value=TRUE,
                            ),
                            sliderInput("cellularity",
                                         "Proportion of tumoral cells:",
                                        min = 0,
                                        max = 100,
                                        value = 100
                            ),
                            sliderInput("minSegLenForFit",
                                        "Minimum length of the segment (in Mb) to be used for fitting the model:",
                                        min = 0,
                                        max = 1000,
                                        value = 1
                            ),
                            actionButton("go", "Run")
                            
                        ),
                        mainPanel(width = "10",
                            splitLayout(
                                verticalLayout(
                                    h2("Segments table"),
                                    DT::dataTableOutput("segTable"),
                                    downloadButton("download_segTable", "Download")
                                ),
                                verticalLayout(
                                    h2("Genes table"),
                                    DT::dataTableOutput("geneTable"),
                                    downloadButton("download_genesTable", "Download")
                                ),
                                # textOutput("debug")
                            )
                        )
                    )
                )




            )),
        tabPanel(title = "ASCAT",
            mainPanel(width = 12,align = "center",
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
            )),
        tabPanel(title = "rCGH",
            mainPanel(width = 12,align = "center",
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
            )),
        tabPanel(title = "Parameters",
            mainPanel(width = 12,align = "center",
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
            )),
        tabPanel(title = "Summary",
            mainPanel(width = 12,align = "center",
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
            )),
        ))





shinyApp(ui=ui, server=server)
