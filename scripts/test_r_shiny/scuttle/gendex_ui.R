library(shiny)
library(shinyBS)
library(shinyjs)
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

ui <- fluidPage(useShinyjs(),
    titlePanel("GenDex - Genomic Index visualization tool"),
    add_busy_spinner(spin = "fading-circle"),
    tabsetPanel(selected="Home",

                
        tabPanel(title = "Home",
            verticalLayout(
                sidebarLayout(
                    sidebarPanel( 
                        fileInput("probeset_txt", "Choose probeset.txt file", accept = ".txt"),
                        checkboxInput("probeSetHeader", "Header", TRUE),
                        fileInput("OSCHP", "Choose .OSCHP file", accept = ".OSCHP"),
                        textInput("prefix", "insert prefix for sample"),
                        
                        ######## App description
                        h3("Quick description"),
                        p("This application aims to help visualize and export the results related to the calculation of the Genomic Index on OncoScan CNV data by three R packages: CGHcall, ASCAT and rCGH."),
                        p("To use it, first load a probeset.txt (CGHcall & rCGH) or .OSCHP file (ASCAT)."),
                        p("Then run CGHcall, ASCAT or rCGH pipelines in their respective panes to compute Genomic Index."),
                        p("Results can be viewed there, along with a recap in the summary pane."),
                        p("The Parameters pane contains explanations on the parameters used by each pipeline."),
                    ),
                    mainPanel(
                        h3("Probeset loaded:"),
                        DT::dataTableOutput("probesetContent"),
                    )
                ),
                
            )
        ),
        
        
        tabPanel(title = "CGHcall",
            mainPanel(width = 12,align = "center",
                verticalLayout(
                    # div(style = "float:left;width:90%;",plotOutput("profilePlot"),plotOutput("probaPlot")),
                    # div(style = "float:right;width:10%;position:relative;top:70px;",
                    #     radioButtons("plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile"),
                    #     selected = "profile")
                    # ),
                    splitLayout(cellWidths = c("80%", "20%"),
                        plotOutput("profilePlot"),
                        radioButtons("plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),
                    ),

                    sidebarLayout(
                        sidebarPanel( width = "2",
                            h3("Parameters"),
                            bsTooltip("undoSD", "This represents the distance under which two segments are fused together.", placement = "bottom", trigger = "hover"),
                            wellPanel(sliderInput("undoSD",
                                                "Undo splits",
                                                min = 0.01,
                                                max = 10,
                                                value = 3
                            )),
                            wellPanel(selectInput("prior",
                                                "How call probabilities should be calculated:", 
                                                choices = c('Per chromosome arm'='not all', 'On whole genome'='all')
                            )),
                            wellPanel(checkboxInput("correctCell", 
                                                "Correct data using proportion of tumoral cells", 
                                                value=TRUE,
                            )),
                            wellPanel(sliderInput("cellularity",
                                                "Proportion of tumoral cells:",
                                                min = 0,
                                                max = 100,
                                                value = 100
                            )),
                            wellPanel(sliderInput("minSegLenForFit",
                                                "Minimum length of the segment (in Mb) to be used for fitting the model:",
                                                min = 0,
                                                max = 1000,
                                                value = 1
                            )),
                            wellPanel(actionButton("go", "Run"))
                        ),
                        mainPanel(width = "10",
                            textOutput("GItext"),
                            splitLayout(
                                verticalLayout(
                                    h2("Segments table"),
                                    DT::dataTableOutput("segTable"),
                                    downloadButton("download_segTable", "Download")
                                ),
                                verticalLayout(
                                    h2("Genes table"),
                                    DT::dataTableOutput("geneTable"),
                                    downloadButton("download_genesTable", "Download"),
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
            mainPanel(width = 12, # align = "center",
                h1("CGHcall"),
                h3("parameter 1: UndoSD"),
                p("After segmenting, the algorithm undoes splits between two consecutive segments if these are less than undoSD (3 by default) sd apart, where sd is the standard deviation of the whole profile. The higher undoSD is, the higher the number of retained splits."),
                h3("parameter 2: prior probability computing"),
                h3("parameter 3: Correct values using cellularity"),
                h3("parameter 4: cellularity value"),
                h3("parameter 5: segment length to include in call computing"),
                
                h1("ASCAT"),
                h3("..."),
                
                h1("rCGH"),
                h3("..."),
            )),
        tabPanel(title = "Summary",
            mainPanel(width = 12, # align = "center",
                h2("GI table"),
                DT::dataTableOutput("GI_table_summary"),
                downloadButton("download_GI_table", "Download"),
                
            )),
        ))





shinyApp(ui=ui, server=server)
