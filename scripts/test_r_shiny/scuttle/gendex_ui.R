library(shiny)
library(tidyverse)
library(DT)

working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
setwd(working_dir)
options("max.print"=100)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

#### load server information
source(file.path(working_dir, "gendex_server.R"))

ui <- fluidPage(
    titlePanel("GenDex - Genomic Index visualization tool"),
    tabsetPanel(selected="Home",

                
        tabPanel(title = "Home",
                 mainPanel(width = 12,align = "center",
                    sidebarLayout(
                        sidebarPanel(
                            fileInput("probeset_txt", "Choose probeset.txt File", accept = ".txt"),
                            checkboxInput("header", "Header", TRUE)
                        ),
                        mainPanel(
                            # tableOutput("contents")
                            DT::dataTableOutput("contents")
                        )
                    )
                           # selectInput("season_year","Select Season",choices=unique(sort(matches$season,
                           #                                                               decreasing=TRUE)), selected = 2019),
                           # submitButton("Go"),
                           # tags$h3("Players table"),
                           # div(style = "border:1px black solid;width:50%",tableOutput("player_table"))
        )),
        
        
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
                    # splitLayout(
                    sidebarLayout(
                        sidebarPanel(width = "2",
                            sliderInput("undoSD",
                                        "Join segments if their difference is greater than this:",
                                        min = 0.1,
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
                            # submitButton("Go"),
                            
                        ),
                        mainPanel(width = "10",
                            splitLayout(
                                verticalLayout(
                                    h2("Segments table"),
                                    DT::dataTableOutput("segTable")
                                ),
                                verticalLayout(
                                    h2("Genes table"),
                                    DT::dataTableOutput("geneTable")
                                ),
                                textOutput("debug")
                            )
                        )
                        ## segs table,
                        ## genes table
                    )
                )
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
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
