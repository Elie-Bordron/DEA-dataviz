library(shiny)
library(tidyverse)

working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
setwd(working_dir)
options("max.print"=100)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

#### load server information
source(file.path(working_dir, "gendex_server.R"))

ui <- fluidPage(
    headerPanel("GenDex - Genomic Index visualization tool"),
    tabsetPanel(

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
                    splitLayout(
                        verticalLayout(
                            sliderInput("undoSD",
                                        "Join segments if their difference is greater than this:",
                                        min = 0.1,
                                        max = 10,
                                        value = 3
                            ),
                            selectInput("prior",
                                        "How call probabilities should be calculated:", 
                                        choices = c('On whole genome'='1','Per chromosome arm'='2')
                            ),
                            checkboxInput ("correctCell", 
                                           "Correct data using proportion of tumoral cells", 
                                           value=FALSE,
                            ),
                            sliderInput("cellularity",
                                        "Proportion of tumoral cells:",
                                        min = 0,
                                        max = 100,
                                        value = 100
                            ),
                            
                        ),
                        tableOutput("segTable"),
                        tableOutput("geneTable")
                        ## segs table,
                        ## genes table
                    )
                )
                      # tags$h3("Team Wins & Points"),
                      # div(style = "float:left;width:36%;",plotOutput("wins_bar_plot")),
                      # div(style = "float:right;width:64%;",plotOutput("points_bar_plot"))
            )),
        tabPanel(title = "Home",
                 mainPanel(width = 12,align = "center",
                           # selectInput("season_year","Select Season",choices=unique(sort(matches$season,
                           #                                                               decreasing=TRUE)), selected = 2019),
                           # submitButton("Go"),
                           # tags$h3("Players table"),
                           # div(style = "border:1px black solid;width:50%",tableOutput("player_table"))
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
