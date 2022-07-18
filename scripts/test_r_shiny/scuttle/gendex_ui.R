library(shiny)
library(shinyBS)
# library(shinyjs)
library(shinybusy)
library(tidyverse)
library(DT)
library(ggplot2)




if(TRUE){ # set this to TRUE on bergo PC 
    print("bergo path")
    working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
} else { # this is used on my own PC
    #S
    # working_dir_shiny = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    # GI_scripts_dir = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
    #C
    working_dir_shiny = "C:/Users/warew/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    GI_scripts_dir = "C:/Users/warew/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
} 


 



#### load server information
source(file.path(working_dir_shiny, "gendex_server.R"))
results_dir = file.path(working_dir_shiny, "gendex_results")
if(!dir.exists(results_dir))dir.create(results_dir)
setwd(working_dir_shiny)
# options("max.print"=10**9)
options("max.print"=10**2)
rstudioapi::filesPaneNavigate(working_dir_shiny)


ui <- fluidPage( #useShinyjs(),
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
                        p("To use it, first load a probeset.txt for CGHcall or rCGH, or load a .OSCHP file for ASCAT."),
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
                    # numericInput("n", "n", 50),
                    conditionalPanel(
                        condition = "output.CGHcall_probesetLoaded == false",  # js
                        wellPanel(
                            h1(span(textOutput("test_warnPanel"), style = 'color:green; font-weight: bold;')),
                        ),
                    ),

                    splitLayout( cellWidths = c("80%", "20%"),
                        verticalLayout(
                            plotOutput("CGHcall_profilePlot"),
                            plotOutput("CGHcall_allDiffPlot"),
                            shinyBS::bsTooltip("CGHcall_profilePlot", "This plot is obtained from the segments calculated by rCGH pipeline and the Log2Ratio values provided in input file probeset.txt", placement = "right", trigger = "hover"),
                            shinyBS::bsTooltip("CGHcall_allDiffPlot", "This plot is obtained from the Allelic Difference values provided in input file probeset.txt", placement = "right", trigger = "hover"),

                        ),
                        radioButtons("CGHcall_plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),

                    ),

                    sidebarLayout(
                        sidebarPanel( width = "2",
                            h3("Parameters"),
                            shinyBS::bsTooltip("undoSD", "This represents the distance under which two segments are fused together.", placement = "bottom", trigger = "hover"),
                            wellPanel(sliderInput("undoSD",
                                                "Undo splits",
                                                min = 0.01,
                                                max = 10,
                                                value = 3)
                            ),
                            shinyBS::bsTooltip("prior", "How call probabilities should be calculated.", placement = "bottom", trigger = "hover"
                            ),
                            wellPanel(selectInput("prior",
                                                "Computing of call probabilities", 
                                                choices = c('Per chromosome arm'='not all', 'On whole genome'='all'))
                            ),
                            shinyBS::bsTooltip("correctCell", "Whether call computing should take the proportion of tumoral cells into account", placement = "bottom", trigger = "hover"
                            ),
                            wellPanel(
                                checkboxInput("correctCell", 
                                                "Correct data using proportion of tumoral cells", 
                                                value=TRUE,
                                ),
                                sliderInput("cellularity",
                                                "Proportion of tumoral cells:",
                                                min = 0,
                                                max = 100,
                                                value = 100
                                )
                            ),
                            wellPanel(sliderInput("minSegLenForFit",
                                                "Minimum length of the segment (in Mb) to be used for fitting the model:",
                                                min = 0,
                                                max = 10,
                                                value = 0.5)
                            ),
                            wellPanel(actionButton("go", "Run")),
                            # wellPanel(actionButton("goPlot", "Run plot")),
                        ),
                        mainPanel(width = "10",
                            verticalLayout(
                                # h3(textOutput("CGHcall_GItext", container=pre)),
                                
                                conditionalPanel(
                                    condition = "output.CGHcall_probesetLoaded == true",  # js
                                    wellPanel(h3(textOutput("CGHcall_GItext")))
                                ),
                                splitLayout(
                                    # verticalLayout(
                                    wellPanel(
                                        h2("Segments table"),
                                        selectInput("selectCN",
                                                "", 
                                                choices = c('Altered (Loss & Gain)'='CN!=2', 'Loss (CN<2)'='CN<2', 'Gain (CN>2)'='CN>2', 'All'='is.numeric(CN)')),
                                        DT::dataTableOutput("CGHcall_segTable"),
                                        downloadButton("CGHcall_download_segTable", "Download")
                                    ),
                                    ### Genes table
                                    # wellPanel(
                                    #     h2("Genes table"),
                                    #     DT::dataTableOutput("geneTable"),
                                    #     downloadButton("download_genesTable", "Download"),
                                    # ),
                                )
                            )
                        )
                    )
                )
            )
        ),

        
        tabPanel(title = "rCGH",
            mainPanel(width = 12,align = "center",
                verticalLayout(
                    conditionalPanel(
                        condition = "output.rCGH_probesetLoaded == false",  # js
                        wellPanel(
                            h1(span(textOutput("rCGH_test_warnPanel"), style = 'color:green; font-weight: bold;')),
                        ),
                    ),

                    splitLayout( cellWidths = c("80%", "20%"),
                        verticalLayout(
                            plotOutput("rCGH_profilePlot"),
                            plotOutput("rCGH_allDiffPlot"),
                            shinyBS::bsTooltip("rCGH_profilePlot", "This plot is obtained from the segments calculated by rCGH pipeline and the Log2Ratio values provided in input file probeset.txt", placement = "right", trigger = "hover"),
                            shinyBS::bsTooltip("rCGH_allDiffPlot", "This plot is obtained from the Allelic Difference values provided in input file probeset.txt", placement = "right", trigger = "hover"),
                        ),
                        # radioButtons("rCGH_plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),

                    ),

                    sidebarLayout(
                        sidebarPanel( width = "2",
                            h3("Parameters"),

                            # ),
                            wellPanel(
                                checkboxInput("useCustomUndoSD", "Override estimated value",value=FALSE,
                                ),
                                sliderInput("undoSD_rCGH",
                                                "Undo splits",
                                                min = 0.01,
                                                max = 10,
                                                value = 3
                                )
                            ),
                            wellPanel(
                                textInput("a_priori_ploidy", "\"A priori\" ploidy", value = "2")
                            ),
                            wellPanel(
                                checkboxInput("scale_adjSignal", "Log2Ratios rescaling",value=TRUE)
                            ),
                            wellPanel(
                                checkboxInput("suppOutliers_adjSignal", "remove outliers",value=TRUE)
                            ),
                            wellPanel(
                                checkboxInput("smooth_segment", "Log2Ratios smoothing",value=TRUE)
                            ),
                            wellPanel(
                                sliderInput("minLen",
                                                "minimal length for a segment (Kb)",
                                                min = 0,
                                                max = 100,
                                                value = 10
                                )
                            ),
                            shinyBS::bsTooltip("useCustomUndoSD", "If this box is not checked, this value will be estimated from the Log2Ratio values", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("a_priori_ploidy", " A priori ploidy value, when known, to adjust the estimation of copy numbers.", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("scale_adjSignal", "Whether Log2Ratios should be rescaled prior to analysis. This step increases significantly the processing time.", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("suppOutliers_adjSignal", "If checked, outliers are removed using an Expectation-Maximization algorithm ", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("smooth_segment", "Should the LRR be smoothed before being segmented", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("undoSD_rCGH", "This represents the distance between means under which two segments are fused together. This distance is equal to this value multiplied by SD, where SD is the standard deviation of the Log2Ratio values.", placement = "bottom", trigger = "hover"),
                            shinyBS::bsTooltip("minLen", " Segments shorter than this value will be merged to the closest adjacent segment", placement = "bottom", trigger = "hover"),

                            wellPanel(actionButton("go_rCGH", "Run")),
                        ),
                        mainPanel(width = "10",
                            verticalLayout(
                                # h3(textOutput("rCGH_GItext", container=pre)),
                                
                                conditionalPanel(
                                    condition = "output.rCGH_probesetLoaded == true",  # js
                                    wellPanel(h3(textOutput("rCGH_GItext")))
                                ),
                                splitLayout(
                                    # verticalLayout(
                                    wellPanel(
                                        h2("Segments table"),
                                        selectInput("selectCN_rCGH",
                                                "", 
                                                choices = c('Altered (Loss & Gain)'='CN!=2', 'Loss (CN<2)'='CN<2', 'Gain (CN>2)'='CN>2', 'All'='is.numeric(CN)')),
                                        DT::dataTableOutput("rCGH_segTable"),
                                        downloadButton("rCGH_download_segTable", "Download")
                                    ),
                                    ### Genes table
                                    # wellPanel(
                                    #     h2("Genes table"),
                                    #     DT::dataTableOutput("geneTable"),
                                    #     downloadButton("download_genesTable", "Download"),
                                    # ),
                                )
                            )
                        )
                    )
                )
            )
        ),




        tabPanel(title = "ASCAT",
                    mainPanel(width = 12,align = "center",
                verticalLayout(
                    conditionalPanel(
                        # condition = "output.rCGH_probesetLoaded == false",  # js
                        condition = "output.rCGH_probesetLoaded == false",  # js
                        wellPanel(
                            # h1(span(textOutput("rCGH_test_warnPanel"), style = 'color:green; font-weight: bold;')),
                            h1(span("Please load a probeset.txt from Home pane", style = 'color:green; font-weight: bold;')),
                        ),
                    ),

                    splitLayout( cellWidths = c("80%", "20%"),
                        verticalLayout(
                            # plotOutput("rCGH_profilePlot"),
                            # plotOutput("rCGH_allDiffPlot"),
                            # shinyBS::bsTooltip("rCGH_profilePlot", "This plot is obtained from the segments calculated by rCGH pipeline and the Log2Ratio values provided in input file probeset.txt", placement = "right", trigger = "hover"),
                            # shinyBS::bsTooltip("rCGH_allDiffPlot", "This plot is obtained from the Allelic Difference values provided in input file probeset.txt", placement = "right", trigger = "hover"),
                        ),
                        # radioButtons("rCGH_plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),

                    ),

                    sidebarLayout(
                        sidebarPanel( width = "2",
                            h3("Parameters"),

                            # # ),
                            # wellPanel(
                            #     checkboxInput("useCustomUndoSD", "Override estimated value",value=FALSE,
                            #     ),
                            #     sliderInput("undoSD_rCGH",
                            #                     "Undo splits",
                            #                     min = 0.01,
                            #                     max = 10,
                            #                     value = 3
                            #     )
                            # ),
                            # wellPanel(
                            #     textInput("a_priori_ploidy", "\"A priori\" ploidy", value = "2")
                            # ),
                            # wellPanel(
                            #     checkboxInput("scale_adjSignal", "Log2Ratios rescaling",value=TRUE)
                            # ),
                            # wellPanel(
                            #     checkboxInput("suppOutliers_adjSignal", "remove outliers",value=TRUE)
                            # ),
                            # wellPanel(
                            #     checkboxInput("smooth_segment", "Log2Ratios smoothing",value=TRUE)
                            # ),
                            # wellPanel(
                            #     sliderInput("minLen",
                            #                     "minimal length for a segment (Kb)",
                            #                     min = 0,
                            #                     max = 100,
                            #                     value = 10
                            #     )
                            # ),
                            # shinyBS::bsTooltip("useCustomUndoSD", "If this box is not checked, this value will be estimated from the Log2Ratio values", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("a_priori_ploidy", " A priori ploidy value, when known, to adjust the estimation of copy numbers.", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("scale_adjSignal", "Whether Log2Ratios should be rescaled prior to analysis. This step increases significantly the processing time.", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("suppOutliers_adjSignal", "If checked, outliers are removed using an Expectation-Maximization algorithm ", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("smooth_segment", "Should the LRR be smoothed before being segmented", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("undoSD_rCGH", "This represents the distance between means under which two segments are fused together. This distance is equal to this value multiplied by SD, where SD is the standard deviation of the Log2Ratio values.", placement = "bottom", trigger = "hover"),
                            # shinyBS::bsTooltip("minLen", " Segments shorter than this value will be merged to the closest adjacent segment", placement = "bottom", trigger = "hover"),

                            wellPanel(actionButton("go_ASCAT", "Run")),
                        ),
                        mainPanel(width = "10",
                            verticalLayout(
                                # h3(textOutput("rCGH_GItext", container=pre)),
                                
                                conditionalPanel(
                                    # condition = "output.rCGH_probesetLoaded == true",  # js
                                    condition = "true",  # js
                                    # wellPanel(h3(textOutput("rCGH_GItext")))
                                    wellPanel(h3())
                                ),
                                splitLayout(
                                    # verticalLayout(
                                    wellPanel(
                                        h2("Segments table"),
                                        # selectInput("selectCN_rCGH",
                                        #         "", 
                                        #         choices = c('Altered (Loss & Gain)'='CN!=2', 'Loss (CN<2)'='CN<2', 'Gain (CN>2)'='CN>2', 'All'='is.numeric(CN)')
                                        # ),
                                        DT::dataTableOutput("ASCAT_segTable"),
                                        downloadButton("ASCAT_download_segTable", "Download")
                                    ),
                                    ### Genes table
                                    # wellPanel(
                                    #     h2("Genes table"),
                                    #     DT::dataTableOutput("geneTable"),
                                    #     downloadButton("download_genesTable", "Download"),
                                    # ),
                                )
                            )
                        )
                    )
                )
            )

        ),

        tabPanel(title = "Summary",
            mainPanel(width = 12,
                verticalLayout(
                    h2("GI table"),
                    DT::dataTableOutput("GI_table_summary"),
                    downloadButton("download_GI_table", "Download"),
                ),
            )
        ), 
        
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
            )
        ),
        ))





shinyApp(ui=ui, server=server)






