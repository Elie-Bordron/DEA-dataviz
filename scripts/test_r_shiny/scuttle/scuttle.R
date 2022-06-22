library(shiny)
# #
# # This is a Shiny web application. You can run the application by clicking
# # the 'Run App' button above.
# #
# # Find out more about building applications with Shiny here:
# #
# #    http://shiny.rstudio.com/
# #
# 
# library(shiny)
# working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts"
# setwd(working_dir)
# options("max.print"=100)
# ## open working directory in Files tab
# rstudioapi::filesPaneNavigate(working_dir)
# 
# # Define UI for application that draws a histogram
# ui <- fluidPage(
#     
#     # Application title
#     titlePanel("Scuttle first shiny app"),
#     
#     # Sidebar with a slider input for number of bins 
#     sidebarLayout(
#         sidebarPanel(
#             sliderInput("bins",
#                         "Number of bins:",
#                         min = 1,
#                         max = 50,
#                         value = 30
#             ),
#             sliderInput("xlim",
#                         "limits of x axis:",
#                         min = 1,
#                         max = 100,
#                         value = 10
#             ),
#             sliderInput("hz_abline",
#                         "abline:",
#                         min = 1,
#                         max = 25,
#                         value = 2
#             )
#         ),
#         
#         # Show a plot of the generated distribution
#         mainPanel(
#             plotOutput("distPlot")
#         )
#     )
# )
# 
# # Define server logic required to draw a histogram
# server <- function(input, output) {
#     
#     output$distPlot <- renderPlot({
#         x    <- faithful[, 2]
#         # x = rnorm(100,10,input$xlim)
#         # generate bins based on input$bins from ui.R
#         bins <- seq(min(x), max(x), length.out = input$bins + 1)
#         
#         # draw the histogram with the specified number of bins
#         hist(x, breaks = bins, col = 'darkgray', border = 'white', xlim = c(-input$xlim+mean(x), input$xlim+mean(x)))
#         abline(h=input$hz_abline)
#     })
# }
# 
# # Run the application 
# shinyApp(ui = ui, server = server)










# 
# ## SERVER.R
# 
# #Initialization
# library(datasets)
# save(iris,file="iris.Rdata")
# save(mtcars,file="m.Rdata")
# 
# server = (function(input, output) {
#     infile <- reactive({
#         infile <- input$datafile
#         if (is.null(infile)) {
#             # User has not uploaded a file yet
#             return(NULL)
#         }
#         objectsLoaded <- load(input$datafile$name) 
#         # the above returns a char vector with names of objects loaded
#         df <- eval(parse(text=objectsLoaded[1])) 
#         # the above finds the first object and returns it
#         return(df)
#     })
#     
#     myData <- reactive({
#         df<-infile()
#         if (is.null(df)) return(NULL)
#         return(df)
#     })
#     output$value1 <- renderPrint({
#         names(iris)
#     })
#     output$value2 <- renderPrint({
#         names(myData())
#     })
#     load("iris.Rdata")   ## data loaded for testing 
# })
# 
# ## UI.R
# ui = (fluidPage(
#     fileInput("datafile", label = h3("File input")),
#     fluidRow(column(4, verbatimTextOutput("value1"))),
#     fluidRow(column(4, verbatimTextOutput("value2")))
# ))
# 
# 


ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose CSV File", accept = ".txt"),
            checkboxInput("header", "Header", TRUE)
        ),
        mainPanel(
            tableOutput("contents")
        )
    )
)

server <- function(input, output) {
    output$contents <- renderTable({
        file <- input$file1
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "txt", "Please upload a txt file"))
        
        read.table(file$datapath, header = input$header, sep="\t")
    })
}

shinyApp(ui, server)


shinyApp(ui=ui, server=server)