
################################ interactive datatable + functional scrollbar in DASHBOARD
# library(shiny)
# library(shinydashboard)
# 
# ## ui.R ##
# 
# ui <- dashboardPage(
# 
#     dashboardHeader(title = "TEST"),
# 
#     dashboardSidebar(
#         # to add scrollbar
#         tags$head(
#             tags$style(HTML("
#                       .sidebar { height: 90vh; overflow-y: auto; }
#                       " )
#             )
#         ),
# 
#         
#         
#         p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), 
#         
#         sidebarMenu(
#             menuItem("See data", tabName = "db")
#         ),
#         radioButtons("rb1", label = "Select data",
#                      choices = list("IRIS" = "iris", "CARS" = "cars"),
#                      selected = "iris")
#     ),
# 
#     dashboardBody(
#         tabItems(
#             tabItem(tabName = "db",
#                     h4("Show selected dataset"),
#                     fluidRow(DT::dataTableOutput('tbl2'))
#             )
#         )
#     )
# )
# 
# 
# ## server.R ##
# library(shiny)
# 
# server <- function(input, output, session) {
# 
#     output$value <- renderPrint({ input$rb1 })
# 
#     data <- reactive({
#         switch(input$rb1,
#                "iris" = iris,
#                cars)
#     })
# 
#     action <- dataTableAjax(session, cars)
#     widget <- datatable(cars,
#                         class = 'display cell-border compact',
#                         filter = 'top',
#                         # server = TRUE,
#                         options = list(ajax = list(url = action))
#     )
# 
#     output$tbl2 <- DT::renderDataTable({
#         DT::datatable(data())
#     })
# }
# 
# shinyApp(ui, server)
################################ interactive datatable








################################ interactive datatable + functional scrollbar in BASE R-SHINY
# library(shiny)

# ## ui.R ##
# ui <- fluidPage(
#     titlePanel("TEST"),
#     verticalLayout(
#         splitLayout(cellWidths = c("80%", "20%"),
#             plotOutput("testPlot"),
#             radioButtons("plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),
#         ),
#         sidebarLayout(
#             sidebarPanel(
#                 # to fill sideBar
#                 p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"),
#                 radioButtons("rb1", label = "Select data",
#                             choices = list("IRIS" = "iris", "CARS" = "cars"),
#                             selected = "iris")
#             ),
#             mainPanel(
#                 DT::dataTableOutput('tbl2')
#             )
#         )
#     )
# )
# ## server.R ##
# server <- function(input, output, session) {
#     output$testPlot = renderPlot({plot(c(8,6,5,3,9,5,2))})
#     output$test_text <- renderPrint({ "sample text" })
#     data <- reactive({
#         switch(input$rb1,
#             "iris" = iris,
#             cars)
#     })
#     action <- dataTableAjax(session, cars)
#     widget <- datatable(cars,
#                         class = 'display cell-border compact',
#                         filter = 'top',
#                         # server = TRUE,
#                         options = list(ajax = list(url = action))
#     )
#     output$tbl2 <- DT::renderDataTable({
#         DT::datatable(data())
#     })
# }
# shinyApp(ui, server)
################################ interactive datatable + functional scrollbar in BASE R-SHINY

################################ find why I don't have a scrollbar -> found
# library(shiny)
# ui <- fluidPage(
#     titlePanel("TEST"),
#     verticalLayout(
#         splitLayout(cellWidths = c("80%", "20%"),
#             plotOutput("profilePlot"),
#             radioButtons("plotChoice","Profile display",choices = c("Call probability"="proba","segmentation results"="profile")),
#         ),
#         sidebarLayout(
#             sidebarPanel(
#                         h3("Parameters"),
#                         wellPanel(sliderInput("undoSD",
#                                                 "Undo splits",
#                                                 min = 0.01,
#                                                 max = 10,
#                                                 value = 3
#                         )),
#                         p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"), p("text"),
#             ),
#             mainPanel(
#                 splitLayout(
#                     verticalLayout(
#                         h2("Segments table"),
#                         DT::dataTableOutput("segTable"),
#                         downloadButton("download_segTable", "Download")
#                     ),
#                     verticalLayout(
#                         h2("Genes table"),
#                         DT::dataTableOutput("geneTable"),
#                         downloadButton("download_genesTable", "Download"),
#                         textOutput("GItext")
#                     ),
#                 )
#             )
#         )
#     )
# )

# server <- function(input, output) {
#     output$profilePlot = renderPlot({plot(c(2,3,6,5,98,56,234))})

# }
# shinyApp(ui, server)





############################### switch between table and plot 
# library(shiny)
# library(shinyjs)
# library(DT)
# library(data.table)
# 
# runApp(list(
#     ui = fluidPage(
#         useShinyjs(),
#         wellPanel(
#             radioButtons("visuBtn", NULL, choices = c(Table = "table", Plot = "plot"))
#         ),
#         wellPanel(
#             # dataTableOutput("mytable"),
#             plotOutput("mytable"),
#             plotOutput("myplot")
#         )
#     ),
#     
#     server = function(input, output, session){
#         dfconc <- data.table(time = c(1,2,3,4,5), concentration = c(0.1, 0.4, 0.5, 0.7, 0.8))
# 
#         # output$mytable <- renderDataTable(
#         output$mytable <- renderPlot(
#             # dfconc, options = list(paging = FALSE, searching = FALSE)
#             plot(c(9,8,5,6))
#         )
# 
#         output$myplot <- renderPlot({
#             plot(dfconc$time, dfconc$concentration, xlab = "Time", ylab = "Concentration")
#         })
# 
#         observeEvent(input$visuBtn,{
#             req(input$visuBtn)
#             if(input$visuBtn == "plot"){
#                 hide("mytable")
#                 show("myplot")
#             }else{
#                 hide("myplot")
#                 show("mytable")
#             }
#         })
#     }
# ))


############################### box CSS
# shinyApp(
#     ui = fluidPage(
#         titlePanel("Some rows with border"),
#         tags$style(HTML("
#       [id*=flower] {
#           border: 4px double red;
#       }
#       #second {
#           border: 2px dashed blue;
#       }
#     ")),
#     fluidRow(id = "firstflower",
#              numericInput("n", "n", 1)
#     ),
#     fluidRow(
#         tags$br()
#     ),
#     fluidRow(id = "secondflower",
#              plotOutput("plot")
#     )
#     ),
#     server = function(input, output, session) {
#         output$plot <- renderPlot( plot(head(cars, input$n)) )
#     }
# )

library(ggplot2)

data <- data.frame(row.names=paste0('row',1:20))
data$x_baseline <- sample(1:20)
data$x_followup <- sample(1:20)
data$y_baseline <- sample(1:20)
data$y_followup <- sample(1:20)


data$lineColor <- factor(ifelse(data$x_followup - data$x_baseline < 0,'increase','decrease'))
# data$lineColor <- factor(ifelse(data$x_followup - data$x_baseline < 0,'yellow','orange'))
segColors = c("yellow", "orange")
names(segColors) = c('increase','decrease')

ggplot(data) + 
geom_point(aes (x = data$x_baseline , y= data$y_baseline), color = "red") + 
geom_point (aes (x = data$x_followup, y = data$y_followup), color = "blue") + 
xlab("X") + 
ylab ("Y") +
geom_segment(data = data , aes(x=data$x_baseline, xend = data$x_followup, y=data$y_baseline, yend = data$y_followup, color=data$lineColor)) + 
scale_color_manual(values=segColors)



