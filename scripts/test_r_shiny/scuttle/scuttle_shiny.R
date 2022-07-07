
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

# ### apply multiple color_scale on same graph
# library(ggplot2)

# data <- data.frame(row.names=paste0('row',1:20))
# data$x_baseline <- sample(1:20)
# data$x_followup <- sample(1:20)
# data$y_baseline <- sample(1:20)
# data$y_followup <- sample(1:20)
# data$lineColor <- factor(ifelse(data$x_followup - data$x_baseline < 0,'increase','decrease'))

# ##colors for segments
# segColors = c("yellow", "orange")
# names(segColors) = c('increase','decrease')
# ## colors for "baseline" points 
# b_pointsColors = c("light green", "dark green")
# names(b_pointsColors) = c('increase','decrease')
# ## colors for "followup" points 
# f_pointsColors = c("light blue", "dark blue")
# names(f_pointsColors) = c('increase','decrease')

# ggplot(data) + 
# # first color scale
# geom_point(aes (x = x_baseline , y= y_baseline, color = lineColor)) +
# scale_color_manual(values=b_pointsColors) +
# # new color scale
# ggnewscale::new_scale_color()+
# geom_point (aes (x = x_followup, y = y_followup, color = lineColor)) + 
# scale_color_manual(values=f_pointsColors) +
# # new color scale
# ggnewscale::new_scale_color()+
# geom_segment(data = data , aes(x=x_baseline, xend = x_followup, y=y_baseline, yend = y_followup, color=lineColor)) + 
# scale_color_manual(values=segColors)

# ####### hide/show an element using renderUI and req()
# library(shiny)
# library(shinydashboard)

# ui <- dashboardPage(
    
#     dashboardHeader(),
#     dashboardSidebar(
#         selectInput(inputId = "select", 
#                     label = "please select an option", 
#                     choices = LETTERS[1:3]),
#         uiOutput("conditional_comment")
#     ),
#     dashboardBody(
#         uiOutput("selection_text"),
#         uiOutput("comment_text")
#     )
# )

# server <- function(input, output) {
    
#     output$selection_text <- renderUI({
#         paste("The selected option is", input$select)
#     })
    
#     output$conditional_comment <- renderUI({
#         req(input$select == "B")
#         textAreaInput(inputId = "comment", 
#                       label = "please add a comment", 
#                       placeholder = "write comment here")
#     })
    
#     output$comment_text <- renderText({
#         input$comment
#     })
# }

# shinyApp(ui = ui, server = server)


# ##### hide/show a an element using conditionalPanel()
# library(shiny)

# ## Module code for 'selectorUI' and 'selector'
# selectorUI <- function(id) {
#   ns <- NS(id)
#   selectizeInput(inputId = ns('select'),
#                  label = 'Make a choice:',
#                  choices = c('Option one', 'Option two'))
# }

# ## Main app
# ui <- shinyUI(fluidPage(
#   selectorUI('id1'),
#   conditionalPanel(condition = "input['id1-select'] == 'Option one'",
#                    p('Option one is selected.'))
# ))

# server <- shinyServer(function(input, output, session) {

# })

# shinyApp(ui = ui, server = server)




# ##### 
# library(shiny)
# ## Module code for 'selectorUI' and 'selector'
# selectorUI <- function(id) {
#   ns <- NS(id)
#   selectizeInput(inputId = ns('select'),
#                  label = 'Make a choice:',
#                  choices = c('Option one', 'Option two'))
# }

# selector <- function(input, output, session) {
#   reactive(input$select)
# }

# ## Main app
# ui <- shinyUI(fluidPage(
#   selectorUI('id1'),
#   uiOutput("dynamic1")
# ))

# server <- shinyServer(function(input, output, session) {


#   output$dynamic1 <- renderUI({
#     condition1 <- callModule(selector, 'id1') # or just callModule(selector, 'id1')()
#     if (condition1() == 'Option one') return(p('Option one is selected.'))
#   })

# })

# shinyApp(ui = ui, server = server)











n <- 200
library("shiny")

# Define the UI
ui <- bootstrapPage(
  numericInput('n', 'Number of obs', n),
  conditionalPanel(condition = "output.cond == true", # here use the condition defined in the server
                   plotOutput('plot') ),
  HTML("Bottom")
)

# Define the server code
server <- function(input, output, session) {
  output$plot <- renderPlot({
    if (input$n > 50) hist(runif(input$n)) else return(NULL)
  })
  # create a condition you use in the ui
  output$cond <- reactive({
    input$n > 50
  })
  outputOptions(output, "cond", suspendWhenHidden = FALSE)
}

# Return a Shiny app object
shinyApp(ui = ui, server = server) 