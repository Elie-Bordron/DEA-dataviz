library(shiny)
library(shinydashboard)

## ui.R ##

ui <- dashboardPage(
    
    dashboardHeader(title = "TEST reactive DT"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("See data", tabName = "db")
        ),
        radioButtons("rb1", label = "Select data", 
                     choices = list("IRIS" = "iris", "CARS" = "cars"),
                     selected = "iris")
    ),
    
    dashboardBody(
        tabItems(
            tabItem(tabName = "db",
                    h4("Show selected dataset"),
                    fluidRow(DT::dataTableOutput('tbl2'))
            )
        )
    )
)  



## server.R ##
library(shiny)

server <- function(input, output, session) {
    
    output$value <- renderPrint({ input$rb1 })
    
    data <- reactive({
        switch(input$rb1,
               "iris" = iris,
               cars)
    })
    
    action <- dataTableAjax(session, cars)
    widget <- datatable(cars, 
                        class = 'display cell-border compact',
                        filter = 'top',
                        # server = TRUE,
                        options = list(ajax = list(url = action))
    )
    
    output$tbl2 <- DT::renderDataTable({
        DT::datatable(data())
    })
}



shinyApp(ui, server)