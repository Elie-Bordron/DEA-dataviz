#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts"
setwd(working_dir)
options("max.print"=100)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Scuttle first shiny app"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30
            ),
            sliderInput("xlim",
                        "limits of x axis:",
                        min = 1,
                        max = 100,
                        value = 10
            ),
            sliderInput("hz_abline",
                        "abline:",
                        min = 1,
                        max = 25,
                        value = 2
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$distPlot <- renderPlot({
        x    <- faithful[, 2]
        # x = rnorm(100,10,input$xlim)
        # generate bins based on input$bins from ui.R
        bins <- seq(min(x), max(x), length.out = input$bins + 1)
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white', xlim = c(-input$xlim+mean(x), input$xlim+mean(x)))
        abline(h=input$hz_abline)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
