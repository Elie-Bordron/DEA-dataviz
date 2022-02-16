library("shiny")

testData = data.frame(c(72, 36, 54, 85, 66), c("a", "b", "c", "d", "e"))
headerPanel()

#define UI
ui = fluidPage(
    
)

#Define server logic required


server = function(input, output){
    
}

#run the app
shinyApp(ui, server)