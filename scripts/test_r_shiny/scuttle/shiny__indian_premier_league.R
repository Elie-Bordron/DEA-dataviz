library(shiny)
library(tidyverse)

working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
dataDir = file.path(working_dir, "IPL_data")
setwd(working_dir)
options("max.print"=100)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)

# Loading Dataset-------------------------------------------------------
deliveries = read.csv(file.path(dataDir, "deliveries.csv"), stringsAsFactors = FALSE)
matches = read.csv(file.path(dataDir, "matches.csv"), stringsAsFactors = FALSE)
# Cleaning Dataset------------------------------------------------------
names(matches)[1] = "match_id"
IPL = dplyr::inner_join(matches,deliveries)


ui <- fluidPage(
    headerPanel("IPL - Indian Premier League"),
    tabsetPanel(
        tabPanel(title = "Season",
                 mainPanel(width = 12,align = "center",
                           selectInput("season_year","Select Season",choices=unique(sort(matches$season,
                                                                                         decreasing=TRUE)), selected = 2019),
                           submitButton("Go"),
                           tags$h3("Players table"),
                           div(style = "border:1px black solid;width:50%",tableOutput("player_table"))
                 )),
        tabPanel(
            title = "Team Wins & Points",
            mainPanel(width = 12,align = "center",
                      tags$h3("Team Wins & Points"),
                      div(style = "float:right;width:36%;",plotOutput("wins_bar_plot")),
                      div(style = "float:left;width:64%;",plotOutput("points_bar_plot"))
            )
        )))

server <- function(input, output) {
    matches_year = reactive({ matches %>% filter(season == input$season_year) })
    playoff = reactive({ nth(sort(matches_year()$match_id,decreasing = TRUE),4) })
    matches_played = reactive({ matches_year() %>% filter(match_id < playoff()) })
    t1 = reactive({ matches_played() %>% group_by(team1) %>% summarise(count = n()) })
    t2 = reactive({ matches_played() %>% group_by(team2) %>% summarise(count = n()) })
    wl = reactive({ matches_played() %>% filter(winner != "") %>% group_by(winner) %>% 
            summarise(no_of_wins = n()) })
    
    wl1=reactive({ matches_played() %>% group_by(winner) %>% summarise(no_of_wins=n()) })
    tied = reactive({ matches_played() %>% filter(winner == "") %>% select(team1,team2) })
    playertable = reactive({data.frame(Teams = t1()$team1,Played=t1()$count+t2()$count,
                                       Wins = wl()$no_of_wins,Points = wl()$no_of_wins*2)})
    
    
    output$player_table = renderTable({ playertable() })
    
    
    
    output$wins_bar_plot = renderPlot({ ggplot(wl1()[2:9,],aes(winner,no_of_wins,fill=winner))+
            geom_bar(stat = "identity")+ theme_classic()+xlab("Teams")+ ylab("Number Of Wins")+theme(
                axis.text.x=element_text(color="white"),legend.position = "none",axis.title=element_text(
                    size=14),plot.background=element_rect(colour="white"))+geom_text(aes(x=winner,(no_of_wins+0.6),
                                                                                         label = no_of_wins,size = 7)) })
    
    output$points_bar_plot = renderPlot({ ggplot(playertable(),aes(Teams,Points,fill=Teams))+
            geom_bar(stat = "identity",size=3)+theme_classic()+theme(axis.text.x=element_text(
                color = "white"),legend.text = element_text(size = 14),axis.title = element_text(size=14))+
            geom_text(aes(Teams,(Points+1),label=Points,size = 7)) })
    
    
    
    
}

shinyApp(ui=ui, server=server)


