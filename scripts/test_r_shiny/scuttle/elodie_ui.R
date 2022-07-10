library(shiny)
library(shinythemes)
library(shinycssloaders)
library(plotly)
library(DT)
library(shinybusy)

row <- function(...) {
  tags$div(class="row", ...)
}

col <- function(width, ...) {
  tags$div(class=paste0("span", width), ...)
}


ui = shinyUI(
  fluidPage(navbarPage(theme = shinytheme("journal"),
                       add_busy_spinner(spin = "fading-circle"),
                       #titlePanel("Cell cycle entry during T-cell activation: a genomic perspective"),
                       tabPanel("Accueil",
                                
                                h1("Cell cycle entry during T-cell activation: a genomic perspective"),
                                verbatimTextOutput("acceuilText"),
                                h3("General"),
                                p("This application aims at help in exploring gene expression, chromatine accessibility and structure during the activation of T-cells."),
                                p("A spinning weel appears in the top right corner when the application is busy computing, loading, displaying data."),
                                p("Sometimes before some actions are made, errors are displayed in bold red because of missing data or something similar. Do not worry, it is normal and once the action is done everything goes to normal."),
                                p("I have made a LOT of tests but I may have missed things, in this case, an error is written in bold red. Usually, the app does not crash and an other action can be done. If nothing work anymore ... just write to me describing exactly what you did and a screeshot showing the error."),
                                h3("t-SNE clustering on SHAPs"),
                                #p("There are two tabs in this section. 't-SNE clustering on SHAPs' aims at exploring t-SNE results made from machine learning analysis of TFBS 
                                #mapping on ATAC peaks associated with differential interactions which involve differentially expressed genes.
                                #'ATAC + CHiC status vs logFC RNAseq' is informative and is inspired from Burren et al. Figure 3c."),
                                p("This section aims at exploring t-SNE results made from machine learning analysis of TFBS 
                                mapping on ATAC peaks associated with differential interactions which involve differentially expressed genes."),
                                #h4("t-SNE clustering on SHAPs"),
                                #"This section is a  bit complex at first seen.",
                                h4("1. The first thing to do is to choose a feature (presence/absence of a given TFBS or TF, ATAC peak annotation)."),
                                p("By default 'db.clust' is selected and represents the t-SNE clusters.
                                A check box allows to display the cluster ID, it needs to be checked before clicking on 'Show features'. It can be selected at any time but a re-click is needed to take it into account."),
                                p("At clicking on 'Show features', the t-SNE clustering with coloured peaks is diplayed. It has to be clicked each time you change of feature."),
                                h4("2. You are interested in a particular set of peaks:"),
                                h5("Select on the graph by click and drag"),
                                p("- This operation is very sensitive, any click on the graph is taken into account."),
                                p("- CAUTION ! If you select a new feature to display, be sure no area is selected in the plot. 
                                If one is selected, just click once in the plot outside of the area."),
                                p("- Bellow the graph, the number of selected peaks is displayed"), 
                                p("Be aware that all the peaks in the area selected are taken into account. 
                                Even the ones not allocated to any cluster and labelled '0'"),
                                h5("Focus on peaks in the area"),
                                p("Above the graph two check boxes are available: 'Focus on cluster' and 'Focus on feature'"),
                                p("- By checking the box 'Focus on cluster' you can select fragments from the most represented cluster in selected area."),
                                p("- By checking the box 'Focus on feature' you can select fragments according to the chosen feature in selected area.
                                In fact, clicking this box will display available value(s) for the displayed feature. 
                                If numeric, a slider appears (often 0 to 1), if categorical, all possible values are displayed and more than one can be checked."),
                                p("- Both Boxes can be checked and the peaks are automatically selected, you can see the number of peaks changing bellow the graph."),
                                p("- You can select the whole graph, focusing or not on features, both are very interesting."),
                                h4("3. Analyse selected peaks."),
                                p("Once you are happy with you selection, you can use the different tabs."),
                                p("Note that the last tab 'Heatmap TFBS/chip frequency clusters' is independent of the selected peaks and shows the complete set."),
                                h5("TF(BS)s frequencies"),
                                p("In this section appear 2 graphs. Left: the frequency of peaks with given TF(BS)s in a decreasing order. Right: heatmap showing the number of occurence of TF(BS)s in selected peaks with some annotations."),
                                p("- You have the possiblity to choose the TF(BS)s according to their frequency in the set of peaks using the slider 'Select x% of peaks'. This selection is applied to all subsections. By default, it is set at 20%."),
                                p("- You can split the peaks according to the clustering showed in the heatmap. It will label the peaks in tables of the other subsections."),
                                h5("Summary table for peaks"),
                                p("Details on selected ATAC peaks are displayed through a dynamic table where elements can be searched and columns can be ordered."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                h5("Summary plots for peaks"),
                                p("Several plots show the distribution of peak annotations. If clusters are asked from the heatmap in 'TF(BS)s frequencies' tab, the plots are splitted accordingly."),
                                h5("Interacting fragments"),
                                p("- Details on selected ATAC peaks and their interacting fragments (with or whitout ATAC peaks) are displayed through a dynamic table where 
                                  elements can be searched and columns can be ordered. The heatmap.info column refers to the heatmap clusters from 'TF(BS)s frequencies'. Cluster1 and 2 columns refer to t-SNE clustering if applicable."),
                                p("- Distribution of ATAC peaks in interacting fragments in t-SNE clusters."),
                                p("- Heatmap showing occurences of TF(BS)s in interacting ATAC peaks. The rows of the heatmap are separated, on the top the TF(BS)s from the selected peaks (labelled in the annotations with suffix '1', e.g. cl1, genom.1 etc.), bottom part concerns the interacting ATAC peaks (labelled in the annotations with suffix '2', e.g. cl2, genom.2 etc.). Only TF(BS)s occurring in at least x% of the peaks ('TF(BS)s frequencies')."),
                                p("As in the 'TF(BS)s frequencies' subsection, you can choose to split the interactions in groups according to the clustering showed in the heatmap."),
                                p("- You can select particular interactions according to the TF(BS)s present in the ATAC peaks of each fragment ('no.ATAC' allows to keep interactions with fragments without ATAC peaks)."),
                                p("Checking the box 'Select on heatmap' will select ATAC peaks and their interaction from the heatmap, meaning interaction with fragment without ATAC peaks are discarded."),
                                p("By default all the interactions from the heatmap are selected. You can select one or more particular clusters from the heatmap by checking the box(es) that appeared in the right. You may deselect the box 'all' if choosing cluster(s)."),
                                p("- TF(BS)s present in at least x% of the ATAC peaks are displayed for both fragments. The OR (default) means that the ATAC peaks carry at least one of the selected TF(BS)s, AND means ATAC peaks have to carry all selected ones."),
                                p("By default, all TF(BS)s are selected (checkbox 'all' that has to be deselected if choosing to select)."),
                                p("- The result of the selection(s) is displayed with a table (similar to the one at the top of the section) with an added column 'interact.heatmap' corresponding to the clustering from the heatmap."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                p("Some details about genes involved in the interactions are displayed. The column selected indicates if the gene is part of the final selection."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                h5("Heatmap bait frequency other clusters"),
                                p("In this tab the results depend on the selection made on the tSNE plot. You can select all fragments for an overall analysis."),
                                p("The heatmap is showing either the count or the frequency (cells) of 'other-fragment' interacting ATAC peaks occurring in each cluster (rows) by bait (columns)."),
                                p("You can select cluster(s) of interest and display the genes involved in interaction with fragments having ATAC in all selected clusters (intersection)."),
                                p("You can download the corresponding matrix as a tab file."),
                                h3("Enjoy !")
                                
                                
                       ),
                       tabPanel("t-SNE clustering on SHAPs",
                                                      uiOutput("TFs.tsne.list"),
                                                      checkboxInput("show.db.clust","Show cluster number"),
                                                      actionButton("go.tSNE.shaps","Show features"),
                                                      h5("To analyse a particular subset of peaks, select an area on the plot and details on selected peaks will appear bellow. By default all peaks in the area are selected."),
                                                      h5("Showing the cluster number is not automatic, you need to check the box and re-click on 'Show features'."),
                                                      h5("CAUTION ! If you select a new feature to display, be sure no area is selected in the plot. If one is selected, just click once in the plot outside of the area."),
                                                      
                                                      h5("---"),
                                                      h5("By checking the box 'Focus on cluster' you can select fragments from the most represented cluster in selected area."),
                                                      h5("By checking the box 'Focus on feature' you can select fragments according to the chosen feature in selected area."),
                                                      splitLayout(checkboxInput("show.all.points.in.area","Focus on cluster"),
                                                       checkboxInput("show.feature.in.area","Focus on feature")),
                                                      uiOutput("select.features"),
                                                      h5("---"),
                                                      h5("NB: You can select the whole picture."),
                                                      plotOutput("tSNE.shaps.icis",width=800,height=750,
                                                                 brush = brushOpts(id = "plot_brush")),
                                                      
                                                      textOutput("nb.frag.selected"),
                                                      h5(""),
                                           tabsetPanel(
                                             tabPanel("TF(BS)s frequencies",
                                                      h5("TF(BS)s present in more than x% of the fragments"),
                                                      sliderInput("min.pc.feature",label = "Select x% of peaks",min = 0,max = 1,value = 0.2,step = 0.1),
                                                      h5("If you want to split the peaks according to the heatmap, choose the number of clusters. This information will appear in the summary table."),
                                                      sliderInput("nb.clust.heatmap",label = "Number of clusters",min = 1,max = 10,value = 1,step = 1),
                                                      splitLayout(plotOutput("feature.shaps.clust",width=500,height=1000),
                                                      
                                                      plotOutput("heatmap.shaps.clust.plot",height=1000,width=500))
                                                      
                                             ),
                                             # tabPanel("Summary table for peaks",
                                             #          dataTableOutput("brush_info"),
                                             #          h5("Prefix for the name of the table to dowload. A column for each TF(BS) present in the heatmap is added."),
                                             #          textInput("prefix.save.shaps.tsne", "", value = ""),
                                             #          actionButton("go.save.shaps.tsne","Dowload table"),
                                             #          textOutput("download.done")
                                        
                                             # ),
                                             tabPanel("Summary plots for peaks",
                                                      plotOutput("tab.shaps.clust.plot",height=800),
                                                      plotOutput("tab.shaps.clust.plot2",height=400),
                                                      plotOutput("tab.shaps.clust.plot3",height=800)
                                                      
                                                      
                                             ),
                                             tabPanel("Interacting fragments",
                                                      #dataTableOutput("tab.interacting.tSNE"),
                                                      plotlyOutput("all.shaps.clust"),
                                                      sliderInput("nb.clust.heatmap.interact",label = "Number of clusters",min = 1,max = 10,value = 1,step = 1),
                                                      plotOutput("interacting.heatmap",height=1000),
                                                      h5("Select interactions according to feature on fragments"),
                                                      splitLayout(checkboxInput("select.on.heatmap.interaction","Select on heatmap"),uiOutput("select.on.heatmap.interaction.cluster")),
                                                      splitLayout(radioButtons("OR.AND.frag1","Fragment with selected ATAC peaks",choices = c("OR","AND"),selected = "OR"),radioButtons("OR.AND.frag2","Interacting fragments",choices = c("OR","AND"),selected = "OR")),
                                                      splitLayout(uiOutput("select.on.frag1"),uiOutput("select.on.frag2")),
                                                      dataTableOutput("selected.interacting.fragments"),
                                                      h5("Prefix for the name of the table to dowload. "),
                                                      textInput("prefix.save.interacting.tab", "", value = ""),
                                                      actionButton("go.save.interacting.tab","Dowload table"),
                                                      textOutput("download.interacting.done"),
                                                      h5("Dowdload the expression table for further functional analyses (e.g. https://maayanlab.cloud/Enrichr/ , http://biit.cs.ut.ee/gprofiler/gost)"),
                                                      dataTableOutput("tab.expression.interacting"),
                                                      textInput("prefix.save.interactingexpr.tab", "", value = ""),
                                                      actionButton("go.save.interacting.expr.tab","Dowload table"),
                                                      textOutput("download.interacting.expr.done")

                                           ),
                                           tabPanel("Genes and cluster",
                                                    selectInput("heatmap.gene.cluster",label = "",choices = c("only.baits","only.other","bait.other")),
                                                    p("Select cluster(s) to ..."),
                                                    selectInput("include.exclude",label = "",choices = c("include","exclude")),
                                                    uiOutput("select.cl.genes.ui"),
                                                    radioButtons("and.or.cl","Intersection or union",choices=c("AND","OR"),inline = T,selected = "AND"),
                                                    p("You can enter one or more gene symbol(s) separated by a white space. A black line appears in the annotation line 'selected' in the main heeatmap. If several genes are searched, it is not possible to know which one it is."),
                                                    textInput("symbol.to.display","",value=""),
                                                    actionButton("plot.genes.and.clusters.heatmap",label = "Plot it !"),
                                                    plotOutput("annotation.gNc",width=1350,brush="zoom_brush",height = 50),
                                                    plotOutput("genes.and.clusters.heatmap",width=1400,height=800),
                                                    plotOutput("zoomplot.gNc",width=1400,height=800),
                                                    textInput("gene_cluster.tab", "Enter a prefix file name for downloading current heatmap data", value = ""),
                                                    actionButton("load.gene_cluster.tab",label = "Download"),
                                                    textOutput("text.gene_cluster.tab")
                                                    
                                                    
                                           )
                                           # tabPanel("Heatmap TFBS/chip frequency clusters",
                                           #          plotOutput("heatmap.freq.TFBS",width=1200,height=600),
                                           #          plotOutput("heatmap.freq.chip",height=600)
                                           #  ),
                                           # tabPanel("Heatmap bait frequency other clusters",
                                           #          selectInput(inputId = "cl.inter",label = "Choose value to plot",choices = c("count","frequency")),
                                           #          p("For visualisation purpose, all the counts over 10 are showed in black."),
                                           #          p("You can add fragments with no ATAC peaks, not in the t-SNE analysis and/or in the t-SNE analysis but not associated with clusters by checking the boxes bellow."),
                                           #          checkboxGroupInput("display.no.clust","Display fragments without ATAC and/or not in tSNE",choices = c("in.tSNE","no.ATAC","no.tSNE","cluster0"),selected = "in.tSNE",inline = T),
                                           #          plotOutput("bait.cluster.interaction.plot",width=800,height=800),
                                           #          p("Selecting 'all' will display and allow downloading the whole matrix. Choosing several clusters end up to select the intersection (AND) or the union (OR)."),
                                           #          radioButtons("AND.OR.cl.choice","",choices=c("OR","AND"),inline = T,selected = "AND"),
                                           #          uiOutput("choose.cls"),
                                           #          actionButton("go.show.me.genes","Show selected genes"),
                                           #          textOutput("show.bait.w.other.cl"),
                                           #          actionButton("download.sel.tab","Download selected genes and cluster"),
                                           #          textOutput("tab.cl.frag.int"),
                                           #          tableOutput("show.bait.w.other.matrix")
                                                    
                                           # )
                                           )

                       ),
                       tabPanel("ChomHMM analysis",
                           
                                
                            sidebarLayout(
                                  sidebarPanel(
                                    h6("Emission states"),
                                    plotOutput("emissions",height=400),
                                    h6("Transitions between states"),
                                    plotOutput("transition",width=300,height=300)
                                  ),
                                  mainPanel(
                                tabsetPanel(
                                  
                                  tabPanel("ChromHMM enrichment peaks",
                                           selectInput("enr.chromHMM","Select group to visualize",choices = c("ATAC","feature")),
                                           h5("G0"),
                                           plotOutput("ATAC.enrichments.G0",height=280),
                                           h5("G1"),
                                           plotOutput("ATAC.enrichments.G1",height=280),
                                           h5("Difference G1/G0"),
                                           plotOutput("ATAC.enrichments",height=280)
                                  ),
                                  tabPanel("ChromHMM enrichment neighbourhood",
                                           uiOutput("select.feature.chromHMM"),
                                           h5("G0"),
                                           plotOutput("ATAC.cl.G0",height=250),
                                           h5("G1"),
                                           plotOutput("ATAC.cl.G1",height=250),
                                           h5("Difference G1/G0"),
                                           plotOutput("ATAC.cl",height=250)
                                  ),
                                  tabPanel("ChromHMM peaks",
                                           selectInput("select.profile.group","",choices=c("ATAC","clusters.tSNE")),
                                           uiOutput("show.profile.chromHMM"),
                                           plotOutput("plot.transition.profile"),
                                           plotOutput("segmentation.neighbours.plot"),
                                           dataTableOutput("dt.transition"),
                                           actionButton("save.chromHMM.tab","Download"),
                                           textOutput("text.save.chromHMM.tab")
                                  )
                                )
                                  )
                       )
                       )
                       
          )      
  )
  
)

source(paste0("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle/", "elodie_server.R"))
shinyApp(ui, server)
