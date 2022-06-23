GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
source(file.path(GI_scripts_dir, "CGHcall.R"))
server <- function(input, output) {
    
    ######## HOME PANE
    output$contents <- DT::renderDataTable({
        file <- input$probeset_txt
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "txt", "Please upload a txt file"))
        
        read.table(file$datapath, header = input$header, sep="\t")
    })
    
    
    
    ######## CGHcall PANE
    output$testPlot <- renderPlot({
        x = rnorm(100,10,5)
        plot(x)

        ###### use this:
        # plotSegTable(segTable,params$sampleNames,savePlot=FALSE)
    })

    # output$segTable = DT::renderDataTable({data.frame(chr = c(1,2,3), nomenclature = c("arr[hg19]1q42.3q44(234,801,409-249,212,365)x1","arr[hg19]2p25.3q36(114,903,549-229,362,995)x4","arr[hg19]2p13.4p12(318,659,327-169,328,446)x3"), CN_state = c(1, 3.5, 3), "Median(Log2)" = c(-0.23, 0.21, 0.18), type=c("Perte heterozygote", "Gain", "Gain"), "size(Kbp)"=c(221405, 196328, 92678), "size(Mb)"=c(200, 196, 92), LOH=c("", "x", "x"), Genes=c("Incluant FH", "LOCUS", "Incluant FGFR3"), "Points de Cassure"=c("", "PID1", "NGEF") )})
    output$segTable = DT::renderDataTable({
        params = getDefParams()
        params$CellularityCorrectSeg = input$correctCell
        params$tumor_prop = input$cellularity
        params$Prior = input$prior
        params$UndoSD = input$undoSD
        params$Minlsforfit = input$minSegLenForFit
        
        file <- input$probeset_txt
        rawProbesData = read.table(file$datapath, header = input$header, sep="\t")
        ncols = dim(rawProbesData)[2]
        params$sampleNames = colnames(rawProbesData)[ncols]
        
        rawProbesData = rawProbesData[2:ncols]
        colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", params$sampleNames)
        resPipeline = pipelineCGHcall(rawProbesData, params)
        CGHcall_segments = getPrbLvSegmentsFromCallObj(resPipeline)
        segTable = getSegTables(CGHcall_segments,params$sampleNames)[[1]]
    })
    output$geneTable = DT::renderDataTable({data.frame(gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2))})
    
}