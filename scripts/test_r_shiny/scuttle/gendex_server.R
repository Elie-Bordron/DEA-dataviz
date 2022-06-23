GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
source(file.path(GI_scripts_dir, "CGHcall.R"))
server <- function(input, output) {
    
    ######## HOME PANE
    output$contents <- DT::renderDataTable({
        file <- input$probeset_txt
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "txt", "Please upload a txt file"))
        
        x = read.table(file$datapath, header = input$header, sep="\t")
        x
    })
    
    ######## CGHcall PANE

    # output$segTable = DT::renderDataTable({data.frame(chr = c(1,2,3), nomenclature = c("arr[hg19]1q42.3q44(234,801,409-249,212,365)x1","arr[hg19]2p25.3q36(114,903,549-229,362,995)x4","arr[hg19]2p13.4p12(318,659,327-169,328,446)x3"), CN_state = c(1, 3.5, 3), "Median(Log2)" = c(-0.23, 0.21, 0.18), type=c("Perte heterozygote", "Gain", "Gain"), "size(Kbp)"=c(221405, 196328, 92678), "size(Mb)"=c(200, 196, 92), LOH=c("", "x", "x"), Genes=c("Incluant FH", "LOCUS", "Incluant FGFR3"), "Points de Cassure"=c("", "PID1", "NGEF") )})
        ### get call results

    rawProbesData = reactive({
        file <- input$probeset_txt
        rawProbesData = read.table(file$datapath, header = input$header, sep="\t")
        rawProbesData = rawProbesData[2:dim(rawProbesData)[2]] 
        # colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", colnames(rawProbesData)[dim(rawProbesData)[2]])
        print(c("rawProbesData: ", rawProbesData))
        rawProbesData
    })
    baseparams = getDefParams()
    params <- reactive({
        baseparams$CellularityCorrectSeg = input$correctCell
        baseparams$tumor_prop = input$cellularity
        baseparams$Prior = input$prior
        baseparams$UndoSD = input$undoSD
        baseparams$Minlsforfit = input$minSegLenForFit
        baseparams$sampleNames = colnames(rawProbesData())[dim(rawProbesData())[2]]
        baseparams$delete_this_value = "PC2979"
        return(baseparams)
    })
    resPipeline = reactive({  
        pipelineCGHcall(rawProbesData(), params())
        # suppressWarnings(pipelineCGHcall(rawProbesData(), params()))
    })
    CGHcall_segments = reactive({
        prbLvSegs = getPrbLvSegmentsFromCallObj(resPipeline())  
        prbLvSegs
    })
    segTable_calculated = reactive({
        params_obj = params()
        segtab = getSegTables(CGHcall_segments(),params_obj$sampleNames) [[1]]
        segtab
    })


    output$segTable = DT::renderDataTable({
        segTab = segTable_calculated()
        # DT::datatable(segTab)
        segTab
    })

    output$geneTable = DT::renderDataTable({data.frame(gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2))})
    
    output$testPlot <- renderPlot({
        params = params()
        plot(y=rawProbesData$Log2Ratio, x=rawProbesData$absPos, pch = 20, cex=0.01, col="dark grey", xlab = "", ylab = "log Ratio", ylim = c(-2,2))
        plotSegTableForWGV(segTable_calculated(), params$sampleNames, savePlot=FALSE, genGrid=FALSE, segColor="dark red", alreadyGoodPos=alreadyGoodPos) # segtables must have these columns: chrom, loc.start, loc.end, CN
        
        # plotSegTable(segTable_calculated(),params$sampleNames,savePlot=FALSE)
    })

    # output$debug = renderText(
    #     "debug text"
    # )
}
