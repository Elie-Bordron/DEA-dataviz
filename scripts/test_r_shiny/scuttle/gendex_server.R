


server <- function(input, output) {
    #### variables
    GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
    working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
    results_dir = file.path(working_dir_shiny, "gendex_results")
    source(file.path(GI_scripts_dir, "CGHcall.R"))
    source(file.path(working_dir_shiny, "gendex_functions.R"))
    
    ######## HOME PANE
    output$probesetContent <- DT::renderDataTable({
        file <- input$probeset_txt
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "txt", "Please upload a txt file"))
        
        x = read.table(file$datapath, header = input$probeSetHeader, sep="\t")
        x
    })
    
    ######## CGHcall PANE
    rawProbesData = reactive({
        print("receiving rawProbesData")
        file <- input$probeset_txt
        rawProbesData = read.table(file$datapath, header = input$probeSetHeader, sep="\t")
        rawProbesData = rawProbesData[2:dim(rawProbesData)[2]] 
        # colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", colnames(rawProbesData)[dim(rawProbesData)[2]])
        # print(c("rawProbesData: ", rawProbesData))
        rawProbesData
    })
    baseparams = getDefParams()
    params <- reactive({
        print("initializing params")
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
        print("calculating pipeline result")
        # pipelineCGHcall(rawProbesData(), params())
        suppressMessages(suppressWarnings(pipelineCGHcall(rawProbesData(), params())))
    })
    CGHcall_segments = reactive({
        print("Extracting probe-level segments from CGHcall result")
        prbLvSegs = getPrbLvSegmentsFromCallObj(resPipeline())  
        prbLvSegs
    })
    segTable_calculated = eventReactive(input$go,{
        print("Extracting segments table from probe-level segments")
        params_obj = params()
        segtab = getSegTables(CGHcall_segments(),params_obj$sampleNames)[[1]]
        segtab
    })
    


    output$segTable = DT::renderDataTable({
        segTab = segTable_calculated()
        # DT::datatable(segTab)
        segTab
    })

    geneTable_obj = data.frame(gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2))
    output$geneTable = DT::renderDataTable({geneTable_obj})
    
    output$testPlot <- renderPlot({
        params = params()
        # plot(y=rawProbesData$Log2Ratio, x=rawProbesData$absPos, pch = 20, cex=0.01, col="dark grey", xlab = "", ylab = "log Ratio", ylim = c(-2,2))
        # plotSegTableForWGV(segTable_calculated(), params$sampleNames, savePlot=FALSE, genGrid=FALSE, segColor="dark red", alreadyGoodPos=alreadyGoodPos) # segtables must have these columns: chrom, loc.start, loc.end, CN
        
        plotSegTable(segTable_calculated(),params$sampleNames,savePlot=FALSE)
    })
    

    
    
    output$download_segTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_segTable")),
        content = function(fileST) {
            print(c("segTable_calculated(): ", segTable_calculated()))
            write.table(segTable_calculated(), fileST, quote=FALSE, row.names=FALSE, sep=";")
        }
    )    
    output$download_genesTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genesTable")),
        content = function(fileGT) {
            print(c("geneTable_obj: ", geneTable_obj))
            write.table(geneTable_obj, fileGT, quote=FALSE, row.names=FALSE, sep=";")
        }
    )    
    
}

