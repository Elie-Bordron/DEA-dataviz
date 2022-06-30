options(shiny.maxRequestSize=102*1024^2)


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
        
        options = list(scrollY = '800px', pageLength = 1000) 

        x = read.table(file$datapath, header = input$probeSetHeader, sep="\t")
        x
    })
    
    ######## CGHcall PANE

    ### results computing
    rawProbesData = reactive({
        print("receiving rawProbesData")
        file <- input$probeset_txt
        rawProbesData = read.table(file$datapath, header = input$probeSetHeader, sep="\t",check.names=FALSE)
        # print(c("rawProbesData: ", rawProbesData))
        ProbeData_filtered = dplyr::select(rawProbesData, c("ProbeSetName", "Chromosome", "Position", starts_with("Log2Ratio")) )
        # print(c("ProbeData_filtered: ", ProbeData_filtered))
        ProbeData_filtered = mutate(ProbeData_filtered, END_POS=Position+20)
        # print(c("ProbeData_filtered after mutating: ", ProbeData_filtered))
        sampleName_within_colName = colnames(dplyr::select(ProbeData_filtered, starts_with("Log2Ratio")))
        split_by_leftBracket = stringr::str_split(sampleName_within_colName, "\\(") [[1]]
        # print(c("split_by_leftBracket: ", split_by_leftBracket))
        split_by_brackets = stringr::str_split(split_by_leftBracket[2], "\\)")[[1]]
        # print(c("split_by_brackets: ", split_by_brackets))
        sampleName = stringr::str_replace(split_by_brackets[1], ".OSCHP", "")
        colnames(ProbeData_filtered) = c("probeID", "CHROMOSOME", "START_POS", sampleName, "END_POS" )
        ProbeData_filtered = dplyr::select(ProbeData_filtered, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleName))
        # print(c("ProbeData_filtered: ", ProbeData_filtered))
        # rawProbesData = rawProbesData[2:dim(rawProbesData)[2]] 
        # colnames(rawProbesData) = c("probeID", "CHROMOSOME", "START_POS", "END_POS", colnames(rawProbesData)[dim(rawProbesData)[2]])
        # print(c("rawProbesData: ", rawProbesData))
        ProbeData_filtered
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
        return(baseparams)
    })
    # resPipeline = reactive({  
    resPipeline = eventReactive(input$go, {  
        print("calculating pipeline result")
        # pipelineCGHcall(rawProbesData(), params())
        suppressMessages(suppressWarnings(pipelineCGHcall(rawProbesData(), params())))
    })
    CGHcall_segments = reactive({
        print("Extracting probe-level segments from CGHcall result")
        # print(c("resPipeline(): ", resPipeline()))
        prbLvSegs = getPrbLvSegmentsFromCallObj(resPipeline())  
        prbLvSegs
    })
    # segTable_calculated = eventReactive(input$go,{
    segTable_calculated = reactive({
        print("Extracting segments table from probe-level segments")
        params_obj = params()
        segtab = getSegTables(CGHcall_segments(),params_obj$sampleNames)[[1]]
        segtab
    })
    



    ### plot
    chosenPlot <- reactive({
    # chosenPlot <- eventReactive(input$goPlot, {
        switch(input$plotChoice,
            profile={
                params = params()
                probeData = rawProbesData()
                colnames(probeData)[c(2:3)] = c("ChrNum", "ChrStart")
                probeData = getAbspos_probeset(probeData)
                colnames(probeData)[c(2:3)] = c("CHROMOSOME", "START_POS")
                currSampleName = params$sampleNames
                colnames(probeData)[which(colnames(probeData)==currSampleName)] <- "Log2Ratio"
                removePoints=10
                plotSegTableForWGV_GG(segTable_calculated(), probeData, removePoints)
                # plotSegTable(segTable_calculated(),params$sampleNames,savePlot=FALSE)
                # plot(c(5,5,5,5,5,5,5,5,6,5))
            },
            proba = {
                CGHcall_obj = resPipeline()
                print(c("CGHcall_obj: ", CGHcall_obj))
                plot(CGHcall_obj)
            }
        )
    })

    output$profilePlot = renderPlot(chosenPlot())
    # output$profilePlot = eventReactive(input$goPlot, {renderPlot(chosenPlot())})

    
    ### Segments table
    output$segTable = DT::renderDataTable({ segTable_calculated() })
    output$download_segTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_segTable")),
        content = function(fileST) {
            print(c("segTable_calculated(): ", segTable_calculated()))
            write.table(segTable_calculated(), fileST, quote=FALSE, row.names=FALSE, sep=";")
        }
    )    
    ### Gene table
    geneTable_obj = data.frame(gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2))
    output$geneTable = DT::renderDataTable({geneTable_obj})
    output$download_genesTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genesTable")),
        content = function(fileGT) {
            print(c("geneTable_obj: ", geneTable_obj))
            write.table(geneTable_obj, fileGT, quote=FALSE, row.names=FALSE, sep=";")
        }
    )
    
    ### GI output as text
    output$GItext <- renderText({
        resGI = calcGI_CGHcall(segTable_calculated()) 
        paste0(resGI[[2]], " alterations were found on ", resGI[[3]], " chromosomes. Genomic Index=", round(resGI[[1]],1) )
    }) 
    
    ######## Summary PANE
    ### GI table
    GI_table = data.frame(sample = c("1-RV","2-AD","3-ES"), GI_CGHcall = c(15,8,22), GI_ASCAT = c(12,18,31), GI_rCGH = c(9,6,13))
    output$GI_table_summary = DT::renderDataTable({GI_table})
    
    output$download_GI_table <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_GI_table")),
        content = function(fileGT) {
            write.table(GI_table, fileGT, quote=FALSE, row.names=FALSE, sep=";")
        }
    )
    ### genes table
    genes_table = data.frame(sample = c("BRCA1","CDK12","p53"), CN_CGHcall = c(2,1,2), CN_ASCAT = c(1,1,2), CN_rCGH = c(2,2,2))
    output$genes_table_summary = DT::renderDataTable({genes_table})
    
    output$download_genes_table_summary <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genes_table")),
        content = function(fileGT) {
            write.table(genes_table, fileGT, quote=FALSE, row.names=FALSE, sep=";")
        }
    )
    
    
}

