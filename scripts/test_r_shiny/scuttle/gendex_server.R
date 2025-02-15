# library(shiny)
# library(shinyBS)
# # library(shinyjs)
# library(shinybusy)
# library(tidyverse)
# library(DT)

options(shiny.maxRequestSize=102*1024^2)



# devtools::install_github('VanLoo-lab/ascat/ASCAT')



server <- function(input, output, session) { ##session is required for conditionalPanel js expression to use output variables
    #### variables
    # working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"


    if(TRUE){ # set this to TRUE on bergo PC 
        print("bergo path")
        working_dir_shiny = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
        GI_scripts_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
    } else { # this is used on my own PC
        #S
        # working_dir_shiny = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
        # GI_scripts_dir = "C:/Users/User/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
        #C
        working_dir_shiny = "C:/Users/warew/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/test_r_shiny/scuttle"
        GI_scripts_dir = "C:/Users/warew/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"    
    }
    results_dir = file.path(working_dir_shiny, "gendex_results")
    ### load functions
    source(file.path(GI_scripts_dir, "CGHcall.R"))
    source(file.path(GI_scripts_dir, "CGHcall_functions.R")) # for several CGHcall.R functions to work
    source(file.path(GI_scripts_dir, "OncoscanR_functions.R")) # calcGI()

    source(file.path(working_dir_shiny, "gendex_functions.R"))
    source(file.path(GI_scripts_dir, "rCGH.R")) # for rCGH pane
    source(file.path(GI_scripts_dir, "rCGH_functions.R")) # for several rCGH.R functions to work
    
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
    

    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ######## CGHcall PANE

    ### results computing
    probesData = reactive({
        print("receiving rawProbesData")
        req(!is.null(input$probeset_txt))
        file <- input$probeset_txt
        rawProbesData = read.table(file$datapath, header = input$probeSetHeader, sep="\t",check.names=FALSE)
        ### remove unused columns
        rawProbesData = dplyr::select(rawProbesData, -contains(c("BAF")))
        rawProbesData = dplyr::select(rawProbesData, -contains(c("WeightedLog2Ratio")))
        rawProbesData = dplyr::select(rawProbesData, -contains(c("NormalDiploid")))
        ### remove NA probes on AllelicDifference and Log2Ratio columns
        ### rename AllDiff column
        name_col_allDiff = colnames(dplyr::select(rawProbesData, contains("AllelicDifference")))
        colnames(rawProbesData)[which(colnames(rawProbesData)==name_col_allDiff)] = "AllelicDifference"
        ### get sampleName
        name_col_LRR = colnames(dplyr::select(rawProbesData, starts_with("Log2Ratio")))
        split_by_leftBracket = stringr::str_split(name_col_LRR, "\\(") [[1]]
        split_by_brackets = stringr::str_split(split_by_leftBracket[2], "\\)")[[1]]
        sampleName = stringr::str_replace(split_by_brackets[1], ".OSCHP", "")
        ###rename LRR column
        colnames(rawProbesData)[which(colnames(rawProbesData)==name_col_LRR)] = sampleName
        ### remove NA probes on AllelicDifference and Log2Ratio columns
        rawProbesData = dplyr::filter(rawProbesData, !is.na(AllelicDifference))
        rawProbesData = dplyr::filter(rawProbesData, !is.na(get(sampleName)))

        ProbeData_filtered = dplyr::select(rawProbesData, c("ProbeSetName", "Chromosome", "Position", sampleName, "AllelicDifference") )
        
        ### check if there are any NA's in the df
        # where_NA = colSums(is.na(ProbeData_filtered))
        # print(c("where_NA: ", where_NA))

        ### add ENDPOS
        ProbeDataWithEndPos = mutate(ProbeData_filtered, END_POS=Position+20)
        colnames(ProbeDataWithEndPos) = c("probeID", "CHROMOSOME", "START_POS", sampleName, "AllelicDifference", "END_POS" )
        ProbeData_renamed = dplyr::select(ProbeDataWithEndPos, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleName, "AllelicDifference"))
        ProbeData_renamed
    })
    params <- reactive({
        print("initializing params")
        baseparams = getDefParams()
        baseparams$CellularityCorrectSeg = input$correctCell
        baseparams$tumor_prop = input$cellularity
        baseparams$Prior = input$prior
        baseparams$UndoSD = input$undoSD
        baseparams$Minlsforfit = input$minSegLenForFit
        baseparams$sampleNames = colnames(probesData())[5]
        return(baseparams)
    })


    resPipeline_ctrlled_by_button <- eventReactive(input$go, { 
    # resPipeline_ctrlled_by_button = reactive({ 
        print("calculating pipeline result")
        probesData = probesData()
        params = params()
        # print(c("params used in pipeline: ", params))
        # print(c("probesData to be used in pipeline: ", probesData))
        ### removing AllelicDifference since it will be interpreted as another sample's LRR values by CGHcall.
        probesData = dplyr::select(probesData, -"AllelicDifference")
        # print(c("probesData used in pipeline: ", probesData))

        # print(c("probesData() used in pipeline: ", probesData))
        # print(c("colnames(probesData()): ", colnames(probesData())))
        # res = pipelineCGHcall(probesData, params())
        res = suppressMessages(suppressWarnings(pipelineCGHcall(probesData, params)))
        print(c("CGHcall pipeline result: ", res))
        res
    })

    randomVals <- eventReactive(input$go, {
        runif(50)
    })

    resPipeline = reactive({resPipeline_ctrlled_by_button()})
    
    CGHcall_segments = reactive({
        print("Extracting probe-level segments from CGHcall result")
        resPipeline = resPipeline()
        # resPipeline = randomVals()
        # print(c("resPipeline: ", resPipeline))
        # print(c("class(resPipeline): ", class(resPipeline)))
        # prbLvSegs = getPrbLvSegmentsFromCallObj(resPipeline, segsType="both")  

        prbLvSegs = getPrbLvSegments(resPipeline, segsType="both")
        prbLvSegs
    })
    segTable_calculated = reactive({
        print("Extracting segments table from probe-level segments")
        params_obj = params()
        CGHcall_segments = CGHcall_segments()
        # print("======================================================")
        print(c("(CGHcall_segments): ", (CGHcall_segments)))
        # print("======================================================")
        # CGHcall_segments = prepareSegtableByProbe(CGHcall_segments)
        segtab = get_seg_table(CGHcall_segments)
        # print(c("colnames(CGHcall_segments): ", colnames(CGHcall_segments)))
        # print(c("colnames(segtab): ", colnames(segtab)))
        
        # print(c("segtab from get_seg_table: ", segtab))
        # print(c("segtab$CN + 2: ", segtab$CN + 2))
        segtab$CN = segtab$CN + 2
        
        segtab
    })
    

    ### plot
    CGHcall_chosenPlot <- reactive({
    # CGHcall_chosenPlot <- eventReactive(input$go, {
        switch(input$CGHcall_plotChoice,
            profile={
                params = params()
                probeData = probesData()
                # colnames(probeData)[c(2:3)] = c("ChrNum", "ChrStart")
                # print(c("probeData sent to getAbspos_probeset : ", probeData))
                probeData = getAbspos_probeset(probeData)
                # colnames(probeData)[c(2:3)] = c("CHROMOSOME", "START_POS")
                currSampleName = params$sampleNames
                # print(c("probeData before renaming: ", probeData))
                # print(c("currSampleName: ", currSampleName))
                # print(c("colnames(probeData): ", colnames(probeData)))
                # print(c("colname==currSampleName: ", colnames(probeData)==currSampleName))
                colnames(probeData)[which(colnames(probeData)==currSampleName)] <- "Log2Ratio"
                # print(c("probeData after renaming: ", probeData))
                removePoints=10
                plotSegTableForWGV_GG(segTable_calculated(), probeData, removePoints)
                # plotSegTable(segTable_calculated(),params$sampleNames,savePlot=FALSE)
                # plot(c(5,5,5,5,5,5,5,5,6,5))
            },
            proba = {
                CGHcall_obj = resPipeline()
                # print(c("CGHcall_obj: ", CGHcall_obj))
                plot(CGHcall_obj)
            }
        )
    })

    output$CGHcall_profilePlot = renderPlot(CGHcall_chosenPlot())
    # output$CGHcall_profilePlot = eventReactive(input$goPlot, {renderPlot(CGHcall_chosenPlot())})

    output$CGHcall_allDiffPlot = renderPlot({
        CGHcall_allDiffPlot_reactive()
    })
    CGHcall_allDiffPlot_reactive <- eventReactive(input$go, {
        print("Allele Difference Plot")
        probeData = probesData()
        # print(c("probeData: ", probeData))
        probeData = dplyr::filter(probeData, CHROMOSOME<23)
        probeData = getAbspos_probeset(probeData)
        # print(c("probeData for alldiff plot: ", probeData))
        ggAllDiff = ggplot(data=probeData, aes(y=AllelicDifference, x=absPos))
        ggAllDiff = ggAllDiff + geom_point(aes(color=factor(CHROMOSOME)), size=0.01, shape=20, show.legend = FALSE)
        # colrVec = rep(c("pink", "green", "orange"), 7); colrVec = c (colrVec, "pink")
        colrVec = rep(c("#fc99d1", "#88ff88", "#fca355"), 7); colrVec = c (colrVec, "#fc99d1")
        names(colrVec) = unique(probeData$CHROMOSOME)
        ggAllDiff = ggAllDiff + scale_color_manual(values = colrVec)
        # print("Allele Difference Plot DONE")
        ggAllDiff = ggAllDiff + theme_bw()
        ggAllDiff

        ### to see if Alldiff data is in the place where Log2Ratio should be in callRes.
        # resPipe = resPipeline()
        # rawLRR = as.data.frame(resPipe@assayData[["copynumber"]])
        # rowsInfo = as.data.frame(fData(resPipe))
        # rawLRR = cbind(rowsInfo, rawLRR)
        # # print(c("rawLRR to be plotted: ", rawLRR))
        # colnames(rawLRR)[1:2] = c("CHROMOSOME", "START_POS")
        # rawLRR = getAbspos_probeset(rawLRR)
        # # print(c("rawLRR to be plotted, but 2 columns have been renamed: CHROMOSOME and START_POS. + absPos has been added : ", rawLRR))
        # LRRcol = rawLRR[["1-RV"]]
        # absPoscol = rawLRR[["absPos"]]
        # # print(c("class(LRRcol): ", class(LRRcol)))
        # # print(c("class(absPoscol): ", class(absPoscol)))
        # plot(y=LRRcol, x=absPoscol)
        # plot(c(9,2,5,6))
    })


    ### Segments table
    segtabToDisplay = reactive({
        segtabToDisplay = segTable_calculated() 
        # segtabToDisplay = dplyr::select(segtab, -contains("abs"))
        colnames(segtabToDisplay)[colnames(segtabToDisplay)=="Log2Ratio"] = "seg.mean"
        # segtabToDisplay$CN = segtabToDisplay$CN + 2

        parsed.str <- parse(text=paste0("dplyr::filter(segtabToDisplay,", input$selectCN, ")"))
        segtabToDisplay = eval(parsed.str)
        print(c("segtabToDisplay: ", segtabToDisplay))
        segtabToDisplay
    })
    output$CGHcall_segTable <- DT::renderDataTable(
        fillContainer=FALSE,
        expr={
            print("Creating table to display")
            req(!is.null(input$probeset_txt))
            segtabToDisplay()
        }, 
        quoted = FALSE,
    )
    output$CGHcall_download_segTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_segTable")),
        content = function(fileST) {
            # print(c("segTable_calculated(): ", segTable_calculated()))
            parsed.str <- parse(text=paste0("dplyr::filter(segtabToDisplay,", input$selectCN, ")"))
            segtabToSave = eval(parsed.str)
            write.table(segtabToSave, fileST, quote=FALSE, row.names=FALSE, sep=";")
        }
    )    
    ### Gene table
    # geneTable_obj = data.frame(
    #     gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2)
    # )
    # geneTableToDisplay = reactive({
    #     segTable = segTable_calculated()
    #     segTable = dplyr::select(segTable, -c("absStart", "absEnd"))
    #     # print(c("colnames(segTable): ", colnames(segTable)))
    #     segTable = segTable[c("Chromosome", "Start", "End", "nbProbes", "Log2Ratio", "seg.med", "probes.Sd", "CN")]
    #     colnames(segTable) = c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.med", "probes.Sd", "estimCopy")
    #     # print(c("segTable: ", segTable))
    #     geneTable = rCGH::byGeneTable(segTable)
    #     geneTable = as.data.frame(geneTable)
    #     # print(c("geneTable: ", geneTable))
    #     print(c("class(geneTable): ", class(geneTable)))
    #     geneTable
    # })
    # output$geneTable = DT::renderDataTable({
    #     geneTableToDisplay()
    #     # geneTable_obj
    #     })
    # output$download_genesTable <- downloadHandler(
    #     filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genesTable")),
    #     content = function(fileGT) {
    #         # print(c("geneTable_obj: ", geneTable_obj))
    #         write.table(geneTable_obj, fileGT, quote=FALSE, row.names=FALSE, sep=";")
    #     }
    # )
    

    ### GI output as text
    output$CGHcall_GItext <- renderText({
        req(!is.null(input$probeset_txt))
        segTab = segTable_calculated()
        # print(c("segTable_calculated() according to output$CGHcall_GItext", segTab))
        # segTab$CN = segTab$CN-2
        # colnames(segTab$seg.mean) = "CN"
        resGI = calcGI_CGHcall(segTab) 
        paste0(resGI[[2]], " alterations were found on ", resGI[[3]], " chromosomes. Genomic Index=", round(resGI[[1]],1) )
    }) 
    

    ###### panel that says "load a probeset.txt" and disappears when it is done
    # output$test_warnPanel = renderUI({
    output$test_warnPanel = renderText({
        req( is.null(input$probeset_txt$name))
        "Please load a probeset.txt from Home pane"
    })

    output$CGHcall_probesetLoaded = reactive({
        pb_loaded = is.character(input$probeset_txt$name)
        print(c("pb_loaded: ", pb_loaded))
        pb_loaded
    })
    
    
    #### output options for different output objects
    outputOptions(output, 'CGHcall_probesetLoaded', suspendWhenHidden = FALSE)







    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ########
    ######## rCGH PANE

    ### results computing
    # probesData = reactive({
    #     print("receiving rawProbesData")
    #     req(!is.null(input$probeset_txt))
    #     file <- input$probeset_txt
    #     rawProbesData = read.table(file$datapath, header = input$probeSetHeader, sep="\t",check.names=FALSE)
    #     # print(c("colnames(rawProbesData): ", colnames(rawProbesData)))
    #     ProbeData_filtered = dplyr::select(rawProbesData, c("ProbeSetName", "Chromosome", "Position", starts_with("Log2Ratio"), starts_with("AllelicDifference")) )
    #     # print(c("ProbeData_filtered: ", ProbeData_filtered))
    #     ProbeDataWithEndPos = mutate(ProbeData_filtered, END_POS=Position+20)
    #     # print(c("ProbeDataWithEndPos after mutating: ", ProbeDataWithEndPos))
    #     sampleName_within_colName = colnames(dplyr::select(ProbeDataWithEndPos, starts_with("Log2Ratio")))
    #     split_by_leftBracket = stringr::str_split(sampleName_within_colName, "\\(") [[1]]
    #     # print(c("split_by_leftBracket: ", split_by_leftBracket))
    #     split_by_brackets = stringr::str_split(split_by_leftBracket[2], "\\)")[[1]]
    #     # print(c("split_by_brackets: ", split_by_brackets))
    #     sampleName = stringr::str_replace(split_by_brackets[1], ".OSCHP", "")
    #     colnames(ProbeDataWithEndPos) = c("probeID", "CHROMOSOME", "START_POS", sampleName, "AllelicDifference", "END_POS" )
    #     ProbeData_renamed = dplyr::select(ProbeDataWithEndPos, c("probeID", "CHROMOSOME", "START_POS", "END_POS", sampleName, "AllelicDifference"))
    #     # print(c("ProbeData_renamed: ", ProbeData_renamed))
    #     # print(c("colnames(ProbeData_renamed): ", colnames(ProbeData_renamed)))
    #     ProbeData_renamed
    # })
    # baseparams = getDefParams()
    # params <- reactive({
    #     print("initializing params")
    #     baseparams$CellularityCorrectSeg = input$correctCell
    #     baseparams$tumor_prop = input$cellularity
    #     baseparams$Prior = input$prior
    #     baseparams$UndoSD = input$undoSD
    #     baseparams$Minlsforfit = input$minSegLenForFit
    #     baseparams$sampleNames = colnames(probesData())[dim(probesData())[2]]
    #     return(baseparams)
    # })

    

    # resrCGHpipeline = reactive({
    resrCGHpipeline = eventReactive(input$go_rCGH, { 
        print("calculating rCGH pipeline result")
        # print(c("colnames(probesData()): ", colnames(probesData())))
        file <- input$probeset_txt
        pathToProbeset_txt = file$datapath
        # res = pipeline_rCGH("1-RV", silent=TRUE)
        params_rCGH = getDefParamsrCGH()
        res = pipeline_rCGH(pathToProbeset_txt, silent=TRUE, params_rCGH)
        # print(c("res: ", res))
        # print(c("res@cnSet: ", res@cnSet))
        res
    })

    rCGH_segments = reactive({
        print("Extracting probe-level segments from rCGH result")
        # print(c("resrCGHpipeline(): ", resrCGHpipeline()))
        # print(c("class(resrCGHpipeline()): ", class(resrCGHpipeline())))
        # prbLvSegs = getPrbLvSegmentsFromCallObj(resrCGHpipeline(), segsType="both")  
        resPipe = resrCGHpipeline()
        print(c("resPipe used to create rCGH_segments: ", resPipe))
        prbLvSegs = getPrbLvSegments_rCGH(resPipe)
        prbLvSegs
    })
    segTable_rCGH = reactive({
        # print("Extracting segments table from probe-level segments")
        params_obj = params()
        rCGH_segments = rCGH_segments()
        # rCGH_segments = prepareSegtableByProbe(rCGH_segments)
        segtab = get_seg_table(rCGH_segments)
        # print(c("colnames(rCGH_segments): ", colnames(rCGH_segments)))
        # print(c("colnames(segtab): ", colnames(segtab)))
        # print(c("segtab: ", segtab))
        # segtab$CN = segtab$CN + 2 # rCGH segTable doesn't need this line.
        
        segtab
    })
    

    ### plot
    rCGH_chosenPlot <- reactive({
    # rCGH_chosenPlot <- eventReactive(input$goPlot, {
        # switch(input$rCGH_plotChoice,
            # profile={
                params = params()
                probeData = probesData()
                # colnames(probeData)[c(2:3)] = c("ChrNum", "ChrStart")
                probeData = getAbspos_probeset(probeData)
                # colnames(probeData)[c(2:3)] = c("CHROMOSOME", "START_POS")
                # print(c("probeData after getAbspos_probeset : ", probeData))
                currSampleName = params$sampleNames
                colnames(probeData)[which(colnames(probeData)==currSampleName)] <- "Log2Ratio"
                removePoints=10
                # print(c("probeData after renaming : ", probeData))
                plotSegTableForWGV_GG(segTable_rCGH(), probeData, removePoints)
                # plotSegTable(segTable_rCGH(),params$sampleNames,savePlot=FALSE)
                # plot(c(5,5,5,5,5,5,5,5,6,5))
            # },
            # proba = {
            #     # plot(rnorm(20))
            #     # plotProfile(resrCGHpipeline()) # profile plot
            #     plotLOH(resrCGHpipeline()) # AllelicDifference plot
            # }
        # )
    })

    output$rCGH_profilePlot = renderPlot(rCGH_chosenPlot())


    output$rCGH_allDiffPlot = renderPlot({
        rCGH_allDiffPlot_reactive()
    })

    #  = renderPlot({
    rCGH_allDiffPlot_reactive <- eventReactive(input$go_rCGH, {
        print("Allele Difference Plot")
        probeData = probesData()
        # print(c("probeData: ", probeData))
        probeData = dplyr::filter(probeData, CHROMOSOME<23)
        probeData = getAbspos_probeset(probeData)
        # print(c("probeData for alldiff plot: ", probeData))
        ggAllDiff = ggplot(data=probeData, aes(y=AllelicDifference, x=absPos))
        ggAllDiff = ggAllDiff + geom_point(aes(color=factor(CHROMOSOME)), size=0.01, shape=20, show.legend = FALSE)
        # colrVec = rep(c("pink", "green", "orange"), 7); colrVec = c (colrVec, "pink")
        colrVec = rep(c("#fc99d1", "#88ff88", "#fca355"), 7); colrVec = c (colrVec, "#fc99d1")
        names(colrVec) = unique(probeData$CHROMOSOME)
        ggAllDiff = ggAllDiff + scale_color_manual(values = colrVec)
        # print("Allele Difference Plot DONE")
        ggAllDiff = ggAllDiff + theme_bw()
        ggAllDiff
    })


    ### Segments table
    segtabToDisplay_rCGH = reactive({
    # segtabToDisplay_rCGH <- eventReactive(input$go_rCGH, {
        segtabToDisplay_rCGH = segTable_rCGH() 
        # segtabToDisplay_rCGH = dplyr::select(segtab, -contains("abs"))
        colnames(segtabToDisplay_rCGH)[colnames(segtabToDisplay_rCGH)=="Log2Ratio"] = "seg.mean"
        # segtabToDisplay_rCGH$CN = segtabToDisplay_rCGH$CN + 2

        parsed.str <- parse(text=paste0("dplyr::filter(segtabToDisplay_rCGH,", input$selectCN_rCGH, ")"))
        segtabToDisplay_rCGH = eval(parsed.str)
        # print(c("segtabToDisplay_rCGH: ", segtabToDisplay_rCGH))
        segtabToDisplay_rCGH
    })
    output$rCGH_segTable <- DT::renderDataTable(
        fillContainer=FALSE,
        expr={
            print("Creating table to display")
            req(!is.null(input$probeset_txt))
            segtabToDisplay_rCGH()
        }, 
        quoted = FALSE,
    )
    output$rCGH_download_segTable <- downloadHandler(
        filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_segTable")),
        content = function(fileST) {
            print(c("segTable_rCGH(): ", segTable_rCGH()))
            write.table(segTable_rCGH(), fileST, quote=FALSE, row.names=FALSE, sep=";")
        }
    )    
    ### Gene table
    # geneTable_obj = data.frame(
    #     gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2)
    # )
    # geneTableToDisplay = reactive({
    #     segTable = segTable_rCGH()
    #     segTable = dplyr::select(segTable, -c("absStart", "absEnd"))
    #     # print(c("colnames(segTable): ", colnames(segTable)))
    #     segTable = segTable[c("Chromosome", "Start", "End", "nbProbes", "Log2Ratio", "seg.med", "probes.Sd", "CN")]
    #     colnames(segTable) = c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.med", "probes.Sd", "estimCopy")
    #     # print(c("segTable: ", segTable))
    #     geneTable = rCGH::byGeneTable(segTable)
    #     geneTable = as.data.frame(geneTable)
    #     # print(c("geneTable: ", geneTable))
    #     print(c("class(geneTable): ", class(geneTable)))
    #     geneTable
    # })
    # output$geneTable = DT::renderDataTable({
    #     geneTableToDisplay()
    #     # geneTable_obj
    #     })
    # output$download_genesTable <- downloadHandler(
    #     filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genesTable")),
    #     content = function(fileGT) {
    #         # print(c("geneTable_obj: ", geneTable_obj))
    #         write.table(geneTable_obj, fileGT, quote=FALSE, row.names=FALSE, sep=";")
    #     }
    # )
    

    ### GI output as text
    rCGH_GI_res <- reactive({
        req(!is.null(input$probeset_txt))
        # print(c("segTable_rCGH() according to output$rCGH_GItext", segTable_rCGH()))
        segTab = segTable_rCGH()
        resGI = calcGI_rCGH(segTab) 
    }) 
    
    output$rCGH_GItext <- renderText({
        resGI = rCGH_GI_res()
        paste0(resGI[[2]], " alterations were found on ", resGI[[3]], " chromosomes. Genomic Index=", round(resGI[[1]],1) )
    }) 

    ###### panel that says "load a probeset.txt" and disappears when it is done
    # output$rCGH_test_warnPanel = renderUI({
    output$rCGH_test_warnPanel = renderText({
        req( is.null(input$probeset_txt$name))
        # # print(c("input$probeset_txt: ", input$probeset_txt))
        # # print(c("input$probeset_txt$name: ", input$probeset_txt$name))
        # fileName = input$probeset_txt$name
        # fileExt = substr(fileName, nchar(fileName)-12+1, nchar(fileName))
        # # print(c("file extension: ", fileExt)) # 12 = nb of char in "probeset.txt"
        
        # if(exists("fileName")){ print(c("fileName: ", fileName))} else {print("no 'name' attribute to input$probeset_txt")}
        # # paste0("input file given is currently ", input$probeset_txt$name, " .")
        
        "Please load a probeset.txt from Home pane"
    })

    output$rCGH_probesetLoaded = reactive({
        pb_loaded = is.character(input$probeset_txt$name)
        print(c("pb_loaded: ", pb_loaded))
        pb_loaded
    })
    
    
    #### output options for different output objects
    outputOptions(output, 'rCGH_probesetLoaded', suspendWhenHidden = FALSE)









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
    # genes_table = data.frame(sample = c("BRCA1","CDK12","p53"), CN_CGHcall = c(2,1,2), CN_ASCAT = c(1,1,2), CN_rCGH = c(2,2,2))
    # output$genes_table_summary = DT::renderDataTable({genes_table})
    
    # output$download_genes_table_summary <- downloadHandler(
    #     filename = buildFileName(res_dir=results_dir, prefix=paste0(input$prefix,"_genes_table")),
    #     content = function(fileGT) {
    #         write.table(genes_table, fileGT, quote=FALSE, row.names=FALSE, sep=";")
    #     }
    # )
    
    
}

