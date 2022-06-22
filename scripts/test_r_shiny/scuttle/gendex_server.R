server <- function(input, output) {
    output$testPlot <- renderPlot({
        x = rnorm(100,10,5)
        plot(x)
    })
    output$segTable = renderTable({data.frame(chr = c(1,2,3), nomenclature = c("arr[hg19]1q42.3q44(234,801,409-249,212,365)x1","arr[hg19]2p25.3q36(114,903,549-229,362,995)x4","arr[hg19]2p13.4p12(318,659,327-169,328,446)x3"), CN_state = c(1, 3.5, 3), "Median(Log2)" = c(-0.23, 0.21, 0.18), type=c("Perte heterozygote", "Gain", "Gain"), "size(Kbp)"=c(221405, 196328, 92678), "size(Mb)"=c(200, 196, 92), LOH=c("", "x", "x"), Genes=c("Incluant FH", "LOCUS", "Incluant FGFR3"), "Points de Cassure"=c("", "PID1", "NGEF") )})
    output$geneTable = renderTable({data.frame(gene_id = c("BRCA1","CDK12","p53"), nb_breakpoints = c(0, 1, 0), nb_segments = c(1, 2, 1), copynumber = c(2,1,2))})
    output$contents <- renderTable({
        file <- input$file1
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "txt", "Please upload a txt file"))
        
        read.table(file$datapath, header = input$header, sep="\t")
    })
    
}