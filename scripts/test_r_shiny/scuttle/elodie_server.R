main.dir <- ""
#main.dir <- "/Users/elodiedarbo/Documents/projects/C-HiC_cell_T/app_integration"
res.dir <- data.dir <- file.path(main.dir,"data")

options(shiny.maxRequestSize=50*1024^2) 

server = shinyServer(function(input, output) {
  ########################################################################
  ## DATA
  ########################################################################
    withProgress(message = 'loading data', {
      message("Loading data ...")
      source(file.path(main.dir,"scripts","functions_for_report.R"))
      load(file.path(data.dir,"refseq_hg19.RData"))
      load(file.path(data.dir,"chrom_objects.RData")) # transitions,emissions,segments.G1,segments.G0,chromHMM.G0,chromHMM.G1,fantom
      load(file.path(data.dir,"other.frag.tab.RData")) # other.frag.tab
      load(file.path(data.dir,"chromHMM.lines.RData"))
      # homer.tmp,shaps,tsne.group,tsne.center,peaks.ATAC.frag,all.interactions,
      # course.expression,homer.num,external.chip.ATAC,YY1,TF,peaks.ATAC,RNAseq,tad_hg19, 
      # tads_borders
      load(file.path(res.dir,"objects_for_integration.RData")) 
      rm(tad_hg19,tads_borders)
      rm(TF)
      rm(peaks.ATAC.frag)
  
      all.interactions <- all.interactions[extend.classif!="stable"]
      #other.frag.stats <- all.interactions[type2=="o",list(SYMBOLs=paste(name1,collapse=","),nb.bait=length(name1)),by=name2]
      #other.frag.stats <- other.frag.stats[nb.bait>1]
      #write.table(other.frag.stats,file.path(res.dir,"other.frag.stats.tab"),sep="\t",quote=F,row.names=F)
      #rm(other.frag.stats)
      time.course <- data.table(read.table(file.path(data.dir,"time_course_GSE138767.tab"),sep="\t",head=T))
      files.H3 <- system(paste("ls",file.path(data.dir,"ChIP_seq","*summary")),intern=T)
      H3 <- do.call("c",lapply(files.H3,function(f){
        nom <- unlist(sapply(strsplit(basename(f),"-"),"[[",1))
        cond <- unlist(sapply(strsplit(nom,"_"),"[[",1))
        mark <- unlist(sapply(strsplit(nom,"_"),"[[",2))
        tmp <- read.table(f,sep="\t",head=T)
        colnames(tmp) <- c("seqnames","start","end","nb.cond","nb.control","pval","score","FDR")
        tmp$cond <- cond
        tmp$mark <- mark
        tmp <- as(tmp,"GRanges")
        tmp <- tmp[tmp$FDR<0.01]
      }))
      load(file.path(data.dir,"YY1_CNR.RData"))
      peaks.ATAC$H3K27me3 <- peaks.ATAC$H3K4me1 <- peaks.ATAC$H3K4me3 <- peaks.ATAC$H3K27Ac <- NULL
      lapply(c("H3K27me3","H3K4me1","H3K4me3","H3K27Ac"),function(mark,H3){
        G0 <- H3[H3$mark==mark & H3$cond=="G0"]
        G1 <- H3[H3$mark==mark & H3$cond=="G1"]
        ovl.G0 <- findOverlaps(peaks.ATAC,G0,minoverlap = 25)
        ovl.G1 <- findOverlaps(peaks.ATAC,G1,minoverlap = 25)
        pos <- 1:length(peaks.ATAC)
        posG1 <- queryHits(ovl.G1)
        posG0 <- queryHits(ovl.G0)
        peaks.ATAC$presence <<- ifelse(!pos%in%c(posG1,posG0),"00",ifelse(pos%in%posG1[posG1%in%posG0],"11",ifelse(pos%in%posG1[!posG1%in%posG0],"01","10")))
        colnames(mcols(peaks.ATAC))[ncol(mcols(peaks.ATAC))] <<- mark
        mark
      },H3)
      
      load(file.path(data.dir,"external.datasets.RData"))
      extern.chip <- c(GSE62486.peaks,GSE116695.peaks,list(YY1.down=YY1[YY1$diff=="down"],YY1.up=YY1[YY1$diff=="up"]))
      annots.extern.chip <- rbind(data.table(GSE62486.annots)[,list(ID,TF,rep=NA,cell,treatment)],data.table(GSE116695.annots)[,list(ID,TF,rep,cell=NA,treatment=NA)],data.table(ID=c("YY1.down","YY1.up"),TF=c("YY1","YY1"),rep=NA,cell=NA,treatment=NA))
      
      p50 <- reduce(c(GSE116695.peaks$p50_1,GSE116695.peaks$p50_3))
      NFAT1 <- reduce(c(GSE116695.peaks$NFAT1_1,GSE116695.peaks$NFAT1_3))
      NFAT2 <- reduce(c(GSE116695.peaks$NFAT2_1,GSE116695.peaks$NFAT2_3))
      polII_RA_Th1 <- reduce(c(GSE62486.peaks$Th1_reactivated_PolII,GSE62486.peaks$Th1_RA_PolII))
      polII_RA_Th2 <- reduce(c(GSE62486.peaks$Th2_reactivated_PolII,GSE62486.peaks$Th2_RA_PolII))
      PolII_US_Th1 <- reduce(c(GSE62486.peaks$Th1_unstimulated_PolII,GSE62486.peaks$Th1_US_PolII))
      PolII_US_Th2 <- reduce(c(GSE62486.peaks$Th2_unstimulated_PolII,GSE62486.peaks$Th2_US_PolII))
      
      exclude.ds <- c("Th1_Native_H3K4me3","Th1_RA_IgG","Th2_RA_IgG","Th2_Native_H3K4me3","p50_1","p50_3","NFAT1_1","NFAT2_3","NFAT2_1","NFAT1_3","Th1_reactivated_PolII","Th1_RA_PolII","Th2_reactivated_PolII","Th2_RA_PolII","Th1_unstimulated_PolII","Th1_US_PolII","Th2_unstimulated_PolII","Th2_US_PolII")
      
      extern.chip <- extern.chip[!names(extern.chip)%in%exclude.ds]
      annots.extern.chip <- annots.extern.chip[!ID%in%exclude.ds]
      
      extern.chip <- c(extern.chip,p50,NFAT1,NFAT2,polII_RA_Th1,polII_RA_Th2,PolII_US_Th1,PolII_US_Th2)
      annots.extern.chip <- rbind(annots.extern.chip,data.table(ID= c("p50","NFAT1","NFAT2","polII_RA_Th1","polII_RA_Th2","PolII_US_Th1","PolII_US_Th2"),TF=NA,rep=NA, cell=NA,treatment=NA))
      
      external.chip.ATAC <- lapply(extern.chip,function(x,peaks.ATAC){
        ovl <- findOverlaps(peaks.ATAC,x,minoverlap = 25,select="first")
        res <- data.table(chip=ifelse(is.na(ovl),0,1))
        colnames(res) <- ""
        res
      },peaks.ATAC)
      external.chip.ATAC <- do.call("cbind",external.chip.ATAC)
      colnames(external.chip.ATAC) <- annots.extern.chip$ID
      external.chip.ATAC$ATAC.id <- peaks.ATAC$id
      rm(extern.chip)
      rm(H3)
      rm(p50,NFAT1,NFAT2,polII_RA_Th1,polII_RA_Th2,PolII_US_Th1,PolII_US_Th2)
      message("Loading data done")
      
    })
  
  #################################################################################################
  ## accueil
  #################################################################################################
  
  output$acceuilText <- renderText({
    txt <- "Application is ready to use."
  })
  
  
  
  #################################################################################################
  ## whole data 
  #################################################################################################
  
  per.fragment.interactions.total <- reactive({
    message("Preparing table with fragments ...")
    int.fragments.stats <- int.fragments()
    int.fragments.stats <- int.fragments.stats[,list(total=length(diff.int),
                                                     nb.gained=sum(diff.int=="gain"),
                                                     nb.lost=sum(diff.int=="loss"),
                                                     nb.differential=sum(diff.int!="stable"),
                                                     direction=sum(direction),
                                                     direction.diff=sum(direction[diff.int!="stable"]),
                                                     directionality = round(sum(direction)/length(diff.int),2)
    ),by=c("name","type")]
    int.fragments.stats$chr <- "0"
    int.fragments.stats$start <- 0
    int.fragments.stats$pos <- "0"
    int.fragments.stats[type=="b"]$chr <- as.vector(seqnames(all.fragments))[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats[type=="b"]$start <- as.vector(start(all.fragments))[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats[type=="b"]$pos <- as.vector(all.fragments$pos)[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats
  })
  
  
  #################################################################################################
  ## frequence per cluster
  #################################################################################################
  
  shaps.cluster <- reactive({
    tsne.group$db.clust <- as.vector(tsne.group$db.clust)
    tsne.group$Feature <- shaps$Feature
    tsne.group$db.clust[!is.na(tsne.group$cl1.split)] <- tsne.group$cl1.split[!is.na(tsne.group$cl1.split)] 
    shaps.cl <- shaps[,-c(437:442)]
    row.names(shaps.cl) <- row.names(tsne.group)
    shaps.cl <- shaps.cl[tsne.group$db.clust!=0,]
    tsne.cl <- tsne.group[tsne.group$db.clust!=0,]
    
    
    
    shaps.cl.dt <- shaps.cl
    shaps.cl.dt$db.clust <- tsne.cl$db.clust
    shaps.cl.dt <- data.table(shaps.cl.dt)
  
    shaps.stats <- lapply(unique(tsne.cl$db.clust),function(x,shaps){
      s <- shaps[db.clust==x]
      s$db.clust <- NULL
      s <- as.data.frame(s)
      Q05 <- apply(s,2,quantile,probs=0.05)
      
      Q25 <- apply(s,2,quantile,probs=0.25)
      med.s <- apply(s,2,median)
      Q75 <- apply(s,2,quantile,probs=0.75)
      Q95 <- apply(s,2,quantile,probs=0.95)
      iqr.s <- apply(s,2,IQR)
      res <- data.table(TFBS=names(med.s),Q05,Q25,med.s,Q75,Q95,iqr.s,clust=x)
      res <- res[order(abs(med.s),decreasing = T)][1:30]
      
    },shaps.cl.dt)
    shaps.stats <- rbindlist(shaps.stats)
    
    shaps.tab <- as.data.frame.matrix(xtabs(med.s~TFBS+clust,shaps.stats))
    
    pheatmap(shaps.tab[apply(shaps.tab,1,function(x) sum(abs(x)>0.05)>1),],
             breaks=seq(-1,1,0.02),
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             border_color = NA,
             color=colorRampPalette(colours()[c(131,131,124,124,1,34,34,35,35)])(100))
    
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- unique(tsne.cl$db.clust)
    sum.shap <- apply(shaps.cl,2,function(x) sum(x!=0)>0)
    color.clust <- color.clust[!is.na(names(color.clust))]
    cl <- "24"
    sum.shap <- apply(shaps.cl[tsne.cl$db.clust==cl,],2,function(x) sum(abs(x)>0.1)>(sum(tsne.cl$db.clust==cl)/2))
    
    pheatmap(shaps.cl[tsne.cl$db.clust==cl,sum.shap],cluster_rows=T,
             show_rownames=F,show_colnames=T,
             breaks=seq(-2,2,0.04),
             annotation_row = tsne.cl[,c("db.clust","Feature")],
             annotation_colors = list(db.clust=color.clust),
             color=colorRampPalette(colours()[c(132,131,124,430,1,419,34,35,36)])(100))
  })
  
  
  freq.TFBS <- reactive({
    message("Preparing TFBS/ChIP frequency per cluster ...")
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    TFBS <- data.table(homer.tmp[,c(1:416,grep("db.clust",colnames(homer.tmp)))])
    chip <- data.table(homer.tmp[,c(432:466,grep("db.clust",colnames(homer.tmp)))])
    freq.TFBS <- TFBS[,lapply(.SD,function(x) sum(x>0)/length(x)),by=db.clust]
    chip.TFBS <- chip[,lapply(.SD,function(x) sum(x>0)/length(x)),by=db.clust]
    list(TFBS=freq.TFBS,chip=chip.TFBS)
  })
  
  output$heatmap.freq.TFBS <- renderPlot({
    mat <- freq.TFBS()$TFBS
    mat <- as.data.frame(mat)
    row.names(mat)  <- mat$db.clust
    mat <- mat[,-1]
    mat <- mat[,apply(mat,2,function(x){max(x)>0.2})]
    pheatmap(mat,clustering_distance_cols = "binary", clustering_distance_rows = "binary",
             clustering_method = "ward.D",fontsize = 7,color=colorRampPalette(colours()[c(1,142,53,34,35,24)])(100),breaks=seq(0,1,0.01))
  })
  
  output$heatmap.freq.chip <- renderPlot({
    mat <- freq.TFBS()$chip
    mat <- as.data.frame(mat)
    row.names(mat)  <- mat$db.clust
    mat <- mat[,-1]
    mat <- mat[,apply(mat,2,function(x){max(x)>0.1})]
    pheatmap(mat,clustering_distance_cols = "binary", clustering_distance_rows = "binary",
             clustering_method = "ward.D",fontsize = 7,color=colorRampPalette(colours()[c(1,142,53,34,35,24)])(100),breaks=seq(0,1,0.01))
  })
  
  bait.cluster.interaction <- reactive({
    interacting.tSNE <- interacting.tSNE()
    # m[,c("heatmap.info","class","cluster1","ATAC.id1","fragment.id1","ATAC.id2","fragment.id2","cluster2")]
    sense <- interacting.tSNE[grepl("^chr",fragment.id2)]
    antisense <- interacting.tSNE[grepl("^chr",fragment.id1)]
    interacting.tSNE <- rbind(sense[,list(fragment.id1,cluster1,fragment.id2,cluster2,ATAC.id2)],
                              antisense[,list(fragment.id1=fragment.id2,cluster1=cluster2,fragment.id2=fragment.id1,cluster2=cluster1,ATAC.id2=ATAC.id1)])
    interacting.tSNE <- unique(interacting.tSNE)
    if (!"no.ATAC"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="no.ATAC"]
    }
    if (!"not.in.tSNE"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="not.in.tSNE"]
    }
    if (!"cluster0"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="0"]
    }
    clust.assoc <- as.data.frame.matrix(table(interacting.tSNE$fragment.id1,interacting.tSNE$cluster2))
    
    list(clust.assoc=clust.assoc,raw=interacting.tSNE)
  })
  
  output$bait.cluster.interaction.plot <- renderPlot({
    bait.cluster.interaction <- bait.cluster.interaction()$clust.assoc
    if (input$cl.inter=="count"){
      pheatmap(as.data.frame(t(bait.cluster.interaction)),show_rownames = T, 
               show_colnames = F,
               border_color = NA,
               clustering_distance_cols = "binary",
               color = colorRampPalette(c("white","orange","red","brown","black"))(100),
               breaks = seq(0,10,0.1))
    }
    else if (input$cl.inter=="frequency"){
      tot <- apply(bait.cluster.interaction,1,sum)
      bait.cluster.interaction <- bait.cluster.interaction/tot
      pheatmap(as.data.frame(t(bait.cluster.interaction)),show_rownames = T, 
               show_colnames = F,
               border_color = NA,
               clustering_distance_cols = "binary",
               color = colorRampPalette(c("white","orange","red","brown","black"))(100),
               breaks = seq(0,1,0.01))
    }
  })
  
  output$choose.cls <- renderUI({
    mat <- bait.cluster.interaction()$clust.assoc
    checkboxGroupInput("choose.cl","Select clusters having other fragments",choices=c("all",colnames(mat)),inline = T)
  })
  
  show.bait.w.other.cl <- eventReactive(input$go.show.me.genes,{
    mat <- bait.cluster.interaction()$clust.assoc
    print(dim(mat))
    cl <- input$choose.cl
    print(cl)
    if (length(cl)>1){
      select.rows <- apply(mat[,cl],1,function(x,opt){
        if (opt=="AND"){
          f <- sum(x>0) == length(x)
        }
        else {
          f <- sum(x>0) > 0
        }
        
      },input$AND.OR.cl.choice)
    }
    else if (length(cl)==1 & cl!="all"){
      select.rows <- mat[,cl]>0
    }
    else if (length(cl)==1 & cl=="all"){
      select.rows <- rep(T,ncol(mat))
    }
    else {
      select.rows <- F
    }
    text.show <- paste("The cluster(s)",paste(cl,collapse = ", "), "contain(s) fragments interacting with", sum(select.rows), "genes: ",paste(row.names(mat)[select.rows],collapse=", "))
    if (input$cl.inter=="frequency"){
      tot <- apply(mat,1,sum)
      mat <- mat/tot
    }
    if (cl[1]!="all") {
      mat <- mat[select.rows,cl,drop=F]
    }
    mat$SYMBOL <- row.names(mat)
    mat <- mat[,sort(colnames(mat),decreasing = T)]
    print(select.rows)
    list(text.show=text.show,mat=mat)
  })
  
  output$show.bait.w.other.cl <- renderText({
    select.rows <- show.bait.w.other.cl()$text.show
  })
  
  output$show.bait.w.other.matrix <- renderTable({
    select.rows <- show.bait.w.other.cl()$mat
  })
  
  tab.cl.frag.int <- eventReactive(input$download.sel.tab,{
    select.rows <- show.bait.w.other.cl()$mat
    cl <- input$choose.cl
    write.table(select.rows,file.path(data.dir,"tables",paste0("selected_genes_with_int_in_cl_",paste(cl,collapse="_"),".tab")),sep="\t",quote=F)
    text <- paste0("File ",paste0("selected_genes_with_int_in_cl_",paste(cl,collapse="_"),".tab")," has been saved.")
  })  
  
  output$tab.cl.frag.int <- renderText({
    tab.cl.frag.int()
  })
  
  #################################################################################################
  ## CHiC plots 
  #################################################################################################
  
  int.fragments <- reactive({
    int.fragments <- all.interactions[,list(id1,id2,name1,name2,type1,type2,type,diff.int=extend.classif,TAD,distance,direction)]
    
    m <- match(int.fragments$name1,RNAseq$ID)
    int.fragments$expr1 <- RNAseq$status[m]
    int.fragments$logFC1 <- RNAseq$logFC[m]
    int.fragments$sigExpr1 <- -log10(RNAseq$adj.P.Val[m])
    m <- match(int.fragments$name2,RNAseq$ID)
    int.fragments$expr2 <- RNAseq$status[m]
    int.fragments$logFC2 <- RNAseq$logFC[m]
    int.fragments$sigExpr2 <- -log10(RNAseq$adj.P.Val[m])
    
    #int.fragments <- int.fragments[!is.na(expr1)]
    #int.fragments <- int.fragments[type=="bo" | (type=="bb" & !is.na(expr2))]
    
    int.fragments <- rbind(int.fragments[,list(id=id1,name=name1,type=type1,B2B=type,diff.int,TAD,distance,direction,sigExpr=sigExpr1,logFC=logFC1)],int.fragments[,list(id=id2,name=name2,type=type2,B2B=type,diff.int,TAD,distance,direction=-direction,sigExpr=sigExpr2,logFC=logFC2)])
    
  })
  
  
  
  
  #################################################################################################
  ## integration - tSNE
  #################################################################################################
  
  
  prepare.shaps <- reactive({
    withProgress(message="Preparing data for SHAP analysis",{
      diff.fragments <- all.interactions[extend.classif!="stable",list(id1,id2,name1,name2,type1,type2,type,diff.int=extend.classif)]
      
      m <- match(diff.fragments$name1,RNAseq$ID)
      diff.fragments$expr1 <- RNAseq$status[m]
      diff.fragments$logFC1 <- RNAseq$logFC[m]
      
      m <- match(diff.fragments$name2,RNAseq$ID)
      diff.fragments$expr2 <- RNAseq$status[m]
      diff.fragments$logFC2 <- RNAseq$logFC[m]
      
      diff.fragments <- diff.fragments[!is.na(expr1)]
      diff.fragments <- diff.fragments[type=="bo" | (type=="bb" & !is.na(expr2))]
      diff.fragments <- rbind(diff.fragments[,list(id=id1,name=name1,type=type1,B2B=type,diff.int,diff.expr=expr1,logFC.expr=logFC1)],diff.fragments[,list(id=id2,name=name2,type=type2,B2B=type,diff.int,diff.expr=ifelse(type=="bo",expr1,expr2),logFC.expr=ifelse(type=="bo",logFC1,logFC2))])
      diff.fragments[,logFC.expr:=NULL]
  
      diff.fragments <- unique(diff.fragments[,-4,with=F],by=NULL)
      
      diff.fragments$seqnames <- as.vector(seqnames(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      diff.fragments$start <- as.vector(start(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      diff.fragments$end <- as.vector(end(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      
      diff.fragments <- as(diff.fragments[,list(seqnames,start,end,id=id,diff.int,bait.SYMBOL=name,type,bait.status=diff.expr)],"GRanges")
      ovl <- findOverlaps(peaks.ATAC,diff.fragments)
      nb.int.diff <- length(unique(subjectHits(ovl)))
      peaks.ATAC.frag <- peaks.ATAC[queryHits(ovl)]
      peaks.ATAC.frag$in.diff.frag <- "yes" 
      peaks.ATAC.frag$diff.int <- diff.fragments$diff.int[subjectHits(ovl)]
      peaks.ATAC.frag$id.frag <- diff.fragments$id[subjectHits(ovl)]
      peaks.ATAC.frag$bait.symbol<- diff.fragments$bait.SYMBOL[subjectHits(ovl)]
      peaks.ATAC.frag$bait.status<- diff.fragments$bait.status[subjectHits(ovl)]
      peaks.ATAC.frag$frag.type<- diff.fragments$type[subjectHits(ovl)]
      peaks.ATAC.frag$all.marks <- paste(peaks.ATAC.frag$H3K27me3,"/",peaks.ATAC.frag$H3K4me3,peaks.ATAC.frag$H3K4me1,peaks.ATAC.frag$H3K27Ac,sep="/") 
      m <- match(peaks.ATAC.frag$id,homer.num$ATAC.id)
      homer.num <- homer.num[m,]
      homer.num$Feature <- paste(peaks.ATAC.frag$diff.int,peaks.ATAC.frag$bait.status,sep="_")
      homer.num$Interaction_Number <- 1:nrow(homer.num)
      #homer.num <- homer.num[order(homer.num$Interaction_Number),]
      homer.num$ATAC.status <- peaks.ATAC.frag$status[!is.na(m)]
      homer.num$frag.type <- peaks.ATAC.frag$frag.type[!is.na(m)]
      homer.num$genomic <- peaks.ATAC.frag$genomic[!is.na(m)]
      homer.num$CG <- peaks.ATAC.frag$CG[!is.na(m)]
      homer.num$bait.status <- peaks.ATAC.frag$bait.status[!is.na(m)]
      homer.num$H3K4me1 <- peaks.ATAC.frag$H3K4me1[!is.na(m)]
      homer.num$H3K4me3 <- peaks.ATAC.frag$H3K4me3[!is.na(m)]
      homer.num$H3K27Ac <- peaks.ATAC.frag$H3K27Ac[!is.na(m)]
      homer.num$H3K27me3 <- peaks.ATAC.frag$H3K27me3[!is.na(m)]
      
      homer.num$all.marks <- peaks.ATAC.frag$all.marks[!is.na(m)]
      
      homer.num$SYMBOL.ATAC <- peaks.ATAC.frag$SYMBOL[!is.na(m)]
      homer.num$frag.ID <- peaks.ATAC.frag$bait.symbol[!is.na(m)]
      
      homer.num <- cbind(homer.num,external.chip.ATAC[match(homer.num$ATAC.id,external.chip.ATAC$ATAC.id),!colnames(external.chip.ATAC)%in%c("ATAC.id","db.clust"),with=F])
      
      
      shaps <- na.omit(shaps)
      shaps <- shaps[!grepl("stable",shaps$feature),]
      
      m <- match(shaps$ID,homer.num$Interaction_Number)
      
      homer.tmp <- homer.num[na.omit(m),!duplicated(colnames(homer.num))]
      row.names(homer.tmp) <- homer.tmp$Interaction_Number
      homer.tmp$db.clust <- factor(tsne.group$db.clust)
      homer.tmp$cl1.split <- tsne.group$cl1.split
      
      colnames(homer.tmp) <- unlist(sapply(strsplit(colnames(homer.tmp),"[.][.]"),"[[",1))
      to.discard <- which(apply(homer.tmp[,1:436],2,sum)==0)
      homer.tmp <- homer.tmp[,-to.discard]
      return(list(homer.tmp=homer.tmp,shaps=shaps,tsne.group=tsne.group,tsne.center=tsne.center,peaks.ATAC.frag=peaks.ATAC.frag,diff.fragments=diff.fragments))
      
    })
  })
  
  output$TFs.tsne.list <- renderUI({
    homer.tmp <- prepare.shaps()$homer.tmp
    toselect <- colnames(homer.tmp)[!colnames(homer.tmp)%in%c("Interaction_Number","frag.ID","SYMBOL.ATAC","ATAC.id")]
    selectInput("select.TFs.tsne.list","Choose feature",choices=rev(toselect),selected="db.clust")
  })
  
  tSNE.shaps.icis <- eventReactive(input$go.tSNE.shaps ,{
    shaps.data <- prepare.shaps()
    homer.tmp <- shaps.data$homer.tmp
    shaps <- shaps.data$shaps
    tsne.group <- shaps.data$tsne.group
    tsne.center <- shaps.data$tsne.center
    
    tsne.center$to.eval <- tsne.center$db.clust
    tsne.center$label <- tsne.center$db.clust
    
    symbols <- rbindlist(lapply(colnames(homer.tmp[1:416]),function(x){
      x <- sub("T[.]box","Tbox",x)
      tmp <- unlist(strsplit(x,"[.]"))
      TF <- paste(tmp[-length(tmp)],collapse=".")
      type <- tmp[length(tmp)]
      data.table(tf=TF,type)
    }))
    TF.labels <- apply(homer.tmp[,1:416],1,function(x,symbols){
      paste(symbols[x>0],collapse ="\n")
    },symbols$tf)
    tsne.group$label <- paste(homer.tmp$Feature,paste0("genomic:",homer.tmp$genomic,":",homer.tmp$SYMBOL.ATAC),paste0("fragment:",homer.tmp$frag.ID),paste0("H3K4me3:",homer.tmp$H3K4me3),paste0("H3K4me1:",homer.tmp$H3K4me1),paste0("H3K27Ac:",homer.tmp$H3K27Ac),TF.labels,sep="\n")
    tsne.group$to.eval <- homer.tmp[,input$select.TFs.tsne.list]
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- as.character(1:length(color.clust))
    if (is.numeric(tsne.group$to.eval)){
      tsne.group <- tsne.group[order(tsne.group$to.eval,decreasing=F),]
    }
    g <- ggplot(aes(x=tSNE1, y=tSNE2, group=to.eval,label=label), data=tsne.group) + 
      geom_point(aes(x=tSNE1, y=tSNE2, colour=to.eval),size=0.5) +   
      theme(panel.background = element_rect(fill = "white", colour = "black")) +
      geom_hline(yintercept=0,linetype="dashed",size=0.3) +
      geom_vline(xintercept=0,linetype="dashed",size=0.3) 
    if (is.numeric(tsne.group$to.eval)){
      g <- g + scale_colour_gradientn(colours = colours()[c(234,420,34,36,24)])
        #scale_colour_gradient2(low = "blue",mid = "lightgrey",high = "red",midpoint = 0)
      #colours = colours()[c(234,420,34,36,24)]
    } 
    else if (input$select.TFs.tsne.list%in%c("H3K4me1","H3K4me3","H3K27Ac","H3K27me3")) {
      g <- g + scale_colour_manual(values=c("00"="grey","10"="blue","01"="red","11"="green"))
    }
    else if (input$select.TFs.tsne.list=="all.marks") {
      names(color.clust) <- unique(tsne.group$to.eval)
      g <- g + scale_colour_manual(values=color.clust)
    }
    else if (input$select.TFs.tsne.list=="cl1.split") {
      names(color.clust) <- unique(na.omit(unique(tsne.group$cl1.split)))
      g <- g + scale_colour_manual(values=color.clust)
    }
    else {
      g <- g + scale_colour_manual(values=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613],closing=colours()[131],opening="red","stable"="grey",gain_up="red",gain_stable=colours()[235],gain_down=colours()[53],loss_up=colours()[96],loss_stable=colours()[613],loss_down=colours()[124],up="red",down=colours()[124],gain="red",loss=colours()[124],b=colours()[613],o=colours()[34],color.clust)) 
    }
    if (input$show.db.clust){
      g <- g + geom_text(data=tsne.center)
    }
    g
  })
  
  output$tSNE.shaps.icis <- renderPlot({
    tSNE.shaps.icis()
  })
  
  
  brush_info <- reactive({
    data.brush <- input$plot_brush
    tsne.group$db.clust <- as.vector(tsne.group$db.clust)
    tsne.group$db.clust[!is.na(tsne.group$cl1.split)] <- tsne.group$cl1.split[!is.na(tsne.group$cl1.split)] 
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    
    
    select.brush <- tsne.group$tSNE1>=data.brush$xmin & tsne.group$tSNE1<=data.brush$xmax & tsne.group$tSNE2>=data.brush$ymin & tsne.group$tSNE2<=data.brush$ymax
    enriched <- homer.tmp[select.brush,c(1:416,432:466)]
    annots <- homer.tmp[select.brush,c("ATAC.id","db.clust","Feature","ATAC.status","frag.type","genomic","SYMBOL.ATAC","H3K4me1","H3K4me3","H3K27Ac","H3K27me3")]
    #annots$nb.annots <- apply(enriched,1,function(x) sum(x>0))
    annots.rows <- annots
    tsne.group <- tsne.group[select.brush,]
    if (input$show.all.points.in.area){
      nb.f <- table(tsne.group$db.clust)
      cl <- names(nb.f)[which.max(nb.f)]
      enriched <- enriched[annots.rows$db.clust==cl,]
      annots.rows <- annots.rows[annots.rows$db.clust==cl,]
    }
    
    list(enriched=enriched,annots.rows=annots.rows)
  })
  
  select.features <- reactive({
    data <- brush_info()
    enriched <- data$enriched
    annots.rows <- data$annots.rows
    
    if (input$show.feature.in.area){
      enriched <- cbind(enriched,annots.rows[,!colnames(annots.rows)%in%colnames(enriched),])
      
      if(is.numeric(enriched[,input$select.TFs.tsne.list])){
        sliderInput("nb.min.feature","Choose minimal number of feature per peak to select",min = 0,max=max(enriched[,input$select.TFs.tsne.list]),value = 0,step=1)
      }
      else if(!is.numeric(enriched[,input$select.TFs.tsne.list])){
        choices <- unique(as.vector(enriched[,input$select.TFs.tsne.list]))
        checkboxGroupInput("select.features.values","Select feature(s) to keep peaks carrying it/them.",choices = choices,inline = T)
      }
    }
  })
  
  output$select.features <- renderUI({
    select.features()
  })
  
  brush_info2 <- reactive({
    data <- brush_info()
    enriched <- data$enriched
    annots.rows <- data$annots.rows
    
    if (input$show.feature.in.area){
      n <- ncol(enriched)
      enriched <- cbind(enriched,annots.rows[,!colnames(annots.rows)%in%colnames(enriched),])
      if(is.numeric(enriched[,input$select.TFs.tsne.list])){
        annots.rows <- annots.rows[enriched[,input$select.TFs.tsne.list]>=input$nb.min.feature,]
        enriched <- enriched[enriched[,input$select.TFs.tsne.list]>=input$nb.min.feature,]
      }
      if(!is.numeric(enriched[,input$select.TFs.tsne.list])){
        annots.rows <- annots.rows[enriched[,input$select.TFs.tsne.list]%in%input$select.features.values,]
        enriched <- enriched[enriched[,input$select.TFs.tsne.list]%in%input$select.features.values,]
      }
      enriched <- enriched[,1:n]
    }
    list(enriched=enriched,annots.rows=annots.rows)
  })
  
  output$nb.frag.selected <- renderText({
    text <- brush_info2()$annots.rows
    text <- paste(nrow(text),"ATAC peaks are selected.")
  })
  
  heatmap.shaps.clust <- reactive({
    brush.data <- brush_info2()
    enriched <- brush.data$enriched
    annots.rows <- brush.data$annots.rows
    annots.rows <- annots.rows[,c("ATAC.id","db.clust","Feature","ATAC.status","SYMBOL.ATAC","frag.type","genomic","H3K4me1","H3K4me3","H3K27Ac","H3K27me3")]
    enriched <- enriched[!duplicated(annots.rows$ATAC.id),]
    annots.rows <- annots.rows[!duplicated(annots.rows$ATAC.id),]
    row.names(enriched) <- annots.rows$ATAC.id
    row.names(annots.rows) <- annots.rows$ATAC.id
    annots.rows <- annots.rows[,-1]
    enriched <- enriched[,apply(enriched,2,function(x) (sum(x>0)/nrow(enriched))>input$min.pc.feature)]
    
    annots.rows <- annots.rows[row.names(annots.rows)%in%row.names(enriched),]
    enriched <- enriched[,order(apply(enriched,2,function(x) sum(x>0)/nrow(enriched)),decreasing=T)]
    #annots.rows <- annots.rows[,-1]
    clust <-  pheatmap(t(enriched),
                       clustering_distance_cols = "binary",
                       clustering_method = "ward.D",silent=T)
    nb.clust <- input$nb.clust.heatmap
    cluster=factor(cutree(clust$tree_col, k = nb.clust))
    annots.rows$heatmap <- cluster[match(row.names(annots.rows),names(cluster))]
    annots.rows <- annots.rows[match(names(cluster),row.names(annots.rows)),]
    #annots.rows$position <- 1:nrow(annots.rows)
    
    data <- list(enriched=enriched,annots.rows=annots.rows,clust.heatmap=clust)    
  })
  
  heatmap.shaps.clust.plot <- reactive({
    data <- heatmap.shaps.clust()
    data$annots.rows$SYMBOL.ATAC <- NULL
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust))
    color.clust <- color.clust[!is.na(names(color.clust))]
    pheatmap(t(data$enriched),
             clustering_distance_cols = "binary",
             #clustering_distance_rows = "binary",
             cluster_rows = F,
             color=colorRampPalette(colours()[c(1,34,35,24)])(100),
             show_rownames = T,show_colnames = F,
             border_color = NA, annotation_col = data$annots.rows,
             clustering_method = "ward.D",
             annotation_colors = list(db.clust=color.clust[names(color.clust)%in%data$annots.rows$db.clust],genomic=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),frag.type=c(o="black",b="grey"),Feature=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124]),ATAC.status=c(closing=colours()[131],opening="red",stable="grey"),H3K4me1=c("00"="grey","10"="blue","01"="red","11"="green"),H3K4me3=c("00"="grey","10"="blue","01"="red","11"="green"),H3K27me3=c("00"="grey","10"="blue","01"="red","11"="green"),H3K27Ac=c("00"="grey","10"="blue","01"="red","11"="green")))
  })
  
  output$heatmap.shaps.clust.plot <- renderPlot({
    heatmap.shaps.clust.plot()
  })
  
  output$feature.shaps.clust <- renderPlot({
    data <- heatmap.shaps.clust()
    enriched <- data$enriched
    features <- sapply(enriched,function(x) sum(x>0)/nrow(enriched))
    features <- data.frame(features=names(features),freq=as.numeric(features))
    features$features <- factor(as.vector(features$features),levels=features$features[order(features$freq,decreasing = F)])
    ggplot(features,aes(x=features,y=freq)) + geom_col(fill="grey",colour="black") + theme_bw() +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12)) + 
      coord_flip()
  })
  
  
  output$brush_info <- renderDataTable({
    peaks.ATAC.frag <- prepare.shaps()$peaks.ATAC.frag
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows$RNAseq.ATAC <- RNAseq[match(annots.rows$SYMBOL.ATAC,ID)]$status
    annots.rows$cinetic.ATAC <- time.course[match(annots.rows$SYMBOL.ATAC,SYMBOL)]$time.course
    m <- match(row.names(annots.rows),peaks.ATAC.frag$id)
    annots.rows$frag.id <-  peaks.ATAC.frag$bait.symbol[m]
    annots.rows$RNAseq.bait <- RNAseq[match(annots.rows$frag.id,ID)]$status
    annots.rows$cinetic.bait <- time.course[match(annots.rows$frag.id,SYMBOL)]$time.course
    tSNE.info <- interacting.tSNE()
    m <- data.table(as.data.frame(matches(row.names(annots.rows),tSNE.info$ATAC.id1,all.x=T,all.y=F)))
    m$frag.id2 <- tSNE.info$fragment.id2[m$y]
    m <- m[,list(frag.id2=paste(unique(frag.id2),collapse=", ")),by=x]
    annots.rows$int.frag[m$x] <- m$frag.id2
    annots.rows[,!grepl("H3K",colnames(annots.rows))]
  })
  
  test.tab <- eventReactive(input$go.save.shaps.tsne,{
    withProgress(message = 'Dowloading files ...', {
      annots.rows <- heatmap.shaps.clust()$annots.rows
      annots.rows$RNAseq <- RNAseq[match(annots.rows$SYMBOL.ATAC,ID)]$status
      annots.rows$time.course <- time.course[match(annots.rows$SYMBOL.ATAC,SYMBOL)]$time.course
      
      t2 <- heatmap.shaps.clust()$enriched
      tab <- cbind(annots.rows,t2)
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.shaps.tsne,".tab")),sep="\t",quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.shaps.tsne,".tab"),"has been written in the table folder.")
  })
  
  output$download.done <- renderText({
    test.tab()
  })
  
  output$tab.shaps.clust.plot <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    g1 <- ggplot(annots.rows,aes(x=Feature)) + geom_bar(aes(fill=Feature),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g2 <- ggplot(annots.rows,aes(x=genomic)) + geom_bar(aes(fill=genomic),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(intergenic="grey",tss=colours()[614],prom=colours()[77],gene.body=colours()[613])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g3 <- ggplot(annots.rows,aes(x=ATAC.status)) + geom_bar(aes(fill=ATAC.status),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(stable="grey",opening="red",closing=colours()[124])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g4 <- ggplot(annots.rows,aes(x=frag.type)) + geom_bar(aes(fill=frag.type),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(b="grey",o="black")) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    multiplot(g1,g2,g3,g4,cols=2)
  })
  
  output$tab.shaps.clust.plot2 <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows <- data.table(annots.rows)[,list(Feature,ATAC.status, heatmap,frag.type,genomic, histones=paste(H3K27me3,H3K4me3,H3K4me1,H3K27Ac,sep="/" ))]
    ggplot(annots.rows,aes(x=histones)) + geom_bar(fill="grey",color="black") + theme_bw() + facet_wrap(~heatmap,ncol=1,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("H3K27me3/H3K4me3/H3K4me1/H3K27Ac")
  })
  
  output$tab.shaps.clust.plot3 <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows <- data.table(annots.rows)[,list(heatmap,H3K27me3,H3K4me3,H3K4me1,H3K27Ac)]
    annots.rows <- melt(annots.rows,id.vars="heatmap")
    annots.rows$value <- factor(as.vector(annots.rows$value),levels=rev(c("00","10","01","11")))
    ggplot(annots.rows,aes(x=variable)) + geom_bar(aes(fill=value),color="black",alpha=0.5) + theme_bw() + scale_fill_manual(values=c("00"="white","10"=colours()[142],"01"=colours()[131],"11"=colours()[613])) + facet_wrap(~heatmap,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("H3K27me3/H3K4me3/H3K4me1/H3K27Ac")
  })
  
  interacting.tSNE <- reactive({
    heatmap.info <- heatmap.shaps.clust()$annots.rows
    brush.data <- brush_info2()
    enriched <- brush.data$enriched
    annots.rows <- brush.data$annots.rows
    homer.tmp <- prepare.shaps()$homer.tmp
    
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    
    shaps.data <- prepare.shaps()
    #all.interactions
    peaks.ATAC.frag <- shaps.data$peaks.ATAC.frag
    diff.fragments <- shaps.data$diff.fragments
    
    tmp <- all.interactions[extend.classif!="stable"]
    
    tmp <- rbind(tmp[,list(id1,id2,type1,type2,name1,name2,type,CHiC=extend.classif)],
                 tmp[,list(id1=id2,id2=id1,type1=type2,type2=type1,name1=name2,name2=name1,type,CHiC=extend.classif)])
    
    m <- as.data.frame(matches(annots.rows$ATAC.id,peaks.ATAC.frag$id,all.y=F))
    annots.rows <- annots.rows[m$x,]
    annots.rows$id.frag <- peaks.ATAC.frag$id.frag[m$y]
    annots.rows$bait.symbol <- peaks.ATAC.frag$bait.symbol[m$y]
    annots.rows$bait.status <- peaks.ATAC.frag$bait.status[m$y]
    
    m1 <- as.data.frame(matches(annots.rows$id.frag,tmp$id1,all.y=F))
    m1$cl <- annots.rows$db.clust[m1$x]
    
    m1$ATAC1 <- annots.rows$ATAC.id[m1$x]
    m1$SYMBOL <- tmp$name1[m1$y]
    m1$Feature1 <- annots.rows$Feature[m1$x]
    
    m1$CHiC.diff <- tmp$CHiC[m1$y]
    m1$id <- tmp$id2[m1$y]
    m1$name2 <- tmp$name2[m1$y]
    m1$bait.status <- diff.fragments$bait.status[match(annots.rows$id.frag[m1$x],diff.fragments$id)]
    
    m1prime <- as.data.frame(matches(m1$id,peaks.ATAC.frag$id.frag,all.y=F))
    m1 <- m1[m1prime$x,]
    m1$ATAC <- peaks.ATAC.frag$id[m1prime$y]
    m1prime <- as.data.frame(matches(m1$id,diff.fragments$id,all.y=F))
    m1 <- m1[m1prime$x,]
    m1$Feature2 <- paste(diff.fragments$diff.int[m1prime$y],diff.fragments$bait.status[m1prime$y],sep="_")
    
    m1$clust <- homer.tmp$db.clust[match(m1$ATAC,homer.tmp$ATAC.id)]
    m1$clust[is.na(m1$ATAC)] <- "no.ATAC"
    on.bait.num <- (1:nrow(annots.rows))[-unique(m1[m1$Feature1==m1$Feature2,]$x)]
    m <- data.table(m1[m1$Feature1==m1$Feature2,])
    m <- m[,list(ATAC.id=ATAC1,cluster1=cl,SYMBOL,Feature=Feature1,interacting.frag=name2,ATAC,cluster2=clust)]
    on.bait <- annots.rows[on.bait.num,]
    on.bait$CHiC <- unlist(sapply(strsplit(on.bait$Feature,"_"),"[[",1))
    on.bait$CHiC <- ifelse(on.bait$CHiC=="gain","nb.gained","nb.lost")
    int.fragments.stats <- per.fragment.interactions.total()
    int.fragments.stats <- int.fragments.stats[match(on.bait$bait.symbol,name)]
    int.fragments.stats <- melt(int.fragments.stats,id.vars="name")
    on.bait$interacting.frag <- paste(on.bait$CHiC,int.fragments.stats$value[int.fragments.stats$variable==on.bait$CHiC],sep=": ")
    m <- rbind(m,data.table(on.bait)[,list(ATAC.id,cluster1=db.clust,SYMBOL=bait.symbol,Feature=bait.status,interacting.frag,ATAC="on.bait",cluster2=db.clust)])
    
    m[is.na(cluster2)]$cluster2 <- "not.in.tSNE"   
    #m[cluster2=="0"]$cluster2 <- "not.in.cluster"
    m$heatmap.info <- heatmap.info$heatmap[match(m$ATAC.id,row.names(heatmap.info))]
    colnames(m) <- c("ATAC.id1","cluster1","fragment.id1","class","fragment.id2","ATAC.id2","cluster2","heatmap.info")
    all.fragments$name <- ifelse(all.fragments$is.bait=="other",all.fragments$IDs,all.fragments$gene)
    gb1 <- match(m$fragment.id1,all.fragments$name)
    gb2 <- match(m$fragment.id2,all.fragments$name)
    m.toprint <- m
    m.toprint$gb1 <- all.fragments$gene.bodies[gb1]
    m.toprint$gb2 <- all.fragments$gene.bodies[gb2]
    #print(colnames(homer.tmp))
    TFs <- rbindlist(apply(homer.tmp,1,function(x,tf){
      TFBS <- paste(tf[1:416][x[1:416]>0],collapse=", ")
      Histones <- paste(tf[c(425:428,432:433,449:450,454:457)][x[c(425:428,432:433,449:450,454:457)]>0],collapse=", ")
      TF <- paste(tf[c(434:448,451:453,458:466)][x[c(434:448,451:453,458:466)]>0],collapse=", ")
      data.table(Histones,TF,TFBS)
    },tf=colnames(homer.tmp)))
    TF1 <- TFs[match(m$ATAC.id1,homer.tmp$ATAC.id)][,list(Histones.1=Histones,TF.1=TF,TFBS.1=TFBS)]
    TF2 <- TFs[match(m$ATAC.id2,homer.tmp$ATAC.id)][,list(Histones.2=Histones,TF.2=TF,TFBS.2=TFBS)]
    m.toprint <- cbind(m.toprint,TF1,TF2)
    m.toprint <- m.toprint[,c("class","cluster1","ATAC.id1","fragment.id1","gb1","ATAC.id2","fragment.id2","cluster2","gb2","Histones.1","TF.1","TFBS.1","Histones.2","TF.2","TFBS.2")]
    type1 <- ifelse(!grepl("^chr",m.toprint$fragment.id1),"b",ifelse(m.toprint$gb1=="no.gene","o","gb"))
    type2 <- ifelse(!grepl("^chr",m.toprint$fragment.id2),"b",ifelse(m.toprint$gb2=="no.gene","o","gb"))
    m.toprint$type <- paste(type1,type2,sep="_")
    message("write file ",file.path(res.dir,"interacting_tSNE_all_info.tab"))
    write.table(m.toprint,file.path(res.dir,"interacting_tSNE_all_info.tab"),sep="\t",quote=F,row.names=F)
    m <- m[,c("heatmap.info","class","cluster1","ATAC.id1","fragment.id1","ATAC.id2","fragment.id2","cluster2")]
    m <- m[!is.na(m$heatmap.info),]
    m <- m[m$class!="stable",]
    m <- data.table(unique(as.data.frame(m)))
  })
    
  output$tab.interacting.tSNE <- renderDataTable({
    interacting.tSNE()
  })
  
  output$all.shaps.clust <- renderPlotly({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    alluvial <- interacting.tSNE()
    noATAC <- nrow(alluvial[cluster2=="no.ATAC"])
    alluvial <- alluvial[ATAC.id2!="on.bait" & cluster2!="no.ATAC"]
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust),"not.in.tSNE")
    color.clust <- color.clust[!is.na(names(color.clust))]
    ggplot(alluvial,aes(x=cluster2)) + ggtitle(paste(noATAC, "interacting fragments do not carry ATAC peaks")) + geom_bar(aes(fill=cluster2),color="black") + xlab("cluster of interacting fragments") + scale_fill_manual(values=c("0"="white",color.clust)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  interacting.heatmap <- reactive({
    data.heatmap <- heatmap.shaps.clust()
    data.interaction <- interacting.tSNE()
    data.interaction <- data.interaction[data.interaction$cluster2!="no.ATAC" & data.interaction$ATAC.id2!="on.bait" ,]
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    #colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    ids <- paste(id1,id2,sep=":")
  
    homer.num$Interaction_Number <- NULL
    data.chip <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
  
    tab1 <- as.data.frame(data.chip[data.chip$ATAC.id%in%id1,])
    row.names(tab1) <- tab1$ATAC.id
    tab1$ATAC.id <- NULL
  
    filter <- apply(tab1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)
    tab1 <- tab1[,filter]
    tab1 <- tab1[,order(apply(tab1,2,function(x) sum(x>0)/length(x)),decreasing=T)]
    
    m1 <- match(id1,row.names(tab1))
    
    tab2 <- as.data.frame(data.chip[data.chip$ATAC.id%in%id2,])
    row.names(tab2) <- tab2$ATAC.id
    tab2$ATAC.id <- NULL
  
    tab2 <- tab2[,apply(tab2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)]
    tab2 <- tab2[,order(apply(tab2,2,function(x) sum(x>0)/length(x)))]
    m2 <- match(id2,row.names(tab2))
  
    tab.interact <- cbind(tab1[m1,],tab2[m2,])
    ids <-  paste(row.names(tab1)[m1],row.names(tab2)[m2],sep=":")
    tab.interact <- tab.interact[!duplicated(ids),]
    row.names(tab.interact) <- ids[!duplicated(ids)]
    annots <- data.frame(id1=row.names(tab1)[m1][!duplicated(ids)],id2=row.names(tab2)[m2][!duplicated(ids)],row.names=row.names(tab.interact))
    cut <- c(frag1=ncol(tab1),frag2=ncol(tab2))
    colnames(tab.interact) <-  unlist(sapply(strsplit(colnames(tab.interact),"[.][.]"),"[[",1))
    clust <-  pheatmap(t(tab.interact),
                       clustering_distance_cols = "binary",
                       clustering_method = "ward.D",silent=T)
    nb.clust <- input$nb.clust.heatmap.interact
    cluster=factor(cutree(clust$tree_col, k = nb.clust))
    annots$heatmap <- cluster[match(row.names(annots),names(cluster))]
    annots <- annots[match(names(cluster),row.names(annots)),]
    #annots$position <- 1:nrow(annots)
    
    list(tab.interact=tab.interact,annots=annots,cut=cut)
  })
  
  genes.and.clusters <- reactive({
    interacting.tSNE <- interacting.tSNE()
    print(interacting.tSNE[fragment.id1=="BCL6" | fragment.id2=="BCL6"])
    fragment1 <- interacting.tSNE[,list(cluster1,ATAC.id1,fragment.id1,class)]
    fragment2 <- interacting.tSNE[,list(cluster2,ATAC.id2,fragment.id2,class)]
    fragment1$type1 <- ifelse(grepl("chr",fragment1$fragment.id1),"o","b")
    fragment1$SYMBOL <- ifelse(fragment1$type1=="b",fragment1$fragment.id1,fragment2$fragment.id2)
    fragment2$type2 <- ifelse(grepl("chr",fragment2$fragment.id2),"o","b")
    fragment2$SYMBOL <- ifelse(fragment2$type2=="b",fragment2$fragment.id2,fragment1$fragment.id1)
    
    switch(input$heatmap.gene.cluster,
           only.baits={
             res <- rbind(fragment1[type1=="b",list(name=fragment.id1,SYMBOL=SYMBOL,ATAC.id=ATAC.id1,cluster=cluster1)],fragment2[type2=="b",list(name=fragment.id2,SYMBOL=SYMBOL,ATAC.id=ATAC.id2,cluster=cluster2)])
             print(res[order(name)])
             res <- unique(res)
             print(res[order(name)])
             res <- res[!cluster%in%c("no.ATAC","0","not.in.tSNE") & !grepl("^nb",SYMBOL),list(counts=length(ATAC.id)),by=c("SYMBOL","cluster")]
             res$SYMBOL <- as.vector(res$SYMBOL)
             res$cluster <- as.vector(res$cluster)
             print(res[SYMBOL=="BCL6"])
             res <- as.data.frame.matrix(xtabs(counts~cluster+SYMBOL,res))
           },
           only.other={
             res <- rbind(fragment1[type1=="o",list(name=fragment.id1,SYMBOL=SYMBOL,ATAC.id=ATAC.id1,cluster=cluster1)],fragment2[type2=="o",list(name=fragment.id2,SYMBOL=SYMBOL,ATAC.id=ATAC.id2,cluster=cluster2)])
             res <- unique(res)
             res <- res[!cluster%in%c("no.ATAC","0","not.in.tSNE") & !grepl("^nb",SYMBOL),list(counts=length(ATAC.id)),by=c("SYMBOL","cluster")]
             res$SYMBOL <- as.vector(res$SYMBOL)
             res$cluster <- as.vector(res$cluster)
             print(res[SYMBOL=="BCL6"])
             res <- as.data.frame.matrix(xtabs(counts~cluster+SYMBOL,res))
           },
           bait.other={
             fragment1$SYMBOL[fragment1$type1=="b" & fragment2$type2=="b"] <- fragment2$SYMBOL[fragment1$type1=="b" & fragment2$type2=="b"] <- paste(as.vector(fragment2$fragment.id2[fragment1$type1=="b" & fragment2$type2=="b"]),as.vector(fragment1$fragment.id1[fragment1$type1=="b" & fragment2$type2=="b"]),sep="/")
             res <- rbind(fragment1[,list(name=fragment.id1,SYMBOL=SYMBOL,ATAC.id=ATAC.id1,cluster=cluster1)],fragment2[,list(name=fragment.id2,SYMBOL=SYMBOL,ATAC.id=ATAC.id2,cluster=cluster2)])
             res <- res[!cluster%in%c("no.ATAC","0","not.in.tSNE") & !grepl("^nb",SYMBOL)]
             res <- unique(res)
             #test <- res[duplicated(SYMBOL) & duplicated(cluster)]$SYMBOL
             #print(res[order(SYMBOL)][SYMBOL%in%test & cluster=="25"])
             res <- res[,list(counts=length(ATAC.id)),by=c("SYMBOL","cluster")]
             #print(res[order(SYMBOL)][SYMBOL%in%test])
             res$SYMBOL <- as.vector(res$SYMBOL)
             res$cluster <- as.vector(res$cluster)
             print(res[SYMBOL=="BCL6"])
             res <- as.data.frame.matrix(xtabs(counts~cluster+SYMBOL,res))
           }
           )
    res 
  })
  
  select.cl.genes <- reactive({
    genes.and.clusters <- genes.and.clusters()
    choices <- unique(unique(row.names(genes.and.clusters)))
    checkboxGroupInput("select.cl.genes","",choices = c("all",choices),selected="all",inline = T)
  })
  
  output$select.cl.genes.ui <- renderUI({
    select.cl.genes()
  })
  
  gene_cluster <- eventReactive(input$load.gene_cluster.tab,{
    genes.and.clusters <- genes.and.clusters()
    print(input$select.cl.genes)
    if (input$select.cl.genes[1]!="all" & length(input$select.cl.genes)>1 & input$include.exclude=="include"){
      genes.and.clusters <- genes.and.clusters[row.names(genes.and.clusters)%in%input$select.cl.genes,]
      genes.and.clusters <- genes.and.clusters[,colSums(genes.and.clusters)>0]
      if (input$and.or.cl=="AND"){
        genes.and.clusters <- genes.and.clusters[,apply(genes.and.clusters,2,function(x) sum(x>0))==length(input$select.cl.genes)]
      }
    }
    else if (input$select.cl.genes[1]!="all" & input$include.exclude=="exclude"){
      genes.and.clusters <- genes.and.clusters[!row.names(genes.and.clusters)%in%input$select.cl.genes,]
      genes.and.clusters <- genes.and.clusters[,colSums(genes.and.clusters)>0]
      
    }
    
    if (!is.null(ranges$x)) {
      data <- genes.and.clusters.heatmap()
      annots.genes <- data$annots.genes
      clust <- data$clust
      genes.and.clusters <- data$genes.and.clusters
      annots.genes$SYMBOL <- row.names(annots.genes)
      annots.genes <- annots.genes[match(colnames(genes.and.clusters)[clust$tree_col$order],annots.genes$SYMBOL),]
      annots.genes$SYMBOL <- colnames(genes.and.clusters)[clust$tree_col$order]
      #annots.genes$status <- "status"
      #annots.genes$SYMBOL <- factor(as.vector(annots.genes$SYMBOL),levels=annots.genes$SYMBOL)
      
      print(ranges$x)
      print(ranges$y)
      keep.genes <- annots.genes$SYMBOL[max(c(0,ranges$x[1])):min(c(nrow(annots.genes),ranges$x[2]))]
      genes.and.clusters <- genes.and.clusters[,keep.genes]
      genes.and.clusters <- genes.and.clusters[rowSums(genes.and.clusters)>0,]
    }
    message("Writing", file.path(res.dir,paste(input$gene_cluster.tab,input$include.exclude,"gene_cluster.tab",sep="_")))
    write.table(genes.and.clusters,file.path(res.dir,paste(input$gene_cluster.tab,input$include.exclude,"gene_cluster.tab",sep="_")),sep="\t")
    text <- paste(file.path(res.dir,paste(input$gene_cluster.tab,input$include.exclude,"gene_cluster.tab",sep="_")),"has been saved")
  })
  
  output$text.gene_cluster.tab <- renderText({
    gene_cluster()
  })
  
  genes.and.clusters.heatmap <- eventReactive(input$plot.genes.and.clusters.heatmap,{
    genes.and.clusters <- genes.and.clusters()
    print(input$select.cl.genes)
    
    if (input$select.cl.genes[1]!="all" & length(input$select.cl.genes)>1 & input$include.exclude=="include"){
      genes.and.clusters <- genes.and.clusters[row.names(genes.and.clusters)%in%input$select.cl.genes,]
      genes.and.clusters <- genes.and.clusters[,colSums(genes.and.clusters)>0]
      if (input$and.or.cl=="AND"){
        genes.and.clusters <- genes.and.clusters[,apply(genes.and.clusters,2,function(x) sum(x>0))==length(input$select.cl.genes)]
      }
    }
    else if (input$select.cl.genes[1]!="all" & input$include.exclude=="exclude"){
      genes.and.clusters <- genes.and.clusters[!row.names(genes.and.clusters)%in%input$select.cl.genes,]
      genes.and.clusters <- genes.and.clusters[,colSums(genes.and.clusters)>0]
    }
    print(dim(genes.and.clusters))
    annots.genes <- data.frame(expr=RNAseq$status,row.names=RNAseq$ID)
    annots.genes$selected <- "0"
    if (input$symbol.to.display!=""){
      genes <- toupper(unlist(strsplit(as.vector(input$symbol.to.display)," ")))
      annots.genes$selected[row.names(annots.genes)%in%genes] <- "1"
    }
    clust <- pheatmap(genes.and.clusters,
             annotation_col = annots.genes,
             annotation_colors = list(expr=c(up="red",down=colours()[124],stable="lightgrey")),
             border_color = NA,
             clustering_distance_cols = "binary",
             color = colorRampPalette(c("white","orange","red","brown","black"))(100),
             breaks = seq(0,max(genes.and.clusters),max(genes.and.clusters)/100))
    return(list(clust=clust,annots.genes=annots.genes,genes.and.clusters=genes.and.clusters))
  })
  
  ranges <- reactiveValues(x=NULL,y=NULL)
  observe({
    brush <- input$zoom_brush
    if (!is.null(brush)) {
      ranges$x <- c(round(brush$xmin), round(brush$xmax))
      ranges$y <- c(round(brush$ymin), round(brush$ymax))
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$annotation.gNc <- renderPlot({
    
    data <- genes.and.clusters.heatmap()
    annots.genes <- data$annots.genes
    clust <- data$clust
    genes.and.clusters <- data$genes.and.clusters
    annots.genes$SYMBOL <- row.names(annots.genes)
    annots.genes <- annots.genes[match(colnames(genes.and.clusters)[clust$tree_col$order],annots.genes$SYMBOL),]
    annots.genes$SYMBOL <- colnames(genes.and.clusters)[clust$tree_col$order]
    annots.genes$status <- "status"
    annots.genes$SYMBOL <- factor(as.vector(annots.genes$SYMBOL),levels=annots.genes$SYMBOL)
    g <- ggplot(annots.genes, aes(SYMBOL, status, fill = expr)) + 
      geom_tile() + xlab("") + ylab("") + theme_light() + theme(axis.text.x=element_blank()) +
      scale_fill_manual(values=c(up="red",down=colours()[124],stable="lightgrey"))
    g
  },  height = 50)
  
  output$zoomplot.gNc <- renderPlot({
    
    if (is.null(ranges$x)) return(NULL)
    
    data <- genes.and.clusters.heatmap()
    annots.genes <- data$annots.genes
    clust <- data$clust
    genes.and.clusters <- data$genes.and.clusters
    annots.genes$SYMBOL <- row.names(annots.genes)
    annots.genes <- annots.genes[match(colnames(genes.and.clusters)[clust$tree_col$order],annots.genes$SYMBOL),]
    annots.genes$SYMBOL <- colnames(genes.and.clusters)[clust$tree_col$order]
    #annots.genes$status <- "status"
    #annots.genes$SYMBOL <- factor(as.vector(annots.genes$SYMBOL),levels=annots.genes$SYMBOL)
    
    print(ranges$x)
    print(ranges$y)
    keep.genes <- annots.genes$SYMBOL[max(c(0,ranges$x[1])):min(c(nrow(annots.genes),ranges$x[2]))]
    genes.and.clusters <- genes.and.clusters[,keep.genes]
    genes.and.clusters <- genes.and.clusters[rowSums(genes.and.clusters)>0,]
    annots.genes <- data$annots.genes
    pheatmap(genes.and.clusters,
             annotation_col = annots.genes,
             annotation_colors = list(expr=c(up="red",down=colours()[124],stable="lightgrey")),
             border_color = NA,
             clustering_distance_cols = "binary",
             color = colorRampPalette(c("white","orange","red","brown","black"))(100),
             breaks = seq(0,max(genes.and.clusters),max(genes.and.clusters)/100),fontsize_col = 8)
  })
  
  output$genes.and.clusters.heatmap <- renderPlot({
    data <- genes.and.clusters.heatmap()
    if (ncol(data$genes.and.clusters)>50){
      showcol <- F
    }
    else{
      showcol <- T
    }
    
    pheatmap(data$genes.and.clusters,
             annotation_col = data$annots.genes,
             annotation_colors = list(expr=c(up="red",down=colours()[124],stable="lightgrey"),selected=c("0"="white","1"="black")),
             border_color = NA, show_colnames = showcol,
             clustering_distance_cols = "binary",
             color = colorRampPalette(c("white","orange","red","brown","black"))(100),
             breaks = seq(0,max(data$genes.and.clusters),max(data$genes.and.clusters)/100))
  })
  
  
    
  output$interacting.heatmap <- renderPlot({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    data <- interacting.heatmap()
    annots <- data$annots
    peaks.ATAC.frag <- prepare.shaps()$peaks.ATAC.frag
    data.interaction <- interacting.tSNE()
    annots$cl1 <- data.interaction$cluster1[match(annots$id1,data.interaction$ATAC.id1)]
    annots$cl2 <- data.interaction$cluster2[match(annots$id2,data.interaction$ATAC.id2)]
    annots$Feature <- data.interaction$class[match(annots$id1,data.interaction$ATAC.id1)]
    annots$frag1 <- peaks.ATAC.frag$frag.type[match(annots$id1,peaks.ATAC.frag$id)]
    annots$frag2 <- peaks.ATAC.frag$frag.type[match(annots$id2,peaks.ATAC.frag$id)]
    annots$status.1 <- peaks.ATAC.frag$status[match(annots$id1,peaks.ATAC.frag$id)]
    annots$status.2 <- peaks.ATAC.frag$status[match(annots$id2,peaks.ATAC.frag$id)]
    annots$genom.1 <- peaks.ATAC.frag$genomic[match(annots$id1,peaks.ATAC.frag$id)]
    annots$genom.2 <- peaks.ATAC.frag$genomic[match(annots$id2,peaks.ATAC.frag$id)]
    
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust),"not.in.tSNE")
    color.clust <- color.clust[!is.na(names(color.clust))]
    color.clust <- color.clust[names(color.clust)%in%c(annots$cl1,annots$cl2)]
    pheatmap(t(data$tab.interact),show_rownames = T,
             clustering_method = "ward.D",
             annotation_col = annots[,rev(c("heatmap","Feature","cl1","cl2","frag1","frag2","genom.1","genom.2","status.1","status.2"))],
             color=colorRampPalette(colours()[c(235,34,35,24)])(100),
             border_color = NA,clustering_distance_cols = "binary",
             show_colnames=F,cluster_rows = F,gaps_row = data$cut[1],
             annotation_colors = list(cl1=color.clust[names(color.clust)%in%annots$cl1],cl2=color.clust[names(color.clust)%in%annots$cl2],genom.1=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),genom.2=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),frag1=c(o="black",b="white"),frag2=c(o="black",b="white"),Feature=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124]),status.1=c(closing=colours()[131],opening="red",stable="grey"),status.2=c(closing=colours()[131],opening="red",stable="grey")))
  })
  
  output$select.on.heatmap.interaction.cluster <- renderUI({
    if(input$select.on.heatmap.interaction){
      data <- interacting.heatmap()
      annots <- data$annots
      checkboxGroupInput("heatmap.interaction.cluster","",choices = c("all",unique(annots$heatmap)),selected = "all",inline = T)
    }
  })
  
  select.on.frag <- reactive({
    data.interaction <- interacting.tSNE()
    #ids <- paste0(data.interaction$cluster1,"_",data.interaction$ATAC.id1,"_",data.interaction$fragment.id1,":",data.interaction$cluster2,"_",data.interaction$ATAC.id2,"_",data.interaction$fragment.id2)
    data.interaction$ATAC.id2[is.na(data.interaction$ATAC.id2)] <- "no.ATAC"
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    ids <- paste(id1,id2,sep=":")
    if(input$select.on.heatmap.interaction){
      data <- interacting.heatmap()
      annots <- data$annots
      tab.interact <- data$tab.interact
      cut.cl <- data$cut
      if (input$heatmap.interaction.cluster[1]!="all"){
       annots <- annots[annots$heatmap%in%input$heatmap.interaction.cluster,]
       tab.interact <- tab.interact[row.names(tab.interact)%in%row.names(annots),]
      }
      frag1 <- tab.interact[,1:cut.cl[1]]
      frag1 <- colnames(frag1)[apply(frag1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)]
      frag2 <- tab.interact[,(cut.cl[1]+1):ncol(tab.interact)]
      frag2 <- c(colnames(frag2)[apply(frag2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)])
    }
    else {
      homer.num$Interaction_Number <- NULL
      homer.num <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
      colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
      
      frag1 <- homer.num[homer.num$ATAC.id%in%id1,]
      frag1 <- colnames(frag1)[apply(frag1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)][-1]
      frag2 <- homer.num[homer.num$ATAC.id%in%id2,]
      frag2 <- c(colnames(frag2)[apply(frag2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)][-1],"no.ATAC")
    }
  list(frag1=frag1,frag2=frag2)
  
  })
  
  output$select.on.frag1 <- renderUI({
    features <- select.on.frag()$frag1
    checkboxGroupInput("frag1.selected","",choices = c("all",features),selected = "all")
  })
  
  output$select.on.frag2 <- renderUI({
    features <- select.on.frag()$frag2
    checkboxGroupInput("frag2.selected","",choices = c("all",features),selected = "all")
  })
  
  tab.to.display <- reactive({
    heatmap.data <- interacting.heatmap()$annots
    data.interaction <- interacting.tSNE()
    data.interaction$ATAC.id2[is.na(data.interaction$ATAC.id2)] <- "no.ATAC"
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    ids <- paste(id1,id2,sep=":")
    
    frag1.f <- select.on.frag()$frag1
    frag2.f <- select.on.frag()$frag2
    if (input$frag1.selected[1]!="all"){
      frag1.f <- input$frag1.selected
    }
    if (input$frag2.selected[1]!="all"){
      frag2.f <- input$frag2.selected
    }
    homer.num$Interaction_Number <- NULL
    homer.num <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
    colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    homer.num <- as.data.frame(homer.num)
    frag1 <- homer.num[match(id1,homer.num$ATAC.id),]
    frag2 <- homer.num[match(id2,homer.num$ATAC.id),]
    frag2$ATAC.id[is.na(frag2$ATAC.id)] <- "no.ATAC"
    if (length(frag1.f)==1){
      keep.frag1 <- frag1[,frag1.f]>0
    }
    else{
      keep.frag1 <- apply(frag1[,frag1.f],1,function(x) sum(x>0))
      keep.frag1 <- ifelse(keep.frag1==length(frag1.f),TRUE,ifelse(keep.frag1==0,FALSE,ifelse(input$OR.AND.frag1=="OR",TRUE,FALSE)))
    }
    
    if ("no.ATAC"%in%frag2.f){
      frag2.f <- frag2.f[frag2.f!="no.ATAC"]
      frag2[is.na(frag2)] <- 1
    }
    else{
      frag2[is.na(frag2)] <- 0
    }
    if (length(frag2.f)==1){
      if (frag2.f=="no.ATAC"){
        keep.frag2 <- keep.frag2$ATAC.id=="no.ATAC"
      }
      else{
        keep.frag2 <- frag2[,frag2.f]>0
      }
    }
    else{
      keep.frag2 <- apply(frag2[,frag2.f],1,function(x) sum(x>0))
      keep.frag2 <- ifelse(keep.frag2==length(frag2.f),TRUE,ifelse(keep.frag2==0,FALSE,ifelse(input$OR.AND.frag2=="OR",TRUE,FALSE)))
    }
    data.interaction$interact.heatmap <- heatmap.data$heatmap[match(ids,row.names(heatmap.data))]
    
    data.interaction <- data.interaction[keep.frag1 & keep.frag2]
    if(input$select.on.heatmap.interaction){
      if (input$heatmap.interaction.cluster[1]!="all"){
        data.interaction <- data.interaction[interact.heatmap%in%input$heatmap.interaction.cluster]
      }
    }
    data.interaction
    
  })
  
  output$selected.interacting.fragments <- renderDataTable({
    tab.to.display()
  })
  
  interacting.tab <- eventReactive(input$go.save.interacting.tab,{
    withProgress(message = 'Dowloading files ...', {
      tab <- tab.to.display()
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.interacting.tab,".tab")),sep="\t",row.names=F,quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.interacting.tab,".tab"),"has been written in the table folder.")
  })
  
  output$download.interacting.done <- renderText({
    interacting.tab()
  })
  
  tab.expression.interacting <- reactive({
    genes <- interacting.tSNE()
    genes <- unique(c(genes$fragment.id1,genes$fragment.id2))
    genes.all <- genes[!grepl("^chr|^nb.",genes)]
    genes <- tab.to.display()
    genes <- unique(c(genes$fragment.id1,genes$fragment.id2))
    genes.selected <- genes[!grepl("^chr|^nb.",genes)]
    RNAseq <- RNAseq[ID%in%genes.all]
    
    time.course <- time.course[,list(ID=SYMBOL,T0.T20,T20.T1,T1.T2 ,T2.T4,T4.T24)]
    expr <- merge(RNAseq[,list(ID,selected=ifelse(ID%in%genes.selected,T,F),RNAseq=status)],time.course,by="ID",all.x=T)
    expr
  })
  
  output$tab.expression.interacting <- renderDataTable({
    tab.expression.interacting()
  })
  
  interacting.tab.expr <- eventReactive(input$go.save.interacting.expr.tab,{
    withProgress(message = 'Dowloading files ...', {
      tab <- tab.expression.interacting()
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.interactingexpr.tab,".tab")),sep="\t",row.names=F,quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.interactingexpr.tab,".tab"),"has been written in the table folder.")
  })
  output$download.interacting.expr.done <- renderText({
    interacting.tab.expr()
  })
  
  other.frag.stats <- reactive({
    other.frag.stats <- all.interactions
    other.frag.stats$bait.expr <- RNAseq$status[match(other.frag.stats$name1,RNAseq$ID)]
    other.frag.stats <- other.frag.stats[type2=="o" & extend.classif!="loss",list(SYMBOLs=paste(name1,collapse=","),nb.bait=length(name1),nb.up=sum(!is.na(bait.expr) & bait.expr=="up"),nb.down=sum(!is.na(bait.expr) & bait.expr=="down"),nb.stable=sum(!is.na(bait.expr) & bait.expr=="stable"),nb.NA=sum(is.na(bait.expr)),nb.w.ATAC=sum(id1%in%peaks.ATAC.frag$id.frag)),by=name2]
    other.frag.stats <- other.frag.stats[nb.bait>1]
    bait.peaks <- all.interactions[type2=="o" & id2%in%other.frag.stats$name2]
    m1 <- matches(bait.peaks$id2,peaks.ATAC.frag$id.frag,all.y=F,all.x=F)
    m1$int.type <- bait.peaks$extend.classif[m1$x]
    m1$o <- bait.peaks$id2[m1$x]
    m1$b <- bait.peaks$id1[m1$x]
    m1$SYMBOL <- bait.peaks$name1[m1$x]
    m1$ATAC.o <- peaks.ATAC.frag$id[m1$y]
    m1$ATAC.o.status <- peaks.ATAC.frag$status[m1$y]
    m2 <- matches(m1$b,peaks.ATAC.frag$id.frag,all.y=F,all.x=F)
    m2$ATAC.o <- m1$ATAC.o[m2$x]
    m2$ATAC.o.status <- m1$ATAC.o.status[m2$x]
    m2$SYMBOL <- m1$SYMBOL[m2$x]
    m2$ATAC.b <- peaks.ATAC.frag$id[m2$y]
    m2$ATAC.b.status <- peaks.ATAC.frag$status[m2$y]
    m2$int.type <- m1$int.type[m2$x]
    m2$bait.status <- RNAseq[match(m2$SYMBOL,ID)]$status
    
    peaks.ATAC.frag$active.marks <- paste(peaks.ATAC.frag$H3K4me3,peaks.ATAC.frag$H3K4me1,peaks.ATAC.frag$H3K27Ac,sep="/")
    interacting.ATAC <- data.table(other=m2$ATAC.o,bait=m2$ATAC.b,int.type=m2$int.type,H3K27me3.o=peaks.ATAC.frag$H3K27me3[match(m2$ATAC.o,peaks.ATAC.frag$id)],H3K4me1.o=peaks.ATAC.frag$H3K4me1[match(m2$ATAC.o,peaks.ATAC.frag$id)],H3K27Ac.o=peaks.ATAC.frag$H3K27Ac[match(m2$ATAC.o,peaks.ATAC.frag$id)],mark.o=peaks.ATAC.frag$active.marks[match(m2$ATAC.o,peaks.ATAC.frag$id)],bait=m2$ATAC.b,ATAC.b.status=m2$ATAC.b.status,mark.b=peaks.ATAC.frag$active.marks[match(m2$ATAC.b,peaks.ATAC.frag$id)],SYMBOL=m2$SYMBOL,bait.status=m2$bait.status)
    interacting.ATAC <- unique(interacting.ATAC)
    
    interacting.ATAC.bait <- interacting.ATAC[,list(count.b=length(bait),SYMBOL=paste(unique(SYMBOL),collapse="/")),by=c("other","mark.o","mark.b","ATAC.b.status","int.type","bait.status")]
    interacting.ATAC.bait <- interacting.ATAC.bait[,list(count.o=length(unique(other)),count.b=sum(count.b),paste(unique(unlist(strsplit(SYMBOL,"/"))),collapse="/"),nb.bait=length(unique(unlist(strsplit(SYMBOL,"/"))))),by=c("mark.o","mark.b","ATAC.b.status","int.type","bait.status")]
    
    interacting.ATAC.bait.tab <- as.data.frame.matrix(xtabs(count.b~mark.b+mark.o,interacting.ATAC.bait))
    interacting.ATAC.bait.tab <- interacting.ATAC.bait.tab[rowSums(interacting.ATAC.bait.tab)>10, colSums(interacting.ATAC.bait.tab)>10]
    
    interacting.ATAC.other.tab <- as.data.frame.matrix(xtabs(count.o~mark.b+mark.o,interacting.ATAC.bait))
    interacting.ATAC.other.tab <- interacting.ATAC.other.tab[rowSums(interacting.ATAC.other.tab)>10, colSums(interacting.ATAC.other.tab)>10]
    
    pheatmap(interacting.ATAC.other.tab,breaks=seq(0,50,0.5))
    pheatmap(interacting.ATAC.bait.tab,breaks=seq(0,50,0.5))
    
    annots.other <- unique(interacting.ATAC[,list(other,mark.o)])
    annots.other <- annots.other[,list(count=length(other)),by=mark.o]
    annots.other <- rbind(annots.other,interacting.ATAC.bait[,list(mark.o=paste(mark.o,mark.b,sep="_"),count=count.b)])
    #interacting.ATAC.bait[,mark.b:=paste(mark.o,mark.b,sep="_")]
    data <- prepare.shaps()
    homer.tmp <- data$homer.tmp
    peaks.ATAC.frag <- data$peaks.ATAC.frag
    bait.peaks <- peaks.ATAC.frag[peaks.ATAC.frag$id.frag%in%bait.peaks$id1]
    #return(list(homer.tmp=homer.tmp,shaps=shaps,tsne.group=tsne.group,tsne.center=tsne.center,peaks.ATAC.frag=peaks.ATAC.frag,diff.fragments=diff.fragments))
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust))
    color.clust <- color.clust[!is.na(names(color.clust))]
    if (!is.null(input$choose.chrom)){
      chr <- input$choose.chrom
      other.frag.stats <- other.frag.stats[grep(chr,name2)]
    }
    other.frag.select <- other.frag.stats[name2 %in% peaks.ATAC.frag[peaks.ATAC.frag$status=="opening"]$id.frag]
    selected.ATAC <- peaks.ATAC.frag[peaks.ATAC.frag$id.frag %in% other.frag.select$name2 & peaks.ATAC.frag$status == "opening"]
    list(selected.ATAC=selected.ATAC,other.frag.select=other.frag.select,peaks.ATAC.frag=peaks.ATAC.frag,bait.peaks=bait.peaks)
  })
  
  many.baits <- reactive({
    data <- other.frag.stats()
    ATAC.all <- unique(as.data.frame(peaks.ATAC)[,c("id","H3K4me1","H3K4me3","H3K27Ac","status")])
    all.marks.comb <- unique(as.vector(paste(ATAC.all$H3K4me1,ATAC.all$H3K4me3,ATAC.all$H3K27Ac,sep="/")))
    ATAC.all$all.marks <- paste(ATAC.all$H3K4me3,ATAC.all$H3K4me1,ATAC.all$H3K27Ac,sep="/")
    stat.whole.ATAC <- rbindlist(lapply(all.marks.comb,function(m,ATAC.all){
      res <- data.table(mark=m,nb.tot=sum(ATAC.all$all.marks==m),nb.no.tot=sum(ATAC.all$all.marks!=m),nb.o=sum(ATAC.all$all.marks==m & ATAC.all$status=="opening"),nb.no.o=sum(ATAC.all$all.marks!=m & ATAC.all$status=="opening"),nb.c=sum(ATAC.all$all.marks==m & ATAC.all$status=="closing"),nb.no.c=sum(ATAC.all$all.marks!=m & ATAC.all$status=="closing"),nb.s=sum(ATAC.all$all.marks==m & ATAC.all$status=="stable"),nb.no.s=sum(ATAC.all$all.marks!=m & ATAC.all$status=="stable"))
      res$pc.tot <- round(res$nb.tot / (res$nb.tot+res$nb.no.tot),2)
      res$pc.o <- round(res$nb.o / (res$nb.o+res$nb.no.o),2)
      res$pc.c <- round(res$nb.c / (res$nb.c+res$nb.no.c),2)
      res$pc.s <- round(res$nb.s / (res$nb.s+res$nb.no.s),2)
      res
    },ATAC.all))
    
    ATAC <- unique(as.data.frame(data$selected.ATAC)[,c("id","H3K4me1","H3K4me3","H3K27Ac","status","all.marks")])
    ATAC$all.marks <- paste(ATAC$H3K4me3,ATAC$H3K4me1,ATAC$H3K27Ac,sep="/")
    stat.selected.ATAC <- rbindlist(lapply(all.marks.comb,function(m,ATAC.all){
      res <- data.table(mark=m,nb.o=sum(ATAC.all$all.marks==m & ATAC.all$status=="opening"),nb.no.o=sum(ATAC.all$all.marks!=m & ATAC.all$status=="opening"))
      res$pc.o <- round(res$nb.o / (res$nb.o+res$nb.no.o),2)
      res
    },ATAC))
    
    ATAC.baits <- unique(as.data.frame(data$bait.peaks)[,c("id","H3K4me1","H3K4me3","H3K27Ac","status","all.marks")])
    ATAC.baits$all.marks <- paste(ATAC.baits$H3K4me3,ATAC.baits$H3K4me1,ATAC.baits$H3K27Ac,sep="/")
    stat.bait.ATAC <- rbindlist(lapply(all.marks.comb,function(m,ATAC.all){
      res <- data.table(mark=m,nb.tot=sum(ATAC.all$all.marks==m),nb.no.tot=sum(ATAC.all$all.marks!=m),nb.o=sum(ATAC.all$all.marks==m & ATAC.all$status=="opening"),nb.no.o=sum(ATAC.all$all.marks!=m & ATAC.all$status=="opening"),nb.c=sum(ATAC.all$all.marks==m & ATAC.all$status=="closing"),nb.no.c=sum(ATAC.all$all.marks!=m & ATAC.all$status=="closing"),nb.s=sum(ATAC.all$all.marks==m & ATAC.all$status=="stable"),nb.no.s=sum(ATAC.all$all.marks!=m & ATAC.all$status=="stable"))
      res$pc.tot <- round(res$nb.tot / (res$nb.tot+res$nb.no.tot),2)
      res$pc.o <- round(res$nb.o / (res$nb.o+res$nb.no.o),2)
      res$pc.c <- round(res$nb.c / (res$nb.c+res$nb.no.c),2)
      res$pc.s <- round(res$nb.s / (res$nb.s+res$nb.no.s),2)
      res
    },ATAC.baits))
    stats.fisher.bait <- rbindlist(lapply(stat.bait.ATAC$mark,function(m,stat.selected.ATAC,stat.whole.ATAC){
      test.w <- fisher.test(matrix(c(stat.selected.ATAC[mark==m]$nb.tot,stat.selected.ATAC[mark==m]$nb.no.tot,stat.whole.ATAC[mark==m]$nb.tot,stat.whole.ATAC[mark==m]$nb.no.tot),ncol=2),alternative = "two.sided")
      data.table(mark=m,nb.bait=stat.selected.ATAC[mark==m]$nb.tot,nb.no.bait=stat.selected.ATAC[mark==m]$nb.no.tot,nb.tot=stat.whole.ATAC[mark==m]$nb.tot,nb.no.tot=stat.whole.ATAC[mark==m]$nb.no.tot,pc.bait=stat.selected.ATAC[mark==m]$pc.tot,pc.whole=stat.whole.ATAC[mark==m]$pc.tot,test.w=test.w$p.value)
    },stat.bait.ATAC,stat.whole.ATAC))
    
    stats.fisher.bait.other <- rbindlist(lapply(stat.bait.ATAC$mark,function(m,stat.selected.ATAC,stat.whole.ATAC){
      test.w <- fisher.test(matrix(c(stat.selected.ATAC[mark==m]$nb.tot,stat.selected.ATAC[mark==m]$nb.no.tot,stat.whole.ATAC[mark==m]$nb.o,stat.whole.ATAC[mark==m]$nb.no.o),ncol=2),alternative = "two.sided")
      data.table(mark=m,nb.bait=stat.selected.ATAC[mark==m]$nb.tot,nb.no.bait=stat.selected.ATAC[mark==m]$nb.no.tot,nb.o=stat.whole.ATAC[mark==m]$nb.o,nb.no.o=stat.whole.ATAC[mark==m]$nb.no.o,pc.bait=stat.selected.ATAC[mark==m]$pc.tot,pc.o=stat.whole.ATAC[mark==m]$pc.o,test.w=test.w$p.value)
    },stat.bait.ATAC,stat.selected.ATAC))
    
    
    stats.fisher <- rbindlist(lapply(stat.selected.ATAC$mark,function(m,stat.selected.ATAC,stat.whole.ATAC){
      test.w <- fisher.test(matrix(c(stat.selected.ATAC[mark==m]$nb.o,stat.selected.ATAC[mark==m]$nb.no.o,stat.whole.ATAC[mark==m]$nb.tot,stat.whole.ATAC[mark==m]$nb.no.tot),ncol=2),alternative = "two.sided")
      test.o <- fisher.test(matrix(c(stat.selected.ATAC[mark==m]$nb.o,stat.selected.ATAC[mark==m]$nb.no.o,stat.whole.ATAC[mark==m]$nb.o,stat.whole.ATAC[mark==m]$nb.no.o),ncol=2),alternative = "two.sided")
      data.table(mark=m,nb.select=stat.selected.ATAC[mark==m]$nb.o,nb.no.select=stat.selected.ATAC[mark==m]$nb.no.o,nb.tot=stat.whole.ATAC[mark==m]$nb.tot,nb.no.tot=stat.whole.ATAC[mark==m]$nb.no.tot,nb.o=stat.whole.ATAC[mark==m]$nb.o,nb.no.o=stat.whole.ATAC[mark==m]$nb.no.o,pc.selected=stat.selected.ATAC[mark==m]$pc.o,pc.whole=stat.whole.ATAC[mark==m]$pc.tot,pc.o=stat.whole.ATAC[mark==m]$pc.o,test.w=test.w$p.value,test.o=test.o$p.value)
    },stat.selected.ATAC,stat.whole.ATAC))
    
    stats.fisher.whole <- rbindlist(lapply(stat.whole.ATAC$mark,function(m,stat.whole.ATAC){
      tmp <- stat.whole.ATAC[mark==m]
      test.o <- fisher.test(matrix(c(tmp$nb.o,tmp$nb.no.o,tmp$nb.tot,tmp$nb.no.tot),ncol=2),alternative="two.sided")
      test.c <- fisher.test(matrix(c(tmp$nb.c,tmp$nb.no.c,tmp$nb.tot,tmp$nb.no.tot),ncol=2),alternative="two.sided")
      test.s <- fisher.test(matrix(c(tmp$nb.s,tmp$nb.no.s,tmp$nb.tot,tmp$nb.no.tot),ncol=2),alternative="two.sided")
      tmp$test.o <- test.o$p.value
      tmp$test.c <- test.c$p.value
      tmp$test.s <- test.s$p.value
      tmp
    },stat.whole.ATAC))
    
    merged.pc <- merge(stat.bait.ATAC[,list(mark,pc.bait=pc.tot)],stats.fisher[,list(mark,pc.selected,pc.whole)],by="mark",all=T)
    merged.pc <- as.data.frame(merged.pc)
    row.names(merged.pc) <- merged.pc$mark
    merged.pc$mark <- NULL
    stats.fisher.bait[test.w==0]$test.w <- 10^-300
    merged.fisher <- merge(stats.fisher.bait[,list(mark,bait=ifelse(pc.bait<pc.whole,log10(test.w),-log10(test.w)))],stats.fisher[,list(mark,selected=ifelse(pc.selected<pc.whole,log10(test.w),-log10(test.w)))],by="mark",all=T)
    merged.fisher <- as.data.frame(merged.fisher)
    row.names(merged.fisher) <- merged.fisher$mark
    merged.fisher$mark <- NULL
    merged.fisher <- merged.fisher[match(row.names(merged.pc),row.names(merged.fisher)),]
    
    
    pheatmap(merged.fisher,annotation_row = merged.pc,color=colorRampPalette(c(colours()[c(125,124,590)],"white",colours()[c(420,34,36)]))(100),breaks=seq(-40,40,0.8))
    
    other.frag.select <- data$other.frag.select
    other.frag.tab <- as.data.frame(other.frag.select)
    m <- matches(other.frag.tab$name2,data$selected.ATAC$id.frag)
    other.frag.tab <- other.frag.tab[m$x,]
    other.frag.tab$ATAC.id <- data$selected.ATAC$id[m$y]
    other.frag.tab <- unique(other.frag.tab)
    other.frag.tab <- other.frag.tab[!duplicated(other.frag.tab$ATAC.id),]
    row.names(other.frag.tab) <- other.frag.tab$ATAC.id
    
    other.frag.tab <- other.frag.tab[,1:5]
    row.names(ATAC) <- ATAC$id
    ATAC <- ATAC[,-1]
    ATAC <- ATAC[,1:3]
    ATAC[ATAC=="00"] <- 0
    ATAC[ATAC=="10"] <- -1
    ATAC[ATAC=="01"] <- 1
    ATAC[ATAC=="11"] <- 2
    ATAC.num <- sapply(ATAC,as.numeric)
    row.names(ATAC.num) <- row.names(ATAC)
    pheatmap(ATAC.num[,c(2,1,3)],show_rownames = F,clustering_method = "ward.D2",cluster_cols = F,annotation_row = other.frag.tab)
  })
  #################################################################################################
  ## chromatine
  #################################################################################################
  
  output$emissions <- renderPlot({
    pheatmap(emissions,cluster_rows = F,cluster_cols = F,display_numbers=T,color = colorRampPalette(c("white",colours()[124],colours()[131],"black"))(100))
  })
  output$emissions2 <- renderPlot({
    pheatmap(emissions,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white",colours()[124],colours()[131],"black"))(100))
  })
  output$transition <- renderPlot({
    pheatmap(transitions,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white",colours()[124],colours()[131],"black"))(100))
  })
  
  output$ATAC.enrichments.G0 <- renderPlot({
    if (input$enr.chromHMM=="ATAC"){
      pheatmap(chromHMM.G0$ATAC[-16,grep("cl[0-9]",colnames(chromHMM.G0$ATAC))],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100),breaks=seq(0,200,2))
    }
    else if (input$enr.chromHMM=="feature"){
      pheatmap(chromHMM.G0$feature[-16,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100))
    }
  })
  output$ATAC.enrichments.G1 <- renderPlot({
    if (input$enr.chromHMM=="ATAC"){
     pheatmap(chromHMM.G1$ATAC[-16,grep("cl[0-9]",colnames(chromHMM.G1$ATAC))],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100),breaks=seq(0,200,2))
    }
    else if (input$enr.chromHMM=="feature"){
      pheatmap(chromHMM.G1$feature[-16,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100))
    }
  })
  
  output$ATAC.enrichments <- renderPlot({
    if (input$enr.chromHMM=="ATAC"){
      mat <- chromHMM.G1$ATAC[-16,grep("cl[0-9]",colnames(chromHMM.G1$ATAC))]-chromHMM.G0$ATAC[-16,grep("cl[0-9]",colnames(chromHMM.G0$ATAC))]
      pheatmap(mat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c(colours()[131],colours()[124],"white","red","brown"))(100),breaks=seq(-200,200,4))
    }
    else if (input$enr.chromHMM=="feature"){
      mat <- chromHMM.G1$feature[-16,]-chromHMM.G0$feature[-16,]
      pheatmap(mat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c(colours()[131],colours()[124],"white","red","brown"))(100),breaks=seq(-50,50,1))
    }
  })
  
  output$ATAC.cl.G0 <- renderPlot({
    cl <- input$cl.chromHMM
    mat <- chromHMM.G0[[cl]]
    mat.tmp <- mat[-1,-1]
    row.names(mat.tmp) <- paste0("E",1:15)
    colnames(mat.tmp) <- mat[1,-1]
    pheatmap(mat.tmp,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100),breaks=seq(0,80,0.8))
  })
  output$ATAC.cl.G1 <- renderPlot({
    cl <- input$cl.chromHMM
    mat <- chromHMM.G1[[cl]]
    mat.tmp <- mat[-1,-1]
    row.names(mat.tmp) <- paste0("E",1:15)
    colnames(mat.tmp) <- mat[1,-1]
    pheatmap(mat.tmp,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("white","red","brown","black"))(100),breaks=seq(0,80,0.8))
  })
  output$ATAC.cl <- renderPlot({
    cl <- input$cl.chromHMM
    mat <- chromHMM.G1[[cl]]
    mat.tmp <- mat[-1,-1]
    row.names(mat.tmp) <- paste0("E",1:15)
    colnames(mat.tmp) <- mat[1,-1]
    mat2 <- chromHMM.G0[[cl]]
    mat2.tmp <- mat2[-1,-1]
    row.names(mat2.tmp) <- paste0("E",1:15)
    colnames(mat2.tmp) <- mat2[1,-1]
    mat.tmp <- mat.tmp - mat2.tmp
    pheatmap(mat.tmp,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c(colours()[131],colours()[124],"white","red","brown"))(100),breaks=seq(-50,50,1))
  })
  output$select.feature.chromHMM <- renderUI({
    selectInput("cl.chromHMM","Select ATAC peaks to visualize",choices = names(chromHMM.G1)[!names(chromHMM.G1)%in%c("ATAC","feature")])
  })
  
  output$show.profile.chromHMM <- renderUI({
    switch(input$select.profile.group,
           ATAC={
             selectInput("ATAC.chromHMM","",choices=c("opening","closing","stable","multi.baits"))
           },
           clusters.tSNE={
             selectInput("clusters.chromHMM","",choices=c(names(chromHMM.G1)[grep("^cl",names(chromHMM.G1))],"selection.tSNE"))
           }
           )
    
  })
  
  segmentation_neighbours <- reactive({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    #segments.G1,segments.G0
    peaks.ATAC.resized <- resize(peaks.ATAC,width=2000,fix="center")
    group <- input$select.profile.group
    switch(group,
           ATAC={
             mat.stats.G0 <- chromHMM.lines[[paste0(input$ATAC.chromHMM,".G0")]]
             mat.stats.G1 <- chromHMM.lines[[paste0(input$ATAC.chromHMM,".G1")]]
           },
           clusters.tSNE={
             if (input$clusters.chromHMM=="selection.tSNE"){
               heatmap.shaps.clust <- heatmap.shaps.clust()$annots.rows
               chromHMM.lines$selection.tSNE$heatmap <- heatmap.shaps.clust$heatmap[match(chromHMM.lines$selection.tSNE$ATAC.id,row.names(heatmap.shaps.clust))]
               mat.stats.G0 <- chromHMM.lines$selection.tSNE[ATAC.id %in% row.names(heatmap.shaps.clust),list(count=length(status)),by=c("POS","G0","heatmap")]
               mat.stats.G1 <- chromHMM.lines$selection.tSNE[ATAC.id %in% row.names(heatmap.shaps.clust),list(count=length(status)),by=c("POS","G1","heatmap")]
               
             }
             else {
               mat.stats.G0 <- chromHMM.lines[[paste0(input$clusters.chromHMM,".G0")]]
               mat.stats.G1 <- chromHMM.lines[[paste0(input$clusters.chromHMM,".G1")]]
             }
           }
    )
    list(mat.stats.G0=mat.stats.G0,mat.stats.G1=mat.stats.G1)
  })
  
  segmentation_peaks <- reactive({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    #segments.G1,segments.G0
    group <- input$select.profile.group
    switch(group,
           ATAC={
             switch(input$ATAC.chromHMM,
                    stable={p <- peaks.ATAC[peaks.ATAC$status=="stable"]},
                    opening={p <- peaks.ATAC[peaks.ATAC$status=="opening"]},
                    closing={p <- peaks.ATAC[peaks.ATAC$status=="closing"]},
                    multi.baits={p <- peaks.ATAC[peaks.ATAC$id%in%other.frag.tab$ATAC.id]}
             )
           },
           clusters.tSNE={
             if (input$clusters.chromHMM=="selection.tSNE"){
               heatmap.shaps.clust <- heatmap.shaps.clust()$annots.rows
               p <- peaks.ATAC[peaks.ATAC$id%in%row.names(heatmap.shaps.clust)]
             }
             else {
               p <- peaks.ATAC[peaks.ATAC$id%in%homer.tmp$ATAC.id[homer.tmp$db.clust%in%c(input$clusters.chromHMM,sub("cl","",input$clusters.chromHMM))]]
             }
             
           }
    )
    ovl.G0 <- data.table(as.data.frame(findOverlaps(p,segments.G0)))
    ovl.G0$state <- segments.G0$state[ovl.G0$subjectHits]
    p.prime <- p[ovl.G0$queryHits]
    segments.G0.prime <-segments.G0[ovl.G0$subjectHits]
    
    ovl.G0$dist <- ifelse(start(p.prime)>start(segments.G0.prime),
                                end(segments.G0.prime) - start(p.prime),
                                end(p.prime) - start(segments.G0.prime))
    ovl.G0 <- ovl.G0[,list(state=state[which.max(dist)]),by=queryHits]
    
    ovl.G1 <- data.table(as.data.frame(findOverlaps(p,segments.G1)))
    ovl.G1$state <- segments.G1$state[ovl.G1$subjectHits]
    p.prime <- p[ovl.G1$queryHits]
    segments.G1.prime <-segments.G1[ovl.G1$subjectHits]
    
    ovl.G1$dist <- ifelse(start(p.prime)>start(segments.G1.prime),
                          end(segments.G1.prime) - start(p.prime),
                          end(p.prime) - start(segments.G1.prime))
    ovl.G1 <- ovl.G1[,list(state=state[which.max(dist)]),by=queryHits]
    trans <- data.table(ATAC.id=p$id,status=p$status,G0=ovl.G0$state,G1=ovl.G1$state,tSNE=homer.tmp$db.clust[match(p$id,homer.tmp$ATAC.id)],bait.multi=other.frag.tab$SYMBOLs[match(p$id,other.frag.tab$ATAC.id)])
  })
  
  output$dt.transition <- renderDataTable({
    segmentation_peaks()
  })
  
  save.chromHMM.tab <- eventReactive(input$save.chromHMM.tab,{
    trans <- segmentation_peaks()
    
    group <- input$select.profile.group
    switch(group,
           ATAC={ 
             tmp <- input$ATAC.chromHMM
           },
           clusters.tSNE={
             tmp <- input$clusters.chromHMM
           }
    )
    
    f <- paste("save.chromHMM_",tmp,".tab",sep="")
    write.table(trans,file.path(data.dir,f),sep="\t",quote=F,row.names=F)
    f
  })
  
  output$text.save.chromHMM.tab <- renderText({
    f <- save.chromHMM.tab()
    text <- paste(f, "has been saved.")
  })
  
  output$plot.transition.profile <- renderPlot({
    trans <- segmentation_peaks()
    trans$G0 <- factor(as.vector(trans$G0),levels=paste0("E",1:15))
    trans$G1 <- factor(as.vector(trans$G1),levels=paste0("E",1:15))
    trans <- trans[order(G1)][order(G0)]
    E.colors <- colours()[c(35,40,59,34,463,53,535,634,147,144,614,88,124,210,232)]
    names(E.colors) <- paste0("E",1:15)
    g0 <- ggplot(trans,aes(x=G0)) + geom_bar(aes(fill=G1),colour="black",size=0.1) + theme_bw() +
      scale_fill_manual(values=E.colors) 
    g1 <- ggplot(trans,aes(x=G1)) + geom_bar(aes(fill=G0),colour="black",size=0.1) + theme_bw() +
      scale_fill_manual(values=E.colors)
    
    multiplot(g0,g1,cols = 2)
  })
  
  output$segmentation.neighbours.plot <- renderPlot({
    E.colors <- colours()[c(35,40,59,34,463,53,535,634,147,144,614,88,124,210,232)]
    names(E.colors) <- paste0("E",1:15)
    mat.stats.G0 <- segmentation_neighbours()$mat.stats.G0
    mat.stats.G1 <- segmentation_neighbours()$mat.stats.G1
    mat.stats.G0$G0 <- factor(as.vector( mat.stats.G0$G0),levels=paste0("E",1:15))
    mat.stats.G1$G1 <- factor(as.vector( mat.stats.G1$G1),levels=paste0("E",1:15))
    n <- sum(mat.stats.G0[POS=="0"]$count)
    mat.stats.G0$count <- mat.stats.G0$count/n
    mat.stats.G1$count <- mat.stats.G1$count/n
    g1 <- ggplot(mat.stats.G0, aes(x=as.numeric(as.vector(POS)),y=count,group=G0)) + 
      geom_smooth(aes(colour=G0),method="loess",span=0.1,fill=NA,size=0.5) + theme_bw() +
      scale_colour_manual(values=E.colors) + 
      labs(x=paste0("distance to peak center (",n,")"),y="percentage in G0") #+ 
      #ylim(c(0,max(mat.stats.G0$count,mat.stats.G1$count)))
    g2 <- ggplot(mat.stats.G1, aes(x=as.numeric(as.vector(POS)),y=count,group=G1)) + 
      geom_smooth(aes(colour=G1),method="loess",span=0.1,fill=NA,size=0.5) + theme_bw() +
      scale_colour_manual(values=E.colors) + 
      labs(x=paste0("distance to peak center (",n,")"),y="percentage in G1") #+ 
      #ylim(c(0,max(mat.stats.G0$count,mat.stats.G1$count)))
    
    if (input$select.profile.group=="clusters.tSNE" & input$clusters.chromHMM=="selection.tSNE"){
       g1 <- g1 + facet_wrap(~heatmap,ncol = 1,scales="free")
       g2 <- g2 + facet_wrap(~heatmap,ncol = 1,scales="free")
    }
    
    multiplot(g1,g2,cols = 2)
  })
  
  
})

