

go_function <- function(v.gns, ontology = "BP"){
  
  df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$acc. %in% v.gns)
  n.genes <- length(unique(df.GO.annot$acc.))
  
  df.GO.annot
}




go_gsea <- function(v.gns, gn.pop = NULL, th = 0.05, title = "", bg.mode = "",
                    mode = "enrichment", n.genome = 27655){
  
  # F, C, P
  df <- read.table(paste(folder_data, "ATH_GO_GOSLIM.txt", sep = "/"), fill = T, sep = "\t")
  ec <- c("HTP", "HDA", "HMP", "HGI", "HEP","EXP","IDA","IPI","IMP","IGI","IEP")
  df <- subset(df, df$V10 %in% ec)
  df <- subset(df, df$V8 %in% c("F","P"))
  
  df.GO.annot <- df[,c(1,5,8)]
  names(df.GO.annot) <- c("acc.", "Term", "Ontology")
  
  if(is.null(gn.pop)){
    gn.pop <- unique(df.GO.annot$acc.)
  }
  
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$acc. %in% gn.pop)
  n.genes <- length(unique(df.GO.annot$acc.))

  df.GO <- subset(df.GO.annot, df.GO.annot$acc. %in% v.gns)
  
  df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
  df.tb.annot.set <- as.data.frame(table(df.GO$Term[df.GO$Term != ""]), stringsAsFactors = FALSE)
  df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
  names(df.tb.annot.set) <- c("term", "in_sample", "in_population")
  n.genes.sset <- length(unique(df.GO$acc.))
  
  n.genes.sset <- length(v.gns) # use this always
  
  # gibberellin mediated signaling pathway
  if(bg.mode == "genome"){
    n.genes <- n.genome # used for vs genome analysis
  }else{
    n.genes <- length(gn.pop)
  }
  

  ## filter go - evidence codes, biological process
  df.BP_enrichment <- data.frame(BP = character(nrow(df.tb.annot.set)), p.val = numeric(nrow(df.tb.annot.set)),
                                 percent = numeric(nrow(df.tb.annot.set)), genes = numeric(nrow(df.tb.annot.set)), 
                                 foldchange = numeric(nrow(df.tb.annot.set)), stringsAsFactors = FALSE)
  
  
  
  for(j in 1:nrow(df.tb.annot.set)){  
    
    n.inset <- df.tb.annot.set$in_sample[j]
    n.genomewide <- df.tb.annot.set$in_population[j]
    
    ### global enrichment test - gene basis
    hitInSample <- n.inset
    sampleSize <- n.genes.sset
    hitInPop <- n.genomewide #sum(tb.rate_limiting_domains$Freq)
    failInPop <- n.genes - hitInPop #(nrow(df.global.domains) - hitInPop)
    
    if(mode == "enrichment"){ # enrichment
      p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    }else if(mode == "depletion"){
      p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= TRUE)
    }
    
    #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
    foldChange <- (hitInSample / sampleSize) / (hitInPop / failInPop)
    #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
    #   p.val <- fisher.test(mat.count)$p.value
    df.BP_enrichment$BP[j] <- df.tb.annot.set$term[j]
    df.BP_enrichment$p.val[j] <- p.val
    df.BP_enrichment$percent[j] <- (n.inset / n.genes.sset)
    df.BP_enrichment$genes[j] = n.inset
    
    df.BP_enrichment$foldchange[j] = foldChange
    
  }
  
  df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$p.val),]
  # df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$foldchange),]
  df.BP_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= th)
  
  # write.csv(df.BP_enrichment, "output/GO_AtBZR1_ChIP_1101_unique_genes.csv", row.names = FALSE)
  
  ### make plot ### 
  
  # vec.significance <- character()
  # for(i in 1:nrow(df.BP_enrichment)){
  #   if(df.BP_enrichment$p.val[i] <= 0.05 & df.BP_enrichment$p.val[i] > 0.01){
  #     vec.significance <- c(vec.significance, "*")
  #   }else if(df.BP_enrichment$p.val[i] <= 0.01){
  #     vec.significance <- c(vec.significance, "**")
  #   }else{
  #     vec.significance <- c(vec.significance, "")  
  #   }
  # }
  # 
  
  # https://cran.r-project.org/web/packages/pathfindR/vignettes/visualization_vignette.html
  data <- df.BP_enrichment
  
  # Most basic bubble plot
  data <- df.BP_enrichment
  data$p.val <- -log10(data$p.val)
  
  data %>%
    arrange(desc(p.val)) %>%
    mutate(BP = factor(BP, BP)) %>%
    ggplot(aes(x=BP, y=p.val, size=genes)) + ### , color=continent)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(2, 10), name="Number genes") + 
    coord_flip() +
    xlab("") +
    ylab("Adjusted P value, -log10") +
    theme_classic() + 
    theme(axis.text=element_text(size=9), axis.title=element_text(size=9,face="italic"), legend.title=element_text(size=8, face="italic"))
  
  
  # over- and under representation
  

  
  # library(ggplot2)
  # p1 <- (ggplot(df.BP_enrichment, aes(x=df.BP_enrichment$BP, y=(df.BP_enrichment$foldchange))) +
  #          geom_bar(colour="black", fill="steelblue", stat="identity") +
  #          coord_flip() +
  #          #coord_trans(y="sqrt") +
  #          guides(fill=FALSE) +
  #          xlab("") +
  #          geom_text(aes(label=vec.significance), vjust=0.6, colour = "black", face = "bold", size = 9) +
  #          ylab("Fold Change") +
  #          #ylab("# of Regulated Genes (Shared Genes Removed)") +
  #          theme(axis.text=element_text(size=9, colour = "black"), axis.title=element_text(size=9)) + 
  #          ggtitle(title) +
  #          theme(plot.title = element_text(lineheight=1.0, face="bold", size = 9)) +
  #          theme_minimal()
  #        
  # ) 
  # 
  # p1  
  
}




go_gsea_separate <- function(v.gns, th = 0.05, ontology = "BP", title = ""){
  
  df.GO.annot <- readRDS(paste(folder_data, "Athaliana_167.annot.2017.rds", sep = "/"))
  n.genes <- length(unique(df.GO.annot$acc.))
  
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
  n.genes <- length(unique(df.GO.annot$acc.))
  
  # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
  df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% v.gns)
  
  df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
  df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
  
  df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
  
  ## filter go - evidence codes, biological process
  df.BP_enrichment <- data.frame(BP = character(nrow(df.tb.annot.set)), p.val = numeric(nrow(df.tb.annot.set)),
                                 percent = numeric(nrow(df.tb.annot.set)), genes = numeric(nrow(df.tb.annot.set)), 
                                 foldchange = numeric(nrow(df.tb.annot.set)), stringsAsFactors = FALSE)
  n.genes.sset <- length(unique(df.GO.annot.set$acc.))
  
  
  for(j in 1:nrow(df.tb.annot.set)){  
    
    n.inset <- df.tb.annot.set$Freq.x[j]
    n.genomewide <- df.tb.annot.set$Freq.y[j]
    
    #n.genomewide <- df.tb.annot$Freq[df.tb.annot$Var1 == df.tb.annot$Var1[j]]
    
    ### global enrichment test - gene basis
    hitInSample <- n.inset
    sampleSize <- n.genes.sset
    hitInPop <- n.genomewide #sum(tb.rate_limiting_domains$Freq)
    failInPop <- n.genes - hitInPop #(nrow(df.global.domains) - hitInPop)
    p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    
    #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
    foldChange <- (hitInSample / sampleSize) / (hitInPop / failInPop)
    #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
    #   p.val <- fisher.test(mat.count)$p.value
    df.BP_enrichment$BP[j] <- df.tb.annot.set$Var1[j]
    df.BP_enrichment$p.val[j] <- p.val
    df.BP_enrichment$percent[j] <- (n.inset / n.genes.sset)
    df.BP_enrichment$genes[j] = n.inset
    
    df.BP_enrichment$foldchange[j] = foldChange
    
  }
  
  df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$p.val),]
  df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$foldchange),]
  df.BP_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= th)
  
  # write.csv(df.BP_enrichment, "output/GO_AtBZR1_ChIP_1101_unique_genes.csv", row.names = FALSE)
  
  ### make plot ### 
  
  vec.significance <- character()
  for(i in 1:nrow(df.BP_enrichment)){
    if(df.BP_enrichment$p.val[i] <= 0.05 & df.BP_enrichment$p.val[i] > 0.01){
      vec.significance <- c(vec.significance, "*")
    }else if(df.BP_enrichment$p.val[i] <= 0.01){
      vec.significance <- c(vec.significance, "**")
    }else{
      vec.significance <- c(vec.significance, "")  
    }
  }
  
  library(ggplot2)
  p1 <- (ggplot(df.BP_enrichment, aes(x=df.BP_enrichment$BP, y=(df.BP_enrichment$foldchange))) +
           geom_bar(colour="black", fill="steelblue", stat="identity") +
           coord_flip() +
           #coord_trans(y="sqrt") +
           guides(fill=FALSE) +
           xlab("") +
           geom_text(aes(label=vec.significance), vjust=0.6, colour = "black", face = "bold", size = 9) +
           ylab("Fold Change") +
           #ylab("# of Regulated Genes (Shared Genes Removed)") +
           theme(axis.text=element_text(size=9, colour = "black"), axis.title=element_text(size=9)) + 
           ggtitle(title) +
           theme(plot.title = element_text(lineheight=1.0, face="bold", size = 9)) +
           theme_minimal()
         
         ) 
  
  p1  
  
}

