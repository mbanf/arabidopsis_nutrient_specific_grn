# TODO: use go slim 

# df.rsk1 <- subset(df, df$V1 == rsk1)


go_function <- function(v.gns, ontology = "BP"){
  
  df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$acc. %in% v.gns)
  n.genes <- length(unique(df.GO.annot$acc.))
  
  df.GO.annot
}


# TODO: write up the number of annotations - rsk1 not in set 
# p < 0.05 and smaller represented as -log10 
# 
# gn
# gns.DE

go_gsea <- function(v.gns, gn.pop = NULL, th = 0.05, title = "", bg.mode = "all",
                    mode = "enrichment"){
  
  #### TODO: 2 variants - population is complete, population is grn
  
  # F, C, P
  df <- read.table("data/ATH_GO_GOSLIM.txt", fill = T, sep = "\t")
  ec <- c("HTP", "HDA", "HMP", "HGI", "HEP","EXP","IDA","IPI","IMP","IGI","IEP")
  df <- subset(df, df$V10 %in% ec)
  df <- subset(df, df$V8 %in% c("F","P"))
  
  df.GO.annot <- df[,c(1,5,8)]
  names(df.GO.annot) <- c("acc.", "Term", "Ontology")
  
  if(is.null(gn.pop)){
    gn.pop <- unique(df.GO.annot$acc.)
  }
  
  # df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
  # n.genes <- length(unique(df.GO.annot$acc.))
  
  
  # df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
  
  # TODO: log fold change, bubble plot, log p-value - explain large fold change
  df.GO.annot <- subset(df.GO.annot, df.GO.annot$acc. %in% gn.pop)
  n.genes <- length(unique(df.GO.annot$acc.))

  # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
  df.GO <- subset(df.GO.annot, df.GO.annot$acc. %in% v.gns)
  
  df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
  df.tb.annot.set <- as.data.frame(table(df.GO$Term[df.GO$Term != ""]), stringsAsFactors = FALSE)
  df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
  names(df.tb.annot.set) <- c("term", "in_sample", "in_population")
  n.genes.sset <- length(unique(df.GO$acc.))
  
  # gibberellin mediated signaling pathway
  if(bg.mode == "all"){
    n.genes <- length(gn.pop)
    n.genes.sset <- length(v.gns) # use this always
  }
  
  n.genes <- 27655 # used for vs genome analysis
  
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
  
  # Libraries
  library(ggplot2)
  library(dplyr)

  
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
  
  
  
  df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
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




####



# 
# #v.genes <- as.character(read.csv("BR_regulated_genes_talk.csv", header = FALSE, quote = "", stringsAsFactors = FALSE)[,1])
# v.genes <- as.character(read.csv("BR_responsive_BZR1_targets.txt", header = FALSE, quote = "", stringsAsFactors = FALSE)[,1])
# v.genes <- as.character(read.csv("BR_responsive_nonBZR1_targets.txt", header = FALSE, quote = "", stringsAsFactors = FALSE)[,1])
# v.genes <- as.character(read.csv("NonBR_responsive_BZR1_targets.txt", header = FALSE, quote = "", stringsAsFactors = FALSE)[,1])
# 
# 
# 
# 
# 
# df.genes <- read.table("BR_targets/Michael BR responsive BZR1 targets.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
# 
# 
# 
# 
# 
# #df.genes <- read.table("genes_filtered_p01.txt", header = TRUE, sep = "\t", quote = "")
# # 
# for(i in 1:nrow(df.genes)){
#   if(grepl("AC", df.genes$name[i])){
#     df.genes$name[i] <- paste(substr(df.genes$name[i], 1, 13),"T", substr(df.genes$name[i],13, 17), sep = "")
#   }
# }

# df.GO.annot <- readRDS("df.GO_11.rds")
# df.GO.annot$acc. <- gsub("\\|.*", "", df.GO.annot$acc.)
# for(i in 1:nrow(df.GO.annot)){
#   cat("Processing... ", round(i/nrow(df.GO.annot) * 100, digits = 2) , "%", "\r") 
#   flush.console()
#   if(grepl("AC", df.GO.annot$acc.[i])){
#     df.GO.annot$acc.[i] <- paste(substr(df.GO.annot$acc.[i], 1, 12),substr(df.GO.annot$acc.[i],14, 17), sep = "")
#   }else{
#     df.GO.annot$acc.[i] <- gsub("_.*", "", df.GO.annot$acc.[i])
#   }
# }
#saveRDS(df.GO.annot, "df.GO_Maize_Trimmed.rds")

# 
# lst.df.BP_enrichment <- list(ncol(df.genes))
# for(i in 1:ncol(df.genes)){
#   
#   v.genes <- as.character(df.genes[,i])
#   
#   df.GO.annot <- readRDS("df.GO_FULL_Maize_Trimmed.rds")
#   #df.GO.annot <- readRDS("df.GO_Maize_Trimmed.rds")
#   df.GO.annot <- unique(df.GO.annot)
#   df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == "BP")
#   
#   df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% v.genes)
#   df.GO.annot.set <- unique(df.GO.annot.set)
#   
#   #df.GO.annot.set <- subset(df.GO.annot.set, df.GO.annot.set$Ontology == "MF")
#   #write.csv(df.GO.annot.set, "MF.csv", row.names = FALSE)
#   
#   #### GO Analysis and Enrichment Tests
#   #library("biomaRt")
#   #genome_names <- c("ath", "bdi", "cre", "gma", "mtr", "osa",  "ppt", "sit", "sly", "smo", "zma")
#   
#   #genome_size.transcript <-  63391
#   #genome_size.locus <- 39469
#   
#   # 4883 genes out of 5600
#   # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
#   df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
#   df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
#   df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
#   
#   ## filter go - evidence codes, biological process
#   df.BP_enrichment <- data.frame(BP = character(nrow(df.tb.annot.set)), p.val = numeric(nrow(df.tb.annot.set)),
#                                  percent = numeric(nrow(df.tb.annot.set)), genes = numeric(nrow(df.tb.annot.set)), 
#                                  foldchange = numeric(nrow(df.tb.annot.set)), stringsAsFactors = FALSE)
#   n.genes.sset <- length(unique(df.GO.annot.set$acc.))
#   n.genes <- length(unique(df.GO.annot$acc.))
#   
#   for(j in 1:nrow(df.tb.annot.set)){  
#     
#     n.inset <- df.tb.annot.set$Freq.x[j]
#     n.genomewide <- df.tb.annot.set$Freq.y[j]
#     
#     #n.genomewide <- df.tb.annot$Freq[df.tb.annot$Var1 == df.tb.annot$Var1[j]]
#     
#     ### global enrichment test - gene basis
#     hitInSample <- n.inset
#     sampleSize <- n.genes.sset
#     hitInPop <- n.genomewide #sum(tb.rate_limiting_domains$Freq)
#     failInPop <- n.genes - hitInPop #(nrow(df.global.domains) - hitInPop)
#     p.val <- print(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
#     
#     #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
#     foldChange <- (hitInSample / sampleSize) / (hitInPop / failInPop)
#     #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
#     #   p.val <- fisher.test(mat.count)$p.value
#     df.BP_enrichment$BP[j] <- df.tb.annot.set$Var1[j]
#     df.BP_enrichment$p.val[j] <- p.val
#     df.BP_enrichment$percent[j] <- (n.inset / n.genes.sset)
#     df.BP_enrichment$genes[j] = n.inset
#     
#     df.BP_enrichment$foldchange[j] = foldChange
#     
#   }
#   
#   df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$p.val),]
#   df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$foldchange),]
#   #df.BP_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= 0.05)
#   df.BP_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= 0.01)
#   
#   df.BP_enrichment["group"] <- names(df.genes)[i]
#   lst.df.BP_enrichment[[i]] <- df.BP_enrichment
# }
# names(lst.df.BP_enrichment) <- names(df.genes)
# #write.csv(df.BP_enrichment, "Enrichment_BR_responsive_BZR1_targets.csv")
# #write.csv(df.BP_enrichment, "Enrichment_BR_responsive_nonBZR1_targets.csv")
# #write.csv(df.BP_enrichment, "Enrichment_NonBR_responsive_BZR1_targets.csv")
# 
# df.total <- rbind(lst.df.BP_enrichment[[1]],
#                   lst.df.BP_enrichment[[2]],
#                   lst.df.BP_enrichment[[3]],
#                   lst.df.BP_enrichment[[4]],
#                   lst.df.BP_enrichment[[5]],
#                   lst.df.BP_enrichment[[6]])
# 
# 
# # df.1 <- read.csv("Enrichment_BR_responsive_BZR1_targets.csv")
# # df.2 <- read.csv("Enrichment_BR_responsive_nonBZR1_targets.csv")
# # df.3 <- read.csv("Enrichment_NonBR_responsive_BZR1_targets.csv")
# # 
# # df.1["group"] <- "BR_responsive_BZR1_targets"
# # df.2["group"] <- "BR_responsive_nonBZR1_targets"
# # df.3["group"] <- "NonBR_responsive_BZR1_targets"
# # 
# # df.total <- rbind(df.1,df.2,df.3)
# 
# #df.total <- df.total[order(df.total$BP),]
# #write.csv(df.total, "Enrichment_all_NonBR_responsive_BZR1_targets.csv")
# #write.csv(df.total, "Enrichment_all_sixCategories.csv")
# 
# 
# for(j in 1:ncol(df.genes)){
#   
#   df.BP_enrichment <- lst.df.BP_enrichment[[j]]
#   
#   vec.significance <- character()
#   for(i in 1:nrow(df.BP_enrichment)){
#     if(df.BP_enrichment$p.val[i] <= 0.05 & df.BP_enrichment$p.val[i] > 0.01){
#       vec.significance <- c(vec.significance, "*")
#     }else if(df.BP_enrichment$p.val[i] <= 0.01){
#       vec.significance <- c(vec.significance, "**")
#     }else{
#       vec.significance <- c(vec.significance, "")  
#     }
#   }
#   
#   library(ggplot2)
#   p1 <- (ggplot(df.BP_enrichment, aes(x=df.BP_enrichment$BP, y=(df.BP_enrichment$foldchange))) +
#            geom_bar(colour="black", aes(fill="gray"),, stat="identity") +
#            coord_flip() +
#            #coord_trans(y="sqrt") +
#            guides(fill=FALSE) +
#            xlab("") +
#            geom_text(aes(label=vec.significance), vjust=0.6, colour = "black", face = "bold", size = 12) +
#            ylab("Fold Change") +
#            #ylab("# of Regulated Genes (Shared Genes Removed)") +
#            theme(axis.text=element_text(size=20, colour = "black"), axis.title=element_text(size=20)) + 
#            #ggtitle("Percentage of Genes associated to Biological Process") +
#            theme(plot.title = element_text(lineheight=1.0, face="bold", size = 20)))
#   
#   names(lst.df.BP_enrichment)[j]
#   p1  
# }
