
df.grn <- c()











df.grn["coDiffExp_0PF"] <- 0
df.grn["coDiffExp_P0F"] <- 0
df.grn["coDiffExp_0P0F"] <- 0
df.grn["coDiffExp"] <- 0


for(j in 1:length(tfs.DE)){
  
  df.regulatoryNetwork.j <- subset(df.regulatoryNetwork.meta, df.regulatoryNetwork.meta$tf ==  tfs.DE[j])
  
  i.min <- min(which(m.anova.set[tfs.DE[j], 1:3] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)]
    i.set <- which(df.regulatoryNetwork.j$target %in% tgs.de)
    df.regulatoryNetwork.j$coDiffExp_0PF[i.set] <- 1
    
  }
  
  
  i.min <- min(which(m.anova.set[tfs.DE[j], 4:6] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)]
    i.set <- which(df.regulatoryNetwork.j$target %in% tgs.de)
    df.regulatoryNetwork.j$coDiffExp_P0F[i.set] <- 1
    
  }
  
  
  i.min <- min(which(m.anova.set[tfs.DE[j], 7:9] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)]
    i.set <- which(df.regulatoryNetwork.j$target %in% tgs.de)
    df.regulatoryNetwork.j$coDiffExp_0P0F[i.set] <- 1
    
  }
  
  df.regulatoryNetwork.j$coDiffExp <- pmax(pmax(df.regulatoryNetwork.j$coDiffExp_0PF, df.regulatoryNetwork.j$coDiffExp_P0F), df.regulatoryNetwork.j$coDiffExp_0P0F)
  df.regulatoryNetwork <- rbind(df.regulatoryNetwork, df.regulatoryNetwork.j)
  
}

