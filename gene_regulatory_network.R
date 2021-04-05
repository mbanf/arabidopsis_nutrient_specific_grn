tgs <- intersect(rownames(m.de), gns.DE)
tfs <- intersect(names(v.regulators), gns.DE)

df.grn <- c()

for(j in 1:length(tfs)){
  
  tf <- tfs[j]
  df.grn.j <- c()
  
  i.min <- min(which(m.anova.set[tf, 1:3] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, i.min:3]) >= 1)])
    
    df <- data.frame(TF = rep(tf, length(tgs.de)), TG = tgs.de)
    df["0PF"] <- 1
    df["P0F"] <- 0
    df["0P0F"] <- 0
    
    df["mode"] <- 0
    for(i in 1:length(tgs.de)){
      idx <- min(which(m.anova.set[tgs.de[i], 1:3] == 1))
      df$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    }
    
    df.grn.j <- df
  }

    
  i.min <- min(which(m.anova.set[tf, 4:6] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (3 + i.min):6]) >= 1)])
    df <- data.frame(TF = rep(tf, length(tgs.de)), TG = tgs.de)
    df["0PF"] <- 0
    df["P0F"] <- 1
    df["0P0F"] <- 0
    
    df["mode"] <- 0
    for(i in 1:length(tgs.de)){
      idx <- min(which(m.anova.set[tgs.de[i], (3 + i.min):6] == 1))
      df$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    }
    
    df.grn.j <- rbind(df.grn.j, df)
  }
  

  i.min <- min(which(m.anova.set[tf, 7:9] == 1))
  
  if(i.min != Inf){
    
    tgs.de <- intersect(tgs, rownames(m.anova.set)[which(rowSums(m.anova.set[, (6 + i.min):9]) >= 1)])
    df <- data.frame(TF = rep(tf, length(tgs.de)), TG = tgs.de)
    df["0PF"] <- 0
    df["P0F"] <- 0
    df["0P0F"] <- 1
    
    df["mode"] <- 0
    for(i in 1:length(tgs.de)){
      idx <- min(which(m.anova.set[tgs.de[i], (6 + i.min):9] == 1))
      df$mode[i] <- sign(m.de[tf, i.min]) * sign(m.de[tgs.de[i], idx])
    }
    
    df.grn.j <- rbind(df.grn.j, df)
  }
  
  if(length(df.grn.j) > 0){
    # df.grn.j$coDiffExp <- pmax(pmax(df.grn.j$`0PF`, df.grn.j$P0F), df.grn.j$`0P0F`)
    df.grn <- rbind(df.grn, df.grn.j)
  }
}

write.csv(df.grn, "manuscript/codifferential_expression/df.grn.codifferential.csv", row.names = FALSE)

