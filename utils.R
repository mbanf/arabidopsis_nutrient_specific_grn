
compute_randomforest_based_GRN <- function(mat.expression=mat.expression, k="sqrt", nb.trees=1000, set.regulators = NULL, set.genes = NULL, seed=1234, importance.measure = "impurity", n.cpus = 2){
  
  library(ranger)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      genes.undefined <- setdiff(regulators, genes)
      if (length(genes.undefined) != 0) {
        stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  # set mtry
  if (class(k) == "numeric") {
    mtry <- K
  } else if (k == "sqrt") {
    mtry <- round(sqrt(n.regulators))
  } else if (k == "all") {
    mtry <- n.regulators-1
  } else {
    stop("Error: invalid parameter k, options: \"sqrt\", \"all\", or integer")
  }
  
  print(paste("Performing random-forest regression based gene regulatory network inference (# of decision trees per gene is: ", nb.trees, ", # of regulators per decision tree node: ", mtry, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  
  strt<-Sys.time()
  for(i in 1:n.genes){
    cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i, drop=FALSE] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    df.xy <- cbind(as.data.frame(x),y)
    rf.model <- ranger(y ~ ., data = df.xy, mtry=mtry, num.trees=nb.trees, importance = "impurity", num.threads = n.cpus)
    imp_scores <- rf.model$variable.importance  
    imp_scores.names <- names(imp_scores)
    mat.rapid[imp_scores.names, gn.i] <- as.numeric(imp_scores)
  }
  print(Sys.time()-strt)
  print("..finished.")
  
  return(mat.rapid / n.samples)
}       


plot_ratios <- function(cormat=cormat, m.value = m.value, v.title = "", v.name = "", v.max = 1.0, v.min = 0, max.line = 100){
  # Melt the correlation matrix
  melted_cormat <- melt(cormat)
  melted_cormat <- na.omit(melted_cormat)
  
  
  melted_value <- melt(m.value)
  melted_value <- na.omit(melted_value)
  
  
  
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "green", high = "red", mid = "white", 
                         midpoint = 0, limit = c(v.min,v.max), name=v.name) +
    
    #               scale_x_discrete(breaks = rownames(cormat), labels=v.x_axis.ticks) + 
    #               scale_y_discrete(breaks = colnames(cormat), labels=v.y_axis.ticks) + 
    
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1,
                                     colour="black")) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                     size = 10, hjust = 0,
                                     colour="black")) +
    ggtitle(v.title) +
    theme(plot.title = element_text(lineheight=1.0, face="bold", size = 12)) + 
    coord_fixed()
  
  
  my.lines<-data.frame(x=c(3.5,6.5,0.5), y=c(0.5,0.5, max.line - 1), xend=c(3.5,6.5,9.5), yend=c(max.line,max.line, max.line - 1))
  
  ggheatmap <- (ggheatmap + 
                  #geom_text(aes(Var2, Var1, label = ""), color = "black", size = 4) + 
                  geom_text(aes(Var2, Var1, label = melted_value$value), color = "black", size = 6, vjust = 1.2) + 
                  theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks = element_blank()
                    #legend.justification = c(1, 0),
                    #legend.position = c(0.6, 0.7),
                    #legend.direction = "horizontal")+
                    #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
                  ) + 
                  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
  )
  
  
  
  
  return(ggheatmap)
}


