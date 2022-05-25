install_and_load_libraries <- function(){
  
  # CRAN
  list.of.packages <- c("ggplot2", "reshape2","ranger", "curl", "gplots", "dendextend", "colorspace", "ggdendro", "grid", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install()
  # 
  # list.of.packages <- c("Biostrings", "VariantAnnotation","BSgenome")
  # BiocManager::install(list.of.packages)
  # 
  
  
  require(reshape2)
  require(ggplot2)
  require(ranger)
  
  require(curl)       # read file from google drive
  require(gplots)     # heatmap.2
  require(dendextend) # make and color dendrogram
  require(colorspace) # diverge_hcl / rainbow_hcl / heat_hcl color palettes
  require(ggdendro)

  library(dplyr)
  
  require(grid)
  
}




compute_randomforest_based_GRN <- function(mat.expression=mat.expression, k="sqrt", nb.trees=1000, set.regulators = NULL, set.genes = NULL, seed=1234, importance.measure = "impurity", n.cpus = 2){
  
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



plot_regulons_tgs <- function(m.exp=m.exp, m.sig = m.sig, 
                              tfs=tfs, tgs=tgs,
                              title = "", v.name = "", 
                              v.max = 1.0, v.min = 0){
  
  # assumption is that the matrix has already been ordered by dendrogram
  dendro.plot <- ggdendrogram(data = dd.col, rotate = TRUE)
  dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank())
  
  n.specs = dim(m.exp)[2]
  
  # Melt the correlation matrix
  melted_cormat <- melt((m.exp))
  melted_cormat <- na.omit(melted_cormat)
  
  melted_value <- melt((m.sig))
  melted_value <- na.omit(melted_value)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "green", high = "red", mid = "white", 
                         midpoint = 0, limit = c(v.min,v.max), name=v.name) +
    
    #               scale_x_discrete(breaks = rownames(cormat), labels=v.x_axis.ticks) + 
    #               scale_y_discrete(breaks = colnames(cormat), labels=v.y_axis.ticks) + 
    
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 6, hjust = 1,
                                     colour="black")) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0,
                                     size = 6, hjust = 0,
                                     colour="black")) +
    
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=6), #change legend title font size
          legend.text = element_text(size=6)) +  #change legend text font size
  
    # ggtitle(title) +
    # theme(plot.title = element_text(lineheight=1.0,  size = 11)) + 
    coord_fixed()
  
  
  ggheatmap <- (ggheatmap + 
                  
                  # commented out for non significance representation
                  geom_text(aes(Var2, Var1, label = melted_value$value), color = "black", size = 6, vjust = 0.8) + 
                  
                  theme(
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(), # remove identifies
                    panel.grid.major = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "left",
                    # legend.justification = c(-1, 1),
                    #legend.position = c(0.6, 0.7),
                    #legend.direction = "horizontal")+
                    #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
                  ) 
                  # geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
                
  )
  
  pdf("test.pdf", height = 15, width = 8, paper = "letter")
  
  grid.newpage()
  print(ggheatmap, vp = viewport(x = 0.32, y = 0.5, width = 0.2, height = 1.0))
  print(dendro.plot, vp = viewport(x = 0.5, y = 0.512, width = 0.2, height = 1.055))
  
  dev.off()
  
  # 8 15
  
  #export figure
  ggsave(p,filename="test.png",height=8,width=15,units="in",dpi=200)
  
  return(ggheatmap)
}



plot_regulons <- function(m.exp=m.exp, m.sig = m.sig, 
                          tfs=tfs, tgs=tgs,
                          title = "", v.name = "", 
                          v.max = 1.0, v.min = 0){
  
  n.specs = dim(m.exp)[2]
  
  # Melt the correlation matrix
  melted_cormat <- melt(t(m.exp[c(tfs, tgs),]))
  melted_cormat <- na.omit(melted_cormat)
  
  melted_value <- melt(t(m.sig[c(tfs, tgs),]))
  melted_value <- na.omit(melted_value)
  
  color.palette  <- colorRampPalette(c("yellow", "orange", "red"))(n=600)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "green", high = "red", mid = "white", 
                         midpoint = 0, limit = c(v.min,v.max), name=v.name) +
    
    #               scale_x_discrete(breaks = rownames(cormat), labels=v.x_axis.ticks) + 
    #               scale_y_discrete(breaks = colnames(cormat), labels=v.y_axis.ticks) + 
    
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 9, hjust = 1,
                                     colour="black")) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                     size = 9, hjust = 0,
                                     colour="black")) +
    ggtitle(title) +
    theme(plot.title = element_text(lineheight=1.0,  size = 11)) + 
    coord_fixed()
  
  
  # my.lines<-data.frame(x=c(3.5,6.5,0.5), y=c(0.5,0.5, max.line - 1), xend=c(3.5,6.5,9.5), yend=c(max.line,max.line, max.line - 1))
  max.line <- length(tfs) + 0.5
  yend = n.specs + 0.5
  my.lines<-data.frame(x=c(max.line), y=c(0.5), xend=c(max.line), yend=c(yend))
  
  
  
  ggheatmap <- (ggheatmap + 
                  
                  # commented out for non significance representation
                  geom_text(aes(Var2, Var1, label = melted_value$value), color = "black", size = 6, vjust = 1) + 
                  
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
  
  
  #export figure
  # ggsave(p,filename="measles-mod3.png",height=5.5,width=8.8,units="in",dpi=200)
  
  return(ggheatmap)
}



plot_clean_targets <- function(m.exp, dd.col,  c_group = 6){
  
  
  dend1 <- dd.col
  
  library(curl)       # read file from google drive
  library(gplots)     # heatmap.2
  library(dendextend) # make and color dendrogram
  library(colorspace) # diverge_hcl / rainbow_hcl / heat_hcl color palettes
  
  
  dev.off()
  
  # number of clusters
  # hc <- hclust(m.exp)
  # ct <- cutree(hc, k = c_group)
  # 
  # dend1 <- as.dendrogram(hc)
  dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
  dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)   # add color to the labels
  
  # reorder the dendrogram, must incl. `agglo.FUN = mean`
  # rMeans <- rowMeans(m.exp, na.rm = T)
  # dend1 <- reorder(dend1, rowMeans(m.exp, na.rm = T), agglo.FUN = mean)
  
  # get the color of the leaves (labels) for `heatmap.2`
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  
  # if plot the dendrogram alone:
  # the size of the labels:
  dend1 <- set(dend1, "labels_cex", 0.5)
  par(mar = c(1,1,1,1))
  plot_horiz.dendrogram(dend1, side = F) # use side = T to horiz mirror if needed
  
  
  color.palette  <- colorRampPalette(c("green", "white", "red"))(n=50)
  ###
  
  ## plot the heatmap with the dendrogram above ##
  par(cex.main=0.5)                   # adjust font size of titles
  heatmap.2(m.exp, main = '',
            # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
            # order by branch mean so the deepest color is at the top
            dendrogram = "row",        # no dendrogram for columns
            Rowv = dend1,              # * use self-made dendrogram
            Colv = "NA",               # make sure the columns follow data's order
            col = color.palette,         # color pattern of the heatmap
            
            trace="none",              # hide trace
            density.info="none",       # hide histogram
            key = T,
            keysize = 1.5,
            key.title = NA,
            margins = c(5,3),         # margin on top(bottom) and left(right) side.
            cexRow=0.4, cexCol = 0.8,      # size of row / column labels
            xlab = "",
            srtCol=90, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
            # margin for the color key
            # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
            key.par=list(mar=c(10,1,2,1)),
            RowSideColors = col_labels, # to add nice colored strips        
            colRow = col_labels,         # add color to label
            labRow=FALSE
  )
 
  
}




plot_dend <- function(){
  # https://liuyanguu.github.io/post/2018/07/16/how-to-draw-heatmap-with-colorful-dendrogram/
  
  dev.off()
  
  c_group <- 8 # number of clusters
  hc <- hclust(dist(F_m2))
  ct <- cutree(hc, k = c_group)
  
  dend1 <- as.dendrogram(hc)
  dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
  dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)   # add color to the labels
  
  # reorder the dendrogram, must incl. `agglo.FUN = mean`
  rMeans <- rowMeans(F_m2, na.rm = T)
  dend1 <- reorder(dend1, rowMeans(F_m2, na.rm = T), agglo.FUN = mean)
  
  # get the color of the leaves (labels) for `heatmap.2`
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  
  # if plot the dendrogram alone:
  # the size of the labels:
  dend1 <- set(dend1, "labels_cex", 0.5)
  par(mar = c(1,1,1,14))
  plot_horiz.dendrogram(dend1, side = F) # use side = T to horiz mirror if needed
  
  ###
  
  ## plot the heatmap with the dendrogram above ##
  par(cex.main=0.5)                   # adjust font size of titles
  heatmap.2(F_m2, main = '',
            # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
            # order by branch mean so the deepest color is at the top
            dendrogram = "row",        # no dendrogram for columns
            Rowv = dend1,              # * use self-made dendrogram
            Colv = "NA",               # make sure the columns follow data's order
            col = diverge_hcl,         # color pattern of the heatmap
            
            trace="none",              # hide trace
            density.info="none",       # hide histogram
            
            margins = c(5,18),         # margin on top(bottom) and left(right) side.
            cexRow=0.4, cexCol = 0.8,      # size of row / column labels
            xlab = "",
            srtCol=90, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
            # margin for the color key
            # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
            key.par=list(mar=c(10,1,2,1)),
            RowSideColors = col_labels, # to add nice colored strips        
            colRow = col_labels         # add color to label
  )
}




plot_condition_sim <- function(){
  
  message("plot condition similarity")
  
  hc <- hclust(dist(t(data.frame("0PF" = v.0PF, "P0F" = v.P0F, "0P0F" = v.0P0F))), "ave")
  ggdendrogram(hc, rotate = FALSE, size = 4) + labs(title="Dendrogram in ggplot2")
  
  dendr <- dendro_data(hc, type="rectangle")
  #your own labels are supplied in geom_text() and label=label
  ggplot() +
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr), aes(x=x, y=y, label=c("  - Ph / + Fe", "  + Ph / - Fe", "  - Ph / - Fe"), hjust=0), size=4) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
}
