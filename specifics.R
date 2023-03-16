# Specific analysis of relationship between Myb49 and RSK1
# Context-specific  (condition, tissue) cases of negative regulation / correlation 
# Examples and statistical analysis of context-specific inhibitive effect of MYB49 on RSK1

folder_data = "data/inhibitor_analysis/"

df.annotation <- read.csv(paste(folder_data, "experiment_annotation.txt" ,sep = "/"),  header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
df.annotation <- subset(df.annotation, !is.na(df.annotation$unique_ID) & df.annotation$condition_treatment_1 != "")

m.fc = data.matrix(read.table(paste(folder_data, "m.foldChange_differentialExpression.txt" ,sep = "/"), header = F, sep = "\t", stringsAsFactors = F), rownames.force = NA)

genes = read.table(paste(folder_data, "genes.txt" ,sep = "/"), header = F, sep = "\t", stringsAsFactors = F)[,1]
experiment_series_ids = read.table(paste(folder_data, "ids_differentialExpressionSamples.txt" ,sep = "/"), header = F, sep = "\t", stringsAsFactors = F)[,1]
experiment_series_ids = as.character(experiment_series_ids)

rownames(m.fc) = genes
colnames(m.fc) = experiment_series_ids


## Extract condition
require(reshape2)
require(ggplot2)
#require(ranger)

require(curl)       # read file from google drive
require(gplots)     # heatmap.2
require(dendextend) # make and color dendrogram
require(colorspace) # diverge_hcl / rainbow_hcl / heat_hcl color palettes
require(ggdendro)

library(dplyr)

context_gsea <- function(tb.context.sample, tb.context.population, th = 0.05, title = "", context = "", col_grad = c("yellow", "red")){
  
  d.gsea <- data.frame(context = names(tb.context.sample), 
                       p.val = rep(1, length(tb.context.sample)),
                       number = rep(0, length(tb.context.sample)), 
                       percent = rep(0, length(tb.context.sample)), 
                       foldchange = rep(1, length(tb.context.sample)), stringsAsFactors = FALSE)
  
  for(i in 1:length(tb.context.sample)){
    
    hitInSample <- tb.context.sample[i]
    sampleSize <- sum(tb.context.sample)
    hitInPop <- tb.context.population[names(tb.context.sample)[i]]
    popSize <- sum(tb.context.population)
    failInPop <- popSize - hitInPop
    fc <- ((hitInSample / sampleSize) / (hitInPop / popSize))
    p.val <- phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
    
    d.gsea$p.val[i] <- p.val
    d.gsea$number[i] <- hitInSample
    d.gsea$percent[i] <- round(hitInSample / sampleSize * 100, 1)
    d.gsea$foldchange[i] = fc
    
  }
  
  d.gsea <- d.gsea[order(d.gsea$foldchange, decreasing = T),]
  
  rangs <- c(min(d.gsea$number), max(d.gsea$number))
  
  data <- d.gsea
  # data$p.val <- -log10(data$p.val)
  
  data %>%
    arrange(desc(foldchange)) %>%
    mutate(context = factor(context, context)) %>%
    ggplot(aes(x=context, y=foldchange, size=number, color=p.val)) +
    scale_colour_gradient(low = col_grad[1], high = col_grad[2], na.value = NA) +
    ggtitle(title) +
    geom_point(alpha=0.9) +
    scale_size(range = rangs, name=paste("Number")) + 
    coord_flip() +
    xlab("") +
    ylab(paste("Fold change (", context, "over representation)")) +
    geom_hline(yintercept=c(1), linetype="dotted") +
    #ylab("Adjusted P value, -log10") +
    theme_classic() + 
    theme(plot.title = element_text(size = 9), axis.text=element_text(size=9), axis.title=element_text(size=9,face="italic"), legend.title=element_text(size=8, face="italic"))
  
}


extract_context <- function(df.annotation, ids = NULL){
  
  if(!is.null(ids)){
    df.annotation <- subset(df.annotation, df.annotation$unique_ID %in% ids)
  }
  
  v.tissues = unique(df.annotation$condition_tissue)
  v.treatments = unique(c(df.annotation$condition_treatment_1, df.annotation$condition_treatment_2))
  v.treatments = v.treatments[!v.treatments == ""]
  
  tb.tissues = numeric(length(v.tissues))
  names(tb.tissues) = v.tissues
  for(i in 1:length(v.tissues)){
    idx = which(df.annotation$condition_tissue == v.tissues[i])
    tb.tissues[i] = length(idx)
  }
  tb.condition_tissues = tb.tissues
  
  
  tb.condition_treatments = numeric(length(v.treatments))
  names(tb.condition_treatments) = unique(v.treatments)
  for(i in 1:length(tb.condition_treatments)){
    idx_1 = which(df.annotation$condition_treatment_1 %in% names(tb.condition_treatments)[i])
    idx_2 = which(df.annotation$condition_treatment_2 %in% names(tb.condition_treatments)[i])
    tb.condition_treatments[i] = length(idx_1) + length(idx_2)
  }
  
  return(list(tb.condition_treatments=tb.condition_treatments, 
              tb.condition_tissues=tb.condition_tissues))
}
  

MYB49 = "AT5G54230"
RSK1 = "AT2G26290"

idx_sig <- intersect(which(abs(m.fc[MYB49,]) > 0.7), which(abs(m.fc[RSK1,]) > 0.7))

idx <- c(intersect(which(m.fc[MYB49,] > 0.7), which(m.fc[RSK1,] < -0.7)),
         intersect(which(m.fc[MYB49,] < -0.7), which(m.fc[RSK1,] > 0.7)))

cols <- c(rgb(0,0,0,alpha = 0.7), rgb(0,0.5,0,alpha=0.7),rgb(0.5,0,0,alpha=0.7))
idx_cols <- as.numeric(1:dim(m.fc)[2] %in% idx_sig) + 1
idx_cols[which(1:dim(m.fc)[2] %in% idx)] = 3

plot(m.fc[MYB49,], m.fc[RSK1,], xlab = "logFC Myb49", ylab = "logFC Rsk1", pch = 1, col = cols[idx_cols], cex = 1.1)
title("Pairwise differential expression across individual conditions between Myb49 (AT5G54230) and Rsk1 (AT2G26290) (circles represent individual conditions)", cex.main = 0.8)

abline(v=c(-0.7,0.7), h = c(-0.7, 0.7), lty = 2)

legend("topright", inset=.02,
       c("Myb49 or Rsk1 differentially expressed below threshold (logFC < 0.7 or > -0.7)", 
         "Myb49 and Rsk1 differentially expressed with similar directionality",
         "Myb49 and Rsk1 differentially expressed with opposite directionality"), fill=cols, horiz=F, cex = 0.75)

#### extract and visualize (new barplot)

# barplot indicates that (above threshold) there is an almost equal ratio between similar and opposite differential expression
counts <- table(idx_cols)

barplot(counts, horiz=F, col = cols,
        names.arg=c("below threshold", "similar directionality", "opposite directionality"))
title("Distribution of pairwise differential expression across 435 individual conditions between Myb49 (AT5G54230) and Rsk1 (AT2G26290)", cex.main = 0.9)

legend("topright", inset=.02,
       c("Myb49 or Rsk1 differentially expressed below threshold (logFC < 0.7 or > -0.7)", 
         "Myb49 and Rsk1 differentially expressed with similar directionality",
         "Myb49 and Rsk1 differentially expressed with opposite directionality"), fill=cols, horiz=F, cex = 1)


# define condition and tissue population
res <- extract_context(df.annotation)
tb.conditions <- res$tb.condition_treatments
tb.tissues <- res$tb.condition_tissues

# non, similar, opposite
idx_similar <- idx_sig[!idx_sig %in% idx]
res <- extract_context(df.annotation, ids = experiment_series_ids[idx_similar])
tb.conditions.sig <- res$tb.condition_treatments
tb.tissues.sig <- res$tb.condition_tissues

title = "Over-representation of identified conditions with pairwise similar differential expression of Myb49 (AT5G54230) and Rsk1 (AT2G26290)"
context_gsea(tb.conditions.sig, tb.conditions, th = 0.05, title = title, context = "condition", col_grad = c("darkgreen", "gray"))

title = "Over-representation of identified tissues with pairwise similar differential expression of Myb49 (AT5G54230) and Rsk1 (AT2G26290)"
context_gsea(tb.tissues.sig, tb.tissues, th = 0.05, title = title, context = "condition", col_grad = c("darkgreen", "gray"))


idx_opposite <- idx
res <- extract_context(df.annotation, ids = experiment_series_ids[idx_opposite])
tb.conditions.sig <- res$tb.condition_treatments
tb.tissues.sig <- res$tb.condition_tissues

title = "Over-representation of identified conditions with pairwise opposite differential expression of Myb49 (AT5G54230) and Rsk1 (AT2G26290)"
context_gsea(tb.conditions.sig, tb.conditions, th = 0.05, title = title, context = "condition", col_grad = c("darkred", "gray"))

title = "Over-representation of identified tissues with pairwise opposite differential expression of Myb49 (AT5G54230) and Rsk1 (AT2G26290)"
context_gsea(tb.tissues.sig, tb.tissues, th = 0.05, title = title, context = "condition", col_grad = c("darkred", "gray"))



