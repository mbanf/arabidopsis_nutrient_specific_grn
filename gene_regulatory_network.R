
specificity = c("-P+Fe", "+P-Fe", "-P-Fe", 
                "-P+Fe, +P-Fe", "-P+Fe, -P-Fe",
                "+P-Fe, -P-Fe", "-P+Fe, +P-Fe, -P-Fe")

code_specificity = seq(1, length(specificity))
names(code_specificity) <- specificity

cols_specificity = c( "cyan", "tomato", "steelblue", "black", "pink",  "purple", "sienna") 
names(cols_specificity) <- seq(1, length(specificity))


m.spec <- matrix(0, nrow=length(v.regulators), ncol = length(gns.DE))
rownames(m.spec) <- names(v.regulators)
colnames(m.spec) <- gns.DE

m.grn <- matrix(0, nrow=length(v.regulators), ncol = length(gns.DE))
rownames(m.grn) <- names(v.regulators)
colnames(m.grn) <- gns.DE


for(i in 1:nrow(df.code)){
  tf = df.code$TF[i]
  tg = df.code$TG[i]
  if(tf != tg){
    
    m.grn[tf, tg] <- 1
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 1}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 2}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 3}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 0){
      m.spec[tf, tg] <- 4}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 0 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 5}
    
    if(df.code$`0PF`[i] == 0 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 6}
    
    if(df.code$`0PF`[i] == 1 & df.code$P0F[i] == 1 & df.code$`0P0F`[i] == 1){
      m.spec[tf, tg] <- 7}
    
  }
}

for(i in 1:nrow(df.dna.de)){
  tf <- as.character(df.dna.de$TF[i])
  tg <- as.character(df.dna.de$TG[i])
  if(tf != tg){
    m.grn[tf, tg] <- m.grn[tf, tg] + df.dna.de$val[i]
  }
}

for(i in 1:nrow(df.rf.grn)){
  tf <- as.character(df.rf.grn$TF[i])
  tg <- as.character(df.rf.grn$TG[i])
  if(tf != tg){
    m.grn[tf, tg] <- m.grn[tf, tg] + as.numeric(df.rf.grn$val[i])
  }
}

df.grn <- melt(m.grn)
names(df.grn) <- c("TF", "TG", "val")

# df.sign <- melt(m.sign)
# names(df.sign) <- c("TF", "TG", "mode")

df.spec <- melt(m.spec)
names(df.spec) <- c("TF", "TG", "spec")

# df.grn["mode"] <- df.sign["mode"]
df.grn["spec"] <- df.spec["spec"]

df.grn <- df.grn[order(-df.grn$val),]

df.grn <- subset(df.grn, df.grn$spec != 0)
df.grn <- subset(df.grn, df.grn$val > 1.8) # 80 % binding evidence or 80 % by additional root expression 

tfs.final = unique(as.character(df.grn$TF))

idx <- sort(unique(df.grn$spec))
specificity = specificity[idx]
cols_specificity = cols_specificity[idx]


#####  

library(igraph)
g=graph.edgelist(as.matrix(df.grn[,1:2]))


# Node size - accorind to number of targets or strength!1
hubs = table(df.grn$TF)
hubs <- hubs / max(hubs)
# specifitiy of regulation - distribution per hub node! 
# highlight RSK1 in size

V(g)$size <-  ifelse(names(V(g)) == rsk1, 2.5, ifelse(names(V(g)) %in% df.grn$TF, hubs[names(V(g))] * 2 + 1, 2))
E(g)$weight <- df.grn$val / max(df.grn$val)

# up and down regultion --- boring :) 
# E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))
E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))

E(g)$color <- cols_specificity[df.grn$spec]

V(g)$color <- ifelse(names(V(g)) == rsk1, adjustcolor("blue", alpha.f = .9) , ifelse(names(V(g)) %in% tfs.final, adjustcolor("orange", alpha.f = .6), adjustcolor("gray", alpha.f = .6)))
V(g)$type <- ifelse(names(V(g)) %in% tfs.final, 0, 1)

vertex.label <- ifelse(names(V(g)) %in% tfs.final, paste(names(V(g)), "(", v.regulators[names(V(g))], ")", sep="") , "")  # TODO: this might not work?
idx <- which(names(V(g)) == rsk1)
vertex.label[idx] <- paste(rsk1, "(RSK1)")

# set.seed(92235) 
# set.seed(992235) 

set.seed(99999) 
plot(g, edge.arrow.size=.1, 
     edge.curved=seq(-0.5, 0.5, length = ecount(g)),
     vertex.label.dist=0.25,
     vertex.label.cex=0.65,
     vertex.label.color = "black",
     vertex.label=vertex.label,edge.width=E(g)$weight)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))


# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("bottomleft", legend=specificity , col = cols_specificity , bty = "n", pch="-" , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.2, 0.2))



# devtools::install_github("USCCANA/netplot")
library(igraph)
library(netplot)


vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), paste(names(V(g)), "(", v.regulators[names(V(g))], ")") , "")  # TODO: this might not work?
idx <- which(names(V(g)) == rsk1)
vertex.label[idx] <- paste(rsk1, "(RSK1)")

size = ifelse(names(V(g)) %in% df.grn$TF, 2.5, 2)
cols = ifelse(names(V(g)) %in% tfs.final, adjustcolor("brown", alpha.f = .6), adjustcolor("gray", alpha.f = .6))
l <- layout_with_fr(g)
nplot(g, layout = l, vertex.label = vertex.label,
      edge.curvature = 0,
      edge.arrow.size=.01,
      edge.color = ~ego(mix=0) + alter(mix=1),
      edge.width=E(g)$weight,
      vertex.color = cols, vertex.size = size, vertex.label.cex=0.5)


ans <- nplot(
  UKfaculty,
  layout                = l,
  vertex.color          = viridis::plasma(5)[V(UKfaculty)$Group + 1],
  vertex.label          = nam,
  vertex.size.range     = c(.01, .04, 4),
  vertex.label.col      =  "black",
  vertex.label.fontface = "bold",
  bg.col                = "transparent",
  vertex.label.show     = .5,
  vertex.label.range    = c(10, 25),
  edge.width.range      = c(1, 4, 5)
)




nplot(x_network, layout = l)

# Putting two plots in the same page (one using igraph and the other network)
gridExtra::grid.arrange(
  nplot(g, layout = l),
  nplot(g, layout = l), ncol=2, nrow=1
)


set.seed(1)
data("UKfaculty", package = "igraphdata")
l <- layout_with_fr(UKfaculty)

plot(UKfaculty, layout = l) # ala igraph

# Random names
set.seed(1)
nam <- sample(babynames::babynames$name, vcount(UKfaculty))

ans <- nplot(
  UKfaculty,
  layout                = l,
  vertex.color          = viridis::plasma(5)[V(UKfaculty)$Group + 1],
  vertex.label          = nam,
  vertex.size.range     = c(.01, .04, 4),
  vertex.label.col      =  "black",
  vertex.label.fontface = "bold",
  bg.col                = "transparent",
  vertex.label.show     = .5,
  vertex.label.range    = c(10, 25),
  edge.width.range      = c(1, 4, 5)
)


# Plot it!
ans



# w / wo dna binding 
test <- subset(df.grn, df.grn$TF == "AT3G16280")
nrow(test[test$val > 6,]) / nrow(test) # binding ratio

# question: is there a similarily regulated module => similar differentially expressed and dna binding (double significance?)
# putative real causal "effect" of the regulator => activation, repression => time delay

# 5 of the 21 genes with dna binding have a similar differential expression pattern! (the also cluster with rsk1 as one of the genes! for this TF)
# TODO: "correlation between the target genes to identify differential expression similarity rank 
# TODO: BP analysis of the regulated moduled?
# TODO: BP analysis of the 21 genes! 

tgs.db <- as.character(test[test$val > 6,]$TG)
m.peak = m.anova.set
colnames(m.peak) <- exps
colSums(m.peak[tgs.db,]) / length(tgs.db)

tgs.ndb <- as.character(test[test$val < 6,]$TG)
colSums(m.peak[tgs.ndb,]) / length(tgs.ndb)

### cluster targets together with rsk1
rsk1 = "AT2G26290"

library(igraph)
g=graph.edgelist(as.matrix(df.grn[,1:2]))
V(g)$size <-  ifelse(names(V(g)) %in% df.grn$TF, 5, 3)

E(g)$weight <- df.grn$val / max(df.grn$val) * 3

E(g)$color <- ifelse(df.grn$mode == 1, "red", ifelse(df.grn$mode == -1, "green", "gray"))

# idx <- which(names(V(g)) == rsk1)
vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), names(V(g)), "")  # TODO: this might not work?
# vertex.label[idx] <- paste(rsk1, "(RSK1)")


ontology = "MF"
df.GO.annot <- readRDS("data/Athaliana_167.annot.2017.rds")
df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == ontology)
n.genes <- length(unique(df.GO.annot$acc.))

df.GO.annot = go_function(gns.DE, ontology = ontology)

gns <- unique(as.character(df.GO.annot$acc.))
na = "no annotation"
terms <- c(unique(as.character(df.GO.annot$Term)), na)
term_map = seq(1,length(terms))
names(term_map) = terms

values <- list()
nodes <- names(V(g))
for(i in 1:length(nodes)){
  
  vals = rep(0, length(terms))
  df = subset(df.GO.annot, df.GO.annot$acc. == nodes[i])
  if(nrow(df) > 0){
    idx = as.integer(term_map[unique(as.character(df$Term))])
    vals[idx] = 1
  }else{
    vals[31] = 1
  }
  values <- append(values, list(vals))
  
}

# install.packages('colorRamps')
library(colorRamps)
cols = list(primary.colors(length(terms),  steps = 3, no.white = TRUE))
term_color_map = cols[[1]]
names(term_color_map) = terms
term_color_map["no annotation"] = "#808080"
cols[[1]][31] = "#808080"

# vals <- lapply(1:10, function(x) sample(1:10,3))
set.seed(1234) 
plot(g, 
     vertex.shape="pie", 
     edge.arrow.size=.2,
     vertex.pie=values,
     vertex.pie.color=cols,
     vertex.label.dist=1,
     vertex.label.cex=0.6,
     vertex.label.color = "black",
     # vertex.size=seq(10,30,length=10), 
     vertex.label=vertex.label)



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("bottomleft", legend=names(term_color_map)  , col = term_color_map , bty = "n", pch=20 , pt.cex = 2, cex = 1, horiz = FALSE) #, inset = c(0.1, 0.1))


plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))




cols = unique(col_labels)
for(j in 1:c_group){
  v.gns <- names(which(ct == j))
  idx <- which(names(V(g)) %in% v.gns)
  V(g)$color[idx] = rep(cols[j], length(idx))
}
  

vertex.label <- ifelse(names(V(g)) %in% names(v.regulators), names(V(g)), "")

# set RSK1 as node 

plot(g, edge.arrow.size=.2, vertex.color=V(g)$color, vertex.label=vertex.label, vertex.label.dist=0.1, vertex.label.color = "black", vertex.label.size = 0.5)#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))


# Add a legend
legend("bottomleft", legend=levels(as.factor(V(network)$carac))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))



plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))




new_cols <- c("white", "red", "black")[membership(wc)]

V(g)$color <- ifelse(names(V(g)) %in% names(v.tf_families), adjustcolor("Orange", alpha.f = .6), adjustcolor("black", alpha.f = .6))

V(g)$type <- ifelse(names(V(g)) %in% names(v.tf_families), 0, 1)

E(g)$weight <- df.regulatoryNetwork.selection$rank / max(df.regulatoryNetwork.selection$rank)


vertex.label <- ifelse(names(V(g)) %in% df.regulatoryNetwork.selection$tf, names(V(g)), "")

# layout <-layout.fruchterman.reingold(g)
plot(g, layout = layout, edge.arrow.size=.2, vertex.label=vertex.label,edge.width=E(g)$weight, vertex.label.dist=0.1, vertex.label.color = "black")#, main = paste(v.tissues[s], " / ", names(l.grn_treatment[[s]])[i], " / ", v.domains[d], sep =""))
