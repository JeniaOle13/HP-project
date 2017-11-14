#############################################################################################################################################################
################################################################ ef - isolates  #############################################################################
#############################################################################################################################################################

# Description

# 1. Heatmap for ARGs and VFs - 
# 2. Make tree - 
# 3. Make lineplots - 

#############################################################################################################################################################

# Import additional libraries ##
library(ggtree)
library(phangorn)
library(colorspace)
library(pheatmap)

# Set work directory ##
workDir <- "/home/acari/Рабочий стол/HP_project/"
setwd(workDir)

# Import libraries and functions ##
source("src/functions.R")

set.seed(10)

# Set color pallete ##
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)

#############################################################################################################################################################
## 1. Heatmap for ARGs and VFs ##############################################################################################################################
#############################################################################################################################################################

# ARGs ## 

setwd("data/blast_results/BLAST_MEGARes/")

list.f <- list.files(pattern = "out")

full.data <- NULL
for (i in list.f){
      df <- fread(i)
      df <- df[df$V3 > 80,]
      df <- as.data.frame(df)
      df$V15 <- df$V4/df$V14
      df <- df[df$V15 > 0.8,]
      df$V16 <- sapply(str_split(df$V2, "\\|"), function(x) tail(x, 1))
      
      dt <- NULL
      for (j in unique(df$V16)){
            df.sbs <- df[df$V16 == j,]
            df.sbs <- df.sbs[df.sbs$V3 == max(df.sbs$V3),]
            dt <- rbind(df.sbs, dt)
      }
      
      genes.list <- data.frame(names = unique(dt$V16), pr = 1, sample = gsub(".dmd.out", "", i))
      full.data <- rbind(full.data, genes.list) 
      
}

full.data <- full.data[full.data$names != "RequiresSNPConfirmation",]

full.data <- spread(full.data, names, pr)
full.data[is.na(full.data)] <- 0
rownames(full.data) <- full.data$sample
full.data <- full.data[-1]

my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)

point_1 <- c("Hp_21-11", "Hp_21-21", "Hp_23-14", "Hp_23-9", "Hp_5-10", "Hp_5-7")
point_2 <- c("Hp_22-10", "Hp_22-12", "Hp_24-1", "Hp_24-3", "Hp_6-10", "Hp_6-9")
point_3 <- c("Hp_7-8")

samples <- as.data.frame(rbind(cbind(sample = point_1, value = "first_point"), 
                               cbind(sample = point_2, value = "second_point"), 
                               cbind(sample = point_3, value = "third_point")))

rownames(samples) <- samples$sample
samples <- samples[-1]
colnames(samples) <- "Time point"

samples$patient[rownames(samples) %in% c("Hp_5-10", "Hp_5-7", "Hp_6-10", "Hp_6-9", "Hp_7-8")] <- "HP003"
samples$patient[rownames(samples) %in% c("Hp_21-11", "Hp_21-21", "Hp_22-10", "Hp_22-12")] <- "HP009"
samples$patient[rownames(samples) %in% c("Hp_23-9", "Hp_23-14", "Hp_24-1", "Hp_24-3")] <- "HP010"

samples.sbs <- samples[!rownames(samples) %in% c("Hp_6-9", "Hp_21-11", "Hp_21-21", "Hp_22-10", "Hp_22-12", "Hp_24-1"),]
samples.sbs <- samples.sbs[order(rownames(samples.sbs)),]

rownames(full.data) <- sapply(str_split(rownames(full.data) , "\\."), function(x) head(x, 1))
rownames(full.data) <- sapply(lapply(str_split(rownames(full.data), "_"), function(x) head(x, 2)), function(x) paste(x, collapse = "_"))
full.data.sbs <- full.data[!rownames(full.data) %in% c("Hp_6-9", "Hp_21-11", "Hp_21-21", "Hp_22-10", "Hp_22-12", "Hp_24-1"),]
full.data.sbs <- full.data.sbs[,colSums(full.data.sbs) > 0]
lol1 <- full.data.sbs

# pheatmap(t(full.data.sbs[!colnames(full.data.sbs) %in% c("ERM", "ERMC", "ERMR", "ERMT")]), 
#          main="", 
#          color= c(my_palette[80], my_palette[200]),
#          annotation_col = samples.sbs, 
#          cluster_rows = F, 
#          cluster_cols = F, 
#          show_colnames = F,
#          fontsize_row=13, 
#          gaps_col =  c(3,3,3), 
#          legend = F) 
# #cellwidth=16.7)

#############################################################################################################################################################

setwd("../BLAST_VFDB/")

list.f <- list.files(pattern = "out")

full.data <- NULL
for (i in list.f){
      df <- fread(i)
      df <- as.data.frame(df)
      df$V15 <- df$V4/df$V14
      df <- df[df$V15 > 0.8,]
      # df$V16 <- sapply(str_split(df$V2, "\\|"), function(x) tail(x, 1))
      
      dt <- NULL
      for (j in unique(df$V2)){
            df.sbs <- df[df$V2 == j,]
            df.sbs <- df.sbs[df.sbs$V3 == max(df.sbs$V3),]
            dt <- rbind(df.sbs, dt)
      }
      
      genes.list <- data.frame(names = unique(dt$V2), pr = 1, sample = gsub(".blast.out", "", i))
      full.data <- rbind(full.data, genes.list)    
} 

full.data <- spread(full.data, names, pr)
full.data[is.na(full.data)] <- 0
rownames(full.data) <- full.data$sample
full.data <- full.data[-1]

rownames(full.data) <- sapply(lapply(str_split(rownames(full.data), "_"), function(x) head(x, 2)), function(x) paste(x, collapse = "_"))
full.data.sbs <- full.data[!rownames(full.data) %in% c("Hp_6-9", "Hp_21-11", "Hp_21-21", "Hp_22-10", "Hp_22-12", "Hp_24-1"),]
colnames(full.data.sbs) <- sapply(str_split(colnames(full.data.sbs), "\\("), function(x) head(x, 1))
lol2 <- full.data.sbs

full.data.sbs <- cbind(lol1,lol2[,colSums(lol2) > 0])
full.data.sbs <- full.data.sbs[!colnames(full.data.sbs) %in% c("ERM", "ERMC", "ERMR", "ERMT")]

genes.annot <- data.frame(names = colnames(full.data.sbs), annotation = "genes")
genes.annot$annotation <- as.character(genes.annot$annotation)
genes.annot$annotation[1:8] <- "ARGs"
genes.annot$annotation[9:13] <- "VF"
rownames(genes.annot) <- genes.annot$names
genes.annot <- genes.annot[-1]

colnames(genes.annot) <- "Genes groups"

colnames(samples.sbs)[2] <- "Patiend ID"

svg("../../../graphs/heatmaps.isolate.svg")
pheatmap(t(full.data.sbs), 
         main="", 
         color= c(my_palette[80], my_palette[200]),
         annotation_col = samples.sbs,
         annotation_row = genes.annot,
         cluster_rows = F, 
         cluster_cols = F, 
         show_colnames = F,
         fontsize_row=13, 
         gaps_col =  c(3,3,3), 
         legend = F, annotation_names_row = F) 
dev.off()

#############################################################################################################################################################
## 2. Make tree #############################################################################################################################################
#############################################################################################################################################################

tree.sh <- read.tree("../../../data/parsnp.tree")
tree.sh$tip.label <- gsub(".fasta", "", tree.sh$tip.label)
tree.sh$tip.label <- gsub("\\'", "", tree.sh$tip.label)
tree.sh$tip.label <- gsub("\\.ref", "", tree.sh$tip.label)

cls <- list(point_1 = c("Hp_23-9", "Hp_23-14", "Hp_5-10", "Hp_5-7"),
            point_2 = c("Hp_6-10", "Hp_24-3"),
            point_3 = c("Hp_7-8"),
            ref = c("EFE10021"))

lol <- groupOTU(tree.sh, cls)

p1 <- data.frame(p1 = c(0.63, 0.08, 0.09), 
                 p2 = c(0.15, 0.05,0.13), 
                 p3 = c(0.11, 0.82, 0.24),
                 p4 = c(0.11, 0.05, 0.54))

svg("../../../graphs/tree.isolate.svg")
ggtree(lol, aes(color=group, linetype=group), size = 2, branch.length='none')+
      geom_tiplab(size=6, vjust=-0.3, hjust = T)+
      scale_color_manual(values=c(my_palette[c(1,30,170,200)]))+ 
      guides(color=FALSE)
dev.off()

#############################################################################################################################################################
## 3. Make lineplots ########################################################################################################################################
#############################################################################################################################################################

Hp_5 <- c(0.63,0.15,0.11,0.11)
Hp_6 <- c(0.08,0.05,0.82,0.05)
Hp_7 <- c(0.09,0.13,0.24,0.54)

Hp_003 <- as.data.frame(rbind(Hp_5, Hp_6, Hp_7))
colnames(Hp_003) <- c(
      "TP1 (Hp_5-10)",
      "TP1 (Hp_5-7)",
      "TP2 (Hp_6-10)", 
      "TP3 (Hp_7-8)"
)

rownames(Hp_003) <- c(
      "TP1 (Hp_5)",
      "TP2 (Hp_6)",
      "TP3 (Hp_7)"
)

Hp_003 <- melt(t(Hp_003))
Hp_003$X3 <- c(1,1,1,1,2,2,2,2,3,3,3,3)
Hp_003$X3 <- as.factor(Hp_003$X3)

Hp_003$X4 <- c("Hp_5-10", "Hp_5-7", "Hp_6-10", "Hp_7-8", "Hp_5-10", "Hp_5-7", "Hp_6-10", "Hp_7-8", "Hp_5-10", "Hp_5-7", "Hp_6-10", "Hp_7-8")
Hp_003$X4 <- as.factor(Hp_003$X4)
colnames(Hp_003)[5] <- "Strain"

svg("../../../graphs/lineplot.isolates1.svg", width = 4, height = 4)
ggplot(Hp_003, aes(X3, value, group = Strain, col = Strain))+
      #geom_point(size = 2)+
      geom_line(size = 1.2, aes(linetype = Strain))+
      theme_linedraw()+
      theme(legend.position="bottom")+
      xlab("Time point")+
      ylab("Relative abundance, %")+
      scale_color_manual(values = my_palette[c(1,1,170,30)])
dev.off()

#############################################################################################################################################################

Hp_23 <- c(0.590, 0.190, 0.220)
Hp_24 <- c(0.007, 0.003, 0.989)

Hp_010 <- as.data.frame(rbind(Hp_23, Hp_24))

rownames(Hp_010) <- c(
      "TP1 (Hp_23)",
      "TP2 (Hp_24)"
)

colnames(Hp_010) <- c(
      "TP1 (Hp_23-14)",
      "TP1 (Hp_23-9)", 
      "TP2 (Hp_24-3)"
)

Hp_010 <- melt(t(Hp_010))
Hp_010$X3 <- c(1,1,1,2,2,2)
Hp_010$Strain <- c("Hp_23-14", "Hp_23-9", "Hp_24-3", "Hp_23-14", "Hp_23-9", "Hp_24-3")
Hp_010$Strain<- as.factor(Hp_010$Strain)
Hp_010$X3 <- as.factor(Hp_010$X3)

svg("../../../graphs/lineplot.isolates2.svg", width = 4, height = 4)
ggplot(Hp_010, aes(X3, value, group = Strain, col = Strain))+
      #geom_point(size = 2)+
      geom_line(size = 1.2, aes(linetype = Strain))+
      theme_linedraw()+
      theme(legend.position="bottom")+
      xlab("Time point")+
      ylab("Relative abundance, %")+
      scale_color_manual(values = my_palette[c(1,1,170,200)])
dev.off()

##############################################################################################################################################################
######################################################################## THE END #############################################################################
##############################################################################################################################################################