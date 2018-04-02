##############################################################################################################################################################
################################################################ Functional data processing ##################################################################
##############################################################################################################################################################

# Description:

# 1. Import datasets - import metadata and functional data obtained by HUMAnN2 pipeline.
# 2. Make NMDS plot - make non-metric MDS plot using functional data (Bray-Curtis dissimilarity).
# 3. LefSe results processing - processing  tables of LefSe results and make relative abundances boxplots for detected features by LefSe analysis.

#############################################################################################################################################################

# Set work directory ##
workDir <- "/home/acari/Рабочий стол/HP_project/"
setwd(workDir)

# Import libraries and functions ##
source("src/functions.R")

set.seed(10)

# Set color pallete ##
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)

##############################################################################################################################################################
## 1. Import datasets ########################################################################################################################################
##############################################################################################################################################################

# Import metadata #
meta.data <- read.csv("data/meta_case.csv")

# Extract first and second time points sampleID
samples <- as.character(meta.data$patient_id[meta.data$point_number %in% c(1,2)])
samples <- samples[which(duplicated(samples))]
meta.sbs <- meta.data[meta.data$patient_id %in% samples & meta.data$point_number %in% c(1,2),]
names.list <- meta.sbs$sample_ID
meta.sbs <- meta.sbs[meta.sbs$sample_ID %in% names.list,]
meta.sbs$patient_id <- as.character(meta.sbs$patient_id)

Point_1 <- as.character(meta.sbs$sample_ID[meta.sbs$point_number == 1])
Point_2 <- as.character(meta.sbs$sample_ID[meta.sbs$point_number == 2])

##############################################################################################################################################################

# Import case data
df.path <- read.csv("data/case.path.hmn", row.names = 1)
df.path.sbs <- df.path[rownames(df.path) %in% c(Point_1, Point_2),]

df.GO <- read.csv("data/humann2/case.GO.csv", sep = "\t")
colnames(df.GO)[-1] <- sapply(str_split(colnames(df.GO)[-1], "\\_"), function(x) tail(x, 1))
rownames(df.GO) <- df.GO$GO
df.GO <- df.GO[colnames(df.GO) %in% c(Point_1, Point_2)]
df.GO.sbs <- df.GO[!colnames(df.GO) %in% "55HP.1"]
df.GO.sbs <- as.data.frame(t(df.GO.sbs))

df.GO <- read.csv("data/humann2/case.EGG.csv", sep = "\t")
colnames(df.GO)[-1] <- sapply(str_split(colnames(df.GO)[-1], "\\_"), function(x) tail(x, 1))
rownames(df.GO) <- df.GO$GO
df.GO <- df.GO[colnames(df.GO) %in% c(Point_1, Point_2)]
df.GO.sbs <- df.GO[!colnames(df.GO) %in% "55HP.1"]
df.GO.sbs <- as.data.frame(t(df.GO.sbs))

df.GO <- read.csv("data/humann2/case.KO.csv", sep = "\t")
colnames(df.GO)[-1] <- sapply(str_split(colnames(df.GO)[-1], "\\_"), function(x) tail(x, 1))
rownames(df.GO) <- df.GO$GO
df.GO <- df.GO[colnames(df.GO) %in% c(Point_1, Point_2)]
df.GO.sbs <- df.GO[!colnames(df.GO) %in% "55HP.1"]
df.GO.sbs <- as.data.frame(t(df.GO.sbs))
##############################################################################################################################################################
df.GO <- read.csv("data/humann2/case.KO.humann.csv", sep = "\t")
df.GO$SampleID <- sapply(str_split(df.GO$SampleID, "\\_"), function(x) x[length(x)-1])
df.GO <- df.GO[df.GO$SampleID %in% c(Point_1, Point_2),]
# df.GO.sbs <- df.GO[!df.GO$SampleID %in% "55HP.1",]
# df.GO.sbs <- as.data.frame(t(df.GO.sbs))
rownames(df.GO) <- 1:nrow(df.GO)
df.GO <- df.GO[-37,]
rownames(df.GO) <- df.GO$SampleID
df.GO <- df.GO[-1]
df.GO[is.na(df.GO)] <- 0

df.GO <- read.csv("data/humann2/case.KEGGpath.csv", sep = "\t")
df.GO$SampleID <- sapply(str_split(df.GO$SampleID, "\\_"), function(x) x[length(x)-1])
df.GO <- df.GO[df.GO$SampleID %in% c(Point_1, Point_2),]
# df.GO.sbs <- df.GO[!df.GO$SampleID %in% "55HP.1",]
# df.GO.sbs <- as.data.frame(t(df.GO.sbs))
rownames(df.GO) <- 1:nrow(df.GO)
df.GO <- df.GO[-37,]
rownames(df.GO) <- df.GO$SampleID
df.GO <- df.GO[-1]
df.GO[is.na(df.GO)] <- 0
##############################################################################################################################################################
## 1.  lol  ##################################################################################################################################################
##############################################################################################################################################################
mds.GO <- metaMDS(df.GO.sbs, "bray")
mds.GO.p <- as.data.frame(mds.GO$points)
mds.GO.p <- merge(meta.data.sbs[c(1,2,10)], cbind(rownames(mds.GO.p), mds.GO.p), by = 1)
mds.GO.p$Time.point.number <- as.factor(mds.GO.p$Time.point.number)

envfit(mds.GO, df.GO, permutations = 1000)

ggplot(mds.GO.p, aes(MDS1, MDS2, col = Time.point.number))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      stat_ellipse()+
      theme(legend.position = "bottom")+
      guides(col=guide_legend(title="Time point"))+
      stat_ellipse(alpha = 0.8, aes(col = Time.point.number), size = 0.7)+
      scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      geom_line(mds.GO.p, mapping = aes(MDS1, MDS2, group = Patient_ID), col = "black", alpha = 0.5)

# df.GO_f <- df.GO.sbs[,apply(df.GO.sbs, 2, var, na.rm=TRUE) != 0]
# df.path.sbs_f <- df.path.sbs[,apply(df.GO, 2, var, na.rm=TRUE) != 0]
# 
# pca <- prcomp(df.GO.sbs, center = T)
# summary(pca)
# 
# # ggbiplot(pca, obs.scale = 1, var.scale = 1, ellipse = T, circle = TRUE, varname.adjust = 0.01, xlim = c(-4000, 2000))+
# #       theme_linedraw()
# 
# pca.point <- as.data.frame(pca$x)[c(1,2)]
# pca.point <- cbind(rownames(pca.point), pca.point)
# colnames(pca.point)[1] <- "SampleID"
# pca.point <- merge(meta.data.sbs[c(1,2,10)], pca.point, by = 1)
# pca.point$Time.point.number <- as.factor(pca.point$Time.point.number)
# 
# ggplot(pca.point, aes(PC1, PC2, col = Time.point.number))+
#       geom_point()+
#       theme_linedraw()+
#       stat_ellipse()+
#       theme(legend.position = "bottom")+
#       xlab("PC1 (96.4%)")+
#       ylab("PC1 (1.2%)")

mds.GO <- metaMDS(df.path.sbs, "bray")
mds.GO.p <- as.data.frame(mds.GO$points)
mds.GO.p <- merge(meta.data.sbs[c(1,2,10)], cbind(rownames(mds.GO.p), mds.GO.p), by = 1)
mds.GO.p$Time.point.number <- as.factor(mds.GO.p$Time.point.number)

ggplot(mds.GO.p, aes(MDS1, MDS2, col = Time.point.number))+
      geom_point()+
      theme_linedraw()+
      stat_ellipse()

df.path.sbs <- df.path.sbs[order(rownames(df.path.sbs)),]
# df.GO.sbs <- df.GO.sbs[order(rownames(df.GO.sbs)),]
# df.KO <- df.GO[order(rownames(df.GO)),]
df.GO <- df.GO[order(rownames(df.GO)),]

adonis(df.path.sbs ~ ., data=meta.data.sbs[c(3:10)], permutations=9999)
adonis(df.GO ~ ., data=meta.data.sbs[c(3:10)][-1,], permutations=9999)

df.GO.lefse <- as.data.frame(t(merge(meta.data.sbs[c(1,10)], cbind(rownames(df.GO), df.GO), by = 1)[-1]))
df.GO.lefse[1,] <- paste0("p_", df.GO.lefse[1,])

# df.GO.lefse <- as.data.frame(t(merge(meta.data.sbs[c(1,10)], cbind(rownames(df.GO.sbs), df.GO.sbs), by = 1)[-1]))
# df.GO.lefse[1,] <- paste0("p_", df.GO.lefse[1,])

df.GO.lefse <- as.data.frame(t(merge(meta.data.sbs[c(1,10)], cbind(rownames(df.path.sbs), df.path.sbs), by = 1)[-1]))
df.GO.lefse[1,] <- paste0("p_", df.GO.lefse[1,])

write.table(df.GO.lefse, "output/lefse.MCpath.txt", quote = F, col.names = F, sep = "\t")

write.table(t(df.gen.sbs[rownames(df.gen.sbs) %in% Point_1,]), "/home/acari/sparCC/data_genera/p1.genera.sparCC.txt", quote = F, sep = "\t")
write.table(t(df.gen.sbs[rownames(df.gen.sbs) %in% Point_2,]), "/home/acari/sparCC/data_genera/p2.genera.sparCC.txt", quote = F, sep = "\t")
##############################################################################################################################################################
## 2.  Make NMDS plot  #######################################################################################################################################
##############################################################################################################################################################

mds.path <- metaMDS(df.path, "bray", k = 2)
mds.path.points <- as.data.frame(mds.path$points)
mds.path.points <- merge(meta.data[c(1,2,11)], cbind(rownames(mds.path.points), mds.path.points), by = 1)

mds.path.points$point_number <- as.factor(mds.path.points$point_number)

svg("graphs/MDS.pathways.svg")
ggplot(mds.path.points, aes(MDS1, MDS2, col = point_number, shape = point_number))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      scale_color_manual(values = my_palette[c(30,200,150)])+
      scale_fill_manual(values = my_palette[c(30,200,150)])+
      scale_shape_manual(values = c(2,16,0))+
      guides(color=guide_legend(title="Time point"), 
             shape=guide_legend(title="Time point"),
             fill = F)+
      stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = point_number))
dev.off()

##############################################################################################################################################################
## 3. LefSe results processing ###############################################################################################################################
##############################################################################################################################################################

# Import LefSe table ##
lefse.path <- read.table("data/lefse_output/case.path.lefse.out", sep = "\t", stringsAsFactors = F)

colnames(lefse.path) <- c("Feature", "LogMaxMean", "Class", "LDA_score", "p.value")
lefse.path <- lefse.path[lefse.path$Class == "point_1",]
lefse.path$p.value <- as.numeric(lefse.path$p.value)
lefse.path <- lefse.path[lefse.path$p.value < 0.1,]
lefse.path <- lefse.path[order(lefse.path$LDA_score),]

write.table(lefse.path, "output/lefse.pathways.table.txt", quote=F, sep='\t', row.names = F)

# pname1 <- sapply(lapply(str_split(lefse.path$Feature, "_"), function(x) head(x,2)), function(x) paste0(x, collapse = "-"))

# svg("graphs/lefse.path.svg", width = 5, height = 2)
# par(mar = c(4.5,12.5,0,0.5)) 
# barplot(name = pname, 
#         height = lefse.path$LDA_score, 
#         horiz = T, 
#         las = 1,
#         col = my_palette[10],
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()

# # Boxplots
# df.path.sbs <- merge(meta.sbs2[c(1,2,11)], cbind(rownames(df.path), df.path[pname]), by = 1)
# pname2 <- c("PWY.5104..L.isoleucine.biosynthesis.IV", 
#             "PWY.6737..starch.degradation.V", 
#             "GLYCOGENSYNTH.PWY..glycogen.biosynthesis.I..from.ADP.D.Glucose.")
# colnames(df.path.sbs)[4:6] <- pname1
# df.path.sbs$point_number <- as.factor(df.path.sbs$point_number)
# df.path.sbs <- melt(df.path.sbs, id.vars = c("sample_ID", "patient_id", "point_number"))
# 
# svg("graphs/boxplots.pathways.3.svg", width = 6, height = 3)
# ggboxplot(df.path.sbs, x = "point_number", y = "log(value)", fill = "point_number", add = "jitter")+
#       theme_linedraw()+
#       facet_wrap(~variable, nrow = 1)+
#       # stat_compare_means(aes(group = point_number), comparisons = my_comparisons, method = "wilcox.test", paired = F)+ 
#       # stat_compare_means(label.y = 4.5, size = 3, paired = T)+
#       guides(fill=FALSE)+
#       scale_fill_manual(values = my_palette[c(1,50,130)])+
#       xlab("Time point")+
#       ylab("log (RPK, %)")
# dev.off()

##############################################################################################################################################################
######################################################################## THE END #############################################################################
##############################################################################################################################################################