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
meta.data <- read.csv("data/meta.case.txt")

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