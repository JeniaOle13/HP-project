#############################################################################################################################################################
################################################################ AR-genes data processing ###################################################################
#############################################################################################################################################################

# Description: 

# 1. Import datasets - import metadata and data obtained by bowtie and MEGARes database.
# 2. Make NMDS plot - make non-metric MDS plot using class of antibiotics genes data (Bray-Curtis dissimilarity).
# 3. LefSe results processing - processing  tables of LefSe results and make relative abundances boxplots for detected features by LefSe analysis.
# 4. Make boxplots - make relative abundance boxplots for features discriminating time points.

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

# Import metadata ##
meta.data <- read.csv("data/meta.case.txt")

# Extract first and second time points sampleID ##
samples <- as.character(meta.data$patient_id[meta.data$point_number %in% c(1,2)])
samples <- samples[which(duplicated(samples))]
meta.sbs <- meta.data[meta.data$patient_id %in% samples & meta.data$point_number %in% c(1,2),]
names.list <- meta.sbs$sample_ID
meta.sbs <- meta.sbs[meta.sbs$sample_ID %in% names.list,]
meta.sbs$patient_id <- as.character(meta.sbs$patient_id)

Point_1 <- as.character(meta.sbs$sample_ID[meta.sbs$point_number == 1])
Point_2 <- as.character(meta.sbs$sample_ID[meta.sbs$point_number == 2])

##############################################################################################################################################################

# Import case data ##
df.group <- read.csv("data/case.group.mgs", sep = " ")
df.mech <- read.csv("data/case.mech.mgs", row.names = 1)
df.class <- read.csv("data/case.class.mgs", row.names = 1)

# Import control data ##
class.control <- read.csv("data/rus.control.class.txt", row.names = 1)

##############################################################################################################################################################
## 2.  Make NMDS plot  #######################################################################################################################################
##############################################################################################################################################################

# class.all <- as.data.frame(UniteMatrices(as.matrix(df.class), as.matrix(class.control)))
# 
# class.all.sbs <- class.all[,which(apply(class.all, 2, max) > 50000)]
# 
# mds <- metaMDS(class.all.sbs, "bray")
# fit <- envfit(mds, class.all.sbs, permutations = 1000)
# 
# mds.points <- as.data.frame(mds$points)
# 
# mds.points$group[rownames(mds.points) %in% rownames(class.control)] <-"RUS_controls" 
# mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <-"first_point" 
# mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <-"second_point" 
# mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <-"third_point" 
# 
# 
# mds.points$color[rownames(mds.points) %in% rownames(class.control)] <- my_palette[200]
# mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <- my_palette[1]
# mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <- my_palette[30]
# mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <- my_palette[150]
# 
# mds.points$pch[rownames(mds.points) %in% rownames(class.control)] <- 8
# mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <- 2
# mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <- 16
# mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <- 0
# mds.points$group <- as.factor(mds.points$group)
# 
# plotPath <-paste('graphs/MDS.class.svg',sep='')
# 
# svg(plotPath, width = 8, height = 8)
# plot(x = mds.points[, 1], y = mds.points[, 2], pch = mds.points[,5], lwd = 2, col = mds.points[,4], # control samples
#      xlab = 'NMDS Axis 1', ylab = 'NMDS Axis 2',
#      xlim = range(mds.points[,1]),
#      ylim = range(mds.points[,2]))
# 
# plot(fit, p.max = 0.001, col = "black", cex = 0.8)
# 
# legend('topright', legend = c("RUS controls", 'Before eradication', '0-2 day after eradication', '20-120 day after eradication'), bty = 'n',
#        # legend('topleft', legend = c('Patients','Russia control',
#        #                              'Denmark control', 'China control', 'USA control'), bty = 'n',
#        col= c(my_palette[200], my_palette[1], my_palette[30], my_palette[150]), 
#        pch=c(8,2,16,0), lwd=c(2,2,2), 
#        lty =c(NA,NA,NA), cex = .9)
# dev.off()

##############################################################################################################################################################
## 3. LefSe results processing ###############################################################################################################################
##############################################################################################################################################################

# Groups ##

lefse.gr <- read.table("data/lefse_output/case.group.lefse.out", sep = "\t", stringsAsFactors = F)
colnames(lefse.gr) <- c("Feature", "LogMaxMean", "Class", "LDA_score", "p.value")

lefse.gr1 <- lefse.gr[lefse.gr$Class == "point_2",]
lefse.gr2 <- lefse.gr[lefse.gr$Class == "point_1",]
lefse.gr <- rbind(lefse.gr1, lefse.gr2)

lefse.gr$p.value <- as.numeric(lefse.gr$p.value)
lefse.gr <- lefse.gr[lefse.gr$p.value < 0.1,]
lefse.gr <- lefse.gr[order(lefse.gr$LDA_score),]

write.table(lefse.gr, "output/lefse.group.table.txt", quote=F, sep='\t', row.names = F)

# lefse.gr$LDA_score[c(6:7)] <- lefse.gr$LDA_score[c(6:7)]*-1
# 
# svg("graphs/lefse.gr.svg", width = 5, height = 2)
# par(mar = c(4.5,7.5,0,0.5)) 
# barplot(name = lefse.gr$Feature, 
#         height = lefse.gr$LDA_score, 
#         horiz = T, 
#         las = 1,
#         col = c(my_palette[50], my_palette[50], my_palette[50], my_palette[50], my_palette[50], my_palette[10], my_palette[10]),
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()

# Mechanisms ## 

lefse.mech <- read.table("data/lefse_output/case.mech.lefse.out", sep = "\t", stringsAsFactors = F)
colnames(lefse.mech) <- c("Feature", "LogMaxMean", "Class", "LDA_score", "p.value")

lefse.mech1 <- lefse.mech[lefse.mech$Class == "point_2",]
lefse.mech2 <- lefse.mech[lefse.mech$Class == "point_1",]
lefse.mech <- rbind(lefse.mech1, lefse.mech2)

lefse.mech$p.value <- as.numeric(lefse.mech$p.value)
lefse.mech <- lefse.mech[lefse.mech$p.value < 0.1,]
lefse.mech <- lefse.mech[order(lefse.mech$LDA_score),]

write.table(lefse.mech, "output/lefse.mech.table.txt", quote=F, sep='\t', row.names = F)

# lefse.mech$LDA_score[4] <- lefse.mech$LDA_score[4]*-1
# 
# svg("graphs/lefse.mech.svg", width = 8, height = 2)
# par(mar = c(4.5,21.5,0,0.5)) 
# barplot(name = lefse.mech$Feature, 
#         height = lefse.mech$LDA_score, 
#         horiz = T, 
#         las = 1,
#         col = c(my_palette[10], my_palette[10], my_palette[10], my_palette[200]),
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()

# Classes ##

lefse.class <- read.table("data/lefse_output/case.class.lefse.out", sep = "\t", stringsAsFactors = F)
colnames(lefse.class) <- c("Feature", "LogMaxMean", "Class", "LDA_score", "p.value")

lefse.class1 <- lefse.class[lefse.class$Class == "point_2",]
lefse.class2 <- lefse.class[lefse.class$Class == "point_1",]
lefse.class <- rbind(lefse.class1, lefse.class2)

lefse.class$p.value <- as.numeric(lefse.class$p.value)
lefse.class <- lefse.class[lefse.class$p.value < 0.1,]
lefse.class <- lefse.class[order(lefse.class$LDA_score),]

write.table(lefse.class, "output/lefse.class.table.txt", quote=F, sep='\t', row.names = F)

# lefse.class$LDA_score[2] <- lefse.class$LDA_score[2]*-1
# 
# svg("graphs/lefse.class.svg", width = 8, height = 2)
# par(mar = c(4.5,6.5,0,0.5)) 
# barplot(name = lefse.class$Feature, xlim = c(-6,6),
#         height = lefse.class$LDA_score, 
#         horiz = T, 
#         las = 1,
#         col = c(my_palette[10], my_palette[50]),
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()

##############################################################################################################################################################
## 4.  Make Boxplots  ########################################################################################################################################
##############################################################################################################################################################

list.class <- colnames(df.class)
list.class <- list.class[-c(7, 9, 11:16)]
df.class.sbs <- df.class[list.class]
df.class.sbs <- merge(meta.data[c(1,2,11)], cbind(rownames(df.class.sbs), df.class.sbs), by = 1)
df.class.sbs <- df.class.sbs[df.class.sbs$sample_ID != "29HP",]

df.class.sbs <- melt(df.class.sbs, id.vars = c("sample_ID", "patient_id", "point_number"))
df.class.sbs$point_number <- as.factor(df.class.sbs$point_number)
df.class.sbs$variable <- as.character(df.class.sbs$variable)
df.class.sbs$variable[df.class.sbs$variable == "Multi.drug.resistance"] <- "MDR"
df.class.sbs$variable <- as.factor(df.class.sbs$variable)

svg("graphs/boxplots.class1.svg", width = 12, height = 4)
ggboxplot(df.class.sbs, x = "point_number", y = "log(value)", fill = "point_number")+
      theme_linedraw(base_size = 15)+
      facet_wrap(~variable, nrow = 1)+
      # stat_compare_means(aes(group = point_number), comparisons = my_comparisons, method = "wilcox.test", paired = F)+ 
      # stat_compare_means(label.y = 4.5, size = 3, paired = T)+
      guides(fill=FALSE)+
      scale_fill_manual(values = my_palette[c(1,50,130)])+
      xlab("Time point")+
      ylab("log (CPM, %)")
dev.off()

df.class.sbs$sample_ID <- as.character(df.class.sbs$sample_ID)

svg("graphs/leniplot.class.svg", width = 7, height = 4)
ggplot(df.class.sbs[df.class.sbs$variable %in% c("MLS", "betalactams", "Tetracyclines"),], aes(x = point_number, y = log(value), col = patient_id, group = patient_id))+
      geom_line()+
      geom_point()+
      theme_linedraw(base_size = 17)+
      facet_wrap(~variable, nrow = 1)+
      # stat_compare_means(aes(group = point_number), comparisons = my_comparisons, method = "wilcox.test", paired = F)+ 
      # stat_compare_means(label.y = 4.5, size = 3, paired = T)+
      guides(col=FALSE)+
      scale_fill_manual(values = my_palette[c(1,50,130)])+
      xlab("Time point")+
      ylab("log (CPM, %)")
dev.off()

##############################################################################################################################################################
######################################################################## THE END #############################################################################
##############################################################################################################################################################