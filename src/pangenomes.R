###########################################################################################################################################################################
################################################################ Pangenomes analysis ######################################################################################
###########################################################################################################################################################################


############################################################################################################################################################################

# Set work directory ##
workDir <- "/home/acari/Рабочий стол/HP_project/"
setwd(workDir)

# Import libraries and functions ##
source("src/functions.R")

set.seed(10)

# Set color pallete ##
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)

#############################################################################################################################################################################
erectale <- read.table("data/panphlan/panphlan_pcopri16_result.csv")
erectale <- as.data.frame(t(erectale))

pca <- prcomp(erectale, center = T)
mds.erectale16 <- metaMDS(erectale, distance = "bray")
# mds.erectale16.p <- as.data.frame(mds.erectale16$points)
mds.erectale16.p <- as.data.frame(pca$x[,c(1,2)])
mds.erectale16.p <- cbind(rownames(mds.erectale16.p), mds.erectale16.p)
colnames(mds.erectale16.p)[1] <- "sampleID"
rownames(mds.erectale16.p) <- 1:nrow(mds.erectale16.p)

mds.erectale16.hp <- mds.erectale16.p[which(!is.na(str_extract(mds.erectale16.p$sampleID, "HP"))),]

mds.erectale16.hp$sampleID <- sapply(str_split(mds.erectale16.hp$sampleID, "_"), function(x) x[length(x)-1])

mds.erectale16.hp <- merge(meta.data.sbs[c(1,2,11)], mds.erectale16.hp, by = 1)

mds.erectale16.hp$group <- "HP"
mds.erectale16.hp$days[mds.erectale16.hp$point_number == 1] <- 0
mds.erectale16.hp$days[mds.erectale16.hp$point_number == 2] <- 14

mds.erectale16.h <- mds.erectale16.p[which(is.na(str_extract(mds.erectale16.p$sampleID, "HP"))),]
mds.erectale16.h$patientID <- sapply(str_split(mds.erectale16.h$sampleID, "\\."), function(x) head(x, 1))
mds.erectale16.h$group <- "healthy"
mds.erectale16.h$days <- sapply(str_split(mds.erectale16.h$sampleID, "\\."), function(x) x[length(x)-1])
mds.erectale16.hp <- mds.erectale16.hp[-3]
colnames(mds.erectale16.hp) <- colnames(mds.erectale16.h[c(1,4,2,3,5,6)])

mds.erectale16.p <- rbind(mds.erectale16.h[c(1,4,2,3,5,6)], mds.erectale16.hp)
mds.erectale16.p$days <- as.numeric(mds.erectale16.p$days)
mds.erectale16.p <- mds.erectale16.p[order(mds.erectale16.p$days),]
mds.erectale16.p$days <- as.factor(mds.erectale16.p$days)

mds.erectale16.hp$days <- as.factor(mds.erectale16.hp$days)

svg("graphs/MDS.pcopri.svg", width = 5, height = 5)
ggplot(mds.erectale16.hp, aes(PC1, PC2, col = group, shape = days))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      # stat_ellipse(alpha = 0.8, aes(col = group, linetype = days), size = 0.7)+
      scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      scale_fill_manual(values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(values = c(16,2,2,0))+
      geom_line(mds.erectale16.hp, mapping = aes(PC1, PC2, group = patient_id), col = "black", alpha = 0.7)
dev.off()

#####################################################################################################################################################################################################

erectale <- read.table("data/panphlan/panphlan_pcopri16_result.csv")
erectale <- as.data.frame(t(erectale))

mds.erectale16 <- metaMDS(erectale, distance = "bray")
mds.erectale16.p <- as.data.frame(mds.erectale16$points)
mds.erectale16.p <- cbind(rownames(mds.erectale16.p), mds.erectale16.p)
colnames(mds.erectale16.p)[1] <- "sampleID"
rownames(mds.erectale16.p) <- 1:nrow(mds.erectale16.p)
mds.erectale16.p$sampleID <- sapply(str_split(mds.erectale16.p$sampleID, "_"), function(x) x[length(x)-1])
mds.erectale16.p <- merge(meta.data.sbs[c(1,2,11)], mds.erectale16.p, by = 1)
mds.erectale16.p$point_number <- as.factor(mds.erectale16.p$point_number)

svg("graphs/MDS.pcopri,svg", width = 5, height = 5)
ggplot(mds.erectale16.p, aes(MDS1, MDS2, col = point_number, shape = point_number))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      # stat_ellipse(alpha = 0.8, aes(col = group, linetype = days), size = 0.7)+
      scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      scale_fill_manual(values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(values = c(16,2,2,0))+
      geom_line(mds.erectale16.p, mapping = aes(MDS1, MDS2, group = patient_id), col = "black", alpha = 0.7)+
      guides(color=guide_legend(title="Time point"), shape=guide_legend(title="Time point"))
dev.off()
