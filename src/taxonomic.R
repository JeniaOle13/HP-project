#############################################################################################################################################################
################################################################ Taxonomic data processing ##################################################################
#############################################################################################################################################################

# Description:

# 1. Import datasets - import metadata and taxonomic data obtained by MetaPhlAn2 (including in HUMAnN2 pipeline).
# 2. Alpha-diversity - calculation alpha-diversity index (Sahnnon index), Wilcoxon statistics calculation between first and second time points 
# and make boxplot.
# 3. Top taxa calculation - calculation MEANs and SDs value for most prevalent genera and species relative abundance.
# 4. MDS family level - make non-metric MDS plot using famuly taxonomic data (Bray-Curtis dissimilarity). 
# 5. Heatmaps (genera and species levels) - make heatmaps using species and genera taxonomic levels (ward linkage for hierarchical clusterization).
# 6. LefSe results processing - processing  tables of LefSe results and make relative abundances boxplots for detected features by LefSe analysis. 
# 7. Opportunistic and pathogens diff between first and second time points - Wilcoxon testing and ploting data of the relative abundance of 
# opportunistic bacteria.
# 8. MDSs for different taxonomic levels (firs and second time points) - make non-metric MDS (NMDS) plots for different taxonomic levels betweent 
# first and second time points.

##############################################################################################################################################################

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
df.org <- read.csv("data/case.org.mpa", row.names = 1)
df.gen <- read.csv("data/case.gen.mpa", row.names = 1)
df.family <- read.table("data/case.family.mpa")

# Import control data
control.family <- read.table("data/rus.control.family.mpa")

##############################################################################################################################################################
## 2. Alpha-diversity  #######################################################################################################################################
##############################################################################################################################################################

# Shannon index calculation
alpha.div <- data.frame(diversity(df.org, index = "shannon"))
colnames(alpha.div) <- "Shannon_index"
alpha.div <- merge(cbind(rownames(alpha.div), alpha.div), meta.data[c(1,11)], by = 1)
alpha.div$point_number <- as.factor(alpha.div$point_number)
colnames(alpha.div)[1] <- "sampleID"

alpha.div.first <- alpha.div$Shannon_index[as.character(alpha.div$sampleID) %in% Point_1]
alpha.div.second <- alpha.div$Shannon_index[as.character(alpha.div$sampleID) %in% Point_2]

# Wilcoxon signed-rank paired test
w.alpha.div <- wilcox.test(alpha.div.first,
                           alpha.div.second, paired = T,
                           alternative = "greater")
print(w.alpha.div)

# Make boxplot
svg("graphs/alpha-div.svg", width = 3, height = 5)
ggplot(alpha.div, aes(point_number, Shannon_index, fill = point_number))+
      geom_boxplot(alpha=0.7, size = 1.2, outlier.alpha = 0)+
      theme_linedraw()+
      scale_fill_manual(values = my_palette[c(1,50,130)])+
      theme(legend.position="none")+
      xlab("Time point")+
      ylab("Shannon index")
dev.off()

##############################################################################################################################################################
## 3. Top taxa calculation ###################################################################################################################################
##############################################################################################################################################################

# Genus level #

df <- NULL
for (i in 1:3){
      ff <- data.matrix(df.gen[rownames(df.gen) %in% meta.data$sample_ID[meta.data$point_number == i],])
      tops <- ff[,order(colMaxs(ff), decreasing = T)]
      tops <- tops[,which(colMaxs(tops) > 1)]
      top <- data.frame(Genus = colnames(tops), 
                        Means = round(colMeans(tops),2), 
                        Sds = round(apply(tops, 2, sd),2))
      top$point_number <- i
      df <- rbind(df, top)
}
rownames(df) <- 1:length(df$Genus)

write.table(df, 'output/top.gen.txt', quote=F, sep='\t', row.names = F)

# Species level #

df <- NULL
for (i in 1:3){
      ff <- data.matrix(df.org[rownames(df.org) %in% meta.data$sample_ID[meta.data$point_number == i],])
      tops <- ff[,order(colMaxs(ff), decreasing = T)]
      tops <- tops[,which(colMaxs(tops) > 1)]
      top <- data.frame(Species = colnames(tops), 
                        Means = round(colMeans(tops),2), 
                        Sds = round(apply(tops, 2, sd),2))
      top$point_number <- i
      df <- rbind(df, top)
}
rownames(df) <- 1:length(df$Species)

write.table(df, 'output/top.org.txt', quote=F, sep='\t', row.names = F)

#############################################################################################################################################################
## 4. MDS family level ######################################################################################################################################
#############################################################################################################################################################

# Merge case and control tables
df.all <- UniteMatrices(as.matrix(df.family), as.matrix(control.family))
df.all <- as.data.frame(df.all)
df.all <- df.all[!colnames(df.all) %in% "Viruses_noname"]
df.all <- df.all[,which(apply(df.all, 2, max) > 1)]

# building MDS
df.all.sbs <- df.all[,which(apply(df.all , 2, max) > 5)]
mds <- metaMDS(df.all.sbs, "bray")
fit <- envfit(mds, df.all.sbs, permutations = 1000)

mds.points <- as.data.frame(mds$points)
#mds.points <- merge(meta.data[c(1,2,11)], cbind(rownames(mds.points), mds.points), by = 1)

mds.points$group[rownames(mds.points) %in% rownames(control.family)] <-"Controls*" 
mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <-"first_point" 
mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <-"second_point" 
mds.points$group[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <-"third_point" 


mds.points$color[rownames(mds.points) %in% rownames(control.family)] <- my_palette[200]
mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <- my_palette[1]
mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <- my_palette[30]
mds.points$color[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <- my_palette[150]

mds.points$pch[rownames(mds.points) %in% rownames(control.family)] <- 8
mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 1]] <- 2
mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 2]] <- 16
mds.points$pch[rownames(mds.points) %in% meta.data$sample_ID[meta.data$point_number == 3]] <- 0
mds.points$group <- as.factor(mds.points$group)

df.fit <- as.data.frame(fit$vectors$arrows)
df.fit$r <- fit$vectors$r
df.fit$pvals <- fit$vectors$pvals

df.fit <- df.fit[df.fit$r > 0.2 & df.fit$pvals < 0.05,]

# Save MDS family plot
svg("graphs/MDS.family.svg", width = 5.2, height = 5.2)
ggplot(mds.points, aes(MDS1, MDS2, col = group, shape = group))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      scale_color_manual(values = my_palette[c(200,1,30,150)])+
      scale_fill_manual(values = my_palette[c(200,1,30,150)])+
      scale_shape_manual(values = c(8,2,16,0))+
      guides(color=guide_legend(title="Group"), 
             shape=guide_legend(title="Group"),
             col = F)+
      stat_ellipse(alpha = 0.5, aes(col = group))+
      annotate("text", x = 0.05180919, y = 0.99865700, label = "Bacteroidaceae")+
      annotate("text", x = 0.99709076, y = -0.07622349, label = "Enterobacteriaceae")+
      annotate("text", x = 0.99627502, y = -0.11623270, label = "Enterococcaceae")+
      annotate("text", x = -0.97930739, y = -0.20237844, label = "Ruminococcaceae")+
      xlim(-1.4, 2)
dev.off()

#############################################################################################################################################################
## 5. Heatmaps (genera and species levels) ##################################################################################################################
#############################################################################################################################################################

# Genera

df.case.gen <- df.gen[,which(apply(df.gen, 2, max) > 2)]
df.case.gen <- df.case.gen[!colnames(df.case.gen) %in% "C2likevirus"]

time_point_1 <- rownames(df.case.gen[rownames(df.case.gen) %in% meta.data$sample_ID[meta.data$point_number == 1],])
time_point_2 <- rownames(df.case.gen[rownames(df.case.gen) %in% meta.data$sample_ID[meta.data$point_number == 2],])
time_point_3 <- rownames(df.case.gen[rownames(df.case.gen) %in% meta.data$sample_ID[meta.data$point_number == 3],])

d1 <- data.frame(ID = time_point_1, Time_point = "first point")
d2 <- data.frame(ID = time_point_2, Time_point = "second point")
d3 <- data.frame(ID = time_point_3, Time_point = "third point") 

d <- rbind(d1, d2, d3)
rownames(d) <- as.character(d$ID); d <- d[-1]
colnames(d) <- "Time point"
df.case.gen <- df.case.gen[names(sort(colSums(df.case.gen), decreasing = T))]

svg("graphs/heatmap.genera.svg", height = 6, width = 10)
pheatmap(t(df.case.gen), 
         main="", 
         color= my_palette,
         annotation_col = d, 
         cluster_rows = F, 
         cluster_cols = T, 
         show_colnames = F,
         fontsize_row=9)
dev.off()

# Organisms
df.case.org <- df.org[,which(apply(df.org, 2, max) > 2)]
df.case.org <- df.case.org[names(sort(colSums(df.case.org), decreasing = T))]

svg("graphs/heatmap.org.svg", height = 12, width = 17)
pheatmap(t(df.case.org), 
         main="", 
         color= my_palette,
         annotation_col = d, 
         cluster_rows = F, 
         cluster_cols = T, 
         show_colnames = F,
         fontsize_row=9)
dev.off()

##############################################################################################################################################################
# 6. LefSe results processing ################################################################################################################################
##############################################################################################################################################################

# Import taxonomic LefSe results #
lefse.org <- read.table("data/lefse_output/case.org.lefse.out", sep = "\t", stringsAsFactors = F)
lefse.gen <- read.table("data/lefse_output/case.gen.lefse.out", sep = "\t", stringsAsFactors = F)

# Genus level #
colnames(lefse.gen) <- c("Feature", "LogMaxMean", "Class", "LDA_score", "p.value")
lefse.gen <- lefse.gen[lefse.gen$Class == "point_1",]
lefse.gen$p.value <- as.numeric(lefse.gen$p.value)
#lefse.gen <- lefse.gen[lefse.gen$p.value < 0.01,]
lefse.gen <- lefse.gen[order(lefse.gen$LDA_score),]
lefse.gen <- lefse.gen[order(lefse.gen$Feature),]

write.table(lefse.gen, "output/lefse.gen.table.txt", quote = F, row.names = F, sep = ",")

# Make barplot #
# svg("graphs/lefse.gen.svg", width = 5, height = 3)
# par(mar = c(4.5,10.5,0,0.5)) 
# barplot(name = lefse.gen$Feature, 
#         height = lefse.gen$LDA_score, 
#         horiz = T, 
#         las = 1,
#         col = my_palette[10],
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()

# Species level #
colnames(lefse.org) <- colnames(lefse.gen)
lefse.org <- lefse.org[lefse.org$Class == "point_1",]
lefse.org$p.value <- as.numeric(lefse.org$p.value)
#lefse.org <- lefse.org[lefse.org$p.value < 0.01,]
lefse.org <- lefse.org[order(lefse.org$LDA_score),]
lefse.gen <- lefse.gen[order(lefse.gen$Feature),]

write.table(lefse.org, "output/lefse.org.table.txt", quote = F, row.names = F, sep = ",")

# Make barplot #
# svg("graphs/lefse.org.svg", width = 7, height = 5)
# par(mar = c(4.5,16.8,0,2.5)) 
# barplot(name = lefse.org$Feature, 
#         height = lefse.org$LDA_score, 
#         horiz = T,
#         xlim=c(0,4),
#         las = 1,
#         col = my_palette[10],
#         cex.axis = 1.5, 
#         xlab = "LDA score (log 10)")
# dev.off()


# Boxplots for LefSe results 
# Genus level #

df.gen.sbs <- df.gen[,colnames(df.gen) %in% lefse.gen$Feature,]
df.gen.sbs <- merge(meta.data[c(1,2,11)], cbind(rownames(df.gen.sbs), df.gen.sbs), by = 1)
df.gen.sbs <- melt(df.gen.sbs, id.vars = c("sample_ID", "patient_id", "point_number"))
df.gen.sbs$point_number <- as.factor(df.gen.sbs$point_number)

svg("graphs/boxplots.gen.svg", width = 16, height = 4)
ggboxplot(df.gen.sbs[df.gen.sbs$variable != "Sutterella",], x = "point_number", y = "log(value)", group = "point_number", fill = "point_number")+
      #geom_jitter(aes(alpha = 0.5))+
      facet_wrap(~variable, nrow = 1, shrink = T)+
      theme_linedraw()+
      theme(strip.text.x = element_text(size = 10, colour = "white"))+
      scale_fill_manual(values = my_palette[c(1,50,130)])+
      guides(fill=FALSE, alpha=FALSE)+
      xlab("Time points")+
      ylab("log (Relativa abundance, %)")
dev.off()

# Species level #
# df.org.sbs <- df.org[,colnames(df.org) %in% lefse.org$Feature,]
# df.org.sbs <- merge(meta.data[c(1,2,11)], cbind(rownames(df.org.sbs), df.org.sbs), by = 1)
# df.org.sbs <- melt(df.org.sbs, id.vars = c("sample_ID", "patient_id", "point_number"))
# df.org.sbs$point_number <- as.factor(df.org.sbs$point_number)
# 
# svg("graphs/boxplots.org.svg", width = 22, height = 10)
# ggboxplot(df.org.sbs, x = "point_number", y = "log(value)", group = "point_number", fill = "point_number")+
#       geom_jitter(aes(alpha = 0.5))+
#       facet_wrap(~variable, ncol = 8, shrink = T)+
#       theme_linedraw()+
#       theme(strip.text.x = element_text(size = 10, colour = "white"))+
#       scale_fill_manual(values = my_palette[c(1,50,130)])+
#       guides(fill=FALSE, alpha=FALSE)+
#       xlab("Time points")+
#       ylab("log (Relativa abundance, %)")
# dev.off()

##############################################################################################################################################################
## 7. Opportunistic and pathogens diff between first and second time points ##################################################################################
##############################################################################################################################################################

# Input opportunistic-pathogenic bacteria list 

pathogenes.list <- c("Enterococcus_faecium", 
                     "Enterococcus_faecalis", 
                     "Escherichia_coli", 
                     "Streptococcus_parasanguinis", 
                     "Klebsiella_pneumoniae", 
                     "Klebsiella_oxytoca",
                     "Haemophilus_parainfluenzae",
                     "Clostridium_perfringens")

# Discovery different relative abundance of taxons by wilcoxon paired test with FDR 

dt <- NULL
for (i in pathogenes.list){
      p1 <- df.org[,i][rownames(df.org) %in% Point_1]
      p2 <- df.org[,i][rownames(df.org) %in% Point_2]
      
      w <- wilcox.test(p1, p2, 
                       paired = T, 
                       alternative = "less")
      
      w <- data.frame(feature = i, 
                      mean_1 = mean(p1), 
                      sd_1 = sd(p1), 
                      mean_2 = mean(p2), 
                      sd_2 = sd(p2),
                      p.value = w$p.value)
      dt <- rbind(dt, w)
}

dt$adj.p.value <- p.adjust(dt$p.value, method = "fdr")

write.table(dt, 'output/diff_pathogens.txt', quote=F, sep='\t', col.names = T, row.names = F)

# Prepare dataset for boxplot
names.opp <- c("Enterococcus_faecium", "Escherichia_coli")
opp.data <- merge(cbind(rownames(df.org), df.org[,colnames(df.org) %in% names.opp]), meta.sbs, by = 1)
colnames(opp.data)[1] <- "sampleID"
opp.data <- melt(opp.data[c(1,2,3,13)], id.vars = c("sampleID", "point_number"))
opp.data$point_number <- as.factor(opp.data$point_number)

# Make boxplot
svg("graphs/boxplots.ef-ec.svg")
ggboxplot(opp.data, x = "point_number", y = "log(value)", group = "point_number", fill = "point_number", size = 2)+
      #geom_jitter(aes(alpha = 0.5))+
      facet_wrap(~variable, nrow = 1, shrink = T)+
      theme_linedraw()+
      theme(strip.text.x = element_text(size = 10, colour = "white"))+
      scale_fill_manual(values = my_palette[c(1,50,130)])+
      guides(fill=FALSE, alpha=FALSE)+
      xlab("Time points")+
      ylab("log (Relativa abundance, %)")
dev.off()

##############################################################################################################################################################
# 8. MDSs for different taxonomic levels (firs and secon time points) ########################################################################################
##############################################################################################################################################################

## Family level #

df.family.2 <- df.family[rownames(df.family) %in% c(Point_1, Point_2),]
df.family.2 <- df.family.2[,which(apply(df.family.2, 2, max) > 1)]

mds.2 <- metaMDS(df.family.2, "bray")
mds.point.2 <- as.data.frame(mds.2$points)
mds.point.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(mds.point.2), mds.point.2), by = 1)
mds.point.2$point_number <- as.factor(mds.point.2$point_number)

# make non-metric MDS-plot #

svg("graphs/MDS.1.family.svg", width = 5, height = 5)
ggplot_boxplot(mds.point.2)
dev.off()

# Plot for Bray-Curtis dissimilarity values 

df.family.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(df.family.2), df.family.2), by = 1)

dt.l <- extract_bray_value(df.family.2)
dt.l <- melt(dt.l)

svg("graphs/lineplot.1.family.svg")
ggplot_lineplot(dt.l)
dev.off()

#############################################################################################################################################################

## Genera level #

df.gen.2 <- df.gen[rownames(df.gen) %in% c(Point_1, Point_2),]
df.gen.2 <- df.gen.2[,which(apply(df.gen.2, 2, max) > 1)]
mds.2 <- metaMDS(df.gen.2, "bray")
mds.point.2 <- as.data.frame(mds.2$points)

mds.point.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(mds.point.2), mds.point.2), by = 1)
mds.point.2$point_number <- as.factor(mds.point.2$point_number)

fit <- envfit(mds.2, df.gen.2, permutations = 1000)

# make non-metric MDS-plot #

svg("graphs/MDS.1.gen.svg", width = 5, height = 5)
ggplot_boxplot(mds.point.2)
dev.off()

# Plot for Bray-Curtis dissimilarity values 

df.gen.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(df.gen.2), df.gen.2), by = 1)

dt.l <- extract_bray_value(df.gen.2)
dt.l <- melt(dt.l)

svg("graphs/lineplot.1.gen.svg")
ggplot_lineplot(dt.l)
dev.off()

#############################################################################################################################################################

## Species level #

df.org.2 <- df.org[rownames(df.org) %in% c(Point_1, Point_2),]
df.org.2 <- df.org.2[,which(apply(df.org.2, 2, max) > 1)]
mds.2 <- metaMDS(df.org.2, "bray")
mds.point.2 <- as.data.frame(mds.2$points)

mds.point.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(mds.point.2), mds.point.2), by = 1)
mds.point.2$point_number <- as.factor(mds.point.2$point_number)

fit <- envfit(mds.2, df.org.2, permutations = 1000)

# make non-metric MDS-plot ##

svg("graphs/MDS.1.org.svg", width = 5, height = 5)
ggplot_boxplot(mds.point.2)
dev.off()

# Plot for Bray-Curtis dissimilarity values ##

df.org.2 <- merge(meta.sbs[c(1,2,11)], cbind(rownames(df.org.2), df.org.2), by = 1)

dt.l <- extract_bray_value(df.org.2)
dt.l <- melt(dt.l)

svg("graphs/lineplot.1.org.svg")
ggplot_lineplot(dt.l)
dev.off()

##############################################################################################################################################################
######################################################################## THE END #############################################################################
##############################################################################################################################################################
