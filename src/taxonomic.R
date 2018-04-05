#############################################################################################################################################################
################################################################ Taxonomic data processing ##################################################################
#############################################################################################################################################################

# Set work directory ##
workDir <- "/home/acari/Рабочий стол/HP_project/"
setwd(workDir)

# Import libraries and functions ##
source("src/functions.R")
source("src/functions_2.R")
source("src/functions_3.R")

set.seed(10)

# Set color pallete ##
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)
##############################################################################################################################################################
## 1. Import datasets ########################################################################################################################################
##############################################################################################################################################################
df.org.sbs <- df.org[rownames(df.org) %in% c(Point_1, Point_2),]
df.org.sbs$Escherichia_coli
# Import metadata ##
meta.data <- read.csv("data/meta_case.csv")

# Extract first and second time points sampleID ##
samples <- as.character(meta.data$Patient_ID[meta.data$Time.point.number %in% c(1,2)])
samples <- samples[which(duplicated(samples))]
meta.sbs <- meta.data[meta.data$Patient_ID %in% samples & meta.data$Time.point.number %in% c(1,2),]
names.list <- meta.sbs$Sample_ID
meta.sbs <- meta.sbs[meta.sbs$Sample_ID %in% names.list,]
meta.sbs$Patient_ID <- as.character(meta.sbs$Patient_ID)

Point_1 <- as.character(meta.sbs$Sample_ID[meta.sbs$Time.point.number == 1])
Point_2 <- as.character(meta.sbs$Sample_ID[meta.sbs$Time.point.number == 2])

meta.data.sbs <- meta.data[meta.data$Sample_ID %in% c(Point_1, Point_2),]
# write.table(meta.data.sbs, "meta.sbs.tsv" , quote = F)

##############################################################################################################################################################
# Import case data ##
metafast.output <- read.table("data/metafast/HP.metafast.txt", row.names = 1)

df.org <- read.csv("data/case.org.mpa", row.names = 1)
df.gen <- read.csv("data/case.gen.mpa", row.names = 1)
df.family <- read.table("data/case.family.mpa")
df.phyla <- read.table("data/case.phyla.mpa")
df.phyla <- df.phyla[rownames(df.phyla) %in% meta.data$sample_ID,]

# Import control data ##
control.family <- read.table("data/rus.control.family.mpa")
healthy.org <- read.table("data/healthy.org.txt")
healthy.fam <- read.table("data/healthy.fam.txt")

metafast.healthy <- read.table("data/metafast/healthy.metafast.txt", row.names = 1)

# df.org.sbs <- df.org[rownames(df.org) %in% meta.data.sbs$sample_ID,]
# ef.df <- merge(meta.data.sbs[c(1,2,11)], data.frame(sample_ID = row.names(df.org.sbs), df.org.sbs$Enterococcus_faecium), by = 1)
# ef.df[order(ef.df$patient_id),]
# unique(ef.df$patient_id)

##############################################################################################################################################################
## 2. Reference-free analysis ################################################################################################################################
##############################################################################################################################################################
colnames(metafast.output) <- sapply(str_split(rownames(metafast.output), "_"), function(x) x[length(x) - 1])
rownames(metafast.output) <- colnames(metafast.output)
metafast.sbs <- metafast.output[rownames(metafast.output) %in% c(Point_1, Point_2),]
metafast.sbs <- metafast.sbs[,colnames(metafast.sbs) %in% c(Point_1, Point_2)]

metafast.mds <- metaMDS(metafast.sbs)
metafast.mds.p <- as.data.frame(metafast.mds$points)
metafast.mds.p <- cbind(rownames(metafast.mds.p), metafast.mds.p)
colnames(metafast.mds.p)[1] <- "sampleID"
metafast.mds.p <- merge(metafast.mds.p, meta.data.sbs[c(1,2,11)], by = 1)
metafast.mds.p$point_number <- as.factor(metafast.mds.p$point_number)

svg("graphs/metafast.mds.svg", height = 5, width = 5)
ggplot(metafast.mds.p, aes(MDS1, MDS2, color = point_number, shape = point_number))+
      theme_linedraw()+
      geom_point(size = 4)+
      stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = point_number))+
      geom_line(metafast.mds.p, mapping = aes(MDS1, MDS2, group = patient_id), col = "black", alpha = 0.5)+
      theme(legend.position = "bottom")+
      scale_color_manual(values = my_palette[c(1,30,170)])+
      scale_shape_manual(values = c(8,2,16,0))+
      guides(color=guide_legend(title="Time point"), 
             shape=guide_legend(title="Time point"),
             fill = F)
dev.off()

anno <- rbind(data.frame(sampleID = Point_1, Time_point = 1),
      data.frame(sampleID = Point_2, Time_point = 2))
rownames(anno) <- anno$sampleID
anno <- anno[-1]
colnames(anno) <- "Time point"
anno$`Time point` <- as.factor(anno$`Time point`)

svg("graphs/metafast.heatmap.svg", height = 5, width = 6.3)
pheatmap(metafast.sbs, 
         annotation_col = anno, 
         show_rownames = F, 
         show_colnames = F)
dev.off()

metafast.sbs.2 <- cbind(rownames(metafast.sbs), metafast.sbs)
colnames(metafast.sbs.2)[1] <- "sampleID"
metafast.sbs.2 <- merge(meta.data.sbs[c(1,2)], metafast.sbs.2, by = 1)
metafast.sbs.2 <- metafast.sbs.2[order(metafast.sbs.2$patient_id),]

patID <- unique(as.character(metafast.sbs.2$patient_id))

df.metafast <- NULL
for (k in patID){
      sbs <- metafast.sbs.2[metafast.sbs.2$patient_id == k,]
      sbs <- sbs[colnames(sbs) %in% as.character(sbs$sample_ID)]
      sbs <- data.frame(patientID = k, metafast_distance = max(sbs))
      df.metafast <- rbind(sbs, df.metafast)
}

colnames(metafast.healthy) <- rownames(metafast.healthy)
sbs.index <- which(!is.na(str_extract(rownames(metafast.healthy), "11-7-0|11-60-0|11-0-0")))
metafast.healthy <- metafast.healthy[sbs.index]
metafast.healthy <- metafast.healthy[sbs.index,]
metafast.healthy$patientID <- sapply(str_split(rownames(metafast.healthy), "-"), function(x) head(x,1))

patID_healthy <- unique(metafast.healthy$patientID)

healthy.dis <- NULL
for (i in patID_healthy){
      metafast.healthy.sbs <- metafast.healthy[metafast.healthy$patientID == i,]
      metafast.healthy.sbs <- metafast.healthy.sbs[colnames(metafast.healthy.sbs) %in% rownames(metafast.healthy.sbs)]
      metafast.healthy.sbs <- metafast.healthy.sbs[1]
      metafast.healthy.sbs$days <- sapply(str_split(rownames(metafast.healthy.sbs), "-"), function(x) x[length(x)-1])
      metafast.healthy.sbs <- metafast.healthy.sbs[metafast.healthy.sbs$days %in% c(7,60),]
      metafast.healthy.sbs$days <- as.numeric(metafast.healthy.sbs$days)
      metafast.healthy.sbs <- metafast.healthy.sbs[order(metafast.healthy.sbs$days),]
      colnames(metafast.healthy.sbs)[1] <- "metafast_distance"
      metafast.healthy.sbs <- cbind(rownames(metafast.healthy.sbs), metafast.healthy.sbs)
      colnames(metafast.healthy.sbs)[1] <- "sampleID"
      rownames(metafast.healthy.sbs) <- 1:nrow(metafast.healthy.sbs)
      healthy.dis <- rbind(metafast.healthy.sbs, healthy.dis)
} 

healthy.dis$days <- as.factor(healthy.dis$days)
healthy.dis$sampleID <- sapply(str_split(healthy.dis$sampleID, "-"), function(x) head(x,1))
colnames(healthy.dis)[1] <- "patientID"

df.metafast$days <- "14"

metafast.all <- rbind(healthy.dis, df.metafast)
metafast.all$days <- as.numeric(as.character(metafast.all$days))
metafast.all <- metafast.all[order(metafast.all$days),]
metafast.all$days <- as.factor(metafast.all$days)

svg("graphs/metafast.changes.svg", width = 8, height = 4)
ggplot(metafast.all, aes(1:nrow(metafast.all), metafast_distance, col = days))+
      geom_point(shape = 1, size = 4.5)+
      geom_point(size = 2)+
      geom_hline(yintercept = mean(metafast.all$metafast_distance), col = "red", linetype = 1, size = 1.2)+
      geom_hline(yintercept = mean(metafast.all$metafast_distance)+sd(metafast.all$metafast_distance), col = "red", linetype = 2, size = 0.8)+
      geom_hline(yintercept = mean(metafast.all$metafast_distance)-sd(metafast.all$metafast_distance), col = "red", linetype = 2, size = 0.8)+
      geom_vline(xintercept = 13.5, linetype = 3, size = 0.8)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      xlab("Patient's #")+
      ylab("Metafast distance")
dev.off()

##############################################################################################################################################################
## 2. Top taxa calculation ###################################################################################################################################
##############################################################################################################################################################

# Genus level #

df <- NULL
for (i in 1:2){
      ff <- data.matrix(df.gen.sbs[rownames(df.gen.sbs) %in% meta.data.sbs$sample_ID[meta.data.sbs$point_number == i],])
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
for (i in 1:2){
      ff <- data.matrix(df.org.sbs[rownames(df.org.sbs) %in% meta.data.sbs$sample_ID[meta.data.sbs$point_number == i],])
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

hclust(dist(df.gen.sbs), d = 3)
#############################################################################################################################################################
## 3. Heatmaps (genera and species levels) ##################################################################################################################
#############################################################################################################################################################
# Genera

df.case.gen <- df.gen[,which(apply(df.gen, 2, max) > 2)]
df.case.gen <- df.case.gen[!colnames(df.case.gen) %in% "C2likevirus"]
df.case.gen <- df.case.gen[rownames(df.case.gen) %in% c(Point_1, Point_2),]

time_point_1 <- rownames(df.case.gen[rownames(df.case.gen) %in% meta.data$sample_ID[meta.data$point_number == 1],])
time_point_2 <- rownames(df.case.gen[rownames(df.case.gen) %in% meta.data$sample_ID[meta.data$point_number == 2],])

d1 <- data.frame(ID = time_point_1, Time_point = 1)
d2 <- data.frame(ID = time_point_2, Time_point = 2)

d <- rbind(d1, d2)
rownames(d) <- as.character(d$ID); d <- d[-1]
colnames(d) <- "Time point"
d$`Time point` <- as.factor(d$`Time point`)
df.case.gen <- df.case.gen[names(sort(colSums(df.case.gen), decreasing = T))]

svg("graphs/heatmap.genera.svg", height = 6.5, width = 10.5)
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
df.case.org <- df.org[,which(apply(df.org, 2, max) > 12)]
df.case.org <- df.case.org[names(sort(colSums(df.case.org), decreasing = T))]
df.case.org <- df.case.org[rownames(df.case.org) %in% c(Point_1, Point_2),]

svg("graphs/heatmap.org.svg", height = 6.5, width = 10.5)
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

df.case.HP_003 <- df.case.org[rownames(df.case.org) %in% c("5HP", "6HP"),]  
df.case.HP_003 <- melt(df.case.HP_003[c("Enterococcus_faecium", "Enterococcus_faecalis", "Bacteroides_dorei")])
df.case.HP_003$time_point <- c(1,2,1,2,1,2)
df.case.HP_003$time_point <- as.factor(df.case.HP_003$time_point)

svg("graphs/bacs.svg", height = 2.5, width = 5)
ggplot(df.case.HP_003, aes(time_point, value, group = variable))+
      facet_wrap(~variable)+
      geom_point(size = 2)+
      geom_line()+
      theme_linedraw()+
      xlab("Time point")+
      ylab("Relative abundance, %")
dev.off()

##############################################################################################################################################################
## 4. Alpha-diversity  #######################################################################################################################################
##############################################################################################################################################################
# Alpha-diversity calculate ##
spec_num <- specnumber(df.org)
shannon.i <- diversity(df.org, index = "shannon")
simpson.i <- diversity(df.org, index = "simpson")

# Make data.frame
alpha.div <- data.frame(spec_num = spec_num, 
                        shannon = shannon.i, 
                        simpson = simpson.i)

alpha.div <- merge(cbind(rownames(alpha.div), alpha.div), meta.data[c(1,10)], by = 1)
alpha.div$Time.point.number <- as.factor(alpha.div$Time.point.number)
colnames(alpha.div)[1] <- "sampleID"

# Wilcoxon signed-rank paired test ##
index.names <- colnames(alpha.div[-c(1,5)])
alpha.div.sbs <- alpha.div[as.character(alpha.div$sampleID) %in% c(Point_1, Point_2),]
colnames(alpha.div.sbs)[5] <- "point_number"

alpha.p <- NULL
for (k in index.names){
      df.l <- alpha.div.sbs[c("point_number", k)]
      
      first.point <- df.l[df.l$point_number == 1,]
      second.point <- df.l[df.l$point_number == 2,]
      
      w.alpha.div <- wilcox.test(first.point[,2],
                                 second.point[,2], paired = T,
                                 alternative = "greater")
      StDev.f <- sd(first.point[,2])
      StDev.s <- sd(second.point[,2])
      alpha.p <- rbind(alpha.p, data.frame(index = k, 
                                           mean.1 = mean(first.point[,2]),
                                           mean.2 = mean(second.point[,2]),
                                           SD.1 = sd(first.point[,2]),
                                           SD.2 = sd(second.point[,2]),
                                           p.val = w.alpha.div$p.value))
}

write.table(alpha.p, 'output/alpha.diff.txt', quote=F, sep='\t', row.names = F)

colnames(alpha.div)[-c(1,5)] <- c("Observed_species", "Shannon_index", "Simpson_index")
index.names2 <- colnames(alpha.div)[-c(1,5)]

# Make boxplots ##
alpha.div.t <- melt(alpha.div)
colnames(alpha.div.t)[3] <- c("index")
alpha.div.t <- alpha.div.t[alpha.div.t$Time.point.number != 3,]
alpha.div.t2 <- alpha.div.t[alpha.div.t$sampleID %in% c(Point_1, Point_2),]

for (i in index.names2){
      svg(paste0("graphs/alpha-", i,".svg"), width = 3, height = 5)
      figure <- ggboxplot(alpha.div.t2[alpha.div.t2$index == i,], x = "Time.point.number", y = "value",
            color = "Time.point.number", add = "jitter",
            palette = my_palette[c(1,50,130)])+
            stat_compare_means(paired = TRUE, label.x.npc = 0.5, label.y.npc = 1.0)+
            theme_linedraw()+
            theme(legend.position="none")+
            xlab("Time point")+
            ggtitle(i)
      print(figure)
      dev.off()      
}
      
# Prevelance plots ## 
# OTUs <- colnames(df.org)
# OTUs[which(!is.na(str_extract(OTUs, "_unclassified")))] <- sapply(str_split(OTUs[which(!is.na(str_extract(OTUs, "_unclassified")))], "_"), 
#                                                                   function(x) head(x, 1))
# taxa.output <- tax_name(query = c(OTUs), get = c("phylum"), db = "ncbi")
# df.org.phyla <- df.org
# colnames(df.org.phyla) <- taxa.output$phylum
# df.org.phyla <- df.org.phyla[!is.na(colnames(df.org.phyla))]
# df.org.phyla <- as.data.frame(t(df.org.phyla))
# df.org.phyla <- cbind(rownames(df.org.phyla), df.org.phyla)
# colnames(df.org.phyla)[1] <- "phyla"
# rownames(df.org.phyla) <- 1:130
# df.org.phyla$phyla <- sapply(str_split(df.org.phyla$phyla, "\\."), function(x) head(x,1))
# df.org.phyla$phyla <- as.factor(df.org.phyla$phyla)
# df.org.phyla.m <- melt(df.org.phyla)
# df.org.phyla.m <- merge(df.org.phyla.m[c(2,1,3)], meta.data[c(1,11)], by = 1)
# df.org.phyla.m <- df.org.phyla.m[df.org.phyla.m$value > 0,]
# 

#############################################################################################################################################################
## 5. Changes in the taxonomy distance and alpha-diversity (first and second time points) ###################################################################
#############################################################################################################################################################
alpha.div.sbs <- merge(meta.data[c(1,2)], alpha.div.sbs, by = 1)
alpha.div.sbs <- alpha.div.sbs[order(alpha.div.sbs$Patient_ID),]

spec_num.changes <- abs(alpha.div.sbs$spec_num[alpha.div.sbs$point_number == 1] - alpha.div.sbs$spec_num[alpha.div.sbs$point_number == 2])
shannon.changes <- abs(alpha.div.sbs$shannon[alpha.div.sbs$point_number == 1] - alpha.div.sbs$shannon[alpha.div.sbs$point_number == 2])
simpson.changes <- abs(alpha.div.sbs$simpson[alpha.div.sbs$point_number == 1] - alpha.div.sbs$simpson[alpha.div.sbs$point_number == 2])

super_plot <- function(x){
      plot(x)
      abline(h = mean(x, color = "red", lty = 5))
      abline(h = mean(x) - sd(x))
      abline(h = mean(x) + sd(x))
}

df.fam.sbs <- df.org[rownames(df.org) %in% c(Point_1, Point_2),]
df.fam.sbs <- cbind(rownames(df.fam.sbs), df.fam.sbs)

colnames(df.fam.sbs)[1] <- "sampleID"
df.fam.sbs <- merge(meta.data[c(1,2)], df.fam.sbs, by = 1)
df.fam.sbs <- df.fam.sbs[order(df.fam.sbs$Patient_ID),]
patID <- as.character(unique(df.fam.sbs$Patient_ID))

bray.dist <- NULL
for (i in patID){
      di <- metaMDSdist(df.fam.sbs[df.fam.sbs$Patient_ID == i,][-c(1:2)], distance = "bray")
      bray.dist <- rbind(data.frame(patientID = i, bray_dist = as.vector(di)), bray.dist)
}

jaccard.dist <- NULL
for (i in patID){
      di <- metaMDSdist(df.fam.sbs[df.fam.sbs$Patient_ID == i,][-c(1:2)], distance = "jaccard")
      jaccard.dist <- rbind(data.frame(patientID = i, jaccard_dist = as.vector(di)), jaccard.dist)
}

euclid.dist <- NULL
for (i in patID){
      di <- metaMDSdist(df.fam.sbs[df.fam.sbs$Patient_ID == i,][-c(1:2)], distance = "euclidean")
      euclid.dist <- rbind(data.frame(patientID = i, euclid_dist = as.vector(di)), euclid.dist)
}

df.dist <- join_all(list(bray.dist, jaccard.dist, euclid.dist), by = 1, type = 'full')
alpha.changes <- data.frame(patientID = unique(alpha.div.sbs$Patient_ID), spec_num.changes, shannon.changes, simpson.changes)
df.changes <- merge(alpha.changes, df.dist, by = 1)

write.table(df.changes, 'output/ch.2p.txt', quote=F, sep='\t', row.names = F)

# Healthy 
healthy.org <- healthy.org[,which(apply(healthy.org, 2, max) > 1)]

healthy.patientID <- sapply(str_split(rownames(healthy.org), "-"), function(x) head(x, 1))
healthy.time_points <- as.numeric(sapply(str_split(rownames(healthy.org), "-"), function(x) x[3]))
meta.healthy <- data.frame(patientID = healthy.patientID, time_points = healthy.time_points)

healthy.spec_num <- specnumber(healthy.org)
healthy.shannon <- diversity(healthy.org, index = "shannon")
healthy.simpson <- diversity(healthy.org, index = "simpson")

alpha.div.healthy <- data.frame(spec_num = healthy.spec_num, 
                                shannon = healthy.shannon, 
                                simpson = healthy.simpson)
alpha.div.healthy <- cbind(rownames(alpha.div.healthy), alpha.div.healthy)
colnames(alpha.div.healthy)[1] <- "sampleID"
rownames(alpha.div.healthy) <- 1:20
alpha.div.healthy <- cbind(alpha.div.healthy, meta.healthy)[c(1,5,2,3,4,6)]

# 7 days
healthy.spec_num.c <- abs(alpha.div.healthy$spec_num[alpha.div.healthy$time_points == 0] - alpha.div.healthy$spec_num[alpha.div.healthy$time_points == 7])
healthy.shannon.c <- abs(alpha.div.healthy$shannon[alpha.div.healthy$time_points == 0] - alpha.div.healthy$shannon[alpha.div.healthy$time_points == 7])
healthy.simpson.c <- abs(alpha.div.healthy$simpson[alpha.div.healthy$time_points == 0] - alpha.div.healthy$simpson[alpha.div.healthy$time_points == 7])

# 60 days 
alpha.div.healthy60 <- alpha.div.healthy[alpha.div.healthy$patientID != "halbarad",]

healthy.spec_num.c2 <- abs(alpha.div.healthy60$spec_num[alpha.div.healthy60$time_points == 0] - alpha.div.healthy60$spec_num[alpha.div.healthy60$time_points == 60])
healthy.shannon.c2 <- abs(alpha.div.healthy60$shannon[alpha.div.healthy60$time_points == 0] - alpha.div.healthy60$shannon[alpha.div.healthy60$time_points == 60])
healthy.simpson.c2 <- abs(alpha.div.healthy60$simpson[alpha.div.healthy60$time_points == 0] - alpha.div.healthy60$simpson[alpha.div.healthy60$time_points == 60])

## Beta-diversity
heal.org.meta <- cbind(meta.healthy, healthy.org)

# 7 days
heal.org.meta.7 <- heal.org.meta[heal.org.meta$time_points %in% c(0,7),]
patID.7 <- as.character(unique(heal.org.meta.7$patientID))

bray.dist.7 <- NULL
for (i in patID.7){
      di <- metaMDSdist(heal.org.meta.7[heal.org.meta.7$patientID == i,][-c(1:2)], distance = "bray")
      bray.dist.7 <- rbind(data.frame(patientID = i, bray_dist = as.vector(di)), bray.dist.7)
}

jaccard.dist.7 <- NULL
for (i in patID.7){
      di <- metaMDSdist(heal.org.meta.7[heal.org.meta.7$patientID == i,][-c(1:2)], distance = "jaccard")
      jaccard.dist.7 <- rbind(data.frame(patientID = i, jaccard_dist = as.vector(di)), jaccard.dist.7)
}

euclid.dist.7 <- NULL
for (i in patID.7){
      di <- metaMDSdist(heal.org.meta.7[heal.org.meta.7$patientID == i,][-c(1:2)], distance = "euclidean")
      euclid.dist.7 <- rbind(data.frame(patientID = i, euclid_dist  = as.vector(di)), euclid.dist.7)
}

healthy.dist.7 <- join_all(list(bray.dist.7, jaccard.dist.7, euclid.dist.7), by = 1, type = 'full')
alpha.changes.7 <- data.frame(patientID = unique(meta.healthy$patientID), 
                              healthy.spec_num.c, healthy.shannon.c, healthy.simpson.c)

healthy.changes.7 <- merge(alpha.changes.7, healthy.dist.7, by = 1)
colnames(healthy.changes.7) <- colnames(df.changes) 

# 60 days 
heal.org.meta.60 <- heal.org.meta[heal.org.meta$time_points %in% c(0,60),]
heal.org.meta.60 <- heal.org.meta.60[heal.org.meta.60$patientID != "halbarad",] 

patID.60 <- as.character(unique(heal.org.meta.60$patientID))

bray.dist.60 <- NULL
for (i in patID.60){
      di <- metaMDSdist(heal.org.meta.60[heal.org.meta.60$patientID == i,][-c(1:2)], distance = "bray")
      bray.dist.60 <- rbind(data.frame(patientID = i, bray_dist = as.vector(di)), bray.dist.60)
}

jaccard.dist.60 <- NULL
for (i in patID.60){
      di <- metaMDSdist(heal.org.meta.60[heal.org.meta.60$patientID == i,][-c(1:2)], distance = "jaccard")
      jaccard.dist.60 <- rbind(data.frame(patientID = i, jaccard_dist = as.vector(di)), jaccard.dist.60)
}

euclid.dist.60 <- NULL
for (i in patID.60){
      di <- metaMDSdist(heal.org.meta.60[heal.org.meta.60$patientID == i,][-c(1:2)], distance = "euclidean")
      euclid.dist.60 <- rbind(data.frame(patientID = i, euclid_dist = as.vector(di)), euclid.dist.60)
}

healthy.dist.60 <- join_all(list(bray.dist.60, jaccard.dist.60, euclid.dist.60), by = 1, type = 'full')
alpha.changes.60 <- data.frame(patientID = unique(meta.healthy$patientID)[-4], 
                               healthy.spec_num.c2, healthy.shannon.c2, healthy.simpson.c2)

healthy.changes.60 <- merge(alpha.changes.60, healthy.dist.60, by = 1)
colnames(healthy.changes.60) <- colnames(df.changes) 

write.table(healthy.changes.7, 'output/ch.healthy.7.txt', quote=F, sep='\t', row.names = F)
write.table(healthy.changes.60, 'output/ch.healthy.60.txt', quote=F, sep='\t', row.names = F)

# Merge tables ## 
df.changes$disease <- "HP"
healthy.changes.7$disease <- "healthy" 
healthy.changes.60$disease <- "healthy" 

df.changes$time_point <- 14
healthy.changes.7$time_point <- 7 
healthy.changes.60$time_point <- 60 

all.changes <- rbind(df.changes, healthy.changes.7, healthy.changes.60)
wilcox.test(all.changes$bray_dist[all.changes$disease == "HP"], all.changes$bray_dist[all.changes$disease == "healthy"], alternative = "greater")
boxplot(all.changes$bray_dist[all.changes$disease == "HP"], all.changes$bray_dist[all.changes$disease == "healthy"])

# Scatterplots making ##
svg("graphs/bray-shannon.svg", width = 6, height = 6.5)
ggplot(all.changes, aes(bray_dist, shannon.changes, shape = as.factor(time_point), col = disease))+
      geom_point(size = 2.5, stroke = 1.5)+
      annotate("rect", xmin = -Inf, xmax = 0.73, ymin = 0.79, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = 0.73, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.73, ymin = 0.4, ymax = 0.79,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = 0.56, xmax = 0.73, ymin = -Inf, ymax = 0.4,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.56, ymin = -Inf, ymax = 0.4,   fill = "darkgreen", alpha = 0.25)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Bray-Curtis distance")+
      ylab("Shannon index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,20))+
      theme(legend.position = "bottom")
dev.off()

svg("graphs/bray-simpson.svg", width = 6, height = 6.5)
ggplot(all.changes, aes(bray_dist, simpson.changes, shape = as.factor(time_point), col = disease))+
      geom_point(size = 2.5, stroke = 1.5)+
      annotate("rect", xmin = -Inf, xmax = 0.73, ymin = 0.21, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = 0.73, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.73, ymin = 0.09, ymax = 0.21,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = 0.56, xmax = 0.73, ymin = -Inf, ymax = 0.09,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.56, ymin = -Inf, ymax = 0.09,   fill = "darkgreen", alpha = 0.25)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Bray-Curtis distance")+
      ylab("Simpson index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,20))+
      theme(legend.position = "bottom")
dev.off()

svg("graphs/jaccard-shannon.svg", width = 6, height = 6.5)
ggplot(all.changes, aes(jaccard_dist, shannon.changes, shape = as.factor(time_point), col = disease))+
      geom_point(size = 2.5, stroke = 1.5)+
      annotate("rect", xmin = -Inf, xmax = 0.85, ymin = 0.79, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = 0.85, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "darkred", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.85, ymin = 0.4, ymax = 0.79,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = 0.7, xmax = 0.85, ymin = -Inf, ymax = 0.4,   fill = "darkblue", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 0.7, ymin = -Inf, ymax = 0.4,   fill = "darkgreen", alpha = 0.25)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Bray-Curtis distance")+
      ylab("Shannon index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,20))+
      theme(legend.position = "bottom")
dev.off()

svg("graphs/euclid-shannon.svg", width = 6, height = 5)
ggplot(all.changes, aes(euclid_dist, shannon.changes, shape = as.factor(time_point), col = disease))+
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "red", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.19, ymin = -Inf, ymax = 0.79,   fill = "blue", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.14, ymin = -Inf, ymax = 0.4,   fill = "green", alpha = 0.2)+
      geom_point(size = 2)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Euclid distance")+
      ylab("Shannon index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,19))
dev.off()

svg("graphs/jaccard-simpson.svg", width = 6, height = 5)
ggplot(all.changes, aes(jaccard_dist, simpson.changes, shape = as.factor(time_point), col = disease))+
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "red", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.85, ymin = -Inf, ymax = 0.21,   fill = "blue", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.70, ymin = -Inf, ymax = 0.09,   fill = "green", alpha = 0.2)+
      geom_point(size = 2)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Jaccard distance")+
      ylab("Simpson index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,19))
dev.off()

svg("graphs/euclid-simpson.svg", width = 6, height = 5)
ggplot(all.changes, aes(euclid_dist, simpson.changes, shape = as.factor(time_point), col = disease))+
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "red", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.20, ymin = -Inf, ymax = 0.21,   fill = "blue", alpha = 0.2)+
      annotate("rect", xmin = -Inf, xmax = 0.15, ymin = -Inf, ymax = 0.09,   fill = "green", alpha = 0.2)+
      geom_point(size = 2)+
      # geom_text(label = df.changes$patientID)+
      theme_linedraw()+
      xlab("Euclid distance")+
      ylab("Simpson index change")+
      scale_color_manual(name = "Status", values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(name = "Days", values = c(0,17,19))
dev.off()

wilcox.test(df.changes$bray_dist[df.changes$patientID %in% meta.data.sbs$Patient_ID[meta.data.sbs$Duodenal.ulcer == 1]],
      df.changes$bray_dist[df.changes$patientID %in% meta.data.sbs$Patient_ID[meta.data.sbs$Duodenal.ulcer == 0]], alternative = "greater")

wilcox.test(df.changes$bray_dist[df.changes$patientID %in% meta.data.sbs$Patient_ID[meta.data.sbs$Gender == "F"]],
            df.changes$bray_dist[df.changes$patientID %in% meta.data.sbs$Patient_ID[meta.data.sbs$Gender == "M"]], alternative = "greater")

#############################################################################################################################################################
# 6. ANOSIM, ADONIS, Maaslin ################################################################################################################################
#############################################################################################################################################################
# ANOSIM ##
df.gen.meta <- merge(cbind(rownames(df.gen), df.gen), meta.data.sbs[c(1,11)], by = 1)
rownames(df.gen.meta) <- df.gen.meta[,1]
df.gen.meta <- df.gen.meta[-1]

fit.anosim <- anosim(df.gen.meta, grouping = df.gen.meta$point_number, permutations = 10000, distance = "bray")
summary(fit.anosim)
plot(fit.anosim)

# ADONIS ##
df.gen.sbs <- df.gen[rownames(df.gen) %in% c(Point_1, Point_2),]
meta.data.sbs <- meta.data[meta.data$Sample_ID %in% c(Point_1, Point_2),]

df.gen.sbs <- cbind(rownames(df.gen.sbs), df.gen.sbs)
colnames(df.gen.sbs)[1] <- "sampleID"
df.gen.sbs <- df.gen.sbs[order(df.gen.sbs$sampleID),][-1]

df.org.sbs <- df.org.sbs[rownames(df.org.sbs) %in% c(Point_1, Point_2),]
df.org.sbs <- df.org.sbs[order(rownames(df.org.sbs)),]

adonis.fit1 <- adonis(df.gen.sbs ~ ., data=meta.data.sbs[c(3:10)], permutations=9999)
adonis.fit2 <- adonis(df.org.sbs ~ ., data=meta.data.sbs[c(3:10)], permutations=9999)

table(meta.data.sbs$Gender[meta.data.sbs$GERD == 1])

df.DU <- rbind(data.frame(changes = df.changes$bray_dist[df.changes$patientID %in% as.character(meta.data.sbs$Patient_ID[meta.data.sbs$Duodenal.ulcer == 1])], group = "DU"),
      data.frame(changes = df.changes$bray_dist[df.changes$patientID %in% as.character(meta.data.sbs$Patient_ID[meta.data.sbs$Duodenal.ulcer == 0])], group = "Non-DU"))

df.Gender <- rbind(data.frame(changes = df.changes$bray_dist[df.changes$patientID %in% as.character(meta.data.sbs$Patient_ID[meta.data.sbs$Gender == "F"])], group = "Female"),
               data.frame(changes = df.changes$bray_dist[df.changes$patientID %in% as.character(meta.data.sbs$Patient_ID[meta.data.sbs$Gender == "M"])], group = "Male"))


svg("graphs/DU.boxplot.svg", width = 3, height = 5)
ggboxplot(df.DU, x = "group", y = "changes",
      color = "group", add = "jitter",
      palette = my_palette[c(1,50,130)])+
      stat_compare_means(label.x.npc = 0.5, label.y.npc = 1.0)+
      theme_linedraw()+
      theme(legend.position="none")+
      ylab("Bray-Сurtis distance")+
      xlab("Status")+
      ggtitle("Duodenal ulcer")
dev.off()

svg("graphs/Gender.boxplot.svg", width = 3, height = 5)
ggboxplot(df.Gender, x = "group", y = "changes",
      color = "group", add = "jitter",
      palette = my_palette[c(1,50,130)])+
      stat_compare_means(label.x.npc = 0.5, label.y.npc = 1.0)+
      theme_linedraw()+
      theme(legend.position="none")+
      ylab("Bray-Сurtis distance")+
      xlab("Gender")+
      ggtitle("Gender")
dev.off()

# write.table(meta.data.sbs, "output/metaSBS.csv", quote = F, row.names = F)    

# Maaslin
# Genus 
df.gen.Maaslin <- cbind(rownames(df.gen.sbs), df.gen.sbs)
df.gen.Maaslin <- merge(meta.data.sbs[c(1:11)], df.gen.Maaslin, by = 1)
df.gen.Maaslin <- df.gen.Maaslin[-2]

write.table(df.gen.Maaslin, "data/Maaslin/genus.maaslin.txt", quote = F, row.names = F, sep = "\t")

# Org
df.org.Maaslin <- cbind(rownames(df.org.sbs), df.org.sbs)
df.org.Maaslin <- merge(meta.data.sbs[c(1:11)], df.org.Maaslin, by = 1)
df.org.Maaslin <- df.org.Maaslin[-2]

write.table(df.org.Maaslin, "data/Maaslin/org.maaslin.txt", quote = F, row.names = F, sep = "\t")

#############################################################################################################################################################
df.shannon <- as.data.frame(shannon.i)
df.shannon <- cbind(rownames(df.shannon), df.shannon)
df.shannon <- merge(meta.data.sbs[c(1,2,10)], df.shannon, by = 1)
df.shannon.sbs <- df.shannon[df.shannon$Time.point.number == 1,]
df.shannon.sbs <- merge(df.shannon.sbs[-1], df.changes[c(1,5)], by = 1)[-2]
df.shannon.sbs$healthy_dist <- df.healthy.dist$healthy_dist

svg("graphs/Shannon_1st-Bray.svg", width = 5.5, height = 5.5)
cor.data <- cor.test(df.shannon.sbs$shannon.i, df.shannon.sbs$bray_dist, method = "spearman")
ggplot(df.shannon.sbs, aes(shannon.i, bray_dist))+
      geom_point(shape = 1, size = 3)+
      geom_point(size = 1)+
      stat_smooth(method = "lm", se = F, col = "red")+
      xlab("Shannon index 1st point")+
      ylab("Bray-Curtis distance between 1st and 2nd point of patients")+
      annotate(x=1.30, y=0.40, 
               label=paste("R = ", round(cor.data$estimate, 2)), 
               geom="text", size=5)+
      annotate(x=1.35, y=0.37, 
               label=paste("p-value = ", round(cor.data$p.value, 2)), 
               geom="text", size=5)+
      theme_linedraw()
dev.off()

rownames(df.fam.sbs) <- df.fam.sbs$sample_ID
df.fam.sbs <- df.fam.sbs[-c(1:2)]

df.healthy.dist <- NULL
for (i in 1:nrow(df.fam.sbs)){
      df.fam.all <- as.data.frame(UniteMatrices(as.matrix(df.fam.sbs[i,]), as.matrix(healthy.org)))
      df.fam.all <- metaMDSdist(df.fam.all, "bray")
      df.fam.all <- data.frame(sampleID = rownames(df.fam.sbs[i,]), dist = mean(df.fam.all[1:20]))
      df.healthy.dist <- rbind(df.fam.all, df.healthy.dist)
}

df.healthy.dist <- merge(meta.data.sbs[c(1,2,10)], df.healthy.dist, by = 1)
df.healthy.dist <- df.healthy.dist[df.healthy.dist$Time.point.number == 2,]
df.healthy.dist <- merge(df.changes[c(1,5)], df.healthy.dist[c(2,4)], by = 1)
colnames(df.healthy.dist) <- c("patientID", "bray_dist", "healthy_dist")

svg("graphs/Bray-Healthy_1st.svg", width = 5.5, height = 5.5)
cor.data <- cor.test(df.healthy.dist$healthy_dist, df.healthy.dist$bray_dist, method = "spearman")
ggplot(df.healthy.dist, aes(healthy_dist, bray_dist))+
      geom_point(shape = 1, size = 3)+
      geom_point(size = 1)+
      stat_smooth(method = "lm", se = F, col = "red")+
      xlab("Mean Bray-Curtis distance between healthy cohort and patient's 1st point")+
      ylab("Bray-Curtis distance between 1st and 2nd point of patients")+
      annotate(x=0.86, y=0.4, 
               label=paste("R = ", round(cor.data$estimate, 2)), 
               geom="text", size=5)+
      annotate(x=0.86, y=0.37, 
               label=paste("p-value = ", round(cor.data$p.value, 2)), 
               geom="text", size=5)+
      theme_linedraw()
dev.off()

svg("graphs/Bray-Healthy_2nd.svg", width = 5.5, height = 5.5)
cor.data <- cor.test(df.healthy.dist$healthy_dist, df.healthy.dist$bray_dist, method = "spearman")
ggplot(df.healthy.dist, aes(healthy_dist, bray_dist))+
      geom_point(shape = 1, size = 3)+
      geom_point(size = 1)+
      stat_smooth(method = "lm", se = F, col = "red")+
      xlab("Mean Bray-Curtis distance between healthy cohort and patient's 2nd point")+
      ylab("Bray-Curtis distance between 1st and 2nd point of patients")+
      annotate(x=0.92, y=0.4, 
               label=paste("R = ", round(cor.data$estimate, 2)), 
               geom="text", size=5)+
      annotate(x=0.92, y=0.37, 
               label=paste("p-value = ", round(cor.data$p.value, 6)), 
               geom="text", size=5)+
      theme_linedraw()
dev.off()

library("KEGG.db")


#############################################################################################################################################################
## 6. MDS plots #############################################################################################################################################
#############################################################################################################################################################

# Merge case and control tables
df.all <- UniteMatrices(as.matrix(df.family[rownames(df.family) %in% c(Point_1, Point_2),]), as.matrix(healthy.fam))
df.all <- as.data.frame(df.all)
df.all <- df.all[!colnames(df.all) %in% "Viruses_noname"]
df.all <- df.all[,which(apply(df.all, 2, max) > 1)]

# building MDS
df.all.sbs <- df.all[,which(apply(df.all , 2, max) > 5)]
mds <- metaMDS(df.all.sbs, "bray")
fit <- envfit(mds, df.all.sbs, permutations = 10000)

mds.points <- as.data.frame(mds$points)

mds.points.hp <- mds.points[rownames(mds.points) %in% c(Point_1, Point_2),]
mds.points.hp <- cbind(rownames(mds.points.hp), mds.points.hp)
colnames(mds.points.hp)[1] <- "sampleID"
mds.points.hp <- merge(mds.points.hp, meta.data.sbs[c(1,2,11)], by = 1)
mds.points.hp$point_number[mds.points.hp$point_number == 1] <- 0
mds.points.hp$point_number[mds.points.hp$point_number == 2] <- 14
colnames(mds.points.hp)[c(4,5)] <- c("patientID", "days")
mds.points.hp$group <- "HP"

mds.points.healthy <- mds.points[!rownames(mds.points) %in% c(Point_1, Point_2),]
mds.points.healthy <- cbind(mds.points.healthy, meta.healthy)
mds.points.healthy <- cbind(rownames(mds.points.healthy), mds.points.healthy)
colnames(mds.points.healthy)[c(1,5)] <- c("sampleID", "days")
mds.points.healthy$group <- "healthy"

mds.points <- rbind(mds.points.healthy, mds.points.hp)
rownames(mds.points) <- 1:nrow(mds.points)
mds.points$days <- as.factor(mds.points$days)
mds.points$group <- as.factor(mds.points$group)

df.fit <- as.data.frame(fit$vectors$arrows)
df.fit$r <- fit$vectors$r
df.fit$pvals <- fit$vectors$pvals

df.fit <- df.fit[df.fit$r > 0.15 & df.fit$pvals < 0.05,]

# Save MDS family plot
svg("graphs/MDS.family.svg", width = 6.5, height = 6.5)
ggplot(mds.points, aes(MDS1, MDS2, col = group, shape = days))+
      geom_point(size = 2.5)+
      theme_linedraw()+
      theme(legend.position="bottom")+
      stat_ellipse(alpha = 0.8, aes(col = group, linetype = days), size = 0.7)+
      scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      scale_fill_manual(values = c("dodgerblue4","firebrick4"))+
      scale_shape_manual(values = c(16,2,17,0))+
      geom_line(mds.points, mapping = aes(MDS1, MDS2, group = patientID), col = "black", alpha = 0.5)+
      annotate("text", x = 0.3091216, y = 0.977889651, label = "Bacteroidaceae")+
      annotate("text", x = -0.9999979, y = 0.052059975, label = "Enterobacteriaceae")+
      annotate("text", x = -0.9890996, y = 0.207248245, label = "Enterococcaceae")+
      annotate("text", x = -0.0627886, y = -0.998026850, label = "Prevotellaceae")+
      annotate("text", x = 0.5123401, y = -0.858782643, label = "Rikenellaceae")+
      annotate("text", x = 0.03771343, y = -1.05, label = "Veillonellaceae")+
      annotate("text", x = -0.02914262, y = 1.05, label = "Porphyromonadaceae")+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=0.2091216, yend=0.977889651), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=-0.9999979, yend=0.002059975), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=-0.9890996, yend=0.147248245), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=-0.0627886, yend=-0.998026850), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=0.5123401, yend=-0.858782643), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=0.03771343, yend=-1.05), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)+
      geom_segment(aes(x=mean(mds.points$MDS1), y=mean(mds.points$MDS2), xend=-0.02914262, yend=1.05), 
                   arrow = arrow(angle = 12, type = "closed"), col = "black", size = 0.2)
dev.off()

meta.data.sbs
##############################################################################################################################################################
# 7. LefSe results processing ################################################################################################################################
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
# lefse.gen <- lefse.gen[order(lefse.gen$Feature),]

write.table(lefse.gen, "output/lefse.gen.table.txt", quote = F, row.names = F, sep = ",")

# Make barplot #
svg("graphs/lefse.gen.svg", width = 5, height = 3)
par(mar = c(4.5,10.5,0,0.5))
barplot(name = lefse.gen$Feature,
        height = lefse.gen$LDA_score,
        horiz = T,
        las = 1,
        col = "darkred",
        cex.axis = 1.5,
        xlab = "LDA score (log 10)")
dev.off()

# Species level #
colnames(lefse.org) <- colnames(lefse.gen)
lefse.org <- lefse.org[lefse.org$Class == "point_1",]
lefse.org$p.value <- as.numeric(lefse.org$p.value)
#lefse.org <- lefse.org[lefse.org$p.value < 0.01,]
lefse.org <- lefse.org[order(lefse.org$LDA_score),]
# lefse.gen <- lefse.gen[order(lefse.gen$Feature),]

write.table(lefse.org, "output/lefse.org.table.txt", quote = F, row.names = F, sep = ",")

# Make barplot #
svg("graphs/lefse.org.svg", width = 7, height = 5)
par(mar = c(4.5,16.8,0,2.5))
barplot(name = lefse.org$Feature,
        height = lefse.org$LDA_score,
        horiz = T,
        xlim=c(0,4),
        las = 1,
        col = "darkred",
        cex.axis = 1.5,
        xlab = "LDA score (log 10)")
dev.off()


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
## 8. Opportunistic and pathogens diff between first and second time points ##################################################################################
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
i <- "Enterococcus_faecium"
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
# 9. MDSs for different taxonomic levels (firs and secon time points) ########################################################################################
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
