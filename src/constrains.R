strains.case <- read.csv("data/constrains/strains.case.csv")
strains.control.7 <- read.csv("data/constrains/strains.control.7.csv")
strains.control.60 <- read.csv("data/constrains/strains.control.60.csv")

df.all <- spread(strains.case, variable, value)
df.all <- df.all[-1]

df.all$Genera <- sapply(str_split(df.all$Annot, "_"), function(x) head(x, 1))
df.all$Genera <- as.factor(df.all$Genera)
df.all$Days <- as.factor(df.all$Days)

svg("graphs/strain.scatter.case.svg", width = 6.5, height = 7.1)
ggplot(df.all, aes(First_point, Second_point, col = Genera))+
      geom_point(size = 3.0, alpha = 0.85, stroke = 1.5, shape = 20)+
      theme_linedraw()+
      xlab("Before eradication")+
      ylab("After eradication")+
      annotate("rect", xmin = 70, xmax = Inf, ymin = -Inf, ymax = 30,   fill = "blue", alpha = 0.25)+
      annotate("rect", xmin = 25, xmax = Inf, ymin = -Inf, ymax = 10,   fill = "blue", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 10, ymin = 25, ymax = Inf,   fill = "green", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 30, ymin = 70, ymax = Inf,   fill = "green", alpha = 0.25)+
      geom_vline(xintercept = 50, col = "black")+
      geom_hline(yintercept = 50, col = "black")+
      # geom_vline(xintercept = 5, col = "red")+
      # geom_hline(yintercept = 5, col = "red")+
      theme(legend.position = "bottom", legend.box = "vertical")+
      scale_color_manual(values = colorRampPalette(c("darkred", "darksalmon", "coral4", "darkblue", "darkcyan", "darkgreen", "darkslategray"))(14))
dev.off()

###################################################################################################################################################################################

df.all <- spread(strains.control.7, variable, value)
df.all$Days <- 7

df.all2 <- spread(strains.control.60, variable, value)
df.all2$Days <- 60

df.all <- rbind(df.all, df.all2)

df.all <- df.all[-1]

df.all$Genera <- sapply(str_split(df.all$Annot, "_"), function(x) head(x, 1))
df.all$Genera <- as.factor(df.all$Genera)
df.all$Days <- as.factor(df.all$Days)

svg("graphs/strain.scatter.control.svg", width = 6.5, height = 7.7)
ggplot(df.all, aes(First_point, Second_point, col = Genera, shape = Days))+
      geom_point(size = 3.0, alpha = 0.85, stroke = 1.5)+
      theme_linedraw()+
      xlab("Before eradication")+
      ylab("After eradication")+
      annotate("rect", xmin = 70, xmax = Inf, ymin = -Inf, ymax = 30,   fill = "blue", alpha = 0.25)+
      annotate("rect", xmin = 25, xmax = Inf, ymin = -Inf, ymax = 10,   fill = "blue", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 10, ymin = 25, ymax = Inf,   fill = "green", alpha = 0.25)+
      annotate("rect", xmin = -Inf, xmax = 30, ymin = 70, ymax = Inf,   fill = "green", alpha = 0.25)+
      geom_vline(xintercept = 50, col = "black")+
      geom_hline(yintercept = 50, col = "black")+
      # geom_vline(xintercept = 5, col = "red")+
      # geom_hline(yintercept = 5, col = "red")+
      theme(legend.position = "bottom", legend.box = "vertical")+
      scale_shape_manual(values=c(20, 15))+
      scale_color_manual(values = colorRampPalette(c("darkred", "darksalmon", "coral4", "darkblue", "darkcyan", "darkgreen", "darkslategray"))(14))
dev.off()

###################################################################################################################################################################################
df.lol <- read.csv("data/constrains/df.lol")
col.plot <- colorRampPalette(c("darkred", "darksalmon", "coral4", "darkblue", "darkcyan", "darkgreen", "darkslategray"))(25)

svg("graphs/strain.scatter.case0.svg", width = 7.5, height = 8.7)
ggplot(df.lol, aes(log(First_point), log(Second_point), col = Genera, size = N, stroke = 1.5))+
      geom_point(alpha = 0.5)+
      theme_linedraw()+
      theme(legend.position = "bottom", legend.box = "vertical")+
      xlab("Overall relative abundance in 1st time point")+
      ylab("Overall relative abundance in 2nd time point")+
      scale_color_manual(values = col.plot)
dev.off()


