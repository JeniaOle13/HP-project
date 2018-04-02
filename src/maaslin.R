#################################################################################################################################################################################
maaslin.gen <- read.csv("output/maaslin_genus/genus-Time_point_number.txt", sep = "\t")
maaslin.gen <- maaslin.gen[which(!is.na(maaslin.gen$Q.value)),]
maaslin.gen <- maaslin.gen[c(2,4,8)]
maaslin.gen <- maaslin.gen[order(maaslin.gen$Coefficient),]
maaslin.gen <- maaslin.gen[maaslin.gen$Q.value < 0.05,]
maaslin.gen$Feature <- factor(maaslin.gen$Feature,levels = rev(as.character(maaslin.gen$Feature)))
maaslin.gen$Taxa[which(!is.na(str_extract(as.character(maaslin.gen$Feature), "\\_")))] <- "Org"
maaslin.gen$Taxa[is.na(maaslin.gen$Taxa)] <- "Gen"
maaslin.gen$Taxa <- as.factor(maaslin.gen$Taxa)

svg("graphs/Maaslin.output.svg", width = 5, height = 6)
ggplot(maaslin.gen, aes(Feature, -Coefficient, fill = Taxa)) +
      geom_bar(stat="identity")+
      coord_flip()+
      theme_linedraw()+
      theme(legend.position="bottom")+
      # scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      scale_fill_manual(values = c("dodgerblue4","firebrick4"))+
      ggtitle("eradication therapy, q-value<0.05")
dev.off()

#################################################################################################################################################################################
maaslin.gen <- read.csv("output/maaslin_genus/genus-Gender.txt", sep = "\t")
maaslin.gen <- maaslin.gen[which(!is.na(maaslin.gen$Q.value)),]
maaslin.gen <- maaslin.gen[c(2,4,8)]
maaslin.gen <- maaslin.gen[order(maaslin.gen$Coefficient),]
maaslin.gen <- maaslin.gen[maaslin.gen$Q.value < 0.2,]
# maaslin.org <- maaslin.gen
maaslin.gen <- rbind(maaslin.gen, maaslin.org)
maaslin.gen <- maaslin.gen[order(maaslin.gen$Coefficient),]
maaslin.gen$Feature <- factor(maaslin.gen$Feature,levels = rev(as.character(maaslin.gen$Feature)))
maaslin.gen$Taxa[which(!is.na(str_extract(as.character(maaslin.gen$Feature), "\\_")))] <- "Org"
maaslin.gen$Taxa[is.na(maaslin.gen$Taxa)] <- "Gen"
maaslin.gen$Taxa <- as.factor(maaslin.gen$Taxa)

svg("graphs/Maaslin.Gender.svg", width = 5, height = 4)
ggplot(maaslin.gen, aes(Feature, Coefficient, fill = Taxa)) +
      geom_bar(stat="identity")+
      coord_flip()+
      theme_linedraw()+
      theme(legend.position="bottom")+
      # scale_color_manual(values = c("dodgerblue4","firebrick4"))+
      scale_fill_manual(values = c("dodgerblue4","firebrick4"))+
      ggtitle("gender Male, q-value<0.2")
dev.off()

write.csv(meta.data.sbs, "output/metaSBS.csv", quote = F, row.names = F)