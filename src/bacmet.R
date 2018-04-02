df.vf <- read.csv("data/bacmet/bacmet.rpkm.txt", sep = " ", stringsAsFactors = F)
df.vf <- df.vf[df.vf$rpkm !=0,]

df.gsa <- read.csv("data/bacmet/bacmet.dmd.out", sep = "\t", header = F)[c(1,2)]
df.gsa$V2 <- sapply(str_split(df.gsa$V2, "\\|"), function(x) x[2])

df.gsa <- loadGSC(df.gsa)

df.vf <- df.vf[duplicated(df.vf$name), ]

df.vf.s <- spread(df.vf[-c(6178, 6175, 6179, 6184, 6187),], name, rpkm)
df.vf.s <- df.vf.s[df.vf.s$sample %in% c(Point_1, Point_2),]
rownames(df.vf.s) <- df.vf.s$sample
df.vf.s <- df.vf.s[-1]

df.vf.s[is.na(df.vf.s)] <- 0

KoData <- df.vf.s

gr1 <- Point_2
gr2 <- Point_1

pval_gen <- c()
lfc <- c()
for (i in colnames(KoData))
{    
      w1<-wilcox.test(as.matrix(KoData[gr1,i]), as.matrix(KoData[gr2,i]), alternative='less', paired = T)
      w2<-wilcox.test(as.matrix(KoData[gr1,i]), as.matrix(KoData[gr2,i]), alternative='greater', paired = T)
      # stub for NA p-values
      if(is.na(w1$p.value) || is.na(w2$p.value))
      {
            if (is.na(w1$p.value) && !is.na(w2$p.value)) {
                  pval_gen<-append(pval_gen, w2$p.value)
                  lfc<-append(lfc, -1)
            }
            if (is.na(w2$p.value) && !is.na(w1$p.value)) {
                  pval_gen<-append(pval_gen, w1$p.value)
                  lfc<-append(lfc, 1)
            }
            if (is.na(w2$p.value) && is.na(w1$p.value)) {
                  pval_gen<-append(pval_gen, 1)
                  lfc<-append(lfc, 0)
            }
      }
      else # detect direction of change 
      {
            
            if(w1$p.value < w2$p.value)
            {
                  lfc<-append(lfc, 1)
                  pval_gen<-append(pval_gen, w1$p.value)
            }
            else
            {
                  lfc<-append(lfc, -1)
                  pval_gen<-append(pval_gen,  w2$p.value)
            }
            
      }  
}

pval <- pval_gen
pval_gen <- p.adjust(pval, method='fdr')

names(pval) <- colnames(KoData)
names(lfc) <- colnames(KoData)

# pval_ap<-append(pval, rep(1,length(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(pval)),1])))) 
# names(pval_ap)<-append(names(pval), as.vector(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(pval)),1])))
# lfc_ap<-append(lfc, rep(0,length(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(lfc)),1]))))
# names(lfc_ap)<-append(names(lfc), as.vector(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(lfc)),1])))

myPval <- pval
myFC <- lfc  

#run Reporter Feature Algorythm from piano package
gsaRes <- runGSA(myPval, myFC, gsc=df.gsa, adjMethod = "fdr",
                 geneSetStat="reporter",
                 signifMethod="geneSampling",
                 nPerm=21000)                 

lol <- GSAsummaryTable(gsaRes)

lol <- lol[c(1,5,8)]
lol[lol$`p adj (dist.dir.up)` < 0.05,]
lol[lol$`p adj (dist.dir.dn)` < 0.05,]

write.table(lol[lol$`p adj (dist.dir.up)` < 0.05,], "lol.tsv", quote = F)
