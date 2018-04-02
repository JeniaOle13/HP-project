library('ade4')
library('made4')
library('cluster')
library('fpc')
library('ecodist')
library('fossil')
library('ggplot2')
library('edgeR')
library('KEGG.db')
library('piano')
library('rgl')
library('reshape')
library('vegan')
library('vioplot')
library('pvclust')
library('ape')
library('gtools')
library('metagenomeSeq')
library("GO.db")
# get sample tags by nation
nationTags <- function(nat) 
{
  rownames(metaAll[which(metaAll[,'nation'] %in% nat),])
}

# get nation by sample tag
tagNation <- function(tag) 
{
  as.character(metaAll[which(rownames(metaAll) == tag), "nation"])
}

# get region by sample tag
tagRegion <- function(tag) 
{
  as.character(metaAll[which(rownames(metaAll) == tag), "geography"])
}


###make  distance in Spearman correlation  - for internal use
####input: feature vectors (table)
####output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}
#####end

####good colours for heatplot
cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, 'RdBu'))(ncol)
  return(rev(hmcol))
}
#####end

##### Ward hclust - for internal use
####input: dist object 
hclustW<- function(x, method='ward', ...) hclust(x, method=method, ...)
######end


###make JS distance - for internal use
####input: feature vectors (table)
####output: class dist object
distJS<-function(x)
{
	result<-c()
	input<-c()
	covs<-x
	covs<-covs+0.000001 # add pseudocount
	covs_norm<-covs/rowSums(covs) #normalize input by 1
	input = t(covs_norm)
	heads<-list(colnames(input))
	result<-matrix(0,nrow = ncol(input), ncol = ncol(input),dimnames = c(heads, heads)) #make result distance table
	kld1<-0
	kld2<-0
	for (i in 1:ncol(input))
		{
			for (p in 1:ncol(input)) 
				{
					xj<-input[,i]
					yj<-input[,p]
					mj<-(xj+yj)/2
					kld1<-kld1 + xj*(log(xj/mj)) #compute Kullback-Leibler coefficients
					kld2<-kld2 + yj*(log(yj/mj))
					result[i,p]<-sqrt(sum((kld1/2) ,(kld2/2))) #compute Jenson-Shannon distance
					kld1<-0
					kld2<-0
				}  
		}
	js_d<-as.dist(result)
	return (js_d)
}
#######end


###make all distances
####input: feature vectors (table), vector with names of metrics ('JS','Spear','Eu','Man','Can','BC') chosen
####output: a list with class dist objects, names of the objects - like metrics chosen
allDist<-function(x, metric=c('JS','Spear','Eu','Man','Can','BC'))
{
	dists<-list()
	num<-0
	for (i in (metric)) 
		{
			{
				num<-num+1
				if (i == 'JS')
					{
						dis<-distJS(x)
					}
				if (i == 'Spear')
					{
						dis<-distSpear(x)
					}
				if (i == 'Eu')
					{
						dis<-dist(x)
					}
				if (i == 'Man')
					{
						dis<-dist(x, method='manhattan')
					}
				if (i == 'Can')
					{
						dis<-dist(x, method='canberra')
					}
				if (i == 'BC')
					{
						dis<-bcdist(x)
					}
				dists[[num]]<-dis
			}	
		}	
names(dists)<-metric
return (dists)
}
########end

###choose top features with total % of abundance across all samples
####input: feature vectors (table), percantage
####output: table with chosen features
chooseTOPfeature<-function(dat,perc)
{
	summof<-sum(colSums(dat))
	sorted<-sort(colSums(dat), decreasing = TRUE)
	count <-0
	num<-0
	name_list<-c()  
	for (i in sorted) #in sorted by total abundance features, take those making %
		{
			count<-count+((i/summof)*100)
			num <-num+1
			if (count>perc) 
				{
					break
				}
			else 
				{
					name_list<-append(name_list, names(sorted[num]))
				}		
	}	
	y<-dat[,name_list]
}
####end


###make all PAM clusterings and return clustering
####input: list of dist objects (from allDist function), outdir 
####output: table with samples as rows and cluster numbers as columns, number of columns=number of distances, column names=distances names
clusPAM <- function(distList) 
{
	allClus<-c()
	sil<-c()

	for (dis in 1:length(distList)) #for each distance object in input
	{
		clus<-pamk(distList[[dis]], diss=T, criterion='ch') #make PAM clustering
		clus<-clus$pamobject		
		nas<-names(distList)[dis]
		sil<-append(sil, clus$silinfo$avg.width) # get ASW
		clus<-t(t(clus$clustering))
		allClus<-cbind(allClus,clus)
	}
	colnames(allClus)<-names(distList)
	silfo<-cbind(names(distList),sil)
	print(silfo)
	return(allClus)
}
########end


### make all PAM clusterings and return ASW
####input: list of dist objects (from allDist function), outdir 
####output: ASW
clusPAM2 <- function(distList) 
{
  allClus<-c()
  sil<-c()
  
  for (dis in 1:length(distList)) #for each distance object in input
  {
    clus<-pamk(distList[[dis]], diss=T, criterion='ch') #make PAM clustering
    clus<-clus$pamobject		
    nas<-names(distList)[dis]
    sil<-append(sil, clus$silinfo$avg.width) # get ASW
    clus<-t(t(clus$clustering))
    allClus<-cbind(allClus,clus)
  }
  colnames(allClus)<-names(distList)
  silfo<-cbind(names(distList),sil)
  return(silfo)  
}
########end


###########Wilcoxon test for cohorts
####input: feature vectors (table), metadata (table), vector with colnames of chosen cohorts in meta
####output: list with 2 tables (p-value meanings, statistic menaings) each table has as much columns as factors. All the values passed limits
cohWilcox<-function(featVectIn, metaData, cohorts)
{
	pv<-colnames(featVectIn)
	stat<-colnames(featVectIn)
	chosenCoh<-metaData[,cohorts]
	namesForCol<-'feature'
	for (i in colnames(chosenCoh)) #for each factor column in cohorts
		{
			featVect<-featVectIn[rownames(chosenCoh[which(!(chosenCoh[,i]%in% c('na','to_count','no','U'))),]),] # choose those rows in feature table that have no 'na', 'no' or other meaningless values for current factor
			toSplit<-as.character(chosenCoh[rownames(featVect),i])
			toSplit<-factor(toSplit)
			num<-0
			for (k in levels(toSplit)) 
				{
					num<-num+1
					if ((num==2)&&(length(levels(toSplit))==2)) #if the factor has only 2 levels, there is no need to calculate for the second
						{
							break
						}
					else
						{
							pval_gen<-c()
							level<-k
							namesForCol<-append(namesForCol,level)					
							for (j in colnames(featVect))
								{
									wil<-wilcox.test(as.matrix(featVect[toSplit==level,j]),as.matrix(featVect[toSplit!=level,j]))
									pval_gen<-append(pval_gen,wil$p.value) #add pvalue meaning
								}
							pval_gen<-p.adjust(pval_gen, method='fdr')
							fin_pval<-c()
							for (p in 1:length(pval_gen))
								{				
									if (pval_gen[p]<0.01) # if the data is above thresholds
										{
											fin_pval<-append(fin_pval,pval_gen[p])

										}
									else
										{
											fin_pval<-append(fin_pval,'no')						
										}
								}
							pv<-cbind(pv,fin_pval)
					}
				}
		}
		colnames(pv)<-namesForCol
		result<-pv
		return (result)
}
########end






########plot 2D MDS
####input: list of dist objects (from allDist function), meta data, meta colname for colouring, meta col name for making symbols ,outdir
####output: MDS plots
makeMDS<-function(distObj, meta, colFact, symbFact, outdir)
{
	num<-0
	for (dis in distObj)
		{			
			num<-num+1
			myMDS<-isoMDS(dis)

			outdirMDS<-paste(outdir,'/MDS_',names(distObj)[num],'_',as.character(colFact),'.pdf',sep='')
			outdirMDS<-file.path(outdirMDS)
			dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]), cols = meta[,colFact], symb=meta[,symbFact])
			plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, shape = factor(symb), color=factor(cols))) + geom_point(size=2)  + scale_colour_hue(name=eval(colFact)) + scale_shape (name=eval(symbFact)) + theme_bw()		
			ggsave(outdirMDS, plotMDS)
						
		}
		
}
########end

########plot 3D MDS
####input: list of dist objects (from allDist function), meta data, meta colname for colouring
####output: 3D MDS plots
make3DMDS<-function(distObj, meta, colFact, outdir)
{
	for (dis in distObj)
		{
			myMDS<-isoMDS(dis, k=3)
			cols<-c(rainbow(7),colours())
			outdirMDS<-paste(outdir,'/MDS3D_',names(distObj[num]),'.pdf',sep='')
			outdirMDS<-file.path(outdirMDS)
			pdf(outdirMDS, pagecentre = FALSE)
			plot3d(myMDS$points, col=cols, size=.7, type='s')
			texts3d(myMDS$points, adj=c(1,1), texts=rownames(myMDS$points), cex=.1)			
			dev.off()
		}
		
}
########end

#############make heatmap
makeHeat<-function(featVect, meta, colRow,outdir)
{
	cols<-rainbow (length(levels(as.factor(meta[,colRow]))))
	cols<-cbind(levels(as.factor(meta[,colRow])), cols)
	rownames(cols)<-cols[,1]
	colVec<-c()
	for (i in rownames(featVect))
		{
			colVec<-append(colVec, cols[meta[i,colRow],2])
		}
	outdirHeat<-paste(outdir,'/Heatmap_',colRow,'.pdf',sep='')
	outdirHeat<-file.path(outdirHeat)
	pdf(outdirHeat, pagecentre = FALSE)
	heatmap.2(as.matrix(featVect), margin=c(12,12), distfun=distSpear, cexRow=.1,col=cols.gentleman(500), cexCol=.6, hclustfun=hclustW, trace='none',RowSideColors=colVec)
	dev.off()
}
#########end


##########make violin plots, using 1 matrix as input (vioplot function modification)
vio3<-function (x, range = 1.5, h = NULL, ylim = NULL, names = NULL,
          horizontal = FALSE, col = "magenta", border = "black", lty = 1,
          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19,
          at, add = FALSE, wex = 1, drawRect = TRUE, las=0, main=NULL, sub=NULL, ylab=NULL, xlab=NULL)
{
  datas <- as.matrix(x)
  n <- nrow(datas)
  if(length(col)<n && length(col)==1)
  {
    col<-rep(col,n)
  }
  
  if (missing(at))
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h)))
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- as.vector(datas[i,])
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i],
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1)
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <-rownames(x)
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add)
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
              c(base[[i]], rev(base[[i]])), col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label, las=las)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                                              rev(at[i] + height[[i]])), col = col[i], border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
          boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  if (!is.null(main) | !is.null(sub) | !is.null(xlab) | !is.null(ylab)){
  title (main = main, sub=sub)}
  invisible(list(upper = upper, lower = lower, median = med,
                 q1 = q1, q3 = q3))
}
###########end

#################make violin plots using 1 or more vectors as input (vioplot function modification)
vio2<-function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL,
horizontal = FALSE, col = "magenta", border = "black", lty = 1,
lwd = 1, rectCol = "black", colMed = "white", pchMed = 19,
at, add = FALSE, wex = 1, drawRect = TRUE, las=0, main=NULL, sub=NULL, ylab=NULL, xlab=NULL)
{
datas <- list(x, ...)
n <- length(datas)
if(length(col)<n && length(col)==1)
  {
    col<-rep(col,n)
  }
if (missing(at))
at <- 1:n
upper <- vector(mode = "numeric", length = n)
lower <- vector(mode = "numeric", length = n)
q1 <- vector(mode = "numeric", length = n)
q3 <- vector(mode = "numeric", length = n)
med <- vector(mode = "numeric", length = n)
base <- vector(mode = "list", length = n)
height <- vector(mode = "list", length = n)
baserange <- c(Inf, -Inf)
args <- list(display = "none")
if (!(is.null(h)))
args <- c(args, h = h)
for (i in 1:n) {
data <- datas[[i]]
data.min <- min(data)
data.max <- max(data)
q1[i] <- quantile(data, 0.25)
q3[i] <- quantile(data, 0.75)
med[i] <- median(data)
iqd <- q3[i] - q1[i]
upper[i] <- min(q3[i] + range * iqd, data.max)
lower[i] <- max(q1[i] - range * iqd, data.min)
est.xlim <- c(min(lower[i], data.min), max(upper[i],
data.max))
smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
args))
hscale <- 0.4/max(smout$estimate) * wex
base[[i]] <- smout$eval.points
height[[i]] <- smout$estimate * hscale
t <- range(base[[i]])
baserange[1] <- min(baserange[1], t[1])
baserange[2] <- max(baserange[2], t[2])
}
if (!add) {
xlim <- if (n == 1)
at + c(-0.5, 0.5)
else range(at) + min(diff(at))/2 * c(-1, 1)
if (is.null(ylim)) {
ylim <- baserange
}
}
if (is.null(names)) {
label <- 1:n
}
else {
label <- names
}
boxwidth <- 0.05 * wex
if (!add)
plot.new()
if (!horizontal) {
if (!add) {
plot.window(xlim = xlim, ylim = ylim)
axis(2)
axis(1, at = at, label = label, las=las)
}
box()
for (i in 1:n) {
polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
c(base[[i]], rev(base[[i]])), col = col[i], border = border,
lty = lty, lwd = lwd)
if (drawRect) {
lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
lty = lty)
rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
q3[i], col = rectCol)
points(at[i], med[i], pch = pchMed, col = colMed)
}
}
}
else {
if (!add) {
plot.window(xlim = ylim, ylim = xlim)
axis(1)
axis(2, at = at, label = label, las=las)
}
box()
for (i in 1:n) {
polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
rev(at[i] + height[[i]])), col = col[i], border = border,
lty = lty, lwd = lwd)
if (drawRect) {
lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
lty = lty)
rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
boxwidth/2, col = rectCol)
points(med[i], at[i], pch = pchMed, col = colMed)
}
}
}
  if (!is.null(main) | !is.null(sub) | !is.null(xlab) | !is.null(ylab)){
  title (main = main, sub=sub, xlab=xlab, ylab=ylab)}
invisible(list(upper = upper, lower = lower, median = med,
q1 = q1, q3 = q3))
}
###############end

####modification of heatmap.2 in order to enable the use of layout() function
hm3<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
               distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                  "row", "column", "none"), symm = FALSE, scale = c("none", 
                                                                                                                    "row", "column"), na.rm = TRUE, revC = identical(Colv, 
                                                                                                                                                                     "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
                                                                                                                                                                       scale != "none", col = "heat.colors", colsep, rowsep, 
               sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
               notecol = "cyan", na.color = par("bg"), trace = c("column", 
                                                                 "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
               vline = median(breaks), linecol = tracecol, margins = c(5, 
                                                                       5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
               cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
               key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                                           "density", "none"), denscol = tracecol, symkey = min(x < 
                                                             0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
               xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
               ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
      c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
      c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
                             else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
    1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  # if (missing(lhei) || is.null(lhei)) 
  #      lhei <- c(keysize, 4)
  #  if (missing(lwid) || is.null(lwid)) 
  #     lwid <- c(keysize, 4)
  # if (missing(lmat) || is.null(lmat)) {
  #    lmat <- rbind(4:3, 2:1)
  #   if (!missing(ColSideColors)) {
  #      if (!is.character(ColSideColors) || length(ColSideColors) != 
  #         nc) 
  #        stop("'ColSideColors' must be a character vector of length ncol(x)")
  #     lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
  #        1)
  #   lhei <- c(lhei[1], 0.2, lhei[2])
  #  }
  #     if (!missing(RowSideColors)) {
  #        if (!is.character(RowSideColors) || length(RowSideColors) != 
  #           nr) 
  #          stop("'RowSideColors' must be a character vector of length nrow(x)")
  #     lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
  #        1), 1), lmat[, 2] + 1)
  #   lwid <- c(lwid[1], 0.2, lwid[2])
  #  }
  #     lmat[is.na(lmat)] <- 0
  #}
  #if (length(lhei) != nrow(lmat)) 
  #      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  # if (length(lwid) != ncol(lmat)) 
  #     stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  # op <- par(no.readonly = TRUE)
  # on.exit(par(op))
  # layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1], 
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
      1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}


# sort table by row names
sort_tab_by_rownames <- function (tab) 
{
  tab[order(rownames(tab)),] 
}

# sort table by column names
sort_tab_by_colnames <- function (tab) 
{
  tab[,order(colnames(tab))] 
}

settlemRusTags <- function(urb) 
{
  rownames(metaRus[which(metaRus[,'urban'] %in% urb),])
}

#############end



###################################
# Detect differentially abundant pathways using piano and Reporter Feature Algorithm
# Input: 2 lists of tags to compare, KO abundance matrix, adj. p-value threshold
# output: list of differentiating pathways, up and down
# Warning!!! Those global variables should be available:
# kegg_parsed - selected KEGG orthologs and pathways
PathwayDiff <- function(gr1, gr2, KoData, myNPerm = 10000) 
{ 
  # find differentially abundant KOs
  pval_gen<-c() # p-values
  lfc<-c() # sign of change (-1 = gr1 > gr2, 1 = gr1 < gr2, 0 = otherwise)
  for (i in colnames(KoData))
  {    
    w1<-wilcox.test(as.matrix(KoData[gr1,i]), as.matrix(KoData[gr2,i]), alternative='less')
    w2<-wilcox.test(as.matrix(KoData[gr1,i]), as.matrix(KoData[gr2,i]), alternative='greater')
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
  pval_gen <- p.adjust(pval_gen, method='fdr') #FDR adjustment for p-values
  pval <- pval_gen
  names(pval) <- colnames(KoData)
  names(lfc) <- colnames(KoData)
  myGsc <- loadGSC(kegg_parsed_inv) #load pathway to KO connection table
  
  ##ADD value 1 to pvalue vector for those KO, which are not present in reference but present in kegg_parsed_inv, same (but zero) for lfc 
  pval_ap<-append(pval, rep(1,length(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(pval)),1])))) 
  names(pval_ap)<-append(names(pval), as.vector(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(pval)),1])))
  lfc_ap<-append(lfc, rep(0,length(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(lfc)),1]))))
  names(lfc_ap)<-append(names(lfc), as.vector(unique(kegg_parsed_inv[-which(kegg_parsed_inv[,1] %in% names(lfc)),1])))
  
  myPval <- pval_ap
  myFC <- lfc_ap  
  
  #run Reporter Feature Algorythm from piano package
  gsaRes <- runGSA(myPval, myFC, gsc=myGsc,
                   geneSetStat="reporter",
                   signifMethod="geneSampling",
                   nPerm=myNPerm)                 
  gsaRes
}  

# extract result table  + pathway description
PathwayDiff_table <- function(gsaRes, pvThresh = 0.05) 
{   
  gsaResTab <- GSAsummaryTable(gsaRes)
  
  # only take directional pathways
  ii <- which((gsaResTab[,5] < pvThresh) | (gsaResTab[,8] < pvThresh))
  
  sf <- gsaResTab[ii,]
  xx <- as.list(KEGGPATHID2NAME)
  sf2 <- sub("ko", "", sf$Name)
  #sf <- cbind(data.matrix(xx[sf2]), sf)
  sf <- data.frame(data.frame(unlist(xx[sf2])), sf)
  
  # not including Stats
  sf<-sf[,c(1,2,3,6,9,13,17)]  
  colnames(sf) <- c("Description", "ID", "Genes tot", "p adj dist.dir.up", "p adj dist.dir.dn", "Genes up", "Genes down")
  sf <- sf[order(sf[,5]),]
  
  ## including Stats
  #sf<-sf[,c(1,2,3,4,6,7,9,13,17)]
  #colnames(sf) <- c("Description", "ID", "Genes tot", "Stat dist.dir.up", "p adj dist.dir.up", "Stat dist.dir.dn", "p adj dist.dir.dn", #"Genes up", "Genes down")  
  #sf <- sf[order(sf[,4]),]
  
  sf  
}
#################

############# Find differentiating EC between 2 groups
#The function makes u-test for each EC, cuts results by p-value<0.01 and 
#then on chosen features make random Forest
#Input: group 1 tags, group 2 tags, data table, factors for groups
#Output: table with EC, their description, direction (-1 - more in first group, 1 - more in second group), connection to pathways, p-values and MAD
compEC<-function (gr1, gr2, EcData, fac)
{
  pval_gen<-c() # p-values
  lfc<-c() # sign of change (-1 = gr1 > gr2, 1 = gr1 < gr2, 0 = otherwise)
  for (i in colnames(EcData))
  {
    w1<-wilcox.test(as.matrix(EcData[gr1,i]), as.matrix(EcData[gr2,i]), alternative='less')
    w2<-wilcox.test(as.matrix(EcData[gr1,i]), as.matrix(EcData[gr2,i]), alternative='greater')
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
  pval_gen <- p.adjust(pval_gen, method='fdr')
  pval <- pval_gen
  names(pval) <- colnames(EcData)
  names(lfc) <- colnames(EcData)
  res<-cbind(pval,lfc)
  res<-res[which(res[,'pval']<0.01),]
  desc<-c()
  maps<-c()
  for (j in rownames(res)) #append descriptions to EC from ecDesc table
  {
    if (j %in% rownames(ecDesc))
    {
      cur_desc<-as.vector(ecDesc[j,])
    }
    else
    {
      cur_desc<-'Description not available'
    }
    desc<-append(desc,cur_desc)
    map_str<-paste(ec2map[j], collapse='|')
    maps<-append(maps,map_str)
  }
  ec<-sub('ko','',j)
  res<-cbind(res,desc,maps)  
  res<-cbind(rownames(res),res)
  
  RF<-randomForest(EcData[c(gr1,gr2),rownames(res)],fac,importance=T,type='classification', ntree=100000) #make randomForest on those features that passed p-value threshold
  RF<-RF$importance[,3]  
  RF<-t(t(RF))
  
  mad<-c() #extract Mean Accuracy Decrease values from randomForest results and append them to the final table
  for (i in rownames(res))
  {
    if(i %in% rownames(RF))
    {
      mad<-append(mad, RF[i,1])
    }
    else
    {
      mad<-append(mad,'no')
    }
  }
  res<-cbind(res,mad)
  return(res)
}
########################

###Calculate metastatistics
calcStat<-function(field,nation)
{
  calc2<-c()
  calc<-metaAll[nationTags(nation),field]
  for (i in calc){
    if(i=='na'|| i=='no'){
      
    } 
    else{
      calc2<-append(calc2,i)}
  }
  calc2<-as.numeric(calc2)
  print(median(calc2))
  print(mean(calc2))
  print(sd(calc2))
  print(min(calc2))
  print(max(calc2))
}



# pick top genera for each sample
get_top_genera_with_perc <- function(g)
{   
  topg <- matrix(0, nrow(g), 3)
  topg_perc <- matrix(0, nrow(g), 3)
  for(i in 1:nrow(g))
  {
    topg[i,] <- names(sort(g[i,], decreasing=TRUE)[1:3])
    topg_perc[i,] <- sort(g[i,], decreasing=TRUE)[1:3]    
  }
  rownames(topg) <- rownames(g)
  rownames(topg_perc) <- rownames(g)
  topg_perc <- round(topg_perc, 2)  
  cbind(topg, topg_perc)
}


# pick top genera for each sample
get_top_genomes_for_genera <- function(g, genome_genus)
{
  topg <- matrix(0, nrow(g), NUM_TOP_GENERA)
  topg_perc <- matrix(0, nrow(g), NUM_TOP_GENERA)
  topgw <- matrix(0, 1, NUM_TOP_GENOMES*NUM_TOP_GENERA)
  topgw_perc <- matrix(0, 1, NUM_TOP_GENOMES*NUM_TOP_GENERA)
  for(i in 1:nrow(g))
  {
    topg[i,] <- names(sort(g[i,], decreasing=TRUE)[1:NUM_TOP_GENERA])
    topg_perc[i,] <- sort(g[i,], decreasing=TRUE)[1:NUM_TOP_GENERA]  
    tt <- matrix(0, 1, NUM_TOP_GENOMES*NUM_TOP_GENERA)
    tt_perc <- matrix(0, 1, NUM_TOP_GENOMES*NUM_TOP_GENERA)    
    idx_tt <- 1
    for(k in 1:NUM_TOP_GENERA)   # for each top genera, pick its top genomes
    {    
      cur <- gw[i,]
      cur <- cur[which(genome_genus[,1] == topg[i,k])]    
      cn <- colnames(sort(cur, decreasing=TRUE)[1:min(NUM_TOP_GENOMES, length(cur))])    
      cn_perc <- as.matrix(sort(cur, decreasing=TRUE)[1:min(NUM_TOP_GENOMES, length(cur))])
      tt[,idx_tt:(idx_tt+length(cn)-1)] <- cn            
      tt_perc[,idx_tt:(idx_tt+length(cn_perc)-1)] <- cn_perc
      idx_tt <- idx_tt + NUM_TOP_GENOMES        
    }  
    topgw <- rbind(topgw, tt)
    topgw_perc <- rbind(topgw_perc, tt_perc)    
  }
  rownames(topg) <- rownames(g)
  rownames(topg_perc) <- rownames(g)
  topgw <- topgw[-1,]
  topgw_perc <- topgw_perc[-1,]  
  topgw_perc <- round(topgw_perc, 2)
  list(topgw, topgw_perc)  
}



rp.joint.fill<-function (DF, map.var, id.type.rp = "samples") 
{
  if (all(is.na(DF[DF$id.type == id.type.rp, map.var]))) {

      temp.vec <- as(DF[, map.var], "character")
      temp.vec[is.na(temp.vec)] <- id.type.rp
      DF[, map.var] <- factor(temp.vec)
    
  }
  return(DF)
}

summarize_taxa <- function(o, level, taxo)
{
  #level <- "family"
  q <- as.data.frame(o)
  q$sample <- rownames(o)
  q.m <- melt(q, id.vars = 'sample')
  q.mm <- merge(q.m, taxo[, c("org", level)], by.x = 'variable', by.y = 'org')
  head(q.m)
  head(taxo[, c("org", level)])
  head(q.mm)  
  colnames(q.mm)[4] <- "level"  
  q.a <- aggregate(data = q.mm, value ~ level + sample, sum)
  head(q.a)  
  q.c <- cast(q.a, sample ~ level, value = 'value')
  rownames(q.c) <- q.c[,"sample"]
  res <- q.c[,-1]
  res  
}
