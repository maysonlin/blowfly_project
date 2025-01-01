library(mixOmics)
#load sample data
data(srbct)
View(srbct)
##Example
Z <- srbct$gene
W <- srbct$class 
#load your own data 

time <- read.table("~/deg1/time.txt")

#load your own data (from 80~130 hrs, 3 replicates)

colnames(time_reorganize) <- c( "80.1", "80.2","80.3", "90.1", "90.2", "90.3", "90.4", "90.5", "90.6", "100.1", "100.2", "100.3","100.4", "100.5", "100.6","110.1", "110.2", "110.3","110.4", "110.5", "110.6", "120.1", "120.2","120.3", "120.4", "120.5", "120.6", "130.1", "130.2","130.3", "130.4", "130.5", "130.6")



## set X (time) and Y (factors)
X <- t(time_reorganize)

X <- t(time_80_to_110)
timecourse <- c("80","80","80","90","90","90","90","90","90","100","100","100","100","100","100","110","110","110","110","110","110","120","120","120","120","120","120","130","130","130","130","130","130")
Y <- factor(timecourse)

timecourse_80_110 <- c("100","100","100","100","100","100","110","110","110","110","110","110","80","80","80","90","90","90","90","90","90")
Y <- factor(timecourse_80_110)

levels(Y)
dim(X)
summary(Y)
length(Y)
sort(table(Y))

##PCA plot

MyResult.pca <- pca(X)
dev.off()
plotIndiv(MyResult.pca, group = agegroup$Group, legend = TRUE)


plotIndiv(MyResult.pca, group = behaviourgroup$Group, legend = TRUE)

#load your own data

X <- t(changedbehaviour)
behaviors <- c("F", "F" ,"F", "F", "F" ,"F", "F", "F" ,"F",  "F", "F" ,"F",  "F", "F" ,"F",  "F", "F" ,"F", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W", "W")

Y <- factor(behaviors)
##maxmize total space
options(max.print=999999)
##maxmize interger space
options(max.print = .Machine$integer.max)

#PLS-DA analysis
MyResult.splsda <- splsda(X, Y, ncomp=2)
plotIndiv(MyResult.splsda)
#cut off 0.7
plotIndiv(MyResult.splsda) 
plotVar(MyResult.splsda)
MyResult.plsda <- splsda(X,Y) 
plotVar(MyResult.splsda, cex = 2)
plotVar(MyResult.splsda, cutoff = 0.94, cex = 2)
  


#assess PLS-DA performance

MyPerf.plsda <- perf(MyResult.splsda, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 10)
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#variance of important score
VIP <- vip(MyResult.splsda)

##PLS-DA on behaviors
plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'behavior',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')


##clustered Image Maps

cim(MyResult.splsda, comp = 1)



VIP <- vip(MyResult.splsda)

ggtitle(VIP)
xlab(gene)
##switch to data frame

VIP.df <- as.data.frame(VIP)

filter.df <- as.data.frame(filter)


ggplot(filter.df, aes(x=BUSCO, y=Numbers.of.transcripts, group=1)) +geom_line() +geom_point()

library(ggplot2)
names(VIP.df)[-1] <- 'Gene'



##variable selection outputs

MyResult.splsda2 <- splsda(X,Y, ncomp=32, keepX=c(15,10,5))
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean',xlim = c(0.0, 0.8), size.title = rel(1),title ="Contribution on comp 1" )

selecV <- selectVar(MyResult.splsda2, comp=1, method="PC")$value

top15_variables_values <- selecV$index



name.var= ID$ID
names(name.var)=rownames(ID$ID)
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')

plotIndiv(MyResult.splsda2, style="3d")

colnames(selected) <- c("100.1", "100.2", "100.3","100.4", "100.5", "100.6","110.1", "110.2", "110.3","110.4", "110.5", "110.6", "120.1", "120.2","120.3", "120.4", "120.5", "120.6", "130.1", "130.2","130.3", "130.4", "130.5", "130.6", "80.1", "80.2","80.3", "90.1", "90.2", "90.3", "90.4", "90.5", "90.6")

select.transcripts <- as.matrix(selected)

colnames(select.transcipts) <- c("100.1", "100.2", "100.3","100.4", "100.5", "100.6","110.1", "110.2", "110.3","110.4", "110.5", "110.6", "120.1", "120.2","120.3", "120.4", "120.5", "120.6", "130.1", "130.2","130.3", "130.4", "130.5", "130.6", "80.1", "80.2","80.3", "90.1", "90.2", "90.3", "90.4", "90.5", "90.6")

rownames(select.transcripts) <- 

heatmap(selected)
library(pheatmap)



heatmap(VIP)
p <- ggplot (VIP.df, aes(x=Gene, y= comp1))+  scale_fill_gradient(low="white", high="blue")
p
py_config()
conda_list(conda = "auto")
repl_python(VIP,)
attr(MyResult.splsda, which=loadings)
str(X)
dimnames(MyResult.splsda)
sample <- splsda(Z, W)
vipsample <- vip(sample)
boxplot(vipsample)

##test PLS-DA on the 15 candidate genes

colnames(testselected) <- c("100.1", "100.2", "100.3","100.4", "100.5", "100.6","110.1", "110.2", "110.3","110.4", "110.5", "110.6", "120.1", "120.2","120.3", "120.4", "120.5", "120.6", "130.1", "130.2","130.3", "130.4", "130.5", "130.6", "80.1", "80.2","80.3", "90.1", "90.2", "90.3", "90.4", "90.5", "90.6")
X <- t(testselected)
timecourse <- c("100","100","100","100","100","100","110","110","110","110","110","110","120","120","120","120","120","120","130","130","130","130","130","130","80","80","80","90","90","90","90","90","90")
Y <- factor(timecourse)

##asess top 15 classification
plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'time',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

#PCA Plot
MyResult.pca <- pca(X)
#load agegroup data
dev.off()
plotIndiv(MyResult.pca, group = agegroup$Group, legend = TRUE)
