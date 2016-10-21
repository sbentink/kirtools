#' @include kirallele.R

library("ggplot2")
exmat.cols <- paste0(c("utr",rep(c("Exon ","Intron "),8),"Exon ","utr"),
                     c(5,rep(1:8,each=2),9,3))
exmat <- sapply(hlatools::HLARangesList(lapply(result,function(x) x$feat)),function(feat) exmat.cols %in% names(feat))
exmat.num <- apply(exmat,2,as.integer)
fu <- sapply(result,function(x) x$fullSeq=="TRUE")
##image(1-apply(exmat[,fu],2,as.numeric))


gg.convert <- function(my.matrix) {
  x <- as.numeric((1:ncol(my.matrix)))
  y <- as.numeric((1:nrow(my.matrix)))
  my.grid <- expand.grid(x=x,y=y)
  z <- sapply(1:nrow(my.grid),function(i) my.matrix[my.grid[i,2],my.grid[i,1]])
  my.grid$z <- z
  return(my.grid)
}

exmat.full <- exmat.num##[,fu]
exmat.clust <- hclust(dist(t(exmat.full),method="binary"))
##exmat.clust$order

gg.mat <- gg.convert(exmat.full[,exmat.clust$order])
gg.mat[,3] <- as.factor(gg.mat[,3])

colnames(gg.mat) <- c("Sequence","Feature","ExonPresent")##,"Order","Full")

ggplot(gg.mat) + geom_tile(aes(x=Feature,y=Sequence,fill=ExonPresent)) + scale_fill_discrete() +
  scale_x_continuous(expand = c(0,0) , limits = c(0,max(gg.mat$Feature))+0.5,
                     breaks=seq(1,length(exmat.cols)),labels=exmat.cols) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(gg.mat$Sequence))+0.5) +
  theme(axis.text.y=element_blank()) +
  theme(axis.ticks.y=element_blank())

head(gg.mat)
##table(gg.mat[,4])
