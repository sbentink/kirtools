#' @include kirallele.R

library("ggplot2")
curr.gene <- KIRAlleleFromDB()


gg.convert <- function(my.matrix) {
  x <- as.numeric((1:ncol(my.matrix)))
  y <- as.numeric((1:nrow(my.matrix)))
  my.grid <- expand.grid(x=x,y=y)
  z <- sapply(1:nrow(my.grid),function(i) my.matrix[my.grid[i,2],my.grid[i,1]])
  my.grid$z <- z
  return(my.grid)
}

raw_mat <- apply(get_feature_mat(curr.gene),2,as.numeric)
colnames(raw_mat) <- paste0(c("utr",rep(c("Exon ","Intron "),8),"Exon ","utr"),
                           c(5,rep(1:8,each=2),9,3))
rownames(raw_mat) <- names(curr.gene)
raw_clust <- hclust(dist(raw_mat,method="binary"))


gg.mat <- gg.convert(t(raw_mat[raw_clust$order,]))
gg.mat[,3] <- as.factor(gg.mat[,3])

colnames(gg.mat) <- c("Sequence","Feature","ExonPresent")##,"Order","Full")

gp <- ggplot(gg.mat) + geom_tile(aes(x=Feature,y=Sequence,fill=ExonPresent)) + scale_fill_discrete() +
  scale_x_continuous(expand = c(0,0) , limits = c(0,max(gg.mat$Feature))+0.5,
                     breaks=seq(1,length(exmat.cols)),labels=exmat.cols) +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(gg.mat$Sequence))+0.5) +
  theme(axis.text.y=element_blank()) +
  theme(axis.ticks.y=element_blank()) + ggtitle(hlatools::elementMetadata(curr.gene)$locus[1])

##use layout:

##open new page:
grid::grid.newpage()
##define layout (e.g. 2 rows, 2 columns, play with height):
gl <- grid::grid.layout(2,2,heights=c(2,1))
##make viewport with layout
vp <- grid::viewport(layout=gl)
##push layout to page
grid::pushViewport(vp)
##the plot shoud appear in the top left corner
##sub-vp for the new plot
vps <- grid::viewport(layout.pos.row=1,layout.pos.col=1)
print(gp,vp=vps)
vps2 <- grid::viewport(layout.pos.row=1,layout.pos.col=2)
print(gp,vp=vps2)








