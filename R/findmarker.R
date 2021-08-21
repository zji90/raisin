#' Marker Gene for Cell Clusters
#'
#' Identify marker genes for each cell cluster using RAISIN model
#'
#' For each cell cluster, this function identifies differential genes between cells within and outside the cluster. The differential is done for all clusters
#' @param expr A numeric matrix of gene expression. Each row is a gene and each column is a cell.
#' @param individual A character vector indicating which individual each cell comes from, corresponding to the columns of \code{expr}.
#' @param cluster A numeric/character vector indicating the cell type/cluster of each cell, corresponding to the columns of \code{expr}.
#' @export
#' @return A list of elements. Each element correponds to one cluster and is a data.frame with six columns. The first two columns are the averaged expression of cells within or outside the cluster. The last four columns are output of \code{raisinfit}.
#' @author Zhicheng Ji<zhicheng.ji@@duke.edu>, Wenpin Hou, Hongkai Ji
#' @examples 
#' res <- findmarker(expr = matrix(rnorm(100*20),nrow=100,dimnames=list(paste0('gene_',1:100),paste0('cell_',1:20))),individual=rep(c('person1','person2'),each=10),cluster=rep(c('celltype1','celltype2'),10))

findmarker <- function(expr,individual,cluster) {
  cn <- colnames(expr)
  uc <- sort(unique(cluster))
  l <- sapply(uc,function(s) {
    ct <- cluster==s
    sample <- paste0(individual,':',ct)
    samplename <- unique(sample)
    group <- sub('.*:','',samplename)
    
    means <- rowsum(t(expr),sample)
    tab <- table(sample)
    means <- means/as.vector(tab[rownames(means)])
    tab <- table(group)
    means <- rowsum(means[samplename,],group)
    means <- means/as.vector(tab[rownames(means)])
    rownames(means)[rownames(means)=='TRUE'] <- 'Withincluster'
    rownames(means)[rownames(means)=='FALSE'] <- 'Outsidecluster'
    design <- data.frame(sample=unique(sample),individual=sub(':.*','',unique(sample)),feature=sub('.*:','',unique(sample)))
    res <- raisintest(raisinfit(expr,sample,design,testtype='paired',filtergene=T))
    data.frame(t(means[,rownames(res)]),res)
  },simplify = F)
  names(l) <- uc
  l
}
