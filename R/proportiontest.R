#' Cell Type Proportion Testing
#'
#' Statistical testing for cell type proportions
#'
#' This function performs statistical testing that tests if the proportions of cell type are differential across groups, and generates corresponding p-value and FDRs.
#' @param celltype A character vector indicating which cell type each cell is.
#' @param sample A character vector indicating the sample that each cell comes from.
#' @param design A numeric matrix for sample design. Row names must be the unique values in \code{sample}. Only the second column of the design will be tested. The intercept term should be added as the first column.
#' @details Within each sample, the cell type proportion is calculated as the number of cells in a cell type/cluster divided by the total number of cells. The proportion is logit transformed and limma is performed to call differential.
#' @import reshape2 limma
#' @export
#' @return A data.frame with three columns: log fold change, p-value and FDR.
#' @author Zhicheng Ji<zhicheng.ji@@duke.edu>, Wenpin Hou, Hongkai Ji
#' @examples 
#' design <- cbind(1,c(1,1,0,0))
#' rownames(design) <- c('person1','person2','person3','person4')
#' res <- proportiontest(celltype=rep(c('celltype1','celltype2'),each=20),sample=rep(c('person1','person2','person3','person4'),10),design = design)

proportiontest <- function(celltype,sample,design) {
  prop <- dcast(as.data.frame(table(celltype,sample)),celltype~sample,value.var='Freq')
  rownames(prop) <- prop[,1]
  prop <- as.matrix(prop[,-1])
  prop <- sweep(prop,2,colSums(prop),'/')
  k <- min(prop[prop > 0])
  prop[prop==0] <- k
  prop[prop==1] <- 1-k
  prop <- log(prop/(1-prop))
  res <- topTable(eBayes(lmFit(prop,design=design[colnames(prop),])),coef = 2,n=nrow(prop))
  res[,c('logFC','P.Value','adj.P.Val')]
}

