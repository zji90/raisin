#' RAISIN Testing
#'
#' Statistical testing for RAISIN model
#'
#' This function performs statistical testing that tests if a coefficient or a combination of coefficients in the design matrix of fixed effects is 0, and generates corresponding p-value and FDRs.
#' @param fit The output from the \code{raisinfit} function.
#' @param coef Only useful if \code{contrast} is \code{NULL}. A numeric value indicating which coefficient is to be tested. 
#' @param contrast A numeric vector indicating the combination of coefficients. Must have the same length as the number of columns of the fixed effect design matrix \code{X}, and the order must agree with its columns.
#' @param fdrmethod A character value of which method to use for multiple testing. Need to be one of the method available in \code{p.adjust}.
#' @details The test statistic follows a t-distribution. A permutation procedure is used to estimate the degree of freedom of the t-distribution. The p-value is calculated as the tail probability of the t-distribution.
#' @export
#' @return A data.frame with four columns: (log-) fold change, statistics, p-value and FDR.
#' @author Zhicheng Ji<zhicheng.ji@@duke.edu>, Wenpin Hou, Hongkai Ji
#' @examples 
#' fit <- raisinfit(expr = matrix(rnorm(100*8),nrow=100,dimnames=list(paste0('gene_',1:100),paste0('cell_',1:8))),sample=rep(c('person1_celltype1','person2_celltype1','person1_celltype2','person3_celltype2'),each=2),testtype='paired',design = data.frame(sample=c('person1_celltype1','person2_celltype1','person1_celltype2','person3_celltype2'),individual=c('person1','person2','person1','person3'),feature=c('celltype1','celltype1','celltype2','celltype2')))
#' res <- raisintest(fit)

raisintest <- function(fit,coef=2,contrast=NULL,fdrmethod='fdr') {
  X <- fit$X
  means <- fit$mean
  G <- nrow(means)
  group <- fit$group
  Z <- fit$Z
  if (is.null(contrast)) {
    contrast <- rep(0,ncol(X))
    contrast[coef] <- 1
  }
  k <- t(contrast) %*% solve(t(X) %*% X) %*% t(X)
  b <- (means %*% t(k))[,1]
  if (identical(unique(fit$group),fit$failgroup)) {
    warning('Unable to estimate variance for all random effects. Setting FDR to 1.')
    res <- data.frame(Foldchange=b,FDR=1,stringsAsFactors=F)
    res[order(-abs(res[,1])),]
  } else {
    a <- colSums((k %*% Z)[1,]^2 * t(fit$sigma2[,group])) + colSums(k[1,]^2 * t(fit$omega2))
    stat <- b/sqrt(a)
    
    simustat <- unlist(lapply(1:10,function(simuid) {
      perX <- X[sample(1:nrow(X)),]
      k <- t(contrast) %*% solve(t(perX) %*% perX) %*% t(perX)
      a <- colSums((k %*% Z)[1,]^2 * t(fit$sigma2[,group])) + colSums(k[1,]^2 * t(fit$omega2))
      (means %*% t(k))[,1]/sqrt(a)
    }))
    
    pnorm <- sum(dnorm(simustat,log=T))
    pt <- sapply(seq(1,100,0.1),function(dt) {
      sum(dt(simustat,df=dt,log=T))  
    })
    
    if (max(pt) > pnorm) {
      df <- seq(1,100,0.1)[which.max(pt)]
      pval <- pt(abs(stat),df,lower.tail = F) * 2
    } else {
      pval <- pnorm(abs(stat),lower.tail = F) * 2
    }
    
    fdr <- p.adjust(pval,method=fdrmethod)
    res <- data.frame(Foldchange=b,stat=stat,pvalue=pval,FDR=fdr,stringsAsFactors=F)
    res <- res[order(res[,4],-abs(res[,2])),]
    res
  }
}
