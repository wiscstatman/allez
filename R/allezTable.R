## Outputs top GO categories ##

## score >= 0
## z.score is one-sided: z.score < 0 indicate enrichment
## for genes outside of gene set

allezTable <- function(allez.out,          # output from allez
                       n.low=5,            # include only gene sets with at least n.low genes
                       n.upp=500,          # include only gene sets with at most n.upp genes
                       n.cell=2,           # include only gene sets with at least n.cell genes present in the data
                       nominal.alpha=0.01, # type I error rate
                       symbol=FALSE,       # whether to display gene symbols   
                       score.threshold=NA  # display genes within genesets if their score exceeds score_threshold. If NA, display all genes
                       ){
  ## gene list of gene_id, probe_id, or symbol, from set.data ##
  idcol <- ifelse(symbol,3,2)
  ## z.score column ##
  zcol <- grep("z.score",colnames(allez.out$setscores))[1]
  
  threshold.is.na <- is.na(score.threshold)
  score.threshold <- ifelse(threshold.is.na, 0, score.threshold)

   ## Number of genes in list and functional set ##
  nc <- tapply(allez.out$aux$set.data$gscores,
               allez.out$aux$set.data[,1],
               function(x) sum(x > score.threshold))
  G <- length(allez.out$aux$globe)

  ## If set.size==G then z.score=NA ##

 # first, check that the sets are the rights sizes
  ok1 <- (allez.out$setscores$set.size >= n.low) &
    (allez.out$setscores$set.size <= n.upp) & (allez.out$setscores$set.size < G) 
  #
  nset <- sum(ok1) ## how many sets in play

  zthr <- qnorm( 1-nominal.alpha/nset ) ## one-sided Bonferroni corrected

  # second, check extreme Z and interesting action in the set

  ok <-   ok1 & (allez.out$setscores[,zcol] >= zthr) &
         (nc[rownames(allez.out$setscores)] >= n.cell)

  allez.table <- allez.out$setscores[ok,
                 -grep("sd",colnames(allez.out$setscores))]
   
   ## Subset set.data ##
  set.data <- allez.out$aux$set.data[
              allez.out$aux$set.data[,1] %in% rownames(allez.table),]
  set.data <- set.data[order(set.data$gscores,decreasing=TRUE),]

   ## rownames(genes) == rownames(allez.table) ##
  genes <- data.frame(
             pos=tapply(set.data[,idcol],set.data[,1],paste,collapse=";"),
             neg=tapply(set.data[,idcol],set.data[,1],function(x)
               paste(rev(x),collapse=";")))
  allez.table$genes <- if(nrow(allez.table)>0)
    genes[cbind(rownames(allez.table),
    ifelse(allez.table[,grep("z.score",colnames(allez.table))[1]]>0,
           "pos","neg"))] else character(0)

  if(!threshold.is.na){
    set.data <- set.data[set.data$gscores > score.threshold,]
    genes <- data.frame(
             pos=tapply(set.data[,idcol],set.data[,1],paste,collapse=";"),
             neg=tapply(set.data[,idcol],set.data[,1],function(x)
               paste(rev(x),collapse=";")))
    allez.table <- cbind(allez.table,
                    in.set=nc[rownames(allez.table)],
                    in.genes=if(nrow(allez.table)>0)
                       genes[cbind(rownames(allez.table),
                       ifelse(allez.table[,grep("z.score",
                          colnames(allez.table))[1]]>0,"pos","neg"))] else
                          character(0))
   }
    ##allez.table$in.set <- allez.table$set.mean*allez.table$n.genes
  ord <- order(allez.table$set.mean,decreasing=TRUE)
  allez.table[ord,]
 }
