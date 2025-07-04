
allez <- function (scores,
                   lib,
                   idtype = c("ENTREZID", "SYMBOL"),
                   library.loc=NULL,
                   sets = c("GO"),
                   locallist = NULL,
                   collapse = c("full", "partial", "none"),
                   reduce = NULL,
                   setstat = c("mean", "var"),
                   universe = c("global", "local"),
                   transform = c("none", "binary", "rank", "nscore"),
                   cutoff = NULL,
                   annotate = TRUE)
{

  ## Disallow NA's in scores ##
  if (any(is.na(scores))){
    stop("Scores containing NA are not permitted. Please handle and try again.")
    scores[is.na(scores)] <- 0
  }
  
  scorenames <- names(scores)
  sets <- match.arg(sets)
  universe <- match.arg(universe)
  transform <- match.arg(transform)
  collapse <- match.arg(collapse)
  setstat <- match.arg(setstat)
  idtype <- match.arg(idtype)
  
  # cannot perform local test if input local lists
  if(universe=="local" & !is.null(locallist)){
    universe="global"
    warning("universe='local' mode is not available when including local list! using universe='global' mode")
  }
  
  if (any(duplicated(scorenames))) 
    stop("input IDs must be unique")
  if( !is.numeric(scores) ) 
    stop("scores must be numeric")
  vv <- apply( X=as.matrix(scores), MARGIN=2, FUN=var)
  if( min(vv) == 0 )
    stop("no variance in at least one profile" ) 
  if(universe=="local"){
    if(collapse=="partial")
      stop("partial correction not implemented for 'local'")
  }
  if(transform=="binary" & !is.numeric(cutoff))
    stop("cutoff must be numeric when transform = 'binary'")

  ## Default reduce ##
  uscores <- unique(scores)
  if(is.null(reduce))
   reduce <- if(length(uscores)==2 & all(uscores %in% 0:1)) max else median
  if(length(uscores)==2 & !all(uscores %in% 0:1))
    warning("if scores are binary, please convert to {0,1}")

  message("Loading necessary libraries...")
  fn_loadPlatformLibraries( Libraries=lib, library.loc=library.loc )

  set_id <- "go_id"
  is.org <- substr(lib,1,3)=="org"
  orgpkg <- ifelse(is.org,lib[1],get(paste(lib[1],"ORGPKG",sep="")))

  ## ANNOTATION ##
  message("Converting annotations to data.frames ...")
  if(!is.org & collapse != "full"){
     set2probe <- toTable(getDataEnv(name=ifelse(sets=="GO",
                        "GO2ALLPROBES","PATH2PROBE"),lib=lib[1]))
     probe.symbol <- toTable(getDataEnv(name="SYMBOL",lib=lib[1]))
     set2probe <- cbind(set2probe, probe.symbol[
        match(set2probe$probe_id, probe.symbol$probe_id),"symbol",drop=FALSE])
     ## remove probes not in scores vector ##
     set2probe <- set2probe[set2probe$probe_id %in% names(scores),]
   }
  if(is.org | collapse != "none"){
    ## Use org info for ENTREZ TO GO ID ##
    set2eg <- toTable(getDataEnv(name=ifelse(sets=="GO",
              "GO2ALLEGS","PATH2EG"), lib=orgpkg))
    org.symbol <- toTable(getDataEnv(name="SYMBOL",lib=orgpkg))
    set2eg <- cbind(set2eg, org.symbol[
              match(set2eg$gene_id,org.symbol$gene_id),"symbol",drop=FALSE])
    ## Remove gene_ids not on microarray or scores##
     if(!is.org){
       probe2eg <- switch(idtype,ENTREZID=toTable(getDataEnv(name="ENTREZID",lib=lib[1])),
         SYMBOL=getDataEnv(name="SYMBOL",lib=lib[1]))
       probe2eg <- probe2eg[probe2eg$probe_id %in% names(scores),]
       egs <- switch(idtype, ENTREZID=unique(probe2eg[,"gene_id"]), SYMBOL=unique(probe2eg[,"symbol"]))
       set2eg <- switch(idtype, ENTREZID=set2eg[set2eg$gene_id %in% egs,], SYMBOL=set2eg[set2eg$symbol %in% egs,])
     } else {
# local list
		if(!is.null(locallist)){
		if(is.null(names(locallist)))names(locallist) = paste0("L",1:length(locallist))
		localunlist=unlist(locallist)
		locallen=length(locallist)
		localeachlen=sapply(locallist,length)
		newlist <- switch(idtype, ENTREZID=data.frame(cbind(gene_id=localunlist,
			go_id=unlist(sapply(1:locallen,function(k)rep(paste0("Local:",names(locallist)[k]),localeachlen[k]),simplify=F)),
			Evidence="local", Ontology="local",symbol=set2eg$symbol[match(localunlist,set2eg$gene_id)]
		),stringsAsFactors=F),
		SYMBOL=data.frame(cbind(gene_id=set2eg$gene_id[match(localunlist,set2eg$symbol)], 
			go_id=unlist(sapply(1:locallen,function(k)rep(paste0("Local:",names(locallist)[k]),localeachlen[k]),simplify=F)),
			Evidence="local", Ontology="local",symbol=localunlist),stringsAsFactors=F))
		names(newlist)[2] <- set_id
		if(set_id=="path_id")newlist <- newlist[c("gene_id","path_id", "symbol")]
		set2eg <- rbind(newlist,set2eg)
		}
		set2eg <- switch(idtype, ENTREZID=set2eg[set2eg$gene_id %in% names(scores),],
			SYMBOL=set2eg[set2eg$symbol %in% names(scores),])
		}}

  ## SCORES ##
  if(collapse == "full" & !is.org){
      message("Reducing probe data to gene data...")
      pscores <- data.frame(probe2eg,scores=scores[probe2eg$probe_id])
      scores <- switch(idtype, ENTREZID=unlist(tapply(pscores$scores,as.character(pscores$gene_id),
        FUN=reduce,simplify=FALSE)),
        SYMBOL=unlist(tapply(pscores$scores,as.character(pscores$symbol),
        FUN=reduce,simplify=FALSE)))
}
  gscores <- switch(transform,
             none = scores,
             rank = rank(scores),
             nscore = qnorm( (rank(scores))/(length(scores)+1) ),
             binary = {
               warning("cutoff used at collapsed gene level not probe level")
               1 * (scores >= cutoff) })

  ## Gene Scores and Annotation ##
  set.data <- if(!is.org & collapse != "full"){
    set2probe <- unique(set2probe[,c(set_id,"probe_id","symbol")])
    data.frame(set2probe,gscores=gscores[set2probe$probe_id])
  } else {
    set2eg <- unique(set2eg[,c(set_id, 'gene_id', 'symbol')])
    switch(idtype, ENTREZID=data.frame(set2eg,gscores=gscores[set2eg$gene_id]),
      SYMBOL=data.frame(set2eg,gscores=gscores[set2eg$symbol]))
  }
  set.mean <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                      mean, simplify=FALSE))
  set.sd <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                    sd, simplify=FALSE))
  set.size <- table(set.data[[set_id]])
  class(set.size) <- "array"

## Globe variable ##
	cl <- switch(idtype, ENTREZID="gene_id",SYMBOL="symbol")
  globe <- if(sets=="GO"){
    ## GO:0008150 = bio process ##
    ## GO:0003674 = mol function ##
    ## GO:0005575 = cel component ##
    bpmfcc <- c("GO:0008150","GO:0003674","GO:0005575")
    gscores[unique(set.data[set.data[,1] %in% bpmfcc,cl])]
  } else gscores[unique(set.data[,cl])]


  mu.globe <- mean(globe)
  sigma.globe <- sd(globe)
  G <- length(globe)
  E.globe <- fn_getE.Globe(globe=globe)

 if (universe == "global"){
    if (setstat == "mean") {
      dd <- sigma.globe * fact(G=G, m=set.size)
      z.score <- (set.mean - mu.globe)/dd
    }
  
    if (setstat == "var") {
      ok <- set.size > 3
      sigma1 <- sapply(X = set.size[ok], FUN = sigma.fun, 
                    E = E.globe, G = G)
      z.score <- (set.sd[ok]^2 - (sigma.globe^2))/sigma1
    }
  
    if (!is.org & collapse == "partial") {
      set.ng <-  switch(idtype, ENTREZID=table(unique(set2eg[,c(set_id,"gene_id")])[,1]),
          SYMBOL=table(unique(set2eg[,c(set_id,"symbol")])[,1]))
      class(set.ng) <- "array"
      z.score <- z.score * sqrt(set.ng[names(z.score)]/
                                set.size[names(z.score)])
    }
  
    res <- cbind(data.frame(set.mean=set.mean, set.sd=set.sd,
                      set.size=set.size)[names(z.score),],
                      z.score=z.score)
    if(!is.org & collapse=="partial") names(res)[4] <- "adj.z.score"
  }

  if(universe=="local"){
    set.names <- if(collapse=="none" & !is.org)
      tapply(set2probe$probe_id, set2probe$go_id,c) else
      switch(idtype, ENTREZID=tapply(set2eg$gene_id,set2eg$go_id,c),
        SYMBOL=tapply(set2eg$symbol,set2eg$go_id,c))

    go.id <- unique(set.data$go_id)

    message("Checking parent/child relationships in GO...")
    ## parent info
    p0 <- rbind(toTable(GOBPPARENTS),
              toTable(GOMFPARENTS),
              toTable(GOCCPARENTS))
    ## only parents annotated on chip ##
    parents <-  p0[p0[,1] %in% go.id & p0[,2] %in% go.id,]
    names(parents)[2] <- "go_parent"
    ## proper subset ##
    in.subset <- apply(parents,1,function(x)
                     all(set.names[[x[1]]] %in% set.names[[x[2]]]))

    parents <- cbind(parents,
        data.frame(set.mean=set.mean[parents[,1]], #v.stat if setstat="mean"
            set.sd=set.sd[parents[,1]], #v.stat^2 if setstat="var"
            set.size=set.size[parents[,1]], #v.nset
            parent.mean=set.mean[parents[,2]], #parent.means | v.means
            parent.sd=set.sd[parents[,2]], #parent.sds | v.sds
            parent.size=set.size[parents[,2]]))[in.subset,]
                 #parent.np | n.par | p.np | v.np

    message("Computing enrichment ..." )
    if( setstat == "mean" ){
      den <- parents$parent.sd*fact(parents$parent.size,
                                   parents$set.size)
      den.globe <- sigma.globe*fact(G, parents$set.size) 
      z1 <- (parents$set.mean - parents$parent.mean)/den   # local z score
      z0 <- (parents$set.mean - mu.globe )/den.globe   # global z score
      res.full <- data.frame(parents, local.zscore=z1, global.z.score=z0 )
    }

    if( setstat=="var" ){
      E.set <- do.call(rbind,tapply(set.data$gscores,
                                  set.data$go_id, fn_getE.Globe))
      colnames(E.set) <- paste("M",1:5,sep="")
      E.par <- E.set[parents$go_parent,] ## EE
    ## remove some iffy cases
      ok <- (parents$set.size > 3) & (parents$parent.sd > 0) &
        (parents$set.size < parents$parent.size-2 )  ## a bit of room
    ## denominator, local scoring
      den <- apply(X = cbind(parents$set.size, parents$parent.size, E.par)[ok,],
                   MARGIN = 1, FUN = function(x) sigma.fun(m=x[1],E=x[3:7],G=x[2]))
      z1 <- (parents$set.sd^2 - parents$parent.sd^2)[ok]/den  # local z score

    ## denominator, global scoring
      den.globe <- sapply(X = parents$set.size[ok], FUN = sigma.fun,
                          E = E.globe, G = G)
      z0 <- (parents$set.sd^2 - sigma.globe^2)[ok]/den.globe # global z score
    ## Return full results, use local.max() to pull out "best" local.zscore ##
      res.full <- data.frame(parents[ok,], local.zscore=z1,
                                 global.zscore=z0 )
    }
    res <- local.max(res.full)
  }

  aux <- list(set.data = set.data, globe = globe)
  if(universe=="local") aux$res.full <- res.full
  
  if (!annotate) {
    main <- res
  }
  if (annotate) {
    message("Labeling output ...")
    fn_loadSetLibraries( sets=sets )
    if (sets == "GO") {
      gterms <- toTable(GOTERM)
      main <- data.frame(gterms[match(rownames(res),gterms$go_id),
                                c("Term","Ontology")],res)
      rownames(main) <- rownames(res)
    }
  }
  out <- list(setscores = main, aux = aux, call = match.call())
  return(out)
}
