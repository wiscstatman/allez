
## Order design matrix ##
## 1. Column order by GO sum
## 2. Row order by GO category, gene reduction value
## 3. Until remaining gene scores sum to 0 (or until at most maxSets sets)

ordMat <- function(aMat,allez.out,maxSets){
  rind <- cind <- character(0)
  zc <- grep("z.score",colnames(allez.out$setscores))
  for(i in 1:ncol(aMat)){
    mat <- if(i==1) aMat*c(allez.out$aux$globe) else
           aMat[-match(rind,rownames(aMat)),
                -match(cind,colnames(aMat)),drop=FALSE]*
           allez.out$aux$globe[-match(rind,names(allez.out$aux$globe))]
    s <- apply(mat,2,sum)
    ok <- ( any(s>0) & i <= maxSets)
    if(ok){    
        ## smax <- which.max(s) ## first set with highest sum
        ## break ties using z.score ##
        smax <- order(s,allez.out$setscores[names(s),zc],decreasing=TRUE)[1]
        cind <- c(cind,colnames(mat)[smax])
        rord <- order(mat[,smax],decreasing=TRUE)
        rind <- c(rind,rownames(mat)[rord][mat[rord,smax]>0])
    } else break
  }
  apply(aMat[rind,cind,drop=FALSE],2,"*",allez.out$aux$globe[rind])
}

allezplot <- function(aOrd, allez.out,
              glab=c("none","gene_id","symbol"),
              slab=c("none","z.score","set.means"),
                      ...){
  
  require(GO.db)

  goterm <- toTable(GOTERM)[,c("go_id","Term")]
  names(goterm) <- c("id","term")
  allterm <- rbind(goterm)
  term <- allterm[match(colnames(aOrd),allterm$id),]
  
  nsets <- ncol(aOrd)
  set.names <- colnames(aOrd)
  ngenes <- nrow(aOrd)
  gene.names <- rownames(aOrd)
  

  glab <- match.arg(glab)
  slab <- match.arg(slab)

  # p <- par(no.readonly=TRUE)
  # rightCol <- 3 ## fix column call for symbols
  # xlabs <- switch(glab,
  #         "none" = 1:ngenes,
  #         "gene_id" = gene.names,
  #         "symbol" = allez.out$aux$set.data[
  #           match(gene.names,allez.out$aux$set.data[,rightcol]),"symbol"])
  # zlabs <- if(slab=="none") 1:nsets else
  #          round(allez.out$setscore[set.names,slab],2)

  ## Space for text in xlim ##
  xpos <- apply(aOrd,2,function(x)
       (1:length(x))[x>0][which.max((1:length(x))[x>0])])
  # mwidth1 <- strwidth(xlabs,units="inches")
  # mwidth2 <- strwidth(term$id,units="inches")
  # mwidth4 <- strwidth(zlabs,units="inches")
  # 
  # ## Check twidth <= par("pin")[1] ##
  # twidth <- strwidth(term$term,units="inches")+0.5*par("cin")[1]
  # pwidth <- (xpos+1)*par("cin")[2]+twidth ## plot width, inches ##
  # mspace <- p$mgp[2]*par("cin")[2]+0.5*p$mai[4]
  # 
  # par(pin=c(if(glab=="none") p$pin[1] else max(pwidth),p$pin[2]),
  #     mai=c(if(glab=="none") p$mai[1] else max(mwidth1)+mspace,
  #       max(mwidth2)+mspace,p$mai[3],
  #       if(slab=="none") p$mai[4] else max(mwidth4)+mspace))
  # 
  # if(glab=="none" & max(twidth)>=par("pin")[1])
  #   warning(paste("Increase width of figure to at least",
  #           round(max(twidth)/diff(par("plt")[1:2]),2),
  #           "inches and re-run allezPlot"))
  # if(glab!="none" & max(pwidth)/diff(par("plt")[1:2])>=p$fin[1])
  #   warning(paste("Increase width of figure to at least",
  #           round(max(pwidth)/diff(par("plt")[1:2]),2),
  #           "inches and re-run allezPlot"))

  # image(1:nrow(aOrd),1:ncol(aOrd),aOrd[,ncol(aOrd):1,drop=FALSE],
  #       xlab="gene",ylab="", yaxt="n",
  #       xaxt=ifelse(glab=="none","s","n"),
  #       col=gray(seq(1,0,length=64)),
  #       xlim=c(0.5,max((xpos+0.5)/(1-twidth/par("pin")[1]))))
  # 
  # if(glab!="none") axis(side=1,at=1:nrow(aOrd),labels=xlabs,las=2, ...)
  # axis(side=2,at=1:ncol(aOrd),labels=colnames(aOrd)[ncol(aOrd):1],las=1, ...)
  # if(slab!="none"){
  #   axis(side=4,at=1:ncol(aOrd),labels=zlabs[ncol(aOrd):1],
  #        las=1,mgp=c(0,max(mwidth4)/par("cin")[2],0)+par("mgp"),hadj=1, ...)
  #   axis(side=4,at=ncol(aOrd)+0.5,tick=FALSE,
  #        labels=slab,las=1,mgp=c(0,max(mwidth4)/par("cin")[2],0)+par("mgp"),
  #        hadj=1,padj=0,xpd=TRUE, ...)
  # }
  # text(xpos+0.5,ncol(aOrd):1,term$term,pos=4, ...)
  # par(p)

  df <- data.frame(
    gene = rep(1:ngenes, nsets),
    gene.name = rep(gene.names, nsets),
    y = rep(nsets:1, each = ngenes),
    set.name = rep(set.names, each = ngenes),
    score = as.vector(aOrd)
  )
  
  

  ggplot(df, aes(gene, y)) + 
    geom_rect(
      aes(
        xmin = gene-1, xmax = gene, 
	#        ymin = y-0.5, ymax = y+0.5, 
	ymin = y-1, ymax = y, 
        colour = score,
        fill = score
        )
      ) + 
#    scale_y_discrete(limits = colnames(aOrd), position="left") + 
    scale_y_discrete(limits = rev(colnames(aOrd)), position="left") + 
    theme( legend.position="none", base_size=11 ) +
    annotate(geom="text", x = xpos + 2, y = (nsets:1) - (1/2), label = term$term, hjust = 0, size=11/(2*.pt) )+ 
	coord_cartesian(clip = 'off') +
    scale_color_gradient(low = "#FFFFFF00", high = "#000000FF")+ 
    scale_fill_gradient(low = "#FFFFFF00", high = "#000000FF") +
    ylab("") +
    xlim(0, 1.4*max(df$gene) ) 
  
}

allezPlot <- function(allez.out,
                     n.low=5,
                     n.upp=500,
                     n.cell=2,
                     nominal.alpha=0.01,
                     gmax=20,
                     maxSets=10,
                     glab=c("none","gene_id","symbol"),
                     slab=c("none","z.score","set.means"),
                     ...){
aMat <- allezMat(allez.out,n.low,n.upp,n.cell,nominal.alpha)
aOrd <- ordMat(aMat,allez.out,maxSets)
allezplot(aOrd, allez.out,
          glab=ifelse(nrow(aOrd)<=gmax,glab,"none"),
          slab=slab, ...)
}
            
            
