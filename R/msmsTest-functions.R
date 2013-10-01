##################################################
###  PROTEOMICS DIFFERENTIAL EXPRESSION TESTS  ###
###         ON SPECTRAL COUNTS MATRIX          ###
###                by J. Gregori               ###
###               September, 2013              ###
##################################################

###  GLM Poisson test, given formula of Ha model 
###   as: y~Treat  or  y~Treat+Batch   etc.
###  against Ho model: y~1   or   y~Batch   etc.
###  Returns a vector with estimated log fold change, residual 
###   deviance, and unadjusted p-value.
################################################################
glm.pois <- function(count,form1,form0,facs,div)
{ off=log(div)
  data <- data.frame(y=round(count),facs,off,stringsAsFactors=TRUE)
  gpois1a <- glm(as.formula(form1),offset=off,data=data,
                  family=poisson(link=log))
  logFC <- coef(gpois1a)[2]*log2(exp(1))
  names(logFC) <- NULL
  gpois1 <- glm(as.formula(form0),offset=off,data=data,
                  family=poisson(link=log))
  anovaPq <- data.frame(anova(gpois1, gpois1a, test="Chisq"))
  p.val <- ifelse(anovaPq[2,4] < 0.1e-15, 1, anovaPq[2,5])
  D <- anovaPq$Deviance[2]
  c(logFC=logFC,D=D,p.val=p.val)
}

###  Model multifactor, regressio de Poisson
##############################################
msms.glm.pois <- function(msnset,form1,form0,facs=NULL,div=NULL)
{ data <- exprs(msnset)
  if(is.null(facs)) facs <- as.matrix(pData(msnset))
  if(is.null(div)) div <- rep(1,ncol(data))
  Np <- nrow(data)
  result <- data.frame(LogFC=numeric(Np),D=numeric(Np),
                       p.value= numeric(Np))
  for(i in 1:Np) 
    result[i,] <- glm.pois(data[i,],form1,form0,facs,div)
  rownames(result) <- rownames(data)
  result
}


###  Quasi-likelihood test, given formula of Ha model 
###   as: y~Treat  or  y~Treat+Batch   etc.
###  against Ho model: y~1   or   y~Batch   etc.
###  Returns a vector with estimated log fold change, residual 
###   deviance, and unadjusted p-value.
################################################################
qp.test <- function(count,form1,form0,facs,div)
{ off=log(div)
  data <- data.frame(y=round(count),facs,off)
  gquasi1a <- glm(as.formula(form1),offset=off,data=data,
                  family=quasi(link=log, variance=mu))
  logFC <- coef(gquasi1a)[2]*log2(exp(1))
  names(logFC) <- NULL
  gquasi1 <- glm(as.formula(form0),offset=off,data=data,
                  family=quasi(link=log, variance=mu))
  anovaPq <- data.frame(anova(gquasi1, gquasi1a, test="F"))
  p.val <- ifelse(anovaPq[2,4] < 0.1e-15, 1, anovaPq[2,6])
  qD <- anovaPq$Deviance[2]
  c(logFC=logFC,D=qD,p.val=p.val)
}

###  Model multifactor, test de quasi-likelihood
###################################################
msms.glm.qlll <- function(msnset,form1,form0,facs=NULL,div=NULL)
{ data <- exprs(msnset)
  if(is.null(facs)) facs <- as.matrix(pData(msnset))
  if(is.null(div)) div <- rep(1,ncol(data))
  Np <- nrow(data)
  result <- data.frame(LogFC=numeric(Np),D=numeric(Np),
                       p.value= numeric(Np))
  for(i in 1:Np) 
    result[i,] <- qp.test(data[i,],form1,form0,facs,div)
  rownames(result) <- rownames(data)
  result
}


###  Multifactor edgeR with an offset, given formula of Ha model 
###   as: y~Treat  or  y~Treat+Batch   etc.
###  against Ho model: y~1   or   y~Batch   etc.
###  Returns a vector with estimated log fold change, likelihood 
###   ratio, and unadjusted p-value.
################################################################  
msms.edgeR <- function(msnset,form1,form0,facs=NULL,div=NULL,fnm=NULL)
{ data <- data.matrix(exprs(msnset))
  if(is.null(facs)) facs <- as.matrix(pData(msnset))
  if(is.null(div)) div <- rep(1,ncol(data))
  if(is.null(fnm)) fnm <- colnames(facs)[1]

  M <- model.matrix(as.formula(substring(form1,2)),data.frame(facs))
  y <- DGEList(counts=data,group=facs[,fnm])
  y <- estimateGLMCommonDisp(y, design=M, offset=log(div))
  y <- estimateGLMTagwiseDisp(y, design=M, offset=log(div))
  fit <- glmFit(y, design= M, offset=log(div))
  Mo <- model.matrix(as.formula(substring(form0,2)),data.frame(facs))
  vc <- setdiff(colnames(M),colnames(Mo))
  res <- glmLRT(fit,coef=vc)				
  ###  Adaptar resultats a l'esquema comu
  tres <- res$table[,c(1,3,4)]
  colnames(tres) <- c("LogFC","LR","p.value")
  tres
}


###  Given the results of a statistical test, and the significance
###   level, evaluates results.
###  method: one among "BH","qval" or"none"
####################################################################
test.results <- function(test,msnset,gpf,gp1,gp2,div=NULL,
                    alpha=0.05,minSpC=2,minLFC=1,method="BH")
{ dat <- exprs(msnset)
  if(!is.factor(gpf))
    gpf <- as.factor(gpf)
  if(is.null(div)) div <- rep(1,length(gpf))
  iA <- which(levels(gpf)==gp1)[1]
  iB <- which(levels(gpf)==gp2)[1]
  if(length(iA)==0 | length(iB)==0)
    stop("Fatal Error: Wrong class names or missing data.")

  o <- order(test[,"p.value"])
  nms <- rownames(test)[o]

  ###  Mean expression
  mcpg <- t(apply(dat[o,],1,function(x) tapply(x,gpf,mean)))
  mcpg <- round(mcpg,1)

  ### 'Normalized' counts and FC by average
  msms.nc <- sweep(dat,MARGIN=2,STATS=div,FUN="/")
  avgc <- t(apply(msms.nc,1,function(v) tapply(v,gpf,mean)))
  lFC.av <- log2(avgc[,iA]/avgc[,iB])
  lFC.av <- lFC.av[o]
  
  ###  Adjust p-values  
  if(method=="BH")
  { adjp <- p.adjust(test[o,"p.value"],method="fdr")
  } else if(method=="qval") 
  { qvobj <- qvalue(test[o,"p.value"])
    adjp <- qvobj$qvalues
  } else adjp <- test[o,"p.value"]
  names(adjp) <- nms 
  
  rm <- signif(data.frame(mcpg[,c(iA,iB)],lFC.Av=lFC.av,test[nms,],adjp=adjp),4)
  
  ###  Post-test filter
  mx.mean <- apply(mcpg[,c(iA,iB)],1,max)
  flags <- (adjp < alpha) & (mx.mean >= minSpC) & 
             abs(test[nms,"LogFC"]) >= minLFC
  rm$DEP <- flags
  
  n <- sum(flags)
  nms <- nms[flags]

  conds <- c(alpha.cut=alpha,SpC.cut=minSpC,LogFC.cut=minLFC)
  if(method=="qval") conds <- c(conds,pi0=qvobj$pi0)
  list(tres=rm,conds=conds)
}


####################################################################
###  Given the LogFC and the p-values builds a table of cumulative
###    values of logFC by q.values cut-offs.
####################################################################
pval.by.fc <- function(pvals,lfc)
{ fc.cuts <- c(0,1/1024,1/4,1/2,1/1.8,1/1.5,1,1.5,1.8,2,4,1024,Inf)
  lfc.cuts <- log2(fc.cuts)
  p.cuts <- c(1,0.2,0.1,0.05,0.01,0.005,0.001,0)
  cl.nms <- c("<=0.001","<=0.005","<=0.01","<=0.05","<=0.1","<=0.2","<=1")
  nc <- length(p.cuts)
  nr <- length(fc.cuts)
  tblt <- table(cut(lfc,breaks=lfc.cuts))
  tbl <- table(cut(pvals,breaks=p.cuts))
  res <- matrix(0,nrow=nr,ncol=nc-1)
  for(i in 1:(length(fc.cuts)-1))
  { flags <- lfc > lfc.cuts[i] & lfc <= lfc.cuts[i+1] 
    tbl <- table(cut(pvals[flags],breaks=p.cuts))
    res[i,] <- as.numeric(tbl)
  }
  res <- t(apply(res,1,cumsum))
  res[nr,] <- apply(res,2,sum)
  dimnames(res) <- list(LogFC=c(names(tblt),"Tot"),p.vals=cl.nms)
  res
}


####################
###  Volcano plot
####################
res.volcanoplot <- function(tres,max.pval=0.05,min.LFC=1,maxx=3,maxy=10,
                            ylbls=20)
{ 
  if(is.null(max.pval)) p.val <- 0.05
  if(is.null(min.LFC)) min.LFC <- 1  
  if(is.null(maxx)) maxx <- max(tres$LogFC)*1.05
  if(is.null(maxy)) maxy <- max(-log10(tres$adjp))*1.05
  if(is.null(ylbls)) ylbls <- 1000
  plot(tres$LogFC,-log10(tres$adjp),pch=19,cex=0.6,xlab="log2(FC)",
       ylab="-log10(adj p.value)",ylim=c(0,maxy),xlim=c(-maxx,maxx))
  grid(); abline(v=0,lty=4,col="gray")
  abline(h=-log10(max.pval),v=c(-min.LFC,min.LFC),lty=4,col="red",lwd=1)
  pflags <- tres$adjp <= max.pval & abs(tres$LogFC) >= min.LFC
  points(tres$LogFC[pflags],-log10(tres$adjp[pflags]),pch=21,bg="blue",
         cex=0.7)
  fflags <- pflags & !tres$DEP
  points(tres$LogFC[fflags],-log10(tres$adjp[fflags]),pch=21,bg="red",
         cex=0.7)
  fflags <- pflags & -log10(tres$adjp) > ylbls
  if( sum(fflags) > 0 )
    text(tres$LogFC[fflags],-log10(tres$adjp[fflags])+maxy/25,cex=0.7,
         font=2,which(fflags),col="blue")
}
