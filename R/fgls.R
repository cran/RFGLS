##Function written by Xiang Li (last update 5/18/11),  Saonli Basu (last update 7/18/12), Rob Kirkpatrick (last update 11/20/12)
###############################################################################
fgls <- function(fixed, data=parent.frame(), tlist=tlist, sizelist=sizelist, 
			 sizeLab="OOPP",Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,indobs=TRUE, get.hessian=FALSE, #RMK 9/8/11
			 vmat=NULL, subset=NULL, weights=NULL, na.action=NULL) {
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    temp <- c("", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]

    temp.fixed <- fixed

    m$formula <- temp.fixed
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    id <-data$ID
    Terms <- terms(fixed)
    X <- model.matrix(Terms, m)
    dimnames(X) <- NULL
    Y <- model.extract(m, "response")
    Y <- as.vector(Y)
    n <- length(Y)
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)

    logfun <- function(theta, X, Y, center,tlist=tlist,sizelist=sizelist,id=id,sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs) {
	if(sizeLab=="OP"){
		if(sum(Mz,Bo,Ad,Mix)==1 | (sum(Mz+Bo)==2 & sum(Ad+Mix)==0)) blocks=.getblocks.A1(theta,tlist,sizelist,parlength=4,indobs=indobs) else
		blocks=.getblocks.A2(theta,tlist,sizelist,parlength=5,indobs=indobs)
	}else
	if(sizeLab=="PP") blocks=.getblocks.E1(theta,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & (sum(Mz,Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & sum(Mz+Bo)==0))) blocks=.getblocks.D1(theta,tlist,sizelist,parlength=3,indobs=indobs) else
	if(sizeLab=="OO" & Mz & (sum(Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & prod(Ad,Mix,Bo)==0))) blocks=.getblocks.D2(theta,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & !Mz & Bo & sum(Ad,Mix)>=1 ) blocks=.getblocks.D3(theta,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & Mz & Bo & sum(Ad,Mix)>=1 ) blocks=.getblocks.D4(theta,tlist,sizelist,parlength=5,indobs=indobs) else
	if(sizeLab=="OPP"){
		if(sum(Mz,Bo,Ad)==1 | (sum(Mz,Bo)==2 & !Ad)) blocks= .getblocks.B1(theta,tlist,sizelist,parlength=7,indobs=indobs) else
		blocks=.getblocks.B2(theta,tlist,sizelist,parlength=9,indobs=indobs)
	}else
	if(sizeLab=="OOPP" & (sum(Mz,Bo,Ad)==1 & !Mix)) blocks=.getblocks.Z1(theta,tlist,sizelist,parlength=8,indobs=indobs) else
	if(sizeLab=="OOPP" & ((sum(Mz,Bo,Ad)==0 & Mix) | (sum(Mz,Bo)==0 & sum(Ad,Mix)==2))) blocks=.getblocks.Z2(theta,tlist,sizelist,parlength=10,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)==0)) blocks=.getblocks.Z3(theta,tlist,sizelist,parlength=9,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==1 & sum(Ad,Mix)>=1)) blocks=.getblocks.Z4(theta,tlist,sizelist,parlength=11,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)>=1)) blocks=.getblocks.Z5(theta,tlist,sizelist,parlength=12,indobs=indobs)


#print(blocks$itheta)
#tkmat = bdiag(blocks)
 tkmat=bdsmatrix(sizelist,blocks$blocks)     	
list.vmat<-listbdsmatrix(tkmat,diag=T,id=F)

vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)
##vmat1<-vmat1+Diagonal(nrow(tkmat),diag(tkmat))
vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
vmat.Inv<-forceSymmetric(vmat.Inv)
gkmat<-as(chol(vmat.Inv),"sparseMatrix")


Lambda<-1/diag(gkmat)
NewLambda<-Diagonal(length(Lambda),Lambda)

  
newz <-gkmat%*%as.matrix(X)
newy <- gkmat%*%Y
lvd<-sum(log(Lambda))
lfit <- lm(newy[,1]~0+as.matrix(newz),subset=subset,weights=weights,na.action=na.action)
       n <- length(newy[,1])
        loglik <- sum(lfit$residuals^2)/2 + lvd   

#print(c(blocks$itheta,loglik))   
  return(loglik)	
}# end of logfunc





    

################################################
### optimze step ###############################
     center <- log(mean((Y-mean(Y))^2)); inivar <- log(var(Y,na.rm=T))
    if(is.null(vmat)){
	########## get the length of parameters to be estimated #################
	if(sizeLab=="OP"){
		if(sum(Mz,Bo,Ad,Mix)==1 | (sum(Mz+Bo)==2 & sum(Ad+Mix)==0)) pars=c(rep(-1,1),rep(inivar,2+indobs)) else
		pars=c(rep(-1,2),rep(inivar,2+indobs))
	}else
	if(sizeLab=="PP") pars=c(rep(-1,1),rep(inivar,2+indobs)) else
	if(sizeLab=="OO" & (sum(Mz,Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & sum(Mz+Bo)==0))) pars=c(rep(-1,1),rep(inivar,1+indobs)) else
	if(sizeLab=="OO" & Mz & (sum(Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & prod(Ad,Mix,Bo)==0))) pars=c(rep(-1,2),rep(inivar,1+indobs)) else
	if(sizeLab=="OO" & !Mz & Bo & sum(Ad,Mix)>=1 ) pars=c(rep(-1,2),rep(inivar,1+indobs)) else
	if(sizeLab=="OO" & Mz & Bo & sum(Ad,Mix)>=1 ) pars=c(rep(-1,3),rep(inivar,1+indobs)) else #PROBLEM! Changed 2+indobs to 1+indobs (RMK 9/8/11)
	if(sizeLab=="OPP"){
		if(sum(Mz,Bo,Ad)==1 | (sum(Mz,Bo)==2 & !Ad)) pars=c(rep(-1,3),rep(inivar,3+indobs)) else
		pars=c(rep(-1,5),rep(inivar,3+indobs))
	}else
	if(sizeLab=="OOPP" & (sum(Mz,Bo,Ad)==1 & !Mix)) pars=c(rep(-1,4),rep(inivar,3+indobs)) else
	if(sizeLab=="OOPP" & ((sum(Mz,Bo,Ad)==0 & Mix) | (sum(Mz,Bo)==0 & sum(Ad,Mix)==2))) pars=c(rep(-1,6),rep(inivar,3+indobs)) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)==0)) pars=c(rep(-1,5),rep(inivar,3+indobs)) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==1 & sum(Ad,Mix)>=1)) pars=c(rep(-1,7),rep(inivar,3+indobs)) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)>=1)) pars=c(rep(-1,8),rep(inivar,3+indobs))
	########################### Optimization #################################	
      nfit <- optim(par=pars,fn=logfun,hessian=get.hessian, #RMK 9/8/11
                      X=X, Y=Y,
                      center=center,tlist=tlist,sizelist=sizelist,id=id,sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs)
 iter <- nfit$counts[1]
########################### construct tkmat matrix from blocks ############
	if(sizeLab=="OP"){
		if(sum(Mz,Bo,Ad,Mix)==1 | (sum(Mz+Bo)==2 & sum(Ad+Mix)==0)) blocks=.getblocks.A1(nfit$par,tlist,sizelist,parlength=4,indobs=indobs) else
		blocks=.getblocks.A2(nfit$par,tlist,sizelist,parlength=5,indobs=indobs)
	}else
	if(sizeLab=="PP") blocks=.getblocks.E1(nfit$par,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & (sum(Mz,Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & sum(Mz+Bo)==0))) blocks=.getblocks.D1(nfit$par,tlist,sizelist,parlength=3,indobs=indobs) else
	if(sizeLab=="OO" & Mz & (sum(Bo,Ad,Mix)==1 | (sum(Ad+Mix)==2 & prod(Ad,Mix,Bo)==0))) blocks=.getblocks.D2(nfit$par,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & !Mz & Bo & sum(Ad,Mix)>=1 ) blocks=.getblocks.D3(nfit$par,tlist,sizelist,parlength=4,indobs=indobs) else
	if(sizeLab=="OO" & Mz & Bo & sum(Ad,Mix)>=1 ) blocks=.getblocks.D4(nfit$par,tlist,sizelist,parlength=5,indobs=indobs) else #indobs=indobs?
	if(sizeLab=="OPP"){
		if(sum(Mz,Bo,Ad)==1 | (sum(Mz,Bo)==2 & !Ad)) blocks= .getblocks.B1(nfit$par,tlist,sizelist,parlength=7,indobs=indobs) else
		blocks=.getblocks.B2(nfit$par,tlist,sizelist,parlength=9,indobs=indobs)
	}else
	if(sizeLab=="OOPP" & (sum(Mz,Bo,Ad)==1 & !Mix)) blocks=.getblocks.Z1(nfit$par,tlist,sizelist,parlength=8,indobs=indobs) else
	if(sizeLab=="OOPP" & ((sum(Mz,Bo,Ad)==0 & Mix) | (sum(Mz,Bo)==0 & sum(Ad,Mix)==2))) blocks=.getblocks.Z2(nfit$par,tlist,sizelist,parlength=10,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)==0)) blocks=.getblocks.Z3(nfit$par,tlist,sizelist,parlength=9,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==1 & sum(Ad,Mix)>=1)) blocks=.getblocks.Z4(nfit$par,tlist,sizelist,parlength=11,indobs=indobs) else
	if(sizeLab=="OOPP" & (sum(Mz,Bo)==2 & sum(Ad,Mix)>=1)) blocks=.getblocks.Z5(nfit$par,tlist,sizelist,parlength=12,indobs=indobs)
    tkmat <- bdsmatrix(sizelist,blocks=blocks$blocks)
       estimates <- blocks$itheta; names(estimates) <- blocks$parname[1:length(blocks$itheta)] #RMK 9-12-11   
    print(estimates)
    if(get.hessian==T){
      hessian.out <- nfit$hessian
      dimnames(hessian.out) <- list(names(estimates)[1:nrow(hessian.out)], names(estimates)[1:nrow(hessian.out)])   
       print(hessian.out)
    } else{hessian.out <- NULL}
#	rm(nfit);gc()
    }
    else{                     #Saonli's suggested change, 11/19/12
      if(is.character(vmat)){
        vmat.temp <- read.csv(vmat,header=T,colClasses="numeric")
        tkmat <- bdsmatrix(sizelist,vmat.temp[,1])
        rm(vmat.temp)
        estimates <- NA
        hessian.out <- NULL
      }
      else{
        tkmat <- vmat
        estimates <- NA
        hessian.out <- NULL
    }}

### end optimize step ########################################################
##############################################################################

list.vmat<-listbdsmatrix(tkmat,diag=T,id=F)
vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)

##vmat1<-vmat1+Diagonal(nrow(tkmat),diag(tkmat))
vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
vmat.Inv<-forceSymmetric(vmat.Inv)
gkmat<-as(chol(vmat.Inv),"sparseMatrix")


Lambda<-1/diag(gkmat)
NewLambda<-Diagonal(length(Lambda),Lambda)

  
newz <-gkmat%*%as.matrix(X)
newy <- gkmat%*%Y
lvd<-sum(log(Lambda))
lfit <- lm(newy[,1]~0+as.matrix(newz),subset=subset,weights=weights,na.action=na.action)
       n <- length(newy[,1])
        loglik <- sum(lfit$residuals^2)/2 + lvd   

names(lfit$coefficients) <- dimnames(X)[[2]]
    ls <- summary(lfit)
    resid.var <- mean(lfit$residuals^2)   #differs from ls$sigma, division by N

    fitted <- c(X %*% lfit$coef)  #fitted, on the original scale
    residuals <- Y - fitted
    #evresiduals <- resid.var*lfit$residuals
    if(is.null(vmat)){
	 iter=iter
    } else {
       iter=0
    }

    fcoef <- lfit$coef
    call$fixed <- fixed
    fit <- list(coefficients=list(fixed=fcoef),
		    estimates = estimates,
		    sigma = tkmat,
                variance= ls$cov.unscaled * ls$sigma^2,
                ctable = ls$coefficients,
                residuals= residuals,
		    #evresiduals=evresiduals,
                fitted.values= fitted,
                #effects=lfit$effects,
                #rank = lfit$rank,
                #assign=lfit$assign,
                df.residual = lfit$df.residual,
                loglik = loglik,
                iter = iter, n=n,
                call = call,
                #method='ML',
                hessian = hessian.out)
    dimnames(fit$ctable) <- dimnames(ls$coef)
    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    oldClass(fit) <- c('fgls')
    return(fit)
}

