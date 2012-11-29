## This function was developed based on the lme.batch function in GWAF package
## by  Qiong Yang and Ming-Huei Chen.
## Function written by Xiang Li (last update 5/18/11) and Saonli Basu (last update 7/18/12), Rob Kirkpatrick (last update 11/20/12)
###############################################################################
gls.batch <- 
function(phenfile,genfile,pedifile,outfile,covmtxfile.in=NULL,
  covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),phen,covars=NULL,med="rfgls",
	sizeLab="OOPP",Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,indobs=TRUE,
	col.names=TRUE,pediheader=FALSE,
  pedicolname=c("FAMID","ID","PID","MID","SEX"),sep.phe=" ",sep.gen=" ",sep.ped=" "){
###########################################################
  med <- "rfgls"
  if(missing(outfile) | is.null(outfile) | !is.character(outfile)){
    stop("Error: argument 'outfile' must be a character string, and has no default.")
  }
  #assign("phen", phen, envir = .GlobalEnv, inherits = T)
  print("Reading in data")
  if(is.character(genfile)){
  test.dat <- read.table(genfile,header=F,na.strings="NA",sep=sep.gen)
  test.dat <- data.frame(t(test.dat))
  } else {
    test.dat <- genfile
    rm(genfile)
  }
  if(is.character(phenfile)){phen.dat=read.table(phenfile,header=TRUE,sep=sep.phe,as.is=T)}
  else{phen.dat <- phenfile}
  names(test.dat)=paste("snp",1:dim(test.dat)[2],sep=".")
  snplist <- names(test.dat)
  if(is.character(pedifile)){pedi.dat=read.table(pedifile,header=pediheader,sep=sep.ped)[,1:length(pedicolname)]}
  else{pedi.dat <- pedifile}
  names(pedi.dat)=pedicolname
  print("Done reading in data")
  
  vmat <- NULL
  varfile <- NULL
  if(!is.null(covmtxfile.in)){
    if(is.character(covmtxfile.in)){
      if(file.exists(covmtxfile.in)==T){
        varfile <- read.csv(covmtxfile.in,header=T,colClasses="numeric")
      }
      else{varfile <- NULL}
    }
    else{vmat <- covmtxfile.in}
  }
#  if(file.exists(paste(phen,"_cov_matrix.txt",sep=""))==T){
#   varfile<-read.csv(paste(phen,"_cov_matrix.txt",sep=""),header=T,colClasses="numeric")
#   }else{
#   varfile=NULL
#   }

  test.dat$ID = pedi.dat$ID
  idlab <- "ID" #"iid"
  result <- NULL
  famid <- "FAMID" #"fid"
  famtype <- "FTYPE" #"ftype"
  sid <- "INDIV"
  if (!is.null(covars))
 	phen.dat <- na.omit(phen.dat[,c(idlab,famid,famtype,sid,phen,covars)]) else
  phen.dat <- na.omit(phen.dat[,c(idlab,famid,famtype,sid,phen)])
  test.dat <- merge(phen.dat,test.dat,by="ID",sort=F)
  
  #create famsize column
  test.dat$famsize = 1
  test.dat$famsize[test.dat$FTYPE!=6]=ave(test.dat$FAMID[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],FUN=length)
  #create unisid column, c-mz twin, b-bio-offspring, a-adopted offspring, f-father, m-mother
  test.dat$unisid=NULL
  test.dat$unisid[test.dat$INDIV==4]="f"
  test.dat$unisid[test.dat$INDIV==3]="m"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==1]="c"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==2]="c"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==1]="a"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==2]="a"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==2]="a"
  #create fam labs
  test.dat$famlab="INDPT"
  test.dat$famlab[test.dat$FTYPE!=6] = ave(test.dat$unisid[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],FUN=function(x) do.call("paste",c(data.frame(matrix(x,1,length(x))),sep="")))
  #get tlist and famsize list; tlist is the list of family labels, and famsize is the list of family sizes
  tlist = tapply(test.dat$famlab[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  names=as.character(unique(test.dat$FAMID[test.dat$famlab!="INDPT"]))
  tlist=tlist[names]
  tlist = c(tlist,rep("INDPT",sum(test.dat$famlab=="INDPT")))
  sizelist = tapply(test.dat$famsize[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  sizelist = sizelist[names]
  sizelist = c(sizelist ,rep(1,sum(test.dat$famlab=="INDPT")))
  test.dat <- rbind(test.dat[test.dat$famlab != "INDPT",], test.dat[test.dat$famlab=="INDPT",])
  id<-test.dat[,idlab]
  #vals <- vector("list",2);  names(vals) <- c("estimates","hessian")
    
    
  #if(med=='fgls'){vmat <- NULL}
  #else{
  if(is.null(vmat)){
    if(is.null(varfile)){
	    if (is.null(covars)){
        lme.out<-try(fgls(test.dat[,phen]~1,data=test.dat,tlist=tlist, sizelist=sizelist, 
		      sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs,get.hessian=F,na.action=na.omit))
        vmat<-bdsmatrix(sizelist,lme.out$sigma@blocks,dimnames=list(id,id)) 
        #vals[[1]] <- lme.out$estimates; vals[[2]] <- lme.out$hessian
        if(is.character(covmtxfile.out)){write.csv(vmat@blocks,file=covmtxfile.out,quote=F,row.names=FALSE)}
      } 
      else{
        x.covar<-as.matrix(test.dat[,covars])	
    	  lme.out<-try(fgls(test.dat[,phen]~1+x.covar,data=test.dat,tlist=tlist, sizelist=sizelist, 
		      sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs,get.hessian=F,na.action=na.omit))
        #vals[[1]] <- lme.out$estimates; vals[[2]] <- lme.out$hessian
        vmat<-bdsmatrix(sizelist,lme.out$sigma@blocks,dimnames=list(id,id))
        if(is.character(covmtxfile.out)){write.csv(vmat@blocks,file=covmtxfile.out,quote=F,row.names=FALSE)}
	  }}
    else{vmat<-bdsmatrix(sizelist,c(varfile[,1]),dimnames=list(id,id))}
	} #By this point, vmat should be defined as other than NULL; possibly, varfile might not be.

  ################# begins single snp analysis ###################
  for (i in snplist) {
    #assign("i",i,envir = .GlobalEnv,inherits=T)
    if (is.null(covars)){test2.dat <- na.omit(test.dat[,c(i,phen,idlab,famid)])}
    else{
	    test2.dat <- na.omit(test.dat[,c(i,phen,covars,idlab,famid)]) 
      x.covar<-as.matrix(test2.dat[,covars])
      #assign("x.covar", x.covar, envir = .GlobalEnv,inherits=T)
	  }
	  
    count<-table(test2.dat[,i]) #Check if current SNP is monomorphic.
	  if(length(count)==1){result<-rbind(result, c(phen,i,rep(NA,7)))}
    else{
    
    if(dim(test2.dat)[1]!=dim(test.dat)[1]){ #Check if anyone is missing the current SNP; if so, cut them out of vc matrix.
      vmat0 <- vmat[which(vmat@Dimnames[[1]] %in% test2.dat[,idlab]),which(vmat@Dimnames[[2]] %in% test2.dat[,idlab])]} 
    else{vmat0 <- vmat}  
      
    id <- test2.dat[,idlab];
    #assign("test2.dat", test2.dat, envir = .GlobalEnv,inherits=T)
    #assign("id",id,envir = .GlobalEnv,inherits=T)
    	
		if (is.null(covars))
      lme.out<-try(fgls(test2.dat[,phen]~test2.dat[,i],data=test2.dat,vmat=vmat0,tlist=tlist, sizelist=sizelist, 
			  sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs,na.action=na.omit)) else
			lme.out<-try(fgls(test2.dat[,phen]~test2.dat[,i]+x.covar,data=test2.dat,vmat=vmat0,tlist=tlist, sizelist=sizelist, 
			 sizeLab=sizeLab,Mz=Mz,Bo=Bo,Ad=Ad,Mix=Mix,indobs=indobs,na.action=na.omit))
          	tmp<-c(lme.out$ctable[2,1],lme.out$ctable[2,2],
          	lme.out$ctable[2,3],lme.out$df.residual,"additive",lme.out$ctable[2,4],med)
          	if (class(tmp)!="try-error"){ 
 			result <- rbind(result, c(phen,i,tmp))
          	} else  result<-rbind(result, c(phen,i,rep(NA,7))) #If fgls() fails for some reason.
	}
  gc()
  } #end of snplist loop

  rm(test.dat)
##  rm(cov_matrix.txt)
##  rm(blocksize.txt)
  colnames(result)<-c("phen","snp","beta","se","t-stat","df","model","pval","method")
  write.table(result, outfile, quote=F,row.names=F, col.names=col.names,sep=" ",na="",append=T)
  #if(!is.null(vals[[1]])){return(vals)}
  return()
}
#mapfile="sampdat.map";phenfile="sampdat.phe";genfile="sampdat.ped";mix=NULL;phen="pheno";sep.gen=" ";sep.phe=" ";outfile="out.out"
