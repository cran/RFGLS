#Adapted from Xiang Li's gls.batch(), by Rob Kirkpatrick.
gls.batch.get <-
  function(phenfile,genfile,pedifile,outfile,covmtxfile.in=NULL,
    covmtxfile.out=paste(phen,"_cov_matrix.txt",sep=""),phen,covars=NULL,med="rfgls",
    sizeLab="OOPP",Mz=TRUE,Bo=TRUE,Ad=TRUE,Mix=TRUE,indobs=TRUE,
    col.names=TRUE,pediheader=FALSE,
    pedicolname=c("FAMID","ID","PID","MID","SEX"),sep.phe=" ",sep.gen=" ",sep.ped=" "){
    ###########################################################
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
  
    getout <- list(test.dat, tlist, sizelist)
    names(getout) <- c("test.dat", "tlist", "sizelist")
    return(getout)
}