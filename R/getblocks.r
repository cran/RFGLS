#Functions written by Xiang Li (5/18/11), Saonli Basu (7/18/12), Rob Kirkpatrick (11/20/12).
###############################################################################
#c-Mz,b-Bio,a-Adopted,f-father,m-mother
#
##OOPP,famsize=4,famtype=Mz or BO or AD,parlength=8
#
#lab=c("ccff","ccf","ccm","cmf","cc","cm","cf","mf","c","m","f","INDPT") or 
#c("bbff","bbf","bbm","bmf","bb","bm","bf","mf","b","m","f","INDPT") or c("aaff","aaf","aam","amf","aa","am","af","mf","a","m","f","INDPT")
#theta=c(cor(P,P),cor(O,m),cor(O,f),cor(O,O),var(O),var(m),var(f),var(ind))
.getblocks.Z1 <- function(theta,tlist,sizelist,parlength=8,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])
cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:4]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="ccmf" | lab=="bbmf" | lab=="aamf") {cors=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cmf" | lab=="bmf" | lab=="amf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="ccm" | lab=="bbm" | lab=="aam") {cors=c(1,itheta[c(4,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="ccf" | lab=="bbf" | lab=="aaf") {cors=c(1,itheta[c(4,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm" | lab=="am") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf" | lab=="af") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cc" | lab=="aa" | lab=="bb") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="a" | lab=="c" | lab=="b") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
   for(cols in 1:sizelist[cnt]){
       vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
   }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,P)","cor(O,m)","cor(O,f)","cor(O,O)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OOPP,famsize=4,famtype=Mixed or Ad+Mixed,parlength=10
#
#lab=c("bamf","bam","baf","bmf","amf","ba","bm","bf","am","af","mf","a","b","m","f","INDPT") or
#theta=c(cor(P,P),cor(b,m),cor(b,f),cor(a,m),cor(a,f),cor(b,a),var(O),var(m),var(f),var(ind))
.getblocks.Z2 <- function(theta,tlist,sizelist,parlength=10,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])
cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(6,2,3)],1,itheta[4:5],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:6]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[2,3],cormt[2,4],cormt[1,2])


blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="bamf" | lab=="aamf") {cors=c(1,itheta[c(6,2,3)],1,itheta[4:5],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bmf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="baf") {cors=c(1,itheta[c(6,3)],1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aaf") {cors=c(1,itheta[c(6,5)],1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bam") {cors=c(1,itheta[c(6,2)],1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="aam") {cors=c(1,itheta[c(6,4)],1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="amf") {cors=c(1,itheta[c(4,5)],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="ba" | lab=="aa") {cors=c(1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="bm") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="bf") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="am") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="af") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,P)","cor(b,m)","cor(b,f)","cor(a,m)","cor(a,f)","cor(b,a)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OOPP,famsize=4,famtype=Mz+Bo,parlen=9
#
#lab=c("ccff","ccf","ccm","cmf","cc","cm","cf","mf","c","m","f","INDPT") or 
#c("bbff","bbf","bbm","bmf","bb","bm","bf","mf","b","m","f","INDPT")
#theta=c(cor(P,P),cor(O,m),cor(O,f),cor(c,c),cor(b,b),var(O),var(m),var(f),var(ind))
.getblocks.Z3 <- function(theta,tlist,sizelist,parlength=9,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])
cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:4]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(5,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1:3,5)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])



blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="ccmf") {cors=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])
	} else
if(lab=="bbmf") {cors=c(1,itheta[c(5,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cmf" | lab=="bmf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="ccm") {cors=c(1,itheta[c(4,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="bbm") {cors=c(1,itheta[c(5,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="ccf") {cors=c(1,itheta[c(4,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bbf") {cors=c(1,itheta[c(5,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cc") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="bb") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="c" | lab=="b") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,P)","cor(O,m)","cor(O,f)","cor(c,c)","cor(b,b)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OOPP,famsize=4,famtype=Mz+Ad or Bo+Ad or Mz+Mixed or Bo+Mixed or Bo+Ad+Mixed or Mz+Ad+Mixed, i.e., Mz famtype and Bo famype are not appearing together,
#but at least Ad or Mixed is present, parlength=11
#
#lab=c("ccff","ccf","ccm","cmf","cc","cm","cf","mf","c","m","f","INDPT") or 
#c("bbff","bbf","bbm","bmf","bb","bm","bf","mf","b","m","f","INDPT") or c("aaff","aaf","aam","amf","aa","am","af","mf","a","m","f","INDPT") or
#c("bamf","bam","baf","bmf","amf","ba","bm","bf","am","af","mf","a","b","m","f","INDPT")
#theta=c(cor(m,f),cor(c/b,m),cor(c/b,f),cor(c/b,c/b),cor(a,m),cor(a,f),cor(a,a),var(O),var(m),var(f),var(ind))
.getblocks.Z4 <- function(theta,tlist,sizelist,parlength=11,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:4]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(7,5,6)],1,itheta[5:6],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1,5,6,7)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(7,2,3)],1,itheta[5:6],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1,2,3,5,6,7)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[2,3],cormt[2,4],cormt[1,2])

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="ccmf" | lab=="bbmf") {cors=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aamf") {cors=c(1,itheta[c(7,5,6)],1,itheta[5:6],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bamf") {cors=c(1,itheta[c(7,2,3)],1,itheta[5:6],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cmf" | lab=="bmf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="amf") {cors=c(1,itheta[5:6],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="ccm" | lab=="bbm") {cors=c(1,itheta[c(4,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="aam") {cors=c(1,itheta[c(7,5)],1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="bam") {cors=c(1,itheta[c(7,2)],1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="ccf" | lab=="bbf") {cors=c(1,itheta[c(4,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aaf") {cors=c(1,itheta[c(7,6)],1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="baf") {cors=c(1,itheta[c(7,3)],1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="am") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="af") {cors=c(1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cc" | lab=="bb") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="aa" | lab=="ba") {cors=c(1,itheta[7],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="a" | lab=="c" | lab=="b") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(m,f)","cor(c/b,m)","cor(c/b,f)","cor(c/b,c/b)","cor(a,m)","cor(a,f)","cor(a,a)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OOPP,famsize=4,famtype=Mz+Bo+Ad or Mz+Bo+Mixed or Mz+Bo+Ad+Mixed, i.e., Mz famtype and Bo famype are both present,
#and at least Ad or Mixed is present, parlength=12
#
#lab=c("ccff","ccf","ccm","cmf","cc","cm","cf","mf","c","m","f","INDPT") or 
#c("bbff","bbf","bbm","bmf","bb","bm","bf","mf","b","m","f","INDPT") or c("aaff","aaf","aam","amf","aa","am","af","mf","a","m","f","INDPT") or
#c("bamf","bam","baf","bmf","amf","ba","bm","bf","am","af","mf","a","b","m","f","INDPT")
#theta=c(cor(m,f),cor(c/b,m),cor(c/b,f),cor(c,c),cor(b,b),cor(a,m),cor(a,f),cor(a,a),var(O),var(m),var(f),var(ind))
.getblocks.Z5 <- function(theta,tlist,sizelist,parlength=12,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:4]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(5,2,3)],1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1,2,3,5)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(8,6,7)],1,itheta[6:7],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1,6,7,8)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[1,2])

cormt = matrix(0,4,4)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[c(8,2,3)],1,itheta[6:7],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[c(1,2,3,6,7,8)]=c(cormt[3,4],cormt[1,3],cormt[1,4],cormt[2,3],cormt[2,4],cormt[1,2])

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="ccmf") {cors=c(1,itheta[c(4,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bbmf") {cors=c(1,itheta[c(5,2,3)],1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aamf") {cors=c(1,itheta[c(8,6,7)],1,itheta[6:7],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bamf") {cors=c(1,itheta[c(8,2,3)],1,itheta[6:7],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cmf" | lab=="bmf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="amf") {cors=c(1,itheta[6:7],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="ccm") {cors=c(1,itheta[c(4,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="bbm") {cors=c(1,itheta[c(5,2)],1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="aam") {cors=c(1,itheta[c(8,6)],1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="bam") {cors=c(1,itheta[c(8,2)],1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="ccf") {cors=c(1,itheta[c(4,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bbf") {cors=c(1,itheta[c(5,3)],1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aaf") {cors=c(1,itheta[c(8,7)],1,itheta[7],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="baf") {cors=c(1,itheta[c(8,3)],1,itheta[7],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="am") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="af") {cors=c(1,itheta[6],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cc") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="bb") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="aa" | lab=="ba") {cors=c(1,itheta[8],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-2-subtracter[1]])} else
if(lab=="a" | lab=="c" | lab=="b") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(m,f)","cor(c/b,m)","cor(c/b,f)","cor(c,c)","cor(b,b)","cor(a,m)","cor(a,f)","cor(a,a)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}

#
##OP,famsize=2,famtype=Mz or BO or AD or Mz+BO,parlen=4
#
#lab=c("cf","cm","f","m","c","INDPT") or c("bf","bm","f","m","b","INDPT") or c("af","am","f","m","a","INDPT")
#theta=c(cor(P,O),var(O),var(P),var(ind))
.getblocks.A1 <- function(theta,tlist,sizelist,parlength=4,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,2)} else {parlen=parlength-1;subtracter=c(0,1)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cf" | lab=="cm" | lab=="bf" |  lab=="bm" | lab=="af" |  lab=="am") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]-1],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]-1])} else
if(lab=="f" | lab=="m") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,O)","var(O)","var(P)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OP,famsize=2,famtype=Mixed, Mz+AD or Mz+Mixed or BO+AD or BO+Mixed or AD+Mixed or Mz+BO+AD or Mz+BO+Mixed or Mz+AD+Mixed or BO+AD+Mixed or Mz+BO+AD+Mixed, parlength=5
#
#lab=c("bf","bm","af","am","f","m","b","a","INDPT") or c("cf","cm","af","am","f","m","c","a","INDPT") or c("cf","cm","bf","bm","af","am","f","m","c","b","a","INDPT")
#theta=c(cor(P,(Mz,BO)),cor(P,AD),var(O),var(P),var(ind))
.getblocks.A2 <- function(theta,tlist,sizelist,parlength=5,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,2)} else {parlen=parlength-1;subtracter=c(0,1)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]


blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cf" | lab=="cm" | lab=="bf" | lab=="bm") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]-1],itheta[parlen-subtracter[1]])} else
if(lab=="af" | lab=="am") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-subtracter[1]-1],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]-1])} else
if(lab=="f" | lab=="m") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
corm[upper.tri(corm,diag=T)]=cors
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,c/b)","cor(P,a)","var(O)","var(P)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##PP,famsize=2,famtype=any type,parlength(which includes the var(independent obs))=4
#
#lab=c("mf","f","m","INDPT") 
#theta=c(cor(P,P),var(m),var(f),var(ind))
.getblocks.E1 <- function(theta,tlist,sizelist,parlength=4,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,2)} else {parlen=parlength-1;subtracter=c(0,1)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="mf") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]-1],itheta[parlen-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-subtracter[1]-1])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
corm[upper.tri(corm,diag=T)]=cors
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,P)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OO,famsize=2,famtype=Mz or BO or AD or Mixed or AD+Mixed,parlen=3
#
#lab=c("cc","c","INDPT") or c("bb","b","INDPT") or c("aa","a","INDPT") or c("ca","a","c","INDPT") or c("ba","b","a","INDPT") or c("ca","aa","a","c","INDPT") or c("ba","aa","a","b","INDPT")
#theta=c(cor(O,O),var(O),var(ind))
.getblocks.D1 <- function(theta,tlist,sizelist,parlength=3,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,1)} else {parlen=parlength-1;subtracter=c(0,0)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cc" | lab=="bb" | lab=="aa" | lab=="ca" | lab=="ba") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
corm[upper.tri(corm,diag=T)]=cors
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(O,O)","var(O)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OO,famsize=2,famtype=Mz+BO or MZ+A or MZ+Mixed or  MZ+A+Mixed,parlength=4
#
#lab=c("cc","bb","b","c","INDPT") or c("cc","aa","c","a","INDPT")  or c("cc","ca","ba","a","c","INDPT") or c("bb","ba","aa","a","b","INDPT")
#theta=c(cor(c,c),cor(b/a, b/a),var(O),var(ind))
.getblocks.D2 <- function(theta,tlist,sizelist,parlength=4,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,1)} else {parlen=parlength-1;subtracter=c(0,0)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[2],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[2]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cc") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bb" | lab=="aa" | lab=="ba") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors

for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(c,c)","cor(b/a,b/a)","var(O)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OO,famsize=2,famtype=BO+AD or BO+Mixed or or BO+A+Mixed,parlength=4
#
#lab=c("bb","aa","a","b","INDPT")  or c("ba","aa","a","b","INDPT") or c("bb","ba","aa","a","b","INDPT")
#theta=c(cor(b,b),cor(b,a),var(O),var(ind))
.getblocks.D3 <- function(theta,tlist,sizelist,parlength=4,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,1)} else {parlen=parlength-1;subtracter=c(0,0)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[2],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[2]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="bb") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aa" | lab=="ba") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
corm[upper.tri(corm,diag=T)]=cors
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(b,b)","cor(b,a)","var(O)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OO,famsize=2,famtype=Mz+BO+Mixed or Mz+BO+AD or Mz+BO+AD+Mixed,parlength=5
#
#lab=c("cc","bb","b","c","INDPT") or c("cc","aa","c","a","INDPT") or c("bb","aa","a","b","INDPT")  or c("ca","aa","a","c","INDPT") or c("ba","aa","a","b","INDPT") or c("cc","ca","aa","a","c","INDPT") or c("bb","ba","aa","a","b","INDPT")
#theta=c(cor(c,c),cor(b,b),cor(b,a),var(O),var(ind))
.getblocks.D4 <- function(theta,tlist,sizelist,parlength=5,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,1)} else {parlen=parlength-1;subtracter=c(0,0)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1]=cormt[1,2]

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[2],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[2]=cormt[1,2]

cormt = matrix(0,2,2)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[3],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[3]=cormt[1,2]

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cc") {cors=c(1,itheta[1],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="bb") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="aa" | lab=="ba" | lab=="ca") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
corm[upper.tri(corm,diag=T)]=cors
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(c,c)","cor(b,b)","cor(b,a)","var(O)","var(ind)")
#if(indobs==F){parname <- parname[1:4]} #RMK 9-12-11
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OPP,famsize=3,famtype=Mz or BO or AD or Mz+BO,parlength=7
#
#lab=c("cmf","cm","cf","m","f","c","INDPT") or c("bmf","bm","bf","m","f","b","INDPT") or c("amf","am","af","m","f","a","INDPT") or c("cmf","cm","cf","bmf","bf","bm","m","f","b","c","INDPT")
#theta=c(cor(P,P),cor(c/b,m),cor(c/b,f),var(O),var(m),var(f),var(ind))
.getblocks.B1 <- function(theta,tlist,sizelist,parlength=7,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])

cormt = matrix(0,3,3)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:3]=c(cormt[2,3],cormt[1,2],cormt[1,3])


blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cmf" | lab=="bmf" | lab=="amf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm" | lab=="am") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf" | lab=="af") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
if(nrow(corm)>1){corm[upper.tri(corm,F)] <- corm[lower.tri(corm,F)]}
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(P,P)","cor(c/b/a,m)","cor(c/b/a,f)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}
#
##OPP,famsize=3,famtype=Mixed or Mz+AD or Mz+Mixed or BO+AD or BO+Mixed or AD+Mixed or Mz+BO+AD or MZ+BO+MIXED OR MZ+AD+MIXED OR BO+AD+MIXED OR MZ+BO+AD+MIXED,parlength=9
#
#lab=c("cmf","cm","cf","amf","am","af","m","f","a","c","INDPT") or c("bmf","bm","bf","amf","am","af","a","m","f","b","INDPT") or c("cmf","cm","cf","bmf","bf","bm","amf","am","af","m","f","a","b","c","INDPT")
#theta=c(cor(f,m),cor((c,b),m),cor((c,b),f),cor(a,m),cor(a,f),var(O),var(m),var(f),var(ind))
.getblocks.B2 <- function(theta,tlist,sizelist,parlength=9,indobs=T){
if(indobs) {parlen=parlength;subtracter=c(1,3)} else {parlen=parlength-1;subtracter=c(0,2)}
itheta=NULL
itheta[1:(parlen-subtracter[2]-1)]=(1-exp(theta[1:(parlen-subtracter[2]-1)]))/(1+exp(theta[1:(parlen-subtracter[2]-1)]))
itheta[(parlen-subtracter[2]):parlen]=exp(theta[(parlen-subtracter[2]):parlen])
cormt = matrix(0,3,3)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[2:3],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:3]=c(cormt[2,3],cormt[1,2],cormt[1,3])

cormt = matrix(0,3,3)
cormt[lower.tri(cormt,diag=T)]=c(1,itheta[4:5],1,itheta[1],1)
cormt=cormt+t(cormt)-diag(diag(cormt))
cormt=nearPD(cormt,keepDiag=T)$mat
itheta[1:3]=c(cormt[2,3],cormt[1,2],cormt[1,3])

blocks=NULL;cnt = 1
for(lab in tlist){
if(lab=="cmf" | lab=="bmf") {cors=c(1,itheta[2:3],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="amf") {cors=c(1,itheta[4:5],1,itheta[1],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="cm" | lab=="bm") {cors=c(1,itheta[2],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="cf" | lab=="bf") {cors=c(1,itheta[3],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="am") {cors=c(1,itheta[4],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-1-subtracter[1]])} else
if(lab=="af") {cors=c(1,itheta[5],1); vars=c(itheta[parlen-2-subtracter[1]],itheta[parlen-subtracter[1]])} else
if(lab=="c" | lab=="b" | lab=="a") {cors=c(1); vars=c(itheta[parlen-2-subtracter[1]])} else
if(lab=="m") {cors=c(1); vars=c(itheta[parlen-1-subtracter[1]])} else
if(lab=="f") {cors=c(1); vars=c(itheta[parlen-subtracter[1]])} else
if(lab=="INDPT") {cors=c(1);vars=c(itheta[parlen])}
corm = matrix(NA,sizelist[cnt],sizelist[cnt]); vcm = matrix(NA,sizelist[cnt],sizelist[cnt])
corm[lower.tri(corm,diag=T)]=cors
if(nrow(corm)>1){corm[upper.tri(corm,F)] <- corm[lower.tri(corm,F)]}
#corm=corm+t(corm)-diag(corm)
corm=nearPD(corm,keepDiag=T)$mat
for(rows in 1:sizelist[cnt]){
    for(cols in 1:sizelist[cnt]){
        vcm[rows,cols] = corm[rows,cols]*sqrt(vars[rows]*vars[cols])
    }
}
blocks=c(blocks,vcm[lower.tri(vcm,diag=T)])
cnt = cnt + 1
} #end of lab
parname=c("cor(f,m)","cor(c/b,m)","cor(c/b,f)","cor(a,m)","cor(a,f)","var(O)","var(m)","var(f)","var(ind)")
return(list(blocks=blocks,parname=parname,itheta=itheta))
}

