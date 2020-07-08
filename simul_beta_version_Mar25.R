
library(MASS)
library(mice)
library(mpath)
library(glmnet)
library(ncvreg)
library(grpreg)

rm(list=ls())

source("sim_aoas_data_generate.R")
source("degree_freedom.R")
beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 5))
#beta.vec = c(-1.00, -0.50, 0.25,  0.50,  1.00,  rep(0.00, 15))
pp = length(beta.vec)

corrln = 0.2; sigma = 1;



########################################
########  MCAR DATA    #################
########################################


miss_dat_mcar = sim.data(gen.size = 100, n = 100, rho = corrln, sig = sigma, beta.coef = beta.vec, prop.miss = 75, miss.type = "MCAR", corr.type = "linear")

### Check proportion of missing ###

miss_dat_mcar$missing.proportion

tab = matrix(0,ncol=15)

for(i in 1:100){ 		##for smaller number of iteration change

  tab1 = NULL
  final.dat = miss_dat_mcar$data[[i]]
  yy = final.dat[,1]
  true.index = 1:5  ##which(beta.vec!=0)
  final.dat = final.dat[,-1]   ## drop the response from the data before imputation
  dsgn.mtrx.org<- final.dat

  ### finding columns with missing data and computing number of missing in each row ###
  frac = apply(final.dat, 2, function(x){length(which(is.na(x)))})
  ind = as.numeric(which(frac != 0)); ind1 = as.numeric(which(frac == 0));


  final.data=data.frame(final.dat)
  imp.dat=mice(final.data,m=5,method="pmm",seed=500)  ## "m" may be 8 or 10; discuss w/ Samiran-da

  ## the five data sets can be obtained by complete(imp.dat,1),
  ## complete(imp.dat,2),...,complete(imp.dat,5)

  ### appending with only missing columns ###
  final.dat.stack=cbind(complete(imp.dat,1),complete(imp.dat,2)[,ind],complete(imp.dat,3)[,ind],
                        complete(imp.dat,4)[,ind],complete(imp.dat,5)[,ind])

  dsgn.mtrx<- final.dat.stack
  grp= c(seq(1,pp,by=1),rep(ind,4)) ## 
  set.seed(1)
  samp.size<- dim(dsgn.mtrx)[1]


  ### Lasso fit ###

  fit.lso.mcar <-  grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grLasso")
  lso.bic<- select(fit.lso.mcar,"EBIC")
  hat.beta.lasso= lso.bic$beta[2:(pp+1)]

  var.nonzero.cv.mcar = which(hat.beta.lasso!=0)

  (lasso.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA")))  ## true nonzero variable detection
  (lasso.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-lasso.true.index.cv.mcar)

  mse.lasso = t(beta.vec - hat.beta.lasso)%*%(beta.vec - hat.beta.lasso)
  sense.lasso = lasso.true.index.cv.mcar/ length(true.index)
  spec.lasso = 1 - lasso.false.positive.cv.mcar/(pp - length(true.index))

  tab1 = c(tab1, mse.lasso, sense.lasso, spec.lasso)

  ### SCAD fit ###

fit <-  grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grSCAD")

p.prime.mat=matrix(0,nrow=length(fit$lambda),ncol=length(fit$beta[-1,1])) ## calculate it for every beta except intercept
indctr=matrix(0,nrow=length(fit$lambda),ncol=length(fit$beta[-1,1]))

for(ii in 1:length(fit$lambda)){
for(jj in 1:length(fit$beta[-1,1])){
if(jj<=pp){                          ## since there 10 non missing variables
indctr[ii,jj]= fit$beta[jj+1,ii]  ## this is for nonmissing predictors
}
else
indctr[ii,jj]= sqrt((fit$beta[jj+1,ii]^2)*5)  ## since we have 5 data sets after imputation from "mice"

if(abs(indctr[ii,jj])<fit$lambda[ii]){
p.prime.mat[ii,jj]=fit$lambda[ii]*sign(indctr[ii,jj])
}
if(abs(indctr[ii,jj])>=fit$lambda[ii] && abs(indctr[ii,jj])< 3.7*fit$lambda[ii]){
p.prime.mat[ii,jj]=(fit$lambda[ii]*3.7*sign(indctr[ii,jj])-indctr[ii,jj])/2.7
}

}
}

rm(ii,jj)

#############################################################################
    ######## we need to subset the design matrix appropriately ##########
#############################################################################

##dsgn.mtrx.sub=list()  ## stores the sub design matrix of non-zero predictors. Being unequal we need array not matrix

df.lam.mar=numeric(0)
gbic.mar=numeric(0)

for(ii in 2:length(fit$lambda)){
non.zero.index=which(fit$beta[-1,ii]!=0)
dsgn.mtrx.sub=as.matrix(dsgn.mtrx[,non.zero.index])   ## this should be a matrix of order n by d w/ d_0
p.prime.sub <- as.vector(p.prime.mat[ii,non.zero.index])   ## this is a vector of length d_0

non.zero.index.incpt=which(fit$beta[,ii]!=0)
beta.nonzero.incpt=fit$beta[non.zero.index.incpt,ii]
beta.nonzero=beta.nonzero.incpt[-1]   ## excluding the intercept

beta.nonzero.abs=numeric(0)
for(jj in 1:length(beta.nonzero)){
ifelse(jj<=pp,beta.nonzero.abs[jj]<-abs(beta.nonzero[jj]),beta.nonzero.abs[jj]<- beta.nonzero[jj])
}

sigma.diag=abs(p.prime.sub)/beta.nonzero.abs
sigma.lam=matrix(0,nrow=length(sigma.diag),ncol=length(sigma.diag))
diag(sigma.lam) <- sigma.diag

df.lam.inv <- ginv(t(dsgn.mtrx.sub)%*%dsgn.mtrx.sub + samp.size*sigma.lam)
df.lam.mat <- dsgn.mtrx.sub%*%df.lam.inv%*%t(dsgn.mtrx.sub)

df.lam.mar[ii-1]= sum(diag(df.lam.mat))

gbic.mar[ii-1]=log(sum((yy-dsgn.mtrx.sub%*%beta.nonzero)-beta.nonzero.incpt[1])^2/samp.size) + 
df.lam.mar[ii-1]*log(samp.size)/samp.size

rm(dsgn.mtrx.sub,p.prime.sub,non.zero.index,non.zero.index.incpt,beta.nonzero,beta.nonzero.incpt,
beta.nonzero.abs,sigma.diag,sigma.lam,df.lam.inv,df.lam.mat)
}

#which.min(gbic.mar)
 hat.beta.scad<-fit$beta[-1,(which.min(gbic.mar))][1:pp]

var.nonzero.cv.mcar=which(fit$beta[-1,(which.min(gbic.mar))][1:pp]!=0)  ## since ii loop has started from 2 
grp.scad.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA"))  ## true nonzero variable detection
grp.scad.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-grp.scad.true.index.cv.mcar


  mse.scad = t(beta.vec - hat.beta.scad)%*%(beta.vec - hat.beta.scad)
  sense.scad = grp.scad.true.index.cv.mcar/ length(true.index)
  spec.scad = 1 - grp.scad.false.positive.cv.mcar/(pp - length(true.index))

  tab1 = c(tab1, mse.scad, sense.scad, spec.scad)

############ in built BIC of R-package
  scad.bic<-select(fit,"BIC")


var.nonzero.cv.mcar=which(scad.bic$beta[c(2:(pp+1))]!=0)
scadbic.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA"))  ## true nonzero variable detection
scadbic.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-scadbic.true.index.cv.mcar

hat.beta.scad.bic<-scad.bic$beta[c(2:(pp+1))]


  mse.scad.bic = t(beta.vec - hat.beta.scad.bic)%*%(beta.vec - hat.beta.scad.bic)
  sense.scad.bic = scadbic.true.index.cv.mcar/ length(true.index)
  spec.scad.bic = 1 - scadbic.false.positive.cv.mcar/(pp - length(true.index))

  tab1 = c(tab1, mse.scad.bic, sense.scad.bic, spec.scad.bic)

  ### MCP fit ###

fit.mcar.mcp <- grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grMCP")


mcp.bic<-select(fit.mcar.mcp,"BIC")


var.nonzero.cv.mcar=which(mcp.bic$beta[c(2:(pp+1))]!=0)
mcp.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA"))  ## true nonzero variable detection
mcp.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-mcp.true.index.cv.mcar

hat.beta.mcp<-mcp.bic$beta[c(2:(pp+1))]


  mse.mcp = t(beta.vec - hat.beta.mcp)%*%(beta.vec - hat.beta.mcp)
  sense.mcp = mcp.true.index.cv.mcar/ length(true.index)
  spec.mcp = 1 - mcp.false.positive.cv.mcar/(pp - length(true.index))

  tab1 = c(tab1, mse.mcp, sense.mcp, spec.mcp)

  ### Group Bridge fit ###

fit.mcar.bg <- gBridge(as.matrix(dsgn.mtrx), yy,group=grp,family="gaussian")

bg.bic<-select(fit.mcar.bg,"BIC")
hat.beta.bg<- bg.bic$beta[c(2:(pp+1))]

var.nonzero.cv.mcar=which(hat.beta.bg!=0)
bg.true.index.cv.mcar=length(which(match(true.index,var.nonzero.cv.mcar)!="NA"))  ## true nonzero variable detection
bg.false.positive.cv.mcar=length(var.nonzero.cv.mcar)-bg.true.index.cv.mcar


  mse.bg = t(beta.vec - hat.beta.bg)%*%(beta.vec - hat.beta.bg)
  sense.bg = bg.true.index.cv.mcar/ length(true.index)
  spec.bg = 1 - bg.false.positive.cv.mcar/(pp - length(true.index))

  tab1 = c(tab1, mse.bg, sense.bg, spec.bg)
  tab = rbind(tab,tab1)

}
tab = tab[-1,]
colnames(tab) = c("mse.lasso","sensitivity.lasso","specificity.lasso","mse.scad","sensitivity.scad","specificity.scad",
"mse.scad.bic","sensitivity.scad.bic","specificity.scad.bic", "mse.mcp","sensitivity.mcp","specificity.mcp",
"mse.Bge","sensitivity.Bge","specificity.Bge")
apply(tab, 2, mean)


rm(mse.lasso,sense.lasso,spec.lasso,mse.scad, sense.scad, spec.scad, mse.scad.bic, sense.scad.bic, spec.scad.bic,
                mse.mcp, sense.mcp, spec.mcp,mse.bg, sense.bg, spec.bg)


########################################
########  MAR DATA    #################
########################################


miss_dat_mar = sim.data(gen.size = 100, n = 100, rho = corrln, sig = sigma, beta.coef = beta.vec, prop.miss = 75, miss.type = "MAR", corr.type = "linear")

### Check proportion of missing ###

miss_dat_mar$missing.proportion

tab = matrix(0,ncol=15)

for(i in 1:100){ 		##for smaller number of iteration change

  tab1 = NULL
  final.dat = miss_dat_mar$data[[i]]
  yy = final.dat[,1]
  true.index = 1:5  ##which(beta.vec!=0)
  final.dat = final.dat[,-1]   ## drop the response from the data before imputation
  dsgn.mtrx.org<- final.dat

  ### finding columns with missing data and computing number of missing in each row ###
  frac = apply(final.dat, 2, function(x){length(which(is.na(x)))})
  ind = as.numeric(which(frac != 0)); ind1 = as.numeric(which(frac == 0));


  final.data=data.frame(final.dat)
  imp.dat=mice(final.data,m=5,method="pmm",seed=500)  ## "m" may be 8 or 10; discuss w/ Samiran-da

  ## the five data sets can be obtained by complete(imp.dat,1),
  ## complete(imp.dat,2),...,complete(imp.dat,5)

  ### appending with only missing columns ###
  final.dat.stack=cbind(complete(imp.dat,1),complete(imp.dat,2)[,ind],complete(imp.dat,3)[,ind],
                        complete(imp.dat,4)[,ind],complete(imp.dat,5)[,ind])

  dsgn.mtrx<- final.dat.stack
  grp= c(seq(1,pp,by=1),rep(ind,4)) ## change this 5 if number of imputation is changed according to Samiran-da
  set.seed(1)
  samp.size<- dim(dsgn.mtrx)[1]


  ### Lasso fit ###

  fit.lso.mar <-  grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grLasso")
  lso.bic<- select(fit.lso.mar,"EBIC")
  hat.beta.lasso= lso.bic$beta[2:(pp+1)]

  var.nonzero.cv.mar = which(hat.beta.lasso!=0)

  (lasso.true.index.cv.mar=length(which(match(true.index,var.nonzero.cv.mar)!="NA")))  ## true nonzero variable detection
  (lasso.false.positive.cv.mar=length(var.nonzero.cv.mar)-lasso.true.index.cv.mar)

  mse.lasso = t(beta.vec - hat.beta.lasso)%*%(beta.vec - hat.beta.lasso)
  sense.lasso = lasso.true.index.cv.mar/ length(true.index)
  spec.lasso = 1 - lasso.false.positive.cv.mar/(pp - length(true.index))

  tab1 = c(tab1, mse.lasso, sense.lasso, spec.lasso)

  ### SCAD fit ###

fit <-  grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grSCAD")

p.prime.mat=matrix(0,nrow=length(fit$lambda),ncol=length(fit$beta[-1,1])) ## calculate it for every beta except intercept
indctr=matrix(0,nrow=length(fit$lambda),ncol=length(fit$beta[-1,1]))

for(ii in 1:length(fit$lambda)){
for(jj in 1:length(fit$beta[-1,1])){
if(jj<=pp){                          ## since there 10 non missing variables
indctr[ii,jj]= fit$beta[jj+1,ii]  ## this is for nonmissing predictors
}
else
indctr[ii,jj]= sqrt((fit$beta[jj+1,ii]^2)*5)  ## since we have 5 data sets after imputation from "mice"

if(abs(indctr[ii,jj])<fit$lambda[ii]){
p.prime.mat[ii,jj]=fit$lambda[ii]*sign(indctr[ii,jj])
}
if(abs(indctr[ii,jj])>=fit$lambda[ii] && abs(indctr[ii,jj])< 3.7*fit$lambda[ii]){
p.prime.mat[ii,jj]=(fit$lambda[ii]*3.7*sign(indctr[ii,jj])-indctr[ii,jj])/2.7
}

}
}

rm(ii,jj)

#############################################################################
    ######## we need to subset the design matrix appropriately ##########
#############################################################################

##dsgn.mtrx.sub=list()  ## stores the sub design matrix of non-zero predictors. Being unequal we need array not matrix

df.lam.mar=numeric(0)
gbic.mar=numeric(0)

for(ii in 2:length(fit$lambda)){
non.zero.index=which(fit$beta[-1,ii]!=0)
dsgn.mtrx.sub=as.matrix(dsgn.mtrx[,non.zero.index])   ## this should be a matrix of order n by d w/ d_0
p.prime.sub <- as.vector(p.prime.mat[ii,non.zero.index])   ## this is a vector of length d_0

non.zero.index.incpt=which(fit$beta[,ii]!=0)
beta.nonzero.incpt=fit$beta[non.zero.index.incpt,ii]
beta.nonzero=beta.nonzero.incpt[-1]   ## excluding the intercept

beta.nonzero.abs=numeric(0)
for(jj in 1:length(beta.nonzero)){
ifelse(jj<=pp,beta.nonzero.abs[jj]<-abs(beta.nonzero[jj]),beta.nonzero.abs[jj]<- beta.nonzero[jj])
}

sigma.diag=abs(p.prime.sub)/beta.nonzero.abs
sigma.lam=matrix(0,nrow=length(sigma.diag),ncol=length(sigma.diag))
diag(sigma.lam) <- sigma.diag

df.lam.inv <- ginv(t(dsgn.mtrx.sub)%*%dsgn.mtrx.sub + samp.size*sigma.lam)
df.lam.mat <- dsgn.mtrx.sub%*%df.lam.inv%*%t(dsgn.mtrx.sub)

df.lam.mar[ii-1]= sum(diag(df.lam.mat))

gbic.mar[ii-1]=log(sum((yy-dsgn.mtrx.sub%*%beta.nonzero)-beta.nonzero.incpt[1])^2/samp.size) + 
df.lam.mar[ii-1]*log(samp.size)/samp.size

rm(dsgn.mtrx.sub,p.prime.sub,non.zero.index,non.zero.index.incpt,beta.nonzero,beta.nonzero.incpt,
beta.nonzero.abs,sigma.diag,sigma.lam,df.lam.inv,df.lam.mat)
}

#which.min(gbic.mar)
 hat.beta.scad<-fit$beta[-1,(which.min(gbic.mar))][1:pp]

var.nonzero.cv.mar=which(fit$beta[-1,(which.min(gbic.mar))][1:pp]!=0)  ## since ii loop has started from 2 
grp.scad.true.index.cv.mar=length(which(match(true.index,var.nonzero.cv.mar)!="NA"))  ## true nonzero variable detection
grp.scad.false.positive.cv.mar=length(var.nonzero.cv.mar)-grp.scad.true.index.cv.mar


  mse.scad = t(beta.vec - hat.beta.scad)%*%(beta.vec - hat.beta.scad)
  sense.scad = grp.scad.true.index.cv.mar/ length(true.index)
  spec.scad = 1 - grp.scad.false.positive.cv.mar/(pp - length(true.index))

  tab1 = c(tab1, mse.scad, sense.scad, spec.scad)

### scad fitting by in built BIC
scad.bic<-select(fit,"BIC")


var.nonzero.cv.mar=which(scad.bic$beta[c(2:(pp+1))]!=0)
scad.bic.true.index.cv.mar=length(which(match(true.index,var.nonzero.cv.mar)!="NA"))  ## true nonzero variable detection
scad.bic.false.positive.cv.mar=length(var.nonzero.cv.mar)-scad.bic.true.index.cv.mar

hat.beta.scad.bic<-scad.bic$beta[c(2:(pp+1))]


  mse.scad.bic = t(beta.vec - hat.beta.scad.bic)%*%(beta.vec - hat.beta.scad.bic)
  sense.scad.bic = scad.bic.true.index.cv.mar/ length(true.index)
  spec.scad.bic = 1 - scad.bic.false.positive.cv.mar/(pp - length(true.index))

  tab1 = c(tab1, mse.scad.bic, sense.scad.bic, spec.scad.bic)

  ### MCP fit ###

fit.mar.mcp <- grpreg(as.matrix(dsgn.mtrx), yy,group=grp,penalty="grMCP")


mcp.bic<-select(fit.mar.mcp,"BIC")


var.nonzero.cv.mar=which(mcp.bic$beta[c(2:(pp+1))]!=0)
mcp.true.index.cv.mar=length(which(match(true.index,var.nonzero.cv.mar)!="NA"))  ## true nonzero variable detection
mcp.false.positive.cv.mar=length(var.nonzero.cv.mar)-mcp.true.index.cv.mar

hat.beta.mcp<-mcp.bic$beta[c(2:(pp+1))]


  mse.mcp = t(beta.vec - hat.beta.mcp)%*%(beta.vec - hat.beta.mcp)
  sense.mcp = mcp.true.index.cv.mar/ length(true.index)
  spec.mcp = 1 - mcp.false.positive.cv.mar/(pp - length(true.index))

  tab1 = c(tab1, mse.mcp, sense.mcp, spec.mcp)

  ### Group Bridge fit ###

fit.mar.bg <- gBridge(as.matrix(dsgn.mtrx), yy,group=grp,family="gaussian")

bg.bic<-select(fit.mar.bg,"BIC")
hat.beta.bg<- bg.bic$beta[c(2:(pp+1))]

var.nonzero.cv.mar=which(hat.beta.bg!=0)
bg.true.index.cv.mar=length(which(match(true.index,var.nonzero.cv.mar)!="NA"))  ## true nonzero variable detection
bg.false.positive.cv.mar=length(var.nonzero.cv.mar)-bg.true.index.cv.mar


  mse.bg = t(beta.vec - hat.beta.bg)%*%(beta.vec - hat.beta.bg)
  sense.bg = bg.true.index.cv.mar/ length(true.index)
  spec.bg = 1 - bg.false.positive.cv.mar/(pp - length(true.index))

  tab1 = c(tab1, mse.bg, sense.bg, spec.bg)
  tab = rbind(tab,tab1)

}
tab = tab[-1,]
colnames(tab) = c("mse.lasso","sensitivity.lasso","specificity.lasso","mse.scad","sensitivity.scad","specificity.scad",
"mse.scad.bic","sensitivity.scad.bic","specificity.scad.bic","mse.mcp","sensitivity.mcp","specificity.mcp",
"mse.Bge","sensitivity.Bge","specificity.Bge")
apply(tab, 2, mean)

rm(mse.lasso,sense.lasso,spec.lasso,mse.scad, sense.scad, spec.scad, mse.scad.bic, sense.scad.bic, spec.scad.bic,
                mse.mcp, sense.mcp, spec.mcp,mse.bg, sense.bg, spec.bg)

