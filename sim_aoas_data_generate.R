sim.data = function(gen.size = 100, n = 200, rho = 0.5, sig = 1, beta.coef, prop.miss = 50, miss.type = c("Non", "MCAR", "MAR"), corr.type = c("linear", "auto"))
{
	set.seed(100)
  p = length(beta.coef)
	mu.vec = rep(0,p)
	times = 1:p
	n.miss = n*(prop.miss/100)
	
	if(corr.type == "linear"){
	  H = matrix(1, p, p)
	  diag(H) = 0
	}
	if(corr.type == "auto"){
	  H = abs(outer(times, times, "-"))
	}

	sig.mat = rho^H

	dsgn.mtrx = mvrnorm(n, mu.vec, sig.mat)
	
	gen.var = c(4, 5, 8)  ## complete variable set
	
	if(prop.miss == 50){
	  if(miss.type == "MCAR"){
	    if(p ==10){prob.miss.mcar = 0.9}
	    if(p ==40){prob.miss.mcar = 0.9}
	  }
	  if(miss.type == "MAR"){
	    if(p ==10){prob.miss.mar = 0.96352}
	    if(p ==40){prob.miss.mar = 0.9961}
	  }
	}
	if(prop.miss == 75){
	  if(miss.type == "MCAR"){
	    if(p ==10){prob.miss.mcar = 0.9}
	    if(p ==40){prob.miss.mcar = 0.9}
	  }
	  if(miss.type == "MAR"){
	    if(p ==10){prob.miss.mar = 0.76308}
	    if(p ==40){prob.miss.mar = 0.97451}
	  }
	}
	
	if(p == 10){
	  if(prop.miss == 50){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = 9
	  }
	  if(prop.miss == 75){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = c(9, 10)
	  }
	}
	if(p == 20){
	  if(prop.miss == 50){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = c(9, 10, 19, 20)
	  }
	  if(prop.miss == 75){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = c(9, 10, 14, 15, 19, 20)
	  }
	}
	if(p == 30){
	  if(prop.miss == 50){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = c(9, 10, 14, 15, 19, 20, 29, 30)
	  }
	  if(prop.miss == 75){
	    miss.var.mar = c(1, 2, 6, 7)  ## missing indexes
	    miss.var.mcar = c(9, 10, 14, 15, 19, 20, 24, 25, 28, 29, 30)
	  }
	}
	
	dat = list(); prop = 0;

	if(miss.type == "Non"){
	  for(i in 1:gen.size){
	    res = dsgn.mtrx%*%beta.vec + rnorm(n, 0, sig) ## response
	    final.dat=cbind(res,dsgn.mtrx)   ## full final data
	    colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
	    dat[[i]] = final.dat
	  }
	  return(data = dat)
	}
  if(miss.type == "MCAR"){
    for(i in 1:gen.size){
      res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
      prob.miss.mcar = 0.9
      miss.index = sort(sample(1:n, n.miss))
	    comp.index = setdiff(1:n, miss.index)
	    miss.var = c(miss.var.mcar, miss.var.mar)
	    comp.var = setdiff(1:p, miss.var)
      miss.mat = matrix(0, nrow = n.miss, ncol = p)
      miss.mat[,comp.var] = 1 
	    a1 = matrix(0,ncol = length(miss.var))
	    ii = 1
      while(ii <= n.miss){
	      samp1 = rbinom(length(miss.var), 1, prob.miss.mcar)
	      if(sum(samp1) != length(miss.var)){
	        a1 = rbind(a1, samp1)
	        ii = ii + 1
	      }
      }
	    a1 = a1[-1,]
	    miss.mat[,miss.var] = a1  ##generating approx 50% missing
	    rm(list = c("ii","a1"))

	    dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
	    dsgn.mtrx.miss[comp.index,] = dsgn.mtrx[comp.index,]

	    for(ii in 1:n.miss){
	      for(jj in 1:p){
	        ifelse(miss.mat[ii,jj]==0,dsgn.mtrx.miss[miss.index[ii],jj]<-NA,dsgn.mtrx.miss[miss.index[ii],jj]<-dsgn.mtrx[miss.index[ii],jj])
	      }
	    }
	    #rm(list=c("ii","miss.mat"))

	    final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
	    colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
	    prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations

	    dat[[i]] = final.dat
  	}
    output = list(missing.proportion = prop/gen.size, data = dat)
    return(output)
	}
  if(miss.type == "MAR"){
    if(p == 10){
    	i = 1
      while(i <= gen.size){
        res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
        comp.var = setdiff(1:p, c(miss.var.mar,miss.var.mcar))
	      x = dsgn.mtrx[,gen.var]  ## design matrix for MAR
        pr1 = 1/(1 + exp(-4.5 + x[,2] + 2*x[,3])) ##probability for missing for index 1
        pr2 = 1/(1 + exp(-3.6 - x[,1] - 2*x[,2] - x[,3] + res)) ##probability for missing for index 2
        pr3 = 1/(1 + exp(-2.25 + x[,1] - 1.5*x[,2] + 0.5*x[,3] + res))  ##probability for missing for index 6
        pr4 = 1/(1 + exp(-2 + 0.25*x[,1] + x[,2] + 0.5*x[,3]))  ##probability for missing for index 7

	      prob = cbind(pr1, pr2, pr3, pr4)
	      rm(list= c("pr1", "pr2", "pr3", "pr4"))

	      miss.mat = matrix(0, nrow = n, ncol = length(miss.var.mar))
        for(ii in 1:n){
          for(jj in 1:length(miss.var.mar)){
            miss.mat[ii,jj] = rbinom(1, 1, prob[ii,jj]) ##generating approx 50% missing
            #print(prob[ii,jj])
          }
        }
        rm(ii)
        miss.mat[miss.mat==0] <- NA
        miss.mat[which(miss.mat==1)] <- dsgn.mtrx[,miss.var.mar][which(miss.mat==1)]
        
        n.miss.mar = which(apply(is.na(miss.mat) , 1, any)==TRUE)
        miss1 = length(n.miss.mar)
        if(miss1 < n.miss){
        	n.miss.mcar = n.miss - miss1
       
        	miss.mat1 = matrix(1, nrow = n)
        	indx1 = sample(setdiff(1:n, n.miss.mar), n.miss.mcar)
        	miss.mat1[indx1,] = NA
        
        	miss.mat1[which(miss.mat1==1)] <- dsgn.mtrx[,miss.var.mcar][which(miss.mat1==1)]
        	#length(which(apply(is.na(miss.mat1) , 1, any)==TRUE))
      
        	dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
        	dsgn.mtrx.miss[,comp.var] = dsgn.mtrx[,comp.var]
        	dsgn.mtrx.miss[,miss.var.mar] = miss.mat
        	dsgn.mtrx.miss[,miss.var.mcar] = miss.mat1
  
        	final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
        	colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
        	prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations
  
        	dat[[i]] = final.dat
        	i = i + 1
        }
      }
      output = list(missing.proportion = prop/gen.size, data = dat)
      return(output)
    }
    if(p == 20){
    	i = 1
      while(i <= gen.size){
        res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
        comp.var = setdiff(1:p, c(miss.var.mar))
        x = dsgn.mtrx[,gen.var]  ## design matrix for MAR
        pr1 = 1/(1 + exp(-4.5 + x[,2] + 2*x[,3])) ##probability for missing for index 1
        pr2 = 1/(1 + exp(-3.6 - x[,1] - 2*x[,2] - x[,3] + res)) ##probability for missing for index 2
        pr3 = 1/(1 + exp(-2.25 + x[,1] - 1.5*x[,2] + 0.5*x[,3] + res))  ##probability for missing for index 6
        pr4 = 1/(1 + exp(-2 + 0.25*x[,1] + x[,2] + 0.5*x[,3]))  ##probability for missing for index 7
        
        prob = cbind(pr1, pr2, pr3, pr4)
        rm(list= c("pr1", "pr2", "pr3", "pr4"))
        
        miss.mat = matrix(0, nrow = n, ncol = length(miss.var.mar))
        for(ii in 1:n){
          for(jj in 1:length(miss.var.mar)){
            miss.mat[ii,jj] = rbinom(1, 1, prob[ii,jj]) ##generating approx 50% missing
            #print(prob[ii,jj])
          }
        }
        
        rm(ii)
        miss.mat[miss.mat==0] <- NA
        miss.mat[which(miss.mat==1)] <- dsgn.mtrx[,miss.var.mar][which(miss.mat==1)]
        
        n.miss.mar = which(apply(is.na(miss.mat) , 1, any)==TRUE)
        miss1 = length(n.miss.mar)
        if(miss1 < n.miss){
       		n.miss.mcar = n.miss - miss1
        
        	a1 = matrix(0,ncol = length(miss.var.mcar))
        	ii = 1
        	while(ii <= n.miss.mcar){
          		samp1 = rbinom(length(miss.var.mcar), 1, 0.5)
          	if(sum(samp1) != length(miss.var.mcar)){
            	a1 = rbind(a1, samp1)
            	ii = ii + 1
          		}
        	}
        	a1 = a1[-1,]
        	miss.mat1 = matrix(1, nrow = n, ncol = length(miss.var.mcar))
        	indx1 = sample(setdiff(1:n, n.miss.mar), n.miss.mcar)
        	miss.mat1[indx1,] = a1
        
        	miss.mat1[miss.mat1==0] <- NA
        	miss.mat1[which(miss.mat1==1)] <- dsgn.mtrx[,miss.var.mcar][which(miss.mat1==1)]
        	#length(which(apply(is.na(miss.mat1) , 1, any)==TRUE))
        	
        	dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
        	dsgn.mtrx.miss[,comp.var] = dsgn.mtrx[,comp.var]
        	dsgn.mtrx.miss[,miss.var.mar] = miss.mat
        	dsgn.mtrx.miss[,miss.var.mcar] = miss.mat1
        
        	final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
        	colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
        	prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations
        
        	dat[[i]] = final.dat
        	i =i + 1
        }
      }
      output = list(missing.proportion = prop/gen.size, data = dat)
      return(output)
    }
    if(p == 30){
      i = 1
      while(i <= gen.size){
        res = dsgn.mtrx%*%beta.coef + rnorm(n, 0, sig) ##response
        comp.var = setdiff(1:p, c(miss.var.mar))
        x = dsgn.mtrx[,gen.var]  ## design matrix for MAR
        pr1 = 1/(1 + exp(-4.5 + x[,2] + 2*x[,3])) ##probability for missing for index 1
        pr2 = 1/(1 + exp(-3.6 - x[,1] - 2*x[,2] - x[,3] + res)) ##probability for missing for index 2
        pr3 = 1/(1 + exp(-2.25 + x[,1] - 1.5*x[,2] + 0.5*x[,3] + res))  ##probability for missing for index 6
        pr4 = 1/(1 + exp(-2 + 0.25*x[,1] + x[,2] + 0.5*x[,3]))  ##probability for missing for index 7
        
        prob = cbind(pr1, pr2, pr3, pr4)
        rm(list= c("pr1", "pr2", "pr3", "pr4"))
        
        miss.mat = matrix(0, nrow = n, ncol = length(miss.var.mar))
        for(ii in 1:n){
          for(jj in 1:length(miss.var.mar)){
            miss.mat[ii,jj] = rbinom(1, 1, prob[ii,jj]) ##generating approx 50% missing
            #print(prob[ii,jj])
          }
        }
        
        rm(ii)
        miss.mat[miss.mat==0] <- NA
        miss.mat[which(miss.mat==1)] <- dsgn.mtrx[,miss.var.mar][which(miss.mat==1)]
        
        n.miss.mar = which(apply(is.na(miss.mat) , 1, any)==TRUE)
        miss1 = length(n.miss.mar)
        if(miss1 < n.miss){
          n.miss.mcar = n.miss - miss1
          
          a1 = matrix(0,ncol = length(miss.var.mcar))
          ii = 1
          while(ii <= n.miss.mcar){
            samp1 = rbinom(length(miss.var.mcar), 1, 0.5)
            if(sum(samp1) != length(miss.var.mcar)){
              a1 = rbind(a1, samp1)
              ii = ii + 1
            }
          }
          a1 = a1[-1,]
          miss.mat1 = matrix(1, nrow = n, ncol = length(miss.var.mcar))
          indx1 = sample(setdiff(1:n, n.miss.mar), n.miss.mcar)
          miss.mat1[indx1,] = a1
          
          miss.mat1[miss.mat1==0] <- NA
          miss.mat1[which(miss.mat1==1)] <- dsgn.mtrx[,miss.var.mcar][which(miss.mat1==1)]
          #length(which(apply(is.na(miss.mat1) , 1, any)==TRUE))
          
          dsgn.mtrx.miss = matrix(0, nrow = n, ncol = p)
          dsgn.mtrx.miss[,comp.var] = dsgn.mtrx[,comp.var]
          dsgn.mtrx.miss[,miss.var.mar] = miss.mat
          dsgn.mtrx.miss[,miss.var.mcar] = miss.mat1
          
          final.dat=cbind(res,dsgn.mtrx.miss)   ## final data with missing covariates and response
          colnames(final.dat) = c("y",paste("x", 1:p, sep=""))
          prop = prop + length(which(apply(is.na(final.dat) , 1, any)==TRUE))/n  ## checking proportion of missing observations
          
          dat[[i]] = final.dat
          i =i + 1
        }
      }
      output = list(missing.proportion = prop/gen.size, data = dat)
      return(output)
    }
  }
}
