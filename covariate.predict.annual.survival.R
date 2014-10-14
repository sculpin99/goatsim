
annual.from.mos<-function(mo.ests,pred.vcv,grouping.factor){
	group.sex.inds<-sapply(1:10,function(x){which(grouping.factor==x)})
	survivals<-lapply(1:10,function(x){
		mo.ests[group.sex.inds[,x]]
	})
	#get a vcv matrix for the survivals from each stageXsex
	vcvs<-lapply(1:10,function(x){
		pred.vcv[group.sex.inds[,x],group.sex.inds[,x]]
		})
	#get estimates from the monthly estimates and vcv matrices
	annual.survivals<-sapply(survivals,prod)
	annual.ses<-sapply(1:10,function(x){
			deltamethod.special('prod',survivals[[x]],vcvs[[x]],ses=T)
		})
	names(annual.survivals)<-names(annual.ses)<-paste(rep(1:5,2),rep(c('F','M'),each=5),sep='')
	annual.survival<-list(annual.survivals,annual.ses)
	names(annual.survival)<-c('annual.survivals','annual.ses')
	return(annual.survival)
}

generate.survivals<-function(RMark.mod,betas,env.parms,goat.mod=T){
	dm<-RMark.mod$design.matrix
	for(i in 1:length(env.parms)){
		dm[grep(names(env.parms)[i],dm)]<-env.parms[i]
	}
	dm2<-matrix(as.numeric(dm),dim(dm)[1],dim(dm)[2])
	monthly.survival<-sapply(1:nrow(dm2),function(x){
					e.sp<-exp(sum(betas*dm2[x,]))
					prob<-e.sp/(1 + e.sp)
					return(prob)
				})
	monthly.survival<-matrix(monthly.survival,7,12,byrow=T)
	annual.survival<-apply(monthly.survival,1,prod)
	if(goat.mod){
		annual.survival<-c(annual.survival[1:5],annual.survival[1:3],annual.survival[6:7])
		names(annual.survival)<-paste(rep(1:5,2),rep(c('F','M'),each=5),sep='')
	}
	return(annual.survival)
}



logit.survival<-function(survivals,ses){
	link.values = log(survivals/(1 - survivals))
	#see the utility RMark function 'logitCI' for where I got the following code.
	deriv.link.values = 1/survivals + 1/(1 - survivals)
	deriv.link.matrix = matrix(0, nrow = length(deriv.link.values), 
            	ncol = length(deriv.link.values))
	diag(deriv.link.matrix) = deriv.link.values
	vcv.real = diag(ses^2)
	link.vcv<-deriv.link.matrix %*% vcv.real %*% t(deriv.link.matrix)
	link.ses = sqrt(diag(link.vcv))
	logit.survivals<-list(link.values,link.ses,link.vcv)
	names(logit.survivals)<-c('logit.ests','logit.ses','logit.vcv')
	return(logit.survivals)
}

predict.survival<-function(RMark.model,covariates,n.samples=1,groupsexmo=1:nrow(RMark.model$design.data$S)){
	covs<-as.data.frame(t(covs))
	preds<-covariate.predictions(RMark.model,covs,indices=groupsexmo)
	#set up survivals by stagesex to obtain annual survivals and vcv's
	group.sex<-rep(1:10,each=12)
	annual.survivals<-annual.from.mos(mo.ests=preds$estimates$estimate,pred.vcv=preds$vcv,grouping.factor=group.sex)
	#Now we want to take a random sample within the CIs of the annual survivals
	#first get the back-transformed survival parameters
	logit.survivals<-logit.survival(survivals=annual.survivals$annual.survivals,ses=annual.survivals$annual.ses)
	#then sample from the normal distribution of the logit survivals
	samp<-mvrnorm(n=n.samples,mu=logit.survivals$logit.ests,Sigma=logit.survivals$logit.vcv)
	if(length(dim(samp))==0){samp<-as.matrix(t(samp))}
	#same as the following: samp<-rnorm(n=10000,mu=link.values,sd=sqrt(diag(vcov.mat.backtransformed.derived))
	#transform all of the samples from logit scale back to probabilities.
	sampled.survivals<-t(sapply(1:nrow(samp),function(x){
						sapply(1:ncol(samp),function(y){
							exp(samp[x,y])/(1 + exp(samp[x,y]))
						})
					}))
	colnames(sampled.survivals)<-groupsex.names<-paste(rep(1:5,2),rep(c('F','M'),each=5),sep='')
	annual.survivals<-list(annual.survivals,sampled.survivals)
	names(annual.survivals)<-c('deterministic','sampled.from.CIs')
	return(annual.survivals)
}

sample.from.rmark.betas<-function(RMark.mod,n.samples){
	betas<-goat.mod$results$beta[,1]
	beta.vcv<-goat.mod$results$beta.vcv
	sample<-t(mvrnorm(n=n.samples,mu=betas,Sigma=RMark.mod$results$beta.vcv))
	return(sample)
}

