rm(list=ls())
library(MASS);library(msm)
t.wd<-'your_path'
setwd(t.wd)
load('goat.mod.Rdata')
source('covariate.predict.annual.survival.R')
agest<-read.table("agestruct.txt")
################INPUTS#######################################################
n.ages<-20 # number of ages, currently there is no living past 20
n.years<-100  #years to project
n.runs<-20 #number of simulations to run
init.N<-300 #inital total N, to be split evenly across ageXsexes
min.fecund.stage<-3 #minimum stage of reproductive females
F0=0.9 # max fecundity
F1=-0.000005 #1st order density dependence in fecundity
F2=-2*10^-7 #2nd order density dependence in fecundity
ratio.males<-0.5 #ratio of males born
density.dep<-1 #should reproductive density dependence be implemented?
betas.at.mean.value<-0 # run sim with betas=mean estimated betas
update.betas.yearly<-0 #should betas be sampled once at the beginning (0) or with a new sample
#at each time step?
kid.survivals<-c(0.50,0.50) #survivals for first stage (doesn't come from model)
ylng.survivals<-rep(0.65,2) #survivals for second stage (doesn't come from model)
#all other survivals come from sampling the betas
init.temp<-5.5;end.temp<-5.5;init.snow<-0;end.snow=0 #choose begin and end temperature and snow values
###############################################################################

env.end.points<-matrix(c(init.temp,init.snow,end.temp,end.snow),2,2,byrow=T)
env.parms<-sapply(1:2,function(x){seq(env.end.points[1,x],env.end.points[2,x],length.out=n.years-1)})
colnames(env.parms)<-c('temp','snw')
sex.ratios<-c(ratio.males,1-ratio.males)
stages<-0:5
stage.lengths<-c(1,1,2,2,3,11)
agesexes<-rep(rep(stages,stage.lengths),each=2)
agesex.names<-paste(agesexes,c('M','F'),sep='')
fecund.inds<-intersect((1:n.ages*2)[(1:n.ages*2)%%2==0],which(agesexes>=min.fecund.stage))
fecundities<-vector(length=length(fecund.inds))
F1=F1*density.dep;F2=F2*density.dep #density dependent parameters
lambdas<-matrix(NA,n.years-1,n.runs)
M<-matrix(0,n.ages*2,n.ages*2)#our transition matrix
N<-matrix(0,n.ages*2,n.years)
Nstor<-matrix(0,n.runs,n.years)
init.pop<-(agest/sum(agest))*init.N
N[,1]<-init.pop[,1]

for(j in 1:n.runs){
	if(betas.at.mean.value){
  		beta.samp<-goat.mod$results$beta[,1]
		} else {
  			beta.samp<-sample.from.rmark.betas(goat.mod,1)
		}
	for(t in 1:(n.years-1)){
		if(update.betas.yearly){
  			beta.samp<-sample.from.rmark.betas(goat.mod,1)
		}

		t.env<-env.parms[t,]
		survivals<-generate.survivals(goat.mod,beta.samp,t.env)
		survivals<-c(kid.survivals[1],survivals[1:5],kid.survivals[2],survivals[6:10])
		names(survivals)[c(1,7)]<-c('0F','0M')
		survivals[c(2,8)]<-ylng.survivals
		age.survivals<-sapply(1:length(agesex.names),function(x){
				survivals[which(names(survivals)==agesex.names[x])]
		})
		for(i in 1:(length(age.survivals)-2)){
			M[i+2,i]<-age.survivals[i]
		}
		fecundity<-F0*exp(F1*sum(N[,t]))*exp(F2*sum(N[,t])^2)
		fecundities[]<-fecundity
		fecund.mat=fecundities/2
		M[1:2,fecund.inds]<-fecund.mat*age.survivals[fecund.inds] #Taal Changed this Oct 22nd
		N[,t+1]=M %*% N[,t]
		N[N<0]=0
		lambdas[t,j]<-Re(eigen(M)$values[1])
	}#end one trajectory
	Nstor[j,]=colSums(N)#store population trajectory j
}#end all trajectories
traj.output<-list(env=c(init.temp,end.temp,init.snow,end.snow),fec.parms=c(F0,F1,F2),
			lambdas=signif(lambdas,5),trajectories=Nstor)
x11()
matplot(t(Nstor),type="l")

