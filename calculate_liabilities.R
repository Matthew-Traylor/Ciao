args<-commandArgs()

file.sample<-args[4]
file.covar.m<-args[5]
file.covar.f<-args[6]
output<-args[7]



results<-NULL;results.beta<-NULL
logistic<-NULL;logistic.beta<-NULL
age.all<-NULL;age.beta<-NULL
att<-NULL
reciproot2pi<-1/sqrt(2*pi)
root2<-sqrt(2)
ind<-NULL
denom.case.m<-NULL;denom.control.m<-NULL
denom.case.f<-NULL;denom.control.f<-NULL





covar.m<-read.table(file.covar.m,header=F,fill=T,strings=F)
covar.f<-read.table(file.covar.f,header=F,fill=T,strings=F)
pheno<-read.table(file.sample,header=T,strings=F,as.is=T,na.strings="NA")
if(dim(pheno)[2]>5){include.pc<-T}else{include.pc<-F}
if(include.pc==T){
	evec<-cbind(pheno[,6:dim(pheno)[2]])
}
IS<-pheno[,4];
pheno<-pheno[,c(1:5)];names(pheno)<-c("ID_1","ID_2","sex","IS","age")



#	filter out male cases
cases.m<-subset(pheno,(IS==1 & sex==1))
#	filter out female cases
cases.f<-subset(pheno,(IS==1 & (sex==0 | sex==2)))

#	filter out male controls
controls.m<-subset(pheno,pheno[,4]==0 & sex==1)
#	filter out female controls
controls.f<-subset(pheno,pheno[,4]==0 & (sex==0 | sex==2))

#	calculate mean value and mean-adjust cases and controls (For PC-adjusted pheno)
cases<-subset(pheno,IS==1)
controls<-subset(pheno,(pheno[,4]==0 | pheno[,4]==2))
mean.value<-dim(cases)[1]/dim(pheno)[1]
controls.m$age.mean.pheno<-0-mean.value
controls.f$age.mean.pheno<-0-mean.value
cases.m$age.mean.pheno<-1-mean.value
cases.f$age.mean.pheno<-1-mean.value





###########
###########	calculate posterior mean residual liability for each individual 
###########
#give cases with missing covariate information the median value for casesa
cases.m$age[is.na(cases.m$age)]<-median(cases.m$age,na.rm=T)
cases.m$age.pheno<-(-(covar.m[2,2]*(cases.m$age-covar.m[2,1])+covar.m[1,1]))
cases.f$age[is.na(cases.f$age)]<-median(cases.f$age,na.rm=T)
cases.f$age.pheno<-(-(covar.f[2,2]*(cases.f$age-covar.f[2,1])+covar.f[1,1]))

#if any controls have liability variable information, then set missing values to median for controls, otherwise
#give controls the median value from cases
#males
if(sum(is.na(controls.m$age))==length(controls.m$age)){
	#controls$age<-rep(covar[2,1],dim(controls)[1])
	controls.m$age<-median(cases$age,na.rm=T)		
}else{
	controls.m$age[is.na(controls.m$age)]<-median(controls.m$age,na.rm=T)
}
controls.m$age.pheno<-(-(covar.m[2,2]*(controls.m$age-covar.m[2,1])+covar.m[1,1]))
#females
if(sum(is.na(controls.f$age))==length(controls.f$age)){
	#controls$age<-rep(covar[2,1],dim(controls)[1])
	controls.f$age<-median(cases$age,na.rm=T)		
}else{
	controls.f$age[is.na(controls.f$age)]<-median(controls.f$age,na.rm=T)
}
controls.f$age.pheno<-(-(covar.f[2,2]*(controls.f$age-covar.f[2,1])+covar.f[1,1]))


numer.case.m<-reciproot2pi*exp((-cases.m$age.pheno*cases.m$age.pheno)/2)
numer.case.f<-reciproot2pi*exp((-cases.f$age.pheno*cases.f$age.pheno)/2)
numer.control.m<-(-reciproot2pi*exp((-controls.m$age.pheno*controls.m$age.pheno)/2))
numer.control.f<-(-reciproot2pi*exp((-controls.f$age.pheno*controls.f$age.pheno)/2))

# calculate posterior mean residual liabilities for cases
#male
for(ind in 1:dim(cases.m)[1]){
	if(cases.m$age.pheno[ind]>0){
		z<-(cases.m$age.pheno[ind]/root2)
		denom.case.m[ind]<-1-(2*pnorm(z*sqrt(2))-1)	# use normal distribution equivalent of error function

	}
	else{
		z<-(-cases.m$age.pheno[ind]/root2)
		denom.case.m[ind]<-1+(2*pnorm(z*sqrt(2))-1)
	}
}
cases.m$liability<-numer.case.m/(denom.case.m/2)
#female
for(ind in 1:dim(cases.f)[1]){
	if(cases.f$age.pheno[ind]>0){
		z<-(cases.f$age.pheno[ind]/root2)
		denom.case.f[ind]<-1-(2*pnorm(z*sqrt(2))-1)	# use normal distribution equivalent of error function

	}
	else{
		z<-(-cases.f$age.pheno[ind]/root2)
		denom.case.f[ind]<-1+(2*pnorm(z*sqrt(2))-1)
	}
}
cases.f$liability<-numer.case.f/(denom.case.f/2)


# calculate posterior mean residual liabilities for controls
#male
for(ind in 1:dim(controls.m)[1]){
        if(controls.m$age.pheno[ind]>0){
                z<-(controls.m$age.pheno[ind]/root2)
                denom.control.m[ind]<-1+(2*pnorm(z*sqrt(2))-1)

        }
        else{
                z<-(-controls.m$age.pheno[ind]/root2)
                denom.control.m[ind]<-1-(2*pnorm(z*sqrt(2))-1)
        }
}
controls.m$liability<-numer.control.m/(denom.control.m/2)
#female
for(ind in 1:dim(controls.f)[1]){
        if(controls.f$age.pheno[ind]>0){
                z<-(controls.f$age.pheno[ind]/root2)
                denom.control.f[ind]<-1+(2*pnorm(z*sqrt(2))-1)

        }
        else{
                z<-(-controls.f$age.pheno[ind]/root2)
                denom.control.f[ind]<-1-(2*pnorm(z*sqrt(2))-1)
        }
}
controls.f$liability<-numer.control.f/(denom.control.f/2)


new.pheno<-rbind(controls.m,controls.f,cases.m,cases.f)
pheno<-pheno[,c(1:2)]
pheno.final<-merge(pheno,new.pheno,by=names(pheno)[1:2],all.x=T,sort=F)
pheno.final<-merge(pheno,pheno.final,by=names(pheno)[1:2],all.x=T,sort=F)
output.liabilities<-pheno.final[,c(1,2,3,4,5,6,7,8)];names(output.liabilities)<-c("ID_1","ID_2","sex","phenotype","liability_variable","liability_variable_mean","liability_threshold","liability")
# calculate mean adjusted liabilities if principal components are included
lines<-is.na(pheno.final$IS)
which.lines<-which(!lines)
pheno.final<-pheno.final[which.lines,]
if(include.pc==TRUE){
        evec<-cbind(evec[which.lines,])
	pheno.final$mean.adj.liability<-pheno.final$liability-sum(pheno.final$liability,na.rm=T)/dim(pheno.final)[1]
	print(paste("running LTDOSAGE with",dim(evec)[2],"covariate(s)",sep=" "))
	} else { 
	print(paste("running LTDOSAGE with no covariate(s)",sep=" "))
}
write.table(output.liabilities,output,row.names=F,col.names=T,sep="\t",quote=F)

