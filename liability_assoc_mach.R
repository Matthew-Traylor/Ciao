# Matthew Traylor 

# Last updated: Friday, April 4, 2014
# R code for implementing LTSCORE function of Zaitlen et al, PLoS Genetics 2012
# Adapted from C++ code available from http://www.hsph.harvard.edu/alkes-price/files/2012/10/ltsoft-1.3.tar_.gz
# Implements LTSCORE approach on SNPTEST format imputed data using expected genotype dosages

args<-commandArgs()
posterior.liabilities.file<-args[4]
mldosefile<-args[5]
mlinfofile<-args[6]
skip<-as.numeric(args[7])
chunk.size<-as.numeric(args[8])
output.filename<-args[9]
include.pc<-args[10]
no.snps<-as.numeric(args[11])
mlphenofile<-args[12]
long<-args[13]

data<-NULL
results<-NULL;results.beta<-NULL
logistic<-NULL;logistic.beta<-NULL
age.all<-NULL;age.beta<-NULL
att<-NULL
reciproot2pi<-1/sqrt(2*pi)
root2<-sqrt(2)
ind<-NULL
denom.case<-NULL;denom.control<-NULL
ltscore<-NULL
eigenstrat.solution<-NULL;eigenstrat.lt.solution<-NULL
rho1<-NULL
zg<-NULL
y<-NULL
ap_p<-NULL;ap_g<-NULL;ap_lp<-NULL


file.mldose<-mldosefile
file.mlinfo<-mlinfofile
covariates<-read.table(mlphenofile,header=T)
pheno.final<-read.table(posterior.liabilities.file,header=T,colClasses=c("character","character",rep("numeric",6)))
lines<-is.na(pheno.final$liability)
which.lines<-which(!lines)

info<-read.table(file.mlinfo,header=T,strings=F)

###########
###########read in data in chunks
###########

chunk.size<-min(chunk.size,no.snps)
alleles<-as.matrix(read.table(file.mldose,header=F,nrows=dim(pheno.final)[1],colClasses=c(rep("numeric",chunk.size))))
pheno.final<-pheno.final[which.lines,]
print(max(which.lines))
print(dim(alleles))
alleles<-alleles[which.lines,]
info<-read.table(file.mlinfo,header=T,skip=(skip-2),nrows=chunk.size,colClasses=c(rep("character",3),rep("numeric",4)))
m<-as.vector(colMeans(alleles))
info_score <- signif(colMeans((t(t(alleles) - m))^2) / (m*(1-m/2)),5)  
freq <- signif(colSums(alleles)/(2*nrow(alleles)),5)
N <- rep(length(which.lines),chunk.size)

if(include.pc==TRUE){
        evec<-covariates[,6:dim(covariates)[2]]
        evec<-cbind(evec[which.lines,])
        pheno.final$mean.adj.liability<-pheno.final$liability-sum(pheno.final$liability,na.rm=T)/dim(pheno.final)[1]
}


if(include.pc==FALSE){
covar.output <- rep("none",chunk.size)
    ###########
        ###########     calculate armitage trend test for case of no covariates
        ###########     
        sumx<-colSums(alleles)
        sumy.att<-sum(pheno.final[,4])
	sumxx<-colSums((alleles)*(alleles))
        sumxy.att<-colSums((alleles)*(pheno.final[,4]))
	sumyy.att<-sum(pheno.final[,4]*pheno.final[,4])
        sum1<-length(which.lines)
        print(alleles)
	mx<-sumx/sum1
	my.att<-sumy.att/sum1
        print(sumxx/sum1-mx*mx)
	sx<-sqrt(sumxx/sum1-mx*mx)
        sy.att<-sqrt(sumyy.att/sum1-my.att*my.att)
        rho<-(sumxy.att/sum1 - mx*my.att)/(sx*sy.att)

        att = sum1*rho*rho



###########
###########calculate LTSCORE for case of no covariates
###########
sumy<-sum(pheno.final$liability)
sumxy<-colSums((alleles)*(pheno.final$liability))
sumyy<-sum((pheno.final$liability)*(pheno.final$liability))
sum1<-dim(pheno.final)[1]

my<-sumy/sum1
sy<-sqrt(sumyy/sum1-my*my)
rho<-(sumxy/sum1 - mx*my)/(sx*sy)

ltscore = sum1*rho*rho
}



###########
###########if principal components are included, correct for population structure
###########
if(include.pc=="TRUE"){
covar.output <- rep(toString(names(evec)),chunk.size)
#calculate mean-adjusted genotypes
sum1<-length(which.lines)
sum<-colSums(alleles)/2
diff<-(-sum/sum1)
diff_m<-as.matrix(diff)
diff_rep<-diff_m[,rep(1,sum1)]
gamma<-alleles/2+t(diff_rep)


###########
###########compute eigenstrat solution
###########
#adjust phenotype for PC
for(cov in 1:dim(evec)[2]){
zp<-sum(pheno.final$liability_variable_mean*evec[,cov],na.rm=T)
y[cov]<-sum(evec[,cov]*evec[,cov],na.rm=T)
if(y[cov]>0){
ap_p<-(zp*y[cov])/(y[cov]*y[cov])
}else{
ap_p<-0
}
if(cov==1){
eigenstrat.pheno<-pheno.final$liability_variable_mean-ap_p*(evec[,cov])
} else {
eigenstrat.pheno<-eigenstrat.pheno-ap_p*(evec[,cov])
}
}



#adjust genotype for PC###
###can this part be vectorized???
###
eigenstrat.geno<-matrix(NA,dim(gamma)[1],dim(gamma)[2])
for(cov in 1:dim(evec)[2]){
for(snp in 1:chunk.size){
zg[snp]<-sum(gamma[,snp]*evec[,cov])
}
if(y[cov]>0){
ap_g<-(zg*y[cov])/(y[cov]*y[cov])
}else{
ap_g<-0
}
if(cov==1){
for(snp in 1:chunk.size){
eigenstrat.geno[,snp]<-gamma[,snp]-ap_g[snp]*evec[,cov]
}
} else {
 for(snp in 1:chunk.size){
                                        eigenstrat.geno[,snp]<-eigenstrat.geno[,snp]-ap_g[snp]*evec[,cov]
                                }
}
}


#compute correlation between corrected genotype and phenotypes
corx<-colSums(eigenstrat.geno)
cory<-sum(eigenstrat.pheno)
corxx<-colSums(eigenstrat.geno*eigenstrat.geno)
coryy<-sum(eigenstrat.pheno*eigenstrat.pheno)
corxy<-colSums(eigenstrat.geno*eigenstrat.pheno)

cormx<-corx/sum1
cormy<-cory/sum1
corsx<-sqrt(corxx/sum1-cormx*cormx)
corsy<-sqrt(coryy/sum1-cormy*cormy)

rho<-(corxy/sum1-cormx*cormy)/(corsx*corsy)
s<-sum1-(dim(evec)[2]+1)
eigenstrat.solution<-s*rho*rho



###########
###########compute eigenstrat solution using liabilities as phenotype
###########
for(cov in 1:dim(evec)[2]){
zlp<-sum(pheno.final$mean.adj.liability*evec[,cov],na.rm=T)
if(y[cov]>0){
ap_lp<-(zlp*y[cov])/(y[cov]*y[cov])
}else{
ap_lp<-0
}
if(cov==1){
eig.lt.pheno<-pheno.final$mean.adj.liability-ap_lp*evec[,cov]
} else {
eig.lt.pheno<-eig.lt.pheno-ap_lp*evec[,cov]
}
}
#compute correlation between corrected genotype and mean-adjusted liabilities corrected for PC

cor.lty<-sum(eig.lt.pheno)
                cor.ltyy<-sum(eig.lt.pheno*eig.lt.pheno)
cor.ltxy<-colSums(eigenstrat.geno*eig.lt.pheno)
                cor.mlty<-cor.lty/sum1
                cor.slty<-sqrt(cor.ltyy/sum1-cor.mlty*cor.mlty)


rho<-(cor.ltxy/sum1-cormx*cor.mlty)/(corsx*cor.slty)
s<-sum1-(dim(evec)[2]+1)
eigenstrat.lt.solution<-s*rho*rho
rho1<-rho
}



###########
###########create data table for results
###########
if(include.pc==TRUE){
        if(!file.exists(output.filename)){
                results<-cbind(info[,c(1,2,3)],info_score,freq,signif(rho,5),signif(eigenstrat.solution,5),signif(eigenstrat.lt.solution,5),signif(pchisq(eigenstrat.solution,lower.tail=F,1),5),signif(pchisq(eigenstrat.lt.solution,lower.tail=F,1),5),N,covar.output)
                names(results)<-c("rsid","A1","A2","info","freq_A1","rho","chisq","chisq_LT","p-value","LT_p-value","N","covars")
                write.table(results,output.filename,row.names=F,col.names=T,quote=F,sep="\t")
        }else{ 
                results<-cbind(info[,c(1,2,3)],info_score,freq,signif(rho,5),signif(eigenstrat.solution,5),signif(eigenstrat.lt.solution,5),signif(pchisq(eigenstrat.solution,lower.tail=F,1),5),signif(pchisq(eigenstrat.lt.solution,lower.tail=F,1),5),N,covar.output)
                write.table(results,output.filename,row.names=F,col.names=F,quote=F,sep="\t",append=T)
        }
}else{ 
        if(!file.exists(output.filename)){
                results<-cbind(info[,c(1,2,3)],info_score,freq,signif(rho,5),signif(att,5),signif(ltscore,5),signif(pchisq(att,lower.tail=F,1),5),signif(pchisq(ltscore,lower.tail=F,1),5),N,covar.output)
                 names(results)<-c("rsid","A1","A2","info","freq_A1","rho","chisq","chisq_LT","p-value","p-value_LT","N","covars")
                 write.table(results,output.filename,row.names=F,col.names=T,quote=F,sep="\t")
        }else{ 
                results<-cbind(info[,c(1,2,3)],info_score,freq,signif(rho,5),signif(att,5),signif(ltscore,5),signif(pchisq(att,lower.tail=F,1),5),signif(pchisq(ltscore,lower.tail=F,1),5),N,covar.output)
                write.table(results,output.filename,row.names=F,col.names=F,quote=F,sep="\t",append=T)
        }
}
