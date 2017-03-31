if(any(ls() == "SAsobEN")){}else{
	SAsobEN <-new.env()
	SAsobEN$.conflicts.OK<-c()
}
SAsobEN$distDict<-data.frame("mass"=as.character(c("cauchy","gamma","lognormal","logistic","negative binomial","normal","weibull","uniform")),"stat"=as.character(c("cauchy","gamma","lnorm","logis","nbinom","norm","weibull","unif")))

#"exp","geom","t",
#"exponential","geometric","t"

####read model's parameters


lwrDens<-function(parVal,shapeA1,shapeA2,shapeB1,shapeB2,distrib){
	denA<-get(distrib)(parVal,shapeA1,shapeA2)
	denB<-get(distrib)(parVal,shapeB1,shapeB2)
	pmin(denA,denB)
}
library(MASS)

ddist<-function(dist){
	funa<-paste("d",as.character(SAsobEN$distDict[which(SAsobEN$distDict[,1]==dist),2]),sep="")
	class(funa)<-"function"
	return(funa)
	}
pdist<-function(dist){
	funa<-paste("p",as.character(SAsobEN$distDict[which(SAsobEN$distDict[,1]==dist),2]),sep="")
	class(funa)<-"function"
	return(funa)
	}
qdist<-function(dist){
	funa<-paste("q",as.character(SAsobEN$distDict[which(SAsobEN$distDict[,1]==dist),2]),sep="")
	class(funa)<-"function"
	return(funa)
	}

library(fitdistrplus)
SAaddPara<-function(){

	cat(c("Type the name of the parameter which sensitivity you want to analyse: \n"),fill=TRUE)
	namePara<-scan(,what="text",nmax=1)
	cat(c("Write down the values " , namePara, " may assume.\n (return blank when done)\n"),fill=TRUE)
	fndPara<-scan()

	seekDist<-function(densi){
		return(suppressWarnings(SAssessDis(fndPara,as.character(densi))))
	}

	#candidateDdf<- unlist(lapply(X=SAsobEN$distDict[,1],FUN=seekDist))
	candidateDdf<- data.frame(distribution=SAsobEN$distDict[,1])
	tmpRes<-matrix(rep(0,5),ncol=5)
	for(enne in seq(1,length(candidateDdf$distribution))){
		risultati <- try(seekDist(candidateDdf$distribution[enne]),silent=TRUE)
		if(is.numeric(risultati)){
			tmpRes<-rbind(tmpRes,risultati)
		}else{
			tmpRes<-rbind(tmpRes,c("E","R","R","O","R"))
		}
	}
	tmpRes<-tmpRes[-1,]
	candidateDdf$distPar1<-round(as.numeric(tmpRes[,1]),2)
	candidateDdf$distPar2<-round(as.numeric(tmpRes[,2]),2)
	candidateDdf$GOFks<-round(as.numeric(tmpRes[,3]),2)
	candidateDdf$singleEffMean<-round(as.numeric(tmpRes[,4]),2)
	candidateDdf$singleEffMax<-round(as.numeric(tmpRes[,5]),2)
####		c(fndDist$estimate[1],fndDist$estimate[2],GoFfndDistr$p.value,meansieff,maxsieff)

	cat(c(namePara ,"fits the following distribution (defined with the firsts 2 columns). \n
	 Goodness Of Fit (comparison with Kolmogorov-Smirnov) is shown in the third column. \n
	 Last Columns are filled with the mean effect of one parameter on the overall distribution and the more sigificant one. \n
	 Which distribution do you like more? \n (consider the number on left and look at the plot) \n"))
	print(candidateDdf[order(candidateDdf$singleEffMax),])

#Preparing Plot
	h<-hist(fndPara,main="Distribution",xlab=namePara)
	xfit<-seq(min(fndPara),max(fndPara),length=40)
	brlen<-diff(h$mids[1:2])
	croma<-rainbow(length(candidateDdf$distribution))
	legend("topright",legend=candidateDdf$distribution,fill=rainbow(length(candidateDdf$distribution)))

	denplot<-function(xfit,disdat,ord,brlen){
		yfit<-get(ddist(disdat$distribution[ord]))(xfit,disdat$distPar1[ord],disdat$distPar1[ord])
		yfit <- yfit*brlen*length(fndPara)
		lines(xfit, yfit, col=croma[ord], lwd=2)
}

	denplotBOOT<-function(nume){try(denplot(xfit,candidateDdf,nume,brlen))}
	lapply(X=seq(1,length(candidateDdf$distribution)),FUN=denplotBOOT)

	promptGo<-scan(,nmax=1)
	while(!any(seq(1,length(candidateDdf$distribution))==promptGo)){
		cat("Which one? \n ")
		promptGo<-scan(,nmax=1)}
	if(any(fndPara%%1 != 0)){
		discretBOOL<-"n"}else{
			cat(c("Parameters values provided are all integers. Do you have a discrete distribution? \n
			y \t only integers allowed for this parameter \n
			n \t continuos values are allowed, just a coincidence \n"))
			discretBOOL<-scan(,what="text",nmax=1)
			while(discretBOOL != "y" & discretBOOL != "n"){
				cat("answer y or n")
				discretBOOL<-scan(,what="text",nmax=1)
			}
		}

	npDist<-data.frame(param=namePara,dist=candidateDdf$distribution[promptGo],P1=as.numeric(tmpRes[promptGo,1]),P2=as.numeric(tmpRes[promptGo,2]),disc=discretBOOL)
	if(any(ls(SAsobEN) == "parDists")){
		SAsobEN$parDists<-rbind(SAsobEN$parDists,npDist)}else{
		SAsobEN$parDists<-npDist
	}
}

SAssessDis<-function(fndPara,distrib){
	deltDens<-function(combPar,distriba=as.character(distrib)){
		altDist <- fitdistr(scrPara[,combPar],distriba)
		sharDens<-integrate(lwrDens,-Inf,Inf,fndDist$estimate[1],fndDist$estimate[2],altDist$estimate[1],altDist$estimate[2],ddist(distriba))$value
#		return(as.numeric(substr(sharDens,start=1,stop=4)[1]))
		return(sharDens)
	}

	if(distrib!="uniform"){
		fndDist<-fitdistr(fndPara,distrib)}else{
		fndDist<-data.frame(estimate=c(min(fndPara),max(fndPara)))
	}



	#GoFfndDistr<-chisq.test(fndPara,p=pnorm(fndPara,fndDist$estimate[1],fndDist$estimate[2]),rescale.p=TRUE,simulate.p.value=TRUE)
	GoFfndDistr<-ks.test(fndPara,as.character(pdist(distrib)),fndDist$estimate[1],fndDist$estimate[2])
	#### and sensibility to parameters amount


	scrPara<-combn(fndPara,(length(fndPara)-1))



	#singParEff<-deltDens(1,"normal")
	singParEff<-unlist(lapply(X=seq(1,length(fndPara)),FUN=deltDens))
	#singParEff<-unlist(lapply(X=seq(1,length(fndPara)),FUN=deltDens))


	meansieff<- 1- mean(singParEff)
	maxsieff <- 1-max(singParEff)

	outDf<-c(fndDist$estimate[1],fndDist$estimate[2],GoFfndDistr$p.value,meansieff,maxsieff)
	return(outDf)
}



####accordingly to Confalonieri define the amount of sets required
modelRuns<- function(){
	cu <- 1
	runs<- (2^(cu+3)*2(length(SAsobEN$parDists[,1])+2))/length(SAsobEN$parDists[,1])
	#while(runs >= (2^(cu+3)*2(length(SAsobEN$parDists[,1])+2))/length(SAsobEN$parDists[,1])){
	#	cu <- cu+1
		return(5)
	#}
}
####adopt a distribution-based deformation of parameter space. Such as a quantile: I'm looking for a (0,1) range!
SobSeqNew<-function(){
	SAsobEN$seqSob<-sobol(30000,2*length(SAsobEN$parDists[,1]),init=TRUE,scrambling=3)
}
modPar4run<-function(){
	SobSeqNew()
	for(field in seq(1,length(SAsobEN$parDists[,1]))){
		SAsobEN$parSeq<-SAsobEN$seqSob
		SAsobEN$parSeq[,field]<-get(qdist(SAsobEN$parDists$dist[field]))(SAsobEN$seqSob[,field],as.numeric(SAsobEN$parDists$P1[field]),as.numeric(SAsobEN$parDists$P2[field]))
		SAsobEN$parSeq<-as.data.frame(SAsobEN$parSeq)
	}
	colnames(SAsobEN$parSeq)<-as.character(SAsobEN$parDists$param)
	SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")]<-round(SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")],digits=0)
}


library(randtoolbox)



biblio2sobol<-function(){
	parAddBOOT<-function(){
		suppressWarnings(SAaddPara())
		cat("Do you want to provide another parameter?\n (y|n) \n")
		morPam<-scan(,what="text",nmax=1)
		while(morPam != "y" & morPam != "n"){
			cat("answer y or n")
			morPam<-scan(,what="text",nmax=1)
		}
	return(morPam)
	}
	morPam<-"y"
	while(morPam =="y"){
		morPam<-parAddBOOT()
	}
	modPar4run()
	cat(c("Where do you wato to save the file for batch processing the EXTERNAL MODEL?"))
	write.table(SAsobEN$parSeq,file=file.choose(),eol = "\r\n" ,sep="\t",row.names=FALSE)
}

####read results from a some kind of files (include compatibility in format with SimLab)
####variance studies.... here I'll have to study deeper!
