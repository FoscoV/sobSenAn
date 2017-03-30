if(any(ls() == "SAsobEN")){}else{
	SAsobEN <-new.env()
	SAsobEN$.conflicts.OK<-c()
}
SAsobEN$distDict<-data.frame("mass"=as.character(c("cauchy","gamma","lognormal","logistic","negative binomial","normal","weibull")),"stat"=as.character(c("cauchy","gamma","lnorm","logis","nbinom","norm","weibull")))

#"exp","geom","t",
#"exponential","geometric","t"

####read model's parameters


lwrDens<-function(parVal,shapeA1,shapeA2,shapeB1,shapeB2,distrib){
	denA<-get(distrib)(parVal,shapeA1,shapeA2)
	denB<-get(distrib)(parVal,shapeB1,shapeB2)
	pmin(denA,denB)
}

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

	print(candidateDdf[order(candidateDdf$GOFks),])
#	cat(c(namePara ,"is distributed with a probability aestimed at ",GoFfndDistr$p.value,"% \n each one of the values provided affect the distribution for ", 1- mean(singParEff),"% (mean)\n", 1-max(singParEff),"% (max) \n do you want to save it for further exploration? \n (y|n)"),fill=T)
#	promptGo<-scan(what="text",nmax=1)
#	while(promptGo != "y" & promptGo != "n"){
#		cat("answer y or n \n ")
#		promptGo<-scan(,what="text",nmax=1)}
#	if (promptGo == "y"){
#		npDist<-data.frame(param=namePara,dist="norm",P1=fndDist$estimate[1],P2=fndDist$estimate[2])
#		if(any(ls(SAsobEN) == "parDists")){
#			SAsobEN$parDists<-rbind(SAsobEN$parDists,npDist)}else{
#			SAsobEN$parDists<-npDist
#		}
#	}
}

SAssessDis<-function(fndPara,distrib){
	deltDens<-function(combPar,distriba=as.character(distrib)){
		altDist <- fitdistr(scrPara[,combPar],distriba)
		sharDens<-integrate(lwrDens,-Inf,Inf,fndDist$estimate[1],fndDist$estimate[2],altDist$estimate[1],altDist$estimate[2],ddist(distriba))$value
#		return(as.numeric(substr(sharDens,start=1,stop=4)[1]))
		return(sharDens)
	}


	fndDist<-fitdistr(fndPara,distrib)

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
	#print(outDf)
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
	SAsobEN$seqSob<-sobol(30000,length(SAsobEN$parDists[,1]),init=TRUE,scrambling=3)
}
modPar4run<-function(){
	SobSeqNew()
	for(field in seq(1,length(SAsobEN$parDists[,1]))){
		SAsobEN$parSeq<-SAsobEN$seqSob
		SAsobEN$parSeq[,field]<-get(paste("q",as.character(SAsobEN$parDists$dist[field]),sep=""))(SAsobEN$seqSob[,field],SAsobEN$parDists$P1[field],SAsobEN$parDists$P2[field])
		SAsobEN$parSeq<-as.data.frame(SAsobEN$parSeq)
	}
}


library(randtoolbox)

####print them out to a specific tab separated file  (include compatibility in format with SimLab)
#### write.csv(fndPara, file=file.choose()) funziona! Andrà usato un write.table con tab comeseparatore di campo e qualche altra impestata opzione per avere l'output formato simlab



#####prompt user for run the model in batch using the provided parameters set

####read results from a some kind of files (include compatibility in format with SimLab)
####variance studies.... here I'll have to study deeper!
