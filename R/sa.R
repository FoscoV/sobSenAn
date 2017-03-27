if(any(ls() == "SAsobEN")){}else{
	SAsobEN <-new.env()
	SAsobEN$.conflicts.OK<-c()
}


####read model's parameters


lwrDens<-function(parVal,shapeA1,shapeA2,shapeB1,shapeB2){
	denA<-dnorm(parVal,shapeA1,shapeA2)
	denB<-dnorm(parVal,shapeB1,shapeB2)
	pmin(denA,denB)
}



SAaddPara<-function(){
	deltDens<-function(combPar){
		altDist <- fitdistr(scrPara[,combPar],"normal")
		sharDens<-integrate(lwrDens,-Inf,Inf,fndDist$estimate[1],fndDist$estimate[2],altDist$estimate[1],altDist$estimate[2])$value
#		return(as.numeric(substr(sharDens,start=1,stop=4)[1]))
		return(sharDens)
	}
	cat(c("Type the name of the parameter which sensitivity you want to analyse: \n"),fill=TRUE)
	namePara<-scan(,what="text",nmax=1)
	cat(c("Write down the values " , namePara, " may assume.\n (return blank when done)\n"),fill=TRUE)
	fndPara<-scan()

	####assess their distribution
		#gaussian
	fndDist<-fitdistr(fndPara,"normal")
		#beta
		#gamma
		#### a scelta supportate da un geom_density?

	#GoFfndDistr<-chisq.test(fndPara,p=pnorm(fndPara,fndDist$estimate[1],fndDist$estimate[2]),rescale.p=TRUE,simulate.p.value=TRUE)
	GoFfndDistr<-ks.test(fndPara,"pnorm",fndDist$estimate[1],fndDist$estimate[2])
	#### and sensibility to parameters amount


	scrPara<-combn(fndPara,(length(fndPara)-1))

	singParEff<-unlist(lapply(X=seq(1,length(fndPara)),FUN=deltDens))

	cat(c(namePara ,"is normally distributed with a probability aestimed at ",GoFfndDistr$p.value,"% \n each one of the values provided affect the distribution for ", 1- mean(singParEff),"% (mean)\n", 1-max(singParEff),"% (max) \n do you want to save it for further exploration? \n (y|n)"),fill=T)
	promptGo<-scan(what="text",nmax=1)
	while(promptGo != "y" & promptGo != "n"){
		cat("answer y or n \n ")
		promptGo<-scan(,what="text",nmax=1)}
	if (promptGo == "y"){
		npDist<-data.frame(param=namePara,dist="norm",P1=fndDist$estimate[1],P2=fndDist$estimate[2])
		if(any(ls(SAsobEN) == "parDists")){
			SAsobEN$parDists<-rbind(SAsobEN$parDists,npDist)}else{
			SAsobEN$parDists<-npDist
		}
	}
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
