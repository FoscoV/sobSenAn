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
deltDens<-function(combPar){
	altDist <- fitdistr(scrPara[,combPar],"normal")
	sharDens<-integrate(lwrDens,-Inf,Inf,fndDist$estimate[1],fndDist$estimate[2],altDist$estimate[1],altDist$estimate[2])
	return(as.numeric(substr(sharDens,start=1,stop=4)[1]))
}


SAaddPara<-function(){
	cat(c("Type the name of the parameter which sensitivity you want to analyse: \n"),fill=TRUE)
	namePara<-scan(,what="text")
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

	cat(c(namePara ,"is normally distributed with a probability aestimed at ",GoF2fndDistr$p.value,"% \n each one of the values provided affect the distribution for ", 1- mean(singParEff),"% (mean)\n", 1-max(singParEff),"% (max) \n do you want to save it for further exploration? \n (y|n)"),fill=T)
	promptGo<-scan(what="text",nmax=1)
	while(promptGo != "y" & promptGo != "n"){
		cat("answer y or n \n ")
		promptGo<-scan(,what="text",nmax=1)}
	if (promptGo == "y"){
		npDist<-data.frame(param=namePara,dist="norm",P1=fndDist$estimate[1],P2=fndDist$estimate[2])
		if(any(ls(SAsobEN) == "parDists")){
			rbind(SAsobEN$parDists,npDist)}else{
			SAsobEN$parDists<-npDist
		}
	}
}

####accordingly to Confalonieri define the amount of sets required
modelRuns<- function(){
	cu <- 1
	runs<- (2^(cu+3)*2(length(SAsobEN$parDists[,1])+2))/length(SAsobEN$parDists[,1])
	while(runs >= (2^(cu+3)*2(length(SAsobEN$parDists[,1])+2))/length(SAsobEN$parDists[,1])){
		cu <- cu+1
	}
}
####adopt a distribustion-based deformation of parameter space. Such as a quantile: I'm looking for a (0,1) range!


#--library(randtoolbox)
####generate pseudo-random LHS.... Sobol sequence, instead!
#--seqSob<-sobol(10,paraNumbe,init=TRUE,scrambling=2)
####print them out to a specific tab separated file  (include compatibility in format with SimLab)
#### write.csv(fndPara, file=file.choose()) funziona! Andrà usato un write.table con tab comeseparatore di campo e qualche altra impestata opzione per avere l'output formato simlab



#####prompt user for run the model in batch using the provided parameters set

####read results from a some kind of files (include compatibility in format with SimLab)
####variance studies.... here I'll have to study deeper!
