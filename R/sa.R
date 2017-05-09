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
rdist<-function(dist){
	funa<-paste("r",as.character(SAsobEN$distDict[which(SAsobEN$distDict[,1]==dist),2]),sep="")
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
		cat("Which one? (number on left) \n ")
		promptGo<-scan(,nmax=1)}

	##Check for discrete distribution.
	if(any(fndPara%%1 != 0)){
		discretBOOL<-"n"}else{
			cat(c("Parameters values provided are all integers. Do you have a discrete distribution? \n
			y \t only integers allowed for this parameter \n
			n \t continuos values are allowed, just a coincidence \n"))
			discretBOOL<-scan(,what="text",nmax=1)
			while(discretBOOL != "y" & discretBOOL != "n" | length(discretBOOL)!=1){
				cat("answer y or n")
				discretBOOL<-scan(,what="text",nmax=1)
			}
		}
	#check for truncated distribution....this is a mess...
	cat(c("Does your distribution have a truncation?\n a minimum value and/or a maximum one? \n (y|n) \n"))
	truncit<-scan(,what="text",nmax=1)
	while(truncit != "y" & truncit != "n"){
				cat("answer y or n")
				truncit<-scan(,what="text",nmax=1)
	}
	if(truncit =="y"){
		cat(c("Do you want to provide \n 1. \t a numeric \n 2. \t a cumulative density \n threshold? \n (1 | 2 ) \n"))
		thretru<-scan(,nmax=1)
		while(thretru != 1 & thretru != 2){
				cat("answer 1 or 2")
				thretru<-scan(,what="text",nmax=1)
		}
		cat(c("Digit the minimum. \n -Inf (case sensitive) for have it open on left \n "))
		minthr<-scan(,nmax=1)
		if(thretru == 2 ){
			minthr<-get(qdist(candidateDdf$distribution[promptGo]))(minthr,as.numeric(tmpRes[promptGo,1]),as.numeric(tmpRes[promptGo,2]))
		}
		cat(c("Digit the maximum. \n Inf (case sensitive) for have it open on right \n "))
		maxthr<-scan(,nmax=1)
		if(thretru == 2 ){

			maxthr<-get(qdist(candidateDdf$distribution[promptGo]))(maxthr,as.numeric(tmpRes[promptGo,1]),as.numeric(tmpRes[promptGo,2]))
		}
	}else{
		#setting a range between 0.1 and 0.9 as default
		minthr<- get(qdist(candidateDdf$distribution[promptGo]))(0.1,as.numeric(tmpRes[promptGo,1]),as.numeric(tmpRes[promptGo,2]))
		maxthr<- get(qdist(candidateDdf$distribution[promptGo]))(0.9,as.numeric(tmpRes[promptGo,1]),as.numeric(tmpRes[promptGo,2]))
	}

	npDist<-data.frame(param=namePara,dist=candidateDdf$distribution[promptGo],P1=as.numeric(tmpRes[promptGo,1]),P2=as.numeric(tmpRes[promptGo,2]),disc=discretBOOL,mintrs=minthr,maxtrs=maxthr,origVal=paste(fndPara,collapse=","))
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
		return(30000)
	#}
}
####adopt a distribution-based deformation of parameter space. Such as a quantile: I'm looking for a (0,1) range!

modPar4run<-function(){

	#truBOOT<-function(numero,campo){truDist(SAsobEN$parDists$dist[campo],get(qdist(SAsobEN$parDists$dist[campo]))(SAsobEN$parDists$mintrs[campo],as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo])),get(qdist(SAsobEN$parDists$dist[campo]))(SAsobEN$parDists$mintrs[campo],as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo])),numero)}
	truBOOT<-function(numero){truDist(SAsobEN$parDists$dist[field],SAsobEN$parDists$mintrs[field],SAsobEN$parDists$maxtrs[field],numero,field)}

	truDist<-function(dista,low,hi,ics,campo){
	if(ics < hi && ics > low ){
		return(get(ddist(dista))(ics,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))/(get(pdist(dista))(hi,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))-get(pdist(dista))(low,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))))
	}else{return(0)}
}

	#seqSob<-sobol(30000,2*length(SAsobEN$parDists[,1]),init=TRUE,scrambling=3)
	seqSob<-sobol(30000,length(SAsobEN$parDists[,1]),init=TRUE,scrambling=3)
	SAsobEN$parSeq<-seqSob
	for(field in seq(1,length(SAsobEN$parDists[,1]))){
		if(SAsobEN$parDists$mintrs[field] != -Inf || SAsobEN$parDists$maxtrs[field] != Inf){
			#in this case we have to find out the truncated distribution

			someRandCDF<-get(rdist(SAsobEN$parDists$dist[field]))(50000,as.numeric(SAsobEN$parDists$P1[field]),as.numeric(SAsobEN$parDists$P2[field]))
			someRandCDF<-subset(someRandCDF,someRandCDF >= SAsobEN$parDists$mintrs[field]&someRandCDF <= SAsobEN$parDists$maxtrs[field])
			trudy<-edfun(someRandCDF,support=range(c(SAsobEN$parDists$mintrs[field],SAsobEN$parDists$maxtrs[field])),dfun=truBOOT)
			SAsobEN$parSeq[,field]<-trudy$qfun(SAsobEN$parSeq[,field])
		}else{
			SAsobEN$parSeq[,field]<-get(qdist(SAsobEN$parDists$dist[field]))(seqSob[,field],as.numeric(SAsobEN$parDists$P1[field]),as.numeric(SAsobEN$parDists$P2[field]))
		}
		SAsobEN$parSeq<-as.data.frame(SAsobEN$parSeq)
	}
	colnames(SAsobEN$parSeq)<-as.character(SAsobEN$parDists$param)
	SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")]<-round(SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")],digits=0)
}
library(spartan)
eFap<-function(thickness=65,cuRvESAMPLE=3){

	#truBOOT<-function(numero,campo){truDist(SAsobEN$parDists$dist[campo],get(qdist(SAsobEN$parDists$dist[campo]))(SAsobEN$parDists$mintrs[campo],as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo])),get(qdist(SAsobEN$parDists$dist[campo]))(SAsobEN$parDists$mintrs[campo],as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo])),numero)}
	truBOOT<-function(numero){truDist(SAsobEN$parDists$dist[field],SAsobEN$parDists$mintrs[field],SAsobEN$parDists$maxtrs[field],numero,field)}

	truDist<-function(dista,low,hi,ics,campo){
		if(ics < hi && ics > low ){
			return(get(ddist(dista))(ics,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))/(get(pdist(dista))(hi,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))-get(pdist(dista))(low,as.numeric(SAsobEN$parDists$P1[campo]),as.numeric(SAsobEN$parDists$P2[campo]))))
		}else{return(0)}
	}
	#Generting samples

	#creating a tempdir named... SAfast!
	dir.create("SAfast")
	efast_generate_sample(FILEPATH="SAfast",NUMCURVES=cuRvESAMPLE,NUMSAMPLES=thickness,PARAMETERS=SAsobEN$parDists$param,PMIN=rep(0,length(SAsobEN$parDists$param)),PMAX=rep(1,length(SAsobEN$parDists$param)))
	#I want the samples uniform between 0 and 1 so that they fit CDFs.
	#Now SAfast is filled with bunnches of files... I'm going to merge them...
	for(curNum in seq(1,cuRvESAMPLE)){
		for(turnPara in SAsobEN$parDists$param){
			nomeSegmFile<-file.path("SAfast",paste("Curve",curNum,"_",turnPara,".csv",sep=""))
			print(nomeSegmFile)
			aSegmPara<-read.csv(nomeSegmFile)
			if(any(ls(SAsobEN) == "parSeq")){
				SAsobEN$parSeq<-rbind(SAsobEN$parSeq,aSegmPara)}else{
				SAsobEN$parSeq<-aSegmPara
			}
		}
	}
	#Burning after read
	unlink("SAfast",recursive=T)

	#adapting parameters to their own distribution
	for(field in seq(1,length(SAsobEN$parDists[,1]))){
		if(SAsobEN$parDists$mintrs[field] != -Inf || SAsobEN$parDists$maxtrs[field] != Inf){
			#in this case we have to find out the truncated distribution

			someRandCDF<-get(rdist(SAsobEN$parDists$dist[field]))(thickness*15,as.numeric(SAsobEN$parDists$P1[field]),as.numeric(SAsobEN$parDists$P2[field]))
			someRandCDF<-subset(someRandCDF,someRandCDF >= SAsobEN$parDists$mintrs[field]&someRandCDF <= SAsobEN$parDists$maxtrs[field])
			trudy<-edfun(someRandCDF,support=range(c(SAsobEN$parDists$mintrs[field],SAsobEN$parDists$maxtrs[field])),dfun=truBOOT)
			SAsobEN$parSeq[,field]<-trudy$qfun(SAsobEN$parSeq[,field])
		}else{
			SAsobEN$parSeq[,field]<-get(qdist(SAsobEN$parDists$dist[field]))(SAsobEN$parSeq[,field],as.numeric(SAsobEN$parDists$P1[field]),as.numeric(SAsobEN$parDists$P2[field]))
		}
		SAsobEN$parSeq<-as.data.frame(SAsobEN$parSeq)
	}
	colnames(SAsobEN$parSeq)<-as.character(SAsobEN$parDists$param)
	SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")]<-round(SAsobEN$parSeq[,which(SAsobEN$parDists$disc == "y")],digits=0)

}


library(randtoolbox)

truDist<-function(dista,low,hi,ics){
	if(ics < hi && ics > low ){
		return(get(ddist(dista))(ics)/(get(pdist(dista))(hi)-get(pdist(dista))(low)))
	}else{return(0)}
}

suppressPackageStartupMessages(library(edfun))


biblio2parameter<-function(straight=FALSE){
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
	SAsobEN$parDists$param<-as.character(SAsobEN$parDists$param)
	SAsobEN$parDists<-SAsobEN$parDists[order((SAsobEN$parDists$param)),]
	if(!straight){
		cat(c("Where do you wato to save the file with parameters distribution?"))
		write.table(SAsobEN$parDists,file=file.choose(),eol = "\r\n" ,sep="\t",row.names=FALSE)
	}
	#modPar4run()
	#cat(c("Where do you wato to save the file for batch processing the EXTERNAL MODEL?"))
	#write.table(SAsobEN$parSeq,file=file.choose(),eol = "\r\n" ,sep="\t",row.names=FALSE)
}
biblio2eFast<-function(){
	#acquire distributions
	biblio2parameter(straight=T)
	#run the script for samples genereation
	SAsobEN$sampleXcur<-65
	#but first, adding the dummy parameter!
	if(any(SAsobEN$parDists$param == "Dummy")){
		SAsobEN$parDists<-SAsobEN$parDists[-which(SAsobEN$parDists$param == "Dummy"),]
	}
	Scemo<-data.frame(param="Dummy",dist="uniform",P1=0,P2=1,disc="n",mintrs=0,maxtrs=1,origVal="0,1")
	SAsobEN$parDists<-rbind(SAsobEN$parDists,Scemo)
	SAsobEN$parDists<-SAsobEN$parDists[order((SAsobEN$parDists$param)),]
	#ready to generate!
	eFap(thickness=SAsobEN$sampleXcur,cuRvESAMPLE=3)
	#export the samples for external run
	fileToWrite<-file.choose()
	write.table(SAsobEN$parSeq,sep="\t", file=fileToWrite,col.names=TRUE,row.names=FALSE,quote=FALSE)
	#sussposing system handle this by himself eol="\r\n",
	#save(list=ls(SAsobEN),file=paste(dirname(fileToWrite),strsplit(fileToWrite,".")[[1]][1],".SAd",sep=""),envir=SAsobEN)
	save(list=ls(SAsobEN),file="Hyperspace.SAd",envir=SAsobEN)
	cat(c("Output created!"))
	SAclean()
}

SAmorSam<-function(sammor){
	cat(c("Select the desired previous session saved \n "),fill=TRUE)
	oldSensSession<-file.choose()
	load(file=oldSensSession,envir=SAsobEN)
	eFap(thickness=sammor)
	fileToWrite<-oldSensSession
	write.table(SAsobEN$parSeq,sep="\t", file=fileToWrite,col.names=TRUE,row.names=FALSE,quote=FALSE)
	#sussposing system handle this by himself eol="\r\n",
	SAsobEN$sampleXcur<-sammor
	save(list=ls(SAsobEN),file=paste(dirname(fileToWrite),strsplit(fileToWrite,".")[[1]][1],".SAd",sep=""),envir=SAsobEN)
	cat(c("Output created!"))
	SAclean()
}


output2Sens<-function(resFile,RISULTATO,hyperspace){
	if(missing(resFile)){
		#find out a file format for ermes to give back the results, supposing tsv
		resFile<-read.table(file.choose(),sep="\t",header=TRUE)
	}else{
		resFile<-read.table(resFile,sep="\t",header=TRUE)
		#resFile<-read.csv(resFile)
		}
	if(missing(hyperspace)){
		cat(c("Where is the .SAd file related to the explored hyperspace?"))
		loadSensSession()
	}


	#imposing 3 resempling curves
	if(length(resFile[,1])%%3 != 0 ){stop("supposed wrong file, not able to trace the resempling sequences")}

	#going to split the output file into the results...
	dir.create("SAfast")
	curvRad<-length(resFile[,1])/3
	#for each curve resampled ... 3 is hard-set...
	for(svolta in seq(1,3)){
		#define last item in current curve
		startcurve<-curvRad*(svolta-1)
	#for each parameter
		for(paNum in seq(1,length(SAsobEN$parDists[,1]))){
			#create in SAfast a file "CurveX_ParameterY_Results.csv"
			strtprm<-startcurve+paNum*SAsobEN$sampleXcur-SAsobEN$sampleXcur+1
			ndprm<-startcurve+paNum*SAsobEN$sampleXcur
			suppressWarnings(write.csv(resFile[strtprm:ndprm,],file=file.path("SAfast",paste("Curve",svolta,"_Parameter",paNum,"_Results.csv",sep="")),row.names=FALSE,col.names=TRUE,quote=FALSE))
		}
	}
	# check who is the one which is going to be valued...]
	SIMoutPT<-setdiff(names(resFile),as.character(SAsobEN$parDists$param))
	#starting eFAST result analysis:

	efast_get_overall_medians("SAfast",3,PARAMETERS=as.character(SAsobEN$parDists$param),NUMSAMPLES=SAsobEN$sampleXcur,MEASURES=SIMoutPT)

	efast_run_Analysis("SAfast",MEASURES=as.array(as.character(SIMoutPT)),PARAMETERS=SAsobEN$parDists$param,NUMCURVES=3,NUMSAMPLES=as.numeric(SAsobEN$sampleXcur),OUTPUTMEASURES_TO_TTEST=1,TTEST_CONF_INT=0.95,GRAPH_FLAG=T,EFASTRESULTFILENAME="SAresults.csv")

	if(!missing(RISULTATO)){
		print("Name your Analysis OUTPUT.zip filename")
		zip(RISULTATO,c(file.path("SAfast",paste(as.array(as.character(SIMoutPT)),".pdf",sep="")),file.path("SAfast","SAresults.csv")))
	}

	#unlink("SAfast",recursive=T)
	unlink("SAfast/Curve*")
}
####read results from a some kind of files (include compatibility in format with SimLab)
####variance studies.... here I'll have to study deeper!

saveSensSession<- function(){
	attach(SAsobEN)
	save(list=ls(SAsobEN),file=paste(Sys.Date(),".SAd",sep=""))
	cat(c("Session saved in ",paste(Sys.Date(),".SAd",sep=""),". \n"),fill=TRUE)
}
loadSensSession<-function(){
	cat(c("Select the desired previous session saved \n "),fill=TRUE)
	oldSensSession<-file.choose()
	load(file=oldSensSession,envir=SAsobEN)
	cat(c("Session in loaded \n "),fill=TRUE)
}
SAclean<-function(){
	rm(list=ls(SAsobEN),envir=SAsobEN)
}

