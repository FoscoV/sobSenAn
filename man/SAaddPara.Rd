\name{SAaddPara}
\alias{SAaddPara}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
prompt user throught the phases in assessing variables distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SAaddPara()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function ()
{
    cat(c("Type the name of the parameter which sensitivity you want to analyse: \n"),
        fill = TRUE)
    namePara <- scan(, what = "text", nmax = 1)
    cat(c("Write down the values ", namePara, " may assume.\n (return blank when done)\n"),
        fill = TRUE)
    fndPara <- scan()
    seekDist <- function(densi) {
        return(suppressWarnings(SAssessDis(fndPara, as.character(densi))))
    }
    candidateDdf <- data.frame(distribution = SAsobEN$distDict[,
        1])
    tmpRes <- matrix(rep(0, 5), ncol = 5)
    for (enne in seq(1, length(candidateDdf$distribution))) {
        risultati <- try(seekDist(candidateDdf$distribution[enne]),
            silent = TRUE)
        if (is.numeric(risultati)) {
            tmpRes <- rbind(tmpRes, risultati)
        }
        else {
            tmpRes <- rbind(tmpRes, c("E", "R", "R", "O", "R"))
        }
    }
    tmpRes <- tmpRes[-1, ]
    candidateDdf$distPar1 <- round(as.numeric(tmpRes[, 1]), 2)
    candidateDdf$distPar2 <- round(as.numeric(tmpRes[, 2]), 2)
    candidateDdf$GOFks <- round(as.numeric(tmpRes[, 3]), 2)
    candidateDdf$singleEffMean <- round(as.numeric(tmpRes[, 4]),
        2)
    candidateDdf$singleEffMax <- round(as.numeric(tmpRes[, 5]),
        2)
    cat(c(namePara, "fits the following distribution (defined with the firsts 2 columns). \n\n\t Goodness Of Fit (comparison with Kolmogorov-Smirnov) is shown in the third column. \n\n\t Last Columns are filled with the mean effect of one parameter on the overall distribution and the more sigificant one. \n\n\t Which distribution do you like more? \n (consider the number on left and look at the plot) \n"))
    print(candidateDdf[order(candidateDdf$singleEffMax), ])
    h <- hist(fndPara, main = "Distribution", xlab = namePara)
    xfit <- seq(min(fndPara), max(fndPara), length = 40)
    brlen <- diff(h$mids[1:2])
    croma <- rainbow(length(candidateDdf$distribution))
    legend("topright", legend = candidateDdf$distribution, fill = rainbow(length(candidateDdf$distribution)))
    denplot <- function(xfit, disdat, ord, brlen) {
        yfit <- get(ddist(disdat$distribution[ord]))(xfit, disdat$distPar1[ord],
            disdat$distPar1[ord])
        yfit <- yfit * brlen * length(fndPara)
        lines(xfit, yfit, col = croma[ord], lwd = 2)
    }
    denplotBOOT <- function(nume) {
        try(denplot(xfit, candidateDdf, nume, brlen))
    }
    lapply(X = seq(1, length(candidateDdf$distribution)), FUN = denplotBOOT)
    promptGo <- scan(, nmax = 1)
    while (!any(seq(1, length(candidateDdf$distribution)) ==
        promptGo)) {
        cat("Which one? (number on left) \n ")
        promptGo <- scan(, nmax = 1)
    }
    if (any(fndPara\%\%1 != 0)) {
        discretBOOL <- "n"
    }
    else {
        cat(c("Parameters values provided are all integers. Do you have a discrete distribution? \n\n\t\t\ty \t only integers allowed for this parameter \n\n\t\t\tn \t continuos values are allowed, just a coincidence \n"))
        discretBOOL <- scan(, what = "text", nmax = 1)
        while (discretBOOL != "y" & discretBOOL != "n") {
            cat("answer y or n")
            discretBOOL <- scan(, what = "text", nmax = 1)
        }
    }
    cat(c("Does your distribution have a truncation?\n a minimum value and/or a maximum one? \n (y|n) \n"))
    truncit <- scan(, what = "text", nmax = 1)
    while (truncit != "y" & truncit != "n") {
        cat("answer y or n")
        truncit <- scan(, what = "text", nmax = 1)
    }
    if (truncit == "y") {
        cat(c("Do you want to provide \n 1. \t a numeric \n 2. \t a cumulative density \n threshold? \n (1 | 2 ) \n"))
        thretru <- scan(, nmax = 1)
        while (thretru != 1 & thretru != 2) {
            cat("answer 1 or 2")
            thretru <- scan(, what = "text", nmax = 1)
        }
        cat(c("Digit the minimum. \n -Inf (case sensitive) for have it open on left \n "))
        minthr <- scan(, nmax = 1)
        if (thretru == 2) {
            minthr <- get(qdist(candidateDdf$distribution[promptGo]))(minthr,
                as.numeric(tmpRes[promptGo, 1]), as.numeric(tmpRes[promptGo,
                  2]))
        }
        cat(c("Digit the maximum. \n Inf (case sensitive) for have it open on right \n "))
        maxthr <- scan(, nmax = 1)
        if (thretru == 2) {
            maxthr <- get(qdist(candidateDdf$distribution[promptGo]))(maxthr,
                as.numeric(tmpRes[promptGo, 1]), as.numeric(tmpRes[promptGo,
                  2]))
        }
    }
    else {
        minthr <- -Inf
        maxthr <- Inf
    }
    npDist <- data.frame(param = namePara, dist = candidateDdf$distribution[promptGo],
        P1 = as.numeric(tmpRes[promptGo, 1]), P2 = as.numeric(tmpRes[promptGo,
            2]), disc = discretBOOL, mintrs = minthr, maxtrs = maxthr)
    if (any(ls(SAsobEN) == "parDists")) {
        SAsobEN$parDists <- rbind(SAsobEN$parDists, npDist)
    }
    else {
        SAsobEN$parDists <- npDist
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
