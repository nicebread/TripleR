\name{parCor}
\alias{parCor}

\title{partial correlation}
\description{
Performs partial correlations between x and y, controlled for z.
}
\usage{
parCor(x,y,z)
}


\arguments{
\item{x}{First variable} 
\item{y}{Second variable}
\item{z}{Control variable. This variable is coerced into a factor; in the TripleR context z usually is the group id.}
}
\details{
Performs partial correlations between x and y, controlled for z. The control variable is coerced into a factor; in the TripleR context z usually is the group id. Do \emph{not} use this function with a continuous control variable - results will be wrong! Degrees of freedom for the \emph{t} test are reduced by \emph{g} - 1 (\emph{g} is the number of groups).
}

\value{
\item{par.cor}{partial correlation}
\item{df}{degrees of freedom for the t test}
\item{t.value}{t value}
\item{p}{p value}
}


\seealso{
\code{\link{RR}}, \code{\link{getEffects}}
}

\examples{
data(multiGroup)
data(multiNarc)

# the function 'head' shows the first few lines of a data structure:
head(multiNarc)

# calculate SRA effects for extraversion ratings
RR.style("p")
RR1 <- RR(ex ~ perceiver.id * target.id | group.id, multiGroup, na.rm=TRUE)

# merge variables to one data set
dat <- merge(RR1$effects, multiNarc, by="id")

# We now have a combined data set with SRA effects and external self ratings:
head(dat)

# function parCor(x, y, z) computes partial correlation between x and y, 
# controlled for group membership z
d1 <- parCor(dat$ex.t, dat$narc, dat$group.id)
d1


# disattenuate for target effect reliability
parCor2 <- d1$par.cor * (1/sqrt(attr(RR1$effects$ex.t, "reliability")))
parCor2
}


\keyword{htest}
