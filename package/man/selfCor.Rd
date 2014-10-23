\name{selfCor}
\alias{selfCor}

\title{partial correlation}
\description{
Performs partial correlations between x and y, controlled for z.
}
\usage{
selfCor(x, digits=3, measure=NA)
}


\arguments{
\item{x}{An RR object} 
\item{digits}{Digits to which values are rounded in the output}
\item{measure}{Either "behavior" or "percetion". Only affects labeling of output; if NA, the current setting of RR.style is taken as default.}
}
\details{
Prints correlations between self ratings (if present in the round robin matrices) and SRA effects. In case of multiple groups, partial correlations are printed (controlled for group membership).
The output of selfCor is also printed in the standard RR-ouput.
}

\value{
A data frame with correlation coefficients and p values.
}


\seealso{
\code{\link{RR}}, \code{\link{getEffects}}
}

\examples{
data(multiGroup)

RR.style("p")
# a single group
RR1 <- RR(ex~perceiver.id*target.id, data=multiGroup[multiGroup$group.id=="2", ], na.rm=TRUE)
selfCor(RR1)

# multiple groups
RR2 <- RR(ex~perceiver.id*target.id|group.id, data=multiGroup, na.rm=TRUE)
selfCor(RR2)

}


\keyword{htest}
