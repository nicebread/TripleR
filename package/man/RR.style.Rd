\name{RR.style}
\alias{RR.style}

\title{Set labeling styles for RR analyses}
\description{
This function sets labels for the printing of RR-objects and for plots. All subsequent calls of RR will be produced in this style (until another style is set). That means, usually you only have to define this style once at the start of your project.
}
\usage{
RR.style(style="behavior", suffixes=NA, minVar=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{style}{a string defining the labeling style - either `behavior' or `perception'.} 
\item{suffixes}{Which suffixes should be append to the actor and partner effects, and to the self ratings? Default is \code{.a}, \code{.p}, and \code{.s}, for actor/ partner/ self in the case of behavior, or \code{.p}, \code{.t}, and \code{.s}, for perceiver/ target/ self in the case of perceptions. If no suffixes are provided in the parameters, these defaults are taken depending on the \code{style} parameter} 
\item{minVar}{Set the minVar parameter for all subsequent analyses. See \code{\link{RR}} for details on this parameter.} 
}

\value{
Printing options and naming conventions are set for all subsequent analyses. If you specify other styles in a print.RR call, this setting temporarily overwrites the settings from RR.style (without changing them).
}



\author{Felix D. Schönbrodt}

\examples{
data("likingLong")

RR.style("behavior")
RR(liking_a ~ perceiver.id*target.id, data=likingLong)

RR.style("p")	# a "p" is enough for "perception"
RR(liking_a ~ perceiver.id*target.id, data=likingLong)

}


