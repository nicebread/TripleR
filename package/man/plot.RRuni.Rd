\name{plot.RRuni}
\alias{plot.RRuni}
\alias{plot.RRbi}
\alias{plot.RRmulti}


\title{Plot variance components from SRAs}
\description{
This function plots variance components for all RR objects.
}
\usage{
\method{plot}{RRuni}(x, \dots, measure=NA, geom="bar")
\method{plot}{RRmulti}(x, \dots, measure=NA, geom="scatter", conf.level=0.95, connect=FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{x}{A RR results object}
\item{\dots}{Further arguments passed to the plot function}
\item{measure}{Sets the labels to `behavior' or `perception'; if no parameter is provided the function uses the defaults from RR.style}   
\item{geom}{Style of plot: set `bar' or `pie' for single groups; `scatter' or `bar' for multiple groups}
\item{conf.level}{Confidence level for error bars in the scatter style}
\item{connect}{Should dots for variance components of each group be connected? (Looks usually very cluttered ...)}
}

\value{
A ggplot2 object is returned - that means, you can subsequently adjust the scale labels, etc.
}

\examples{
\dontrun{
	data(likingLong)
	RR1 <- RR(liking_a ~ perceiver.id*target.id, data=likingLong)
	plot(RR1)
	plot(RR1, geom="pie")
	
	RR2 <- RR(liking_a + metaliking_a ~ perceiver.id*target.id, data=likingLong)
	plot(RR2)
	
	
	data("multiLikingLong")
	RR1m <- RR(liking_a ~ perceiver.id*target.id|group.id, data=multiLikingLong)
	plot(RR1m)
	plot(RR1m, measure="perception")
	plot(RR1m, measure="perception", geom="bar")
	plot(RR1m, measure="perception", connect=TRUE)

	RR2m <- RR(liking_a + metaliking_a ~ perceiver.id*target.id|group.id, data=multiLikingLong)
	plot(RR2m)
}
}


