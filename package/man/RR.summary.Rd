\name{RR.summary}
\alias{RR.summary}

\title{Print group descriptives}
\description{
Shows descriptives for multiple groups (how many groups, which groups are excluded, group sizes)
}

\usage{
RR.summary(formule, data)
}

\arguments{
\item{formule}{A formula specifying the variable names. Example usage: \code{DV ~ perceiver.id*target.id | group.id} (group.id only necessary in multi group data sets)}
\item{data}{A data frame in long format)}
}

\value{
Printed output
}


\examples{

data("multiGroup")
RR.summary(ex~perceiver.id*target.id|group.id, data=multiGroup) 
}
