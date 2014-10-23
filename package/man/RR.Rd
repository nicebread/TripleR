\name{RR}
\alias{RR}
\alias{liking_a}
\alias{likingLong}
\alias{multiLikingLong}
\alias{liking_b}
\alias{metaliking_a}
\alias{metaliking_b}
\alias{multiGroup}
\alias{multiNarc}


\title{Triple R: Round-Robin Analyses Using R}
\description{
The function RR performs Social Relation Model (SRM) analyses for single or multiple round-robin groups. 
}
\usage{
RR(formula, data, na.rm=FALSE, minData=1, verbose=TRUE, g.id=NULL, index="", 
	exclude.ids="", varname=NA, minVar=localOptions$minVar, ...)
}


\arguments{
\item{formula}{a formula (see details) consisting of a measure (e.g. a rating) obtained with a round-robin group} 
\item{data}{A data frame in long format}
\item{na.rm}{If missing values are present, you have to set this parameter to TRUE}
\item{minData}{Sets the minimum of data points which have to be present in each row and column}
\item{verbose}{Defines if warnings and additional information should be printed on execution}
\item{g.id}{For internal use only; do not set this parameter}
\item{index}{set \code{index = 'enhance'} for additionally requesting an index for self enhancement (self rating - perceiver effect - target effect - group mean of self ratings; Kwan, John, Kenny, Bond, & Robins, 2004) along with the actor and partner effects.}
\item{exclude.ids}{For internal use only; do not set this parameter}
\item{varname}{The name stem of the effects variables. By default, this is the first variable passed in the formula. In case of latent constructs, however, it might be preferable to set a new name for the latent construct.}
\item{minVar}{Actor and partners effects are only calculated if the respective relative variance component is greater than minVar. Set minVar to NA if this cleaning should not be performed. For small groups, Kenny (1994) suggests a minimum relative variance of 10\% for the interpretation of SRM effects. In any case, actor/ partner effects and correlations with these variables should not be interpreted if these components have negative variance estimates. minVar defaults to zero; with RR.style this default can be changed for all subsequent analyses.}

\item{...}{Further undocumented or internal arguments passed to the function}
}
\details{
Please note: detailed instructions on how to use the TripleR package are provided in the built in pdf document `How to use TripleR'. You can open this document either by on this link: \url{../doc/TripleR.pdf}, or you can browse all included vignettes by opening the index of the package documentation (scroll down to the very end of this page and click on `Index'; than click on `Overview of user guides and package vignettes'). These help files are only for quick references.

The estimation of the parameters of the Social Relation Model (SRM) is based on formulae provided by Kenny (1994; p. 236-244). For tests of significance of a single group, Triple R computes standard errors by using formulae published by Bond and Lashley (1996) for the case of a univariate SRM analysis. The formulae for the standard errors of the bivariate SRM parameters were kindly provided to us by C.F. Bond in personal communication. If multiple groups are provided, by default a between-group t-test is employed to calculate the significance. If you have very few groups with a considerable size (n>15), even in the multiple group scenario it might be preferable to use the within-group test of significance. You can inspect the within-group test of significance for each of the multiple groups in the return value (see \code{groups.univariate}).

The formula to specify the round robin design has following notation:
\code{DV ~ perceiver.id * target.id | group.id} (group.id is only provided if multiple groups are present). If two variables are used to describe one latent construct, both are connected with a slash on the left hand side: \code{DV1a/DV1b ~ perceiver.id * target.id}. If two variables are used to describe two manifest constructs, both are connected with a \code{+} on the left hand side: \code{DV1+DV2 ~ perceiver.id * target.id}. A latent analysis of two constructs would be notated as following: \code{DV1a/DV1b + DV2a/DV2b ~ perceiver.id * target.id}.


Data sets from the Mainz Freshman Study (see Back, Schmukle, & Egloff, in press) are included (liking_a, liking_b, metaliking_a, metaliking_b, likingLong), as well as an additional data set containing missing values (multiGroup)

The handling for missing data (\code{na.rm=TRUE}) is performed in three steps:
\itemize{
\item Rows and columns which have less then \code{minData} data points are removed from the matrix (i.e. both the `missing' row or column and the corresponding column or row, even if they have data in them)
\item For the calculation of actor and partner variances, actor-partner-covariances and the respective effects, missing values are imputed as the average of the respective row and col means. The calculation of relationship variances and covariances as well as relationship effects is also based on the imputed data set; however, single relationship effects which were missing in the original data set are set to missing again.
\item In the case of multiple variables (i.e., latent or bivariate analyses), participants are excluded listwise to ensure that all analyses are based on the same data set.
}
}

\value{
Printed are both unstandardized and standardized SRM estimates with the corresponding standard errors and t-values for actor variance, partner variance, relationship variance, error variance, actor-partner covariance, and relationship covariance. In case of a bivariate analysis values are additionally provided for actor-actor covariance, partner-partner covariance, actor-partner covariance, partner-actor covariance, intrapersonal relationship covariance, and interpersonal relationship covariance.

Reliabilities of actor, partner, and relationship effects (the latter only in the case of latent analyses) are printed according to Bonito & Kenny (2010).

The returned values depend on the kind of analyses that is performed:

Univariate, single group:

\item{effects}{actor and partner effects for all participants; if self ratings are provided they also are returned as group mean centered values}
\item{effects.gm}{actor and partner effects for all participants, group mean added}
\item{effectsRel}{relationship effects for all participants}
\item{varComp}{variance components}

Bivariate, single group:
\item{univariate}{List of results of univariate of SRM analyses of both variables- specify variable in index: univariate[[1]] or univariate[[2]]. That means, each of the both \code{univariate} objects is a complete RR object for the univariate analyses, nested in the results objects. If you want to retrieve the effects for the first variable, you have to type \code{RR2$univariate[[1]]$effects}. If you want to retrieve the variance components, you have to type \code{RR2$univariate[[1]]$varComp}}
\item{bivariate}{Results of bivariate SRM analyses}


In the multiple group case, following values are returned:

\item{univariate}{The weighted average of univariate SRM results}
\item{bivariate}{The weighted average of bivariate SRM results}
\item{groups}{SRM results for each group}
\item{effects}{actor and partner effects for all participants}
\item{effectsRel}{relationship effects for all participants}
\item{varComp.group}{a list of variance components in all single groups}
\item{group.var}{group variance estimate}


If self ratings are present in the data set, the function also prints the correlation between self ratings and actor/partner effects. In case of multiple groups, these are corrected for group membership (partial correlations). These correlations with self-ratings can also directly be computed with the function \code{\link{selfCor}}. Partial correlations to external (non-SRM) variables can be computed with the function \code{\link{parCor}}

}



\note{In case that a behavior was measured, the variances of an SRM analysis are labeled as actor variance, partner variance and relationship variance (default output labels). In case that a perception was measured, perceiver is printed instead of actor and target is printed instead of partner. You can set appropriate output labels by using the function \code{\link{RR.style}}.
These settings from RR.style, however, can be overwritten for each \code{print} call: :

\itemize{
	\item \code{print(RRobject, measure1="behavior")}: prints output for a univariate analysis of behavior data.
	
	\item \code{print(RRobject, measure1="perception")}: prints output for a univariate analysis of perception data.

	\item \code{print(RRobject, measure1="behavior", measure2="behavior")}: prints output for a bivariate analysis of behavior data.
	
	\item \code{print(RRobject, measure1="perception", measure2="perception")}: prints output for a bivariate analysis of perception data.
	
	\item \code{print(RRobject, measure1="behavior", measure2="perception")} or \cr
	\code{print(RRobject, measure1="perception", measure2="behavior")}: prints output for a bivariate analysis of behavior and perception data.
	
	\item \code{print(RRobject, measure1="perception", measure2="metaperception")}: is for the special case that a perception and a corresponding metaperception was measured. In this case, additionally the appropriate output labels for bivariate SRM analyses of other- and metaperceptions (reciprocities, assumed reciprocities, meta-accuracies; see Kenny, 1994) are presented.
}


You can plot any RR object using \code{plot(RR)}. See \code{\link{plot.RRuni}} for possible parameters.
}

\references{

The main reference for the TripleR package is:
\itemize{
\item{Schönbrodt, F. D., Back, M. D., & Schmukle, S. C. (2012). TripleR: An R package for social relations analyses based on round-robin designs. \emph{Behavior Research Methods, 44}, 455-470. doi:10.3758/s13428-011-0150-4}
}

Please cite this paper if you use TripleR in your research.

Further information on SRM and its application:
\itemize{  
\item{Back, M. D., Schmukle, S. C. & Egloff, B. (in press). A closer look at first sight: Social relations lens model analyses of personality and interpersonal attraction at zero acquaintance. \emph{European Journal of Personality}.}

\item{Bonito, J. A., & Kenny, D. A. (2010). The measurement of reliability of social relations components from round-robin designs. \emph{Personal Relationships, 17}, 235 - 251.}

\item{Kenny, D. A. (1994). \emph{Interpersonal perception: A social relations analysis}. New York: Guilford Press.}

\item{Kwan, V. S. Y., John, O. P., Kenny, D. A., Bond, M. H., & Robins, R. W. (2004). Reconceptualizing individual differences in self-enhancement bias: An interpersonal approach. \emph{Psychological Review, 111}, 94–110.}

\item{Lashley, B. R., & Bond, C. F., Jr. (1997). Significance testing for round robin data. \emph{Psychological Methods}, \emph{2}, 278-291.}}
}


\author{Felix D. Schönbrodt, Mitja D. Back, Stefan C. Schmukle}

\seealso{
\code{\link{RR.style}}, \code{\link{getEffects}}, \code{\link{plot_missings}}, \code{\link{long2matrix}}, \code{\link{matrix2long}}, \code{\link{plot.RRuni}}, \code{\link{RR.summary}}, \code{\link{selfCor}}, \code{\link{parCor}}
}

\examples{
# The example data are taken from the "Mainz Freshman Study" and consist 
# of ratings of liking as well as ratings of the metaperception of 
# liking at zero-acquaintance using a Round-Robin group of 54 participants 

#------------------------------------------------------------
# ----  Single group   --------------------------------------
#------------------------------------------------------------

# Load data frame in long format - it contains 4 variables:
#liking ratings indicator a (liking_a, "How likeable do you find this person?")  
#liking ratings indicator b (liking_b, "Would you like to get to know this person?")
#metaliking ratings indicator a (metaliking_a, "How likeable does this person find you?")
#metaliking ratings indicator b (metaliking_b, "Would this person like to get to know you?")
data("likingLong")

#manifest univariate SRM analysis
RR1 <- RR(liking_a ~ perceiver.id*target.id, data=likingLong)

\dontrun{
#manifest bivariate SRM analysis
RR2 <- RR(liking_a + metaliking_a ~ perceiver.id*target.id, data=likingLong)

#latent (construct-level) univariate SRM analysis
RR3 <- RR(liking_a / liking_b ~ perceiver.id*target.id, data=likingLong)

#latent (construct-level) univariate SRM analysis, define new variable name for the latent construct
RR3b <- RR(liking_a / liking_b ~ perceiver.id*target.id, data=likingLong, varname="liking")
# compare:
head(RR3$effects)
head(RR3b$effects)

#latent (construct-level) bivariate SRM analysis
RR4 <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id, data=likingLong)


# prints output of the manifest univariate analysis
# in terms of actor and partner variance (default output labels)
print(RR1, measure1="behavior") 

# prints output of the manifest univariate analysis 
# in terms of perceiver and target variance (appropriate for perception data)
print(RR1, measure1="perception")

# prints output of the manifest bivariate SRM analysis appropriate 
# for perception-metaperception data  
print(RR2, measure1="perception", measure2="metaperception")

#prints output of the latent univariate SRM analysis
# appropriate for perception data  
print(RR3, measure1="perception") 

#prints output of the latent bivariate SRM analysis
# appropriate for perception-perception data  
# Note: you can use abbreviations of the strings "behavior", "perception", and "metaperception"
print(RR4, measure1="p", measure2="p") 
}


#------------------------------------------------------------
# ----  Multiple groups --------------------------------------
#------------------------------------------------------------

# data("multiLikingLong") is a variant of the liking data set (see above) with multiple groups
data("multiLikingLong")

# set RR.style to "perception" (affects subsequent printing of objects)
RR.style("perception")

#manifest univariate SRM analysis
RR1m <- RR(liking_a ~ perceiver.id*target.id|group.id, data=multiLikingLong)

\dontrun{
#manifest bivariate SRM analysis
RR2m <- RR(liking_a + metaliking_a ~ perceiver.id*target.id|group.id, data=multiLikingLong)

#latent (construct-level) univariate SRM analysis
RR3m <- RR(liking_a / liking_b ~ perceiver.id*target.id|group.id, data=multiLikingLong)

#latent (construct-level) bivariate SRM analysis
RR4m <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id|group.id,
		 		data=multiLikingLong)

# prints output of the manifest univariate analysis
# in terms of actor and partner variance (default output labels)
print(RR1m, measure1="behavior") 

# prints output of the manifest univariate analysis 
# in terms of perceiver and target variance (appropriate for perception data)
print(RR1m, measure1="perception")


#------------------------------------------------------------
# ----  Multiple groups with missing values --------------------------------------
#------------------------------------------------------------

# a multi group data set with two variables:
# ex = extraversion ratings, and ne = neurotizism ratings
data("multiGroup")

#manifest univariate SRM analysis, data set with missings
RR1miss <- RR(ex~perceiver.id*target.id|group.id, data=multiGroup, na.rm=TRUE)

#manifest univariate SRM analysis, data set with missings, 
# minimum 10 data points are requested for each participant
RR1miss <- RR(ex~perceiver.id*target.id|group.id, data=multiGroup, na.rm=TRUE, minData=10)
}

}

\keyword{univar}
\keyword{htest}
