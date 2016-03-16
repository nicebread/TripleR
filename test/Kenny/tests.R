#library(rmarkdown); render("tests.R", "pdf_document")

#' ## Kenny data set
#' ### Univariate manifest
library(rio)
dat <- import("roundrobin.sav")

(RR1 <- RR(y ~ actor*partner|group, data=dat, se="SOREMO"))
(RR2 <- RR(y ~ actor*partner|group, data=dat, se="LashleyBond"))

#' ## Stefan data set
dat2 <- import("Raw_data_Ratings_Social_Mimicry.sav")

#' ### Univariate manifest
(SMZ_Mean_GlobalandActions.compounds <- RR(ZSM_Mean_GlobalandActions ~ subject * target | group, 
	data=dat2, se="SOREMO"))
(SMZ_Mean_GlobalandActions.compounds <- RR(ZSM_Mean_GlobalandActions ~ subject * target | group, 
	data=dat2, se="LashleyBond"))

(SMcompounds <- RR(social_mimicry_global_rater1_rater2_rater3 ~ subject * target | group, 
	data=dat2, se="SOREMO"))
(SMcompounds <- RR(social_mimicry_global_rater1_rater2_rater3 ~ subject * target | group, 
	data=dat2, se="LashleyBond"))

(Liking_1compounds <- 
	RR(liking_1 ~ subject * target | group, data=dat2, se="SOREMO"))
(Liking_1compounds <- 
	RR(liking_1 ~ subject * target | group, data=dat2, se="LashleyBond"))


#' ### Univariate manifest



#' ## Builtin data set

#' ### Bivariate manifest
data(multiLikingLong)
#manifest bivariate SRM analysis
(RR2m <- RR(liking_a + metaliking_a ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="SOREMO"))
(RR2m <- RR(liking_a + metaliking_a ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="LashleyBond"))

#' ### Univariate latent
data(multiLikingLong)
(RR2m <- RR(liking_a/liking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="SOREMO"))
str(RR2m)
(RR2m <- RR(liking_a/liking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="LashleyBond"))

#' ### Bivariate latent
#latent (construct-level) bivariate SRM analysis
(RR4m <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="SOREMO"))
(RR4m <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="LashleyBond"))
		
#' ### Bivariate latent, multi group
(RR2m <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="SOREMO"))
	
(RR2m <- RR(liking_a/liking_b + metaliking_a/metaliking_b ~ perceiver.id*target.id|group.id, 
	data=multiLikingLong, se="LashleyBond"))
	
	
	
	
# ---------------------------------------------------------------------
# mutliGroup data set

data(multiGroup)
RR1 <- RR(ex+ex~perceiver.id*target.id|group.id, data=multiGroup, na.rm=TRUE, se="LashleyBond")
print(RR1, digits=6)

data(likingLong)

RR2 <- RR(liking_a+liking_a~perceiver.id*target.id, data=likingLong, na.rm=TRUE, se="LashleyBond")


# single group bivariate manifest:OK
print(RR(liking_a+liking_b~perceiver.id*target.id, data=likingLong, na.rm=TRUE, se="LashleyBond"), digits=5)

# single group bivariate latent:OK
print(RR(liking_a/liking_b~perceiver.id*target.id, data=likingLong, na.rm=TRUE, se="LashleyBond"), digits=5)

data(multiLikingLong)

RR1 <- RR(liking_a~perceiver.id*target.id|group.id, data=multiLikingLong, na.rm=TRUE, se="LashleyBond")
print(RR1, digits=5)

print(RR(liking_a+liking_b~perceiver.id*target.id, data=multiLikingLong, na.rm=TRUE, se="LashleyBond"), digits=5)


RRSM <- import("Raw_data_Ratings_Social_Mimicry.sav")

RR(Zsocial_mimicry_global_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)

RR(Zsocial_mimicry_actions_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)

RR(Zsocial_mimicry_global_rater1_rater2_rater3+Zsocial_mimicry_actions_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)
RR(Zsocial_mimicry_global_rater1_rater2_rater3/Zsocial_mimicry_actions_rater1_rater2_rater3 ~ subject*target|group, data=RRSM)