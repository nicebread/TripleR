

#------------------------------------
#-- Skinning variables
#------------------------------------


localOptions <- new.env(parent=globalenv())
localOptions$suffixes <- c(".a", ".p", ".s")
localOptions$style <- "behavior"
localOptions$minVar <- 0

role <- list()
role$behavior <- c("Actor", "Partner", "Relationship")
role$perception <- c("Perceiver", "Target", "Relationship")
role$metaperception <- c("Perceiver", "Target", "Relationship")

# Labels actor-partner
unilabels_b <- c("actor variance", "partner variance", "relationship variance", "error variance", "actor-partner covariance", "relationship covariance")
# Labels target-perceiver
unilabels_p <- c("perceiver variance", "target variance", "relationship variance", "error variance", "perceiver-target covariance", "relationship covariance")

unilabels2 <- c("estimate", "standardized", "se", "t.value")

# labels for metaperception
unilabels_b_meta1 <- c("perceiver variance otherperception", "target variance otherperception",  "relationship variance otherperception", "error variance otherperception", "generalized reciprocity otherperception", "dyadic reciprocity otherperception")
unilabels_b_meta2 <- c("perceiver variance metaperception", "target variance metaperception", "relationship variance metaperception", "error variance metaperception", "generalized reciprocity metaperception", "dyadic reciprocity metaperception")

# Labels for bivariate analyses
bilabels_bb <- c("actor-actor covariance","partner-partner covariance",
"actor-partner covariance","partner-actor covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_pp <- c("perceiver-perceiver covariance","target-target covariance",
"perceiver-target covariance","target-perceiver covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_bp <- c("actor-perceiver covariance","partner-target covariance",
"actor-target covariance","partner-perceiver covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_pb <- c("perceiver-actor covariance","target-partner covariance",
"perceiver-partner covariance","target-actor covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")


bilabels_meta <- c("Perceiver assumed reciprocity","Generalized assumed reciprocity",
"Perceiver meta-accuracy", "Generalized meta-accuracy", "Dyadic assumed reciprocity", "Dyadic meta-accuracy")


#' @export
RR.summary <- function(formule, data) {
	# save warnings for later output
	l <- long2matrix(formule, data)
	cat(paste("Number of groups:",length(l))); cat("\n")
	cat(paste("Number of valid groups (i.e., more than 3 group members):\n",sum(laply(l, function(x) ncol(x) > 3)))); cat("\n")
	
	# which groups are excluded?
	excl.id <- c()
	for (i in 1: length(l)) {if (nrow(l[[i]]) <=3) excl.id <- c(excl.id, attr(l[[i]], "group.id"))}
	if (length(excl.id>0)) {
		cat("Following groups are excluded because they have 3 or less members:")
		cat(paste("\n\t",excl.id)) 
		cat("\n\n")
	}

	cat("Group sizes:\n")
	for (i in 1:length(l)) {cat(paste(attr(l[[i]], "group.id"),": n=",ncol(l[[i]]),"\n", sep=""))}

	cat("\n\n")
	cat("Group member statistics:"); cat("\n")
	gs <- laply(l, function(x) ncol(x))	# group sizes
	gm <- table(gs)
	names(gm) <- paste("n=",names(gm), sep="")
	print(gm); cat("\n")
	cat(paste("Min:",min(gs))); cat("\n")
	cat(paste("Max:",max(gs))); cat("\n")
	cat(paste("Mean:",round(mean(gs), 2))); cat("\n\n")
}

# set options for printing results etc.
# style = c("behavior", "perception")
#' @export
RR.style <- function(style="behavior", suffixes=NA, minVar=NA) {
	localOptions$style <- style <- match.arg(style, c("behavior", "perception"))
	
	if (is.na(suffixes)) {
		if (style=="behavior") {
			localOptions$suffixes <- c(".a", ".p", ".s")
		} else 
		if (style=="perception") {
			localOptions$suffixes <- c(".p", ".t", ".s")
		}
	} else {
		localOptions$suffixes <- suffixes
	}
	
	if (is.na(minVar)) {
		localOptions$minVar <- 0
	} else {
		localOptions$minVar <- minVar
	}
}


# corrects values < -1 to -1 and values > 1 to 1
clamp <- function(...) {
	x <- c(...)
	x[x < -1] <- -1
	x[x >  1] <-  1
	return(x)
}



clearLongData <- function(formule, data, minData=1) {
	ll1 <- long2matrix(formule, data, reduce=TRUE, minData=minData)
	
	lhs <- strsplit(gsub(" ","",as.character(formule)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(formule)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]
	
	var.id <- lhs
	actor.id <- rhs[1]
	partner.id <- rhs[2]
	if (length(rhs)>=3) {group.id <- rhs[3]} else {group.id="group.id"}
	
	
	
	ll2 <- ldply(ll1, function(x) {
		matrix2long(x, new.ids=FALSE, var.id=var.id)
	})
	colnames(ll2)[1:3] <- c(group.id, actor.id, partner.id)

	return(ll2)
}




# berechnet schnell die Effekte, ohne sonstigen Krimskrams
quickeffects <- function(RRMatrix) {
	RRMatrix <- as.matrix(RRMatrix)
	mip <- rowMeans(RRMatrix, na.rm=TRUE)
	mpj <- colMeans(RRMatrix, na.rm=TRUE)
	mpp <- mean(RRMatrix, na.rm=TRUE)  # grand mean

   n <- length(mip)
   a <- ((n-1)^2)/(n*(n-2))*mip + (n-1)/(n*(n-2))*mpj - (n-1)/(n-2)*mpp #actor effects
   b <- ((n-1)^2)/(n*(n-2))*mpj + (n-1)/(n*(n-2))*mip - (n-1)/(n-2)*mpp #partner effects
   
   am <- matrix(a, nrow(RRMatrix), ncol(RRMatrix))
   bm <- t(matrix(b, nrow(RRMatrix), ncol(RRMatrix)))
   c <- RRMatrix - am - bm - mpp		 # relationship effect

	return(list(a=a,b=b,c=c,m=mpp))
}

impute <- function(RRMatrix, na.rm="meansNA", stress.max = 0.01, maxIt=100) {

	# in Matrix umwandeln, sonst geht's nicht ...
	RRMatrix2 <- as.matrix(RRMatrix)
	rownames(RRMatrix2) <- rownames(RRMatrix)
	colnames(RRMatrix2) <- colnames(RRMatrix)
	RRMatrix <- RRMatrix2
	
	NAs <- is.na(RRMatrix)		# Matrix, die die Position der NAs ausserhalb der Diagonale speichert
	
	save.diag <- diag(RRMatrix)	# self ratings aus der Diagonale abspeichern
	diag(RRMatrix) <- NA

	eff0 <- eff <- quickeffects(RRMatrix)
	stress <- 1
	it <- 0

	# save evolution of parameters
	as <- matrix(eff$a, nrow=1, ncol=ncol(RRMatrix))
	bs <- matrix(eff$b, nrow=1, ncol=ncol(RRMatrix))
	ms <- c()
	
	while (stress > stress.max) {
	
		rM <- rowMeans(RRMatrix, na.rm=TRUE)
		cM <- colMeans(RRMatrix, na.rm=TRUE)
		rM_mean <- mean(rM)
		cM_mean <- mean(cM)
		grandmean <- mean(RRMatrix, na.rm=TRUE)
		
		for (i in 1:ncol(RRMatrix)) {
			for (j in 1:nrow(RRMatrix)) {
				if (NAs[j,i]==TRUE) {
					
					# Ersetzung durch mittleres Zeilen/ Spaltenmittel
					if (grepl("mean", na.rm)) {
						RRMatrix[j,i] <- (rM[j] + cM[i])/2
					}
				
				}
			}
		}

		if (grepl("1", na.rm)) break; # bei means1, chi1: beim ersten Durchgang gleich rausspringen

		eff <- quickeffects(RRMatrix)

		stress <- max(max(abs(eff$a - eff0$a)), max(abs(eff$b - eff0$b)))
		eff0 <- eff
		it <- it + 1
		
		if (it > maxIt) {
			warning("Maximum iterations exceeded; fall back to single imputation.", call.=FALSE)
			return(impute(RRMatrix2, paste(na.rm,"1",sep="")))
		}
		
		as <- rbind(as, eff$a)
		bs <- rbind(bs, eff$b)
		ms <- c(ms, eff$m)
	}	
	
	diag(RRMatrix) <- save.diag
	if (!grepl("NA", na.rm)) {NAs[] <- FALSE}
	
	return(list(RRMatrix=RRMatrix, NAs=NAs, iterations=it, as=as, bs=bs, ms=ms))
}




# calculates Actor-, Partner- and Relationship-Effects from a single RR-Matrix
RR.effects <- function(RRMatrix, name=NA, na.rm=FALSE, index="", varname="NA") {
	if (!is.na(varname)) {name <- varname} else {
		if (!is.null(attr(RRMatrix, "varname"))) name <- attr(RRMatrix, "varname")
	}
	
	RRold <- RRMatrix		
	if (na.rm & sum(is.na(RRMatrix))>nrow(RRMatrix)) {
		imp <- impute(RRMatrix)
		RRMatrix <- imp$RRMatrix
		imputation <- TRUE
	} else {
		imputation <- FALSE
	}
	

	RRMatrix2 <- as.matrix(RRMatrix)
	mip <- rowMeans(RRMatrix2, na.rm=TRUE)
	mpj <- colMeans(RRMatrix2, na.rm=TRUE)
	mpp <- mean(RRMatrix2, na.rm=TRUE)  # grand mean

	n <- length(mip)
	a <- ((n-1)^2)/(n*(n-2))*mip + (n-1)/(n*(n-2))*mpj - (n-1)/(n-2)*mpp #actor effects
	b <- ((n-1)^2)/(n*(n-2))*mpj + (n-1)/(n*(n-2))*mip - (n-1)/(n-2)*mpp #partner effects

	am <- matrix(a, nrow(RRMatrix2), ncol(RRMatrix2))
	bm <- t(matrix(b, nrow(RRMatrix2), ncol(RRMatrix2)))
	c <- RRMatrix2 - am - bm - mpp		 # relationship effect
	rownames(c) <- colnames(c) <- rownames(RRMatrix)
	
	# delete all relationship effects which had a NA in the original matrix
	if (imputation) {
		c[imp$NAs==TRUE] <- NA
	}
	


	# actor and partner effects
	self <- FALSE
	
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		self <- TRUE
		self.centered <- attr(RRMatrix, "self.ratings")-mean(attr(RRMatrix, "self.ratings"), na.rm=TRUE)
		
		eff <- data.frame(id = rownames(RRMatrix), actor=a, partner=b, self=self.centered)
		if (!is.null(name)) {colnames(eff)[2:4] <- paste(name, localOptions$suffixes, sep="")}
		attr(eff[,4], "type") <- "self"
		
		## Add selfenhancement index?
		if (index != "") {
			if (match.arg(index, c("enhance"))=="enhance") {
				kwan <- self.centered - a - b
				eff <- data.frame(eff, enhance=kwan)
				if (!is.null(name)) {colnames(eff)[colnames(eff)=="enhance"] <- paste(name, ".enhance", sep="")}
			}
		}
		
	} else {
		eff <- data.frame(id = rownames(RRMatrix), actor=a, partner=b)
		if (!is.null(name)) {colnames(eff)[2:3] <- paste(name, localOptions$suffixes[1:2], sep="")}
	}
	
	attr(eff[,2], "type") <- "actor"
	attr(eff[,3], "type") <- "partner"
	
	
	
	## Relationship effects
	
	effRel <- reshape2::melt(c)
	effRel <- effRel[apply(effRel, 1, function(x) {x[1] != x[2]}),]
	colnames(effRel) <- c("actor.id", "partner.id", "relationship")
	effRel[,1] <- factor(effRel[,1])
	effRel[,2] <- factor(effRel[,2])
	
	# sort relEffects according to dyad
	digits <- floor(log10(n))+1
	effRel$dyad <- factor(apply(effRel, 1, function(x) paste(sort(x[1:2], decreasing=FALSE), collapse="_")))
	effRel$dyad <- factor(effRel$dyad, labels=paste(attr(RRold, "group.id"), "_", f2(1:length(levels(effRel$dyad)), 0, paste("0",digits,sep="")), sep=""))
	effRel <- effRel[,c(1,2,4,3)]
	effRel <- effRel[order(effRel$dyad, effRel$actor.id),]
	
	
	
	## group means added
	
	eff.gm <- eff
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		eff.gm[,2:3] <- eff.gm[,2:3]+mpp
		eff.gm[,4] <- attr(RRMatrix, "self.ratings")
	} else {
		eff.gm[,2:3] <- eff.gm[,2:3]+mpp
	}
	


	## construct return object
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		res <- list(actor = a, partner = b, relationship = c, eff=eff, effRel=effRel, eff.gm=eff.gm, self=self.centered)
	} else {
		res <- list(actor = a, partner = b, relationship = c, eff=eff, effRel=effRel, eff.gm=eff.gm)
	}
	
	attr(res, "self") <- self
	
	return(res)
}




# calculates variance components from a single RR-Matrix
RR.univariate <- function(RRMatrix, na.rm=FALSE, verbose=TRUE, corr.fac="1", index="", varname=NA) {
	
	if (is.null(RRMatrix)) return();
	
	if (nrow(RRMatrix)<4) {
		warning(paste("WARNING: group",attr(RRMatrix, "group.id"),"has 3 or fewer subjects. For calculation of SRM variables minimum group size is 4."), call.=FALSE);
		return();
	}
	
	
	# emit some warnings about missings if there are NAs outside the diagonale
	if ((sum(is.na(RRMatrix)) > nrow(RRMatrix)) & na.rm==FALSE)
		stop("There are NAs outside the diagonale. Calculations are aborted.")
		
	n <- nrow(RRMatrix)	
	n.NA <- sum(is.na(RRMatrix)) - n
	
	
	## Warning if too many missings are present
	if (n.NA > 0 & na.rm==TRUE & verbose==TRUE) {
		wa <- FALSE
		if (n==4 & n.NA > 1) {wa <- TRUE}
		if (n==5 & n.NA > 2) {wa <- TRUE}
		if (n==6 & n.NA > 4) {wa <- TRUE}
		if (n==7 & n.NA > 6) {wa <- TRUE}
		if (n==8 & n.NA > 8) {wa <- TRUE}
		if (n>=10 & (n.NA/(n^2-n)) > .20) {wa <- TRUE}
		
		if (wa==TRUE) {
			warning(paste(attr(RRMatrix, "varname"),": The number of missing values (n.NA=",n.NA,"; ",round((n.NA/(n^2-n))*100, 1),"%) in group ",attr(RRMatrix, "group.id")," exceeds the recommended maximum number of missings according to Schoenbrodt, Back, & Schmukle (in prep.). Estimates might be biased.", sep=""), call.=FALSE)
		}
	}
	
	eff <- RR.effects(RRMatrix, name=attr(RRMatrix, "varname"), na.rm=na.rm, index=index, varname=varname)
	
	A <- sum(eff$actor^2)/(n-1)
	B <- sum(eff$partner^2)/(n-1)
	C <- sum(eff$actor*eff$partner)/(n-1)
	e <- 0.5*(eff$relationship + t(eff$relationship))
	d <- eff$relationship - t(eff$relationship)
	
	if (na.rm==TRUE) {
		
		if (is.na(corr.fac)) {
			corr.fac <- (n*(n-1)) / (n*(n-1) - sum(is.na(e)) + n)
		} else {
			corr.fac <- eval(parse(text=corr.fac))
		}
		
		D <- (sum(e^2, na.rm=TRUE) * corr.fac)  / (((n-1)*(n-2)/2)-1)
		E <- ((sum(d^2,na.rm=TRUE)/2)  * corr.fac) / ((n-1)*(n-2))
	} else {
		D <- sum(e^2, na.rm=TRUE) / (((n-1)*(n-2)/2)-1)
		E <- (sum(d^2,na.rm=TRUE)/2) / ((n-1)*(n-2))
	}
	
	
	scc <- (D+E)/2 #relationship variance
	sccs <- (D-E)/2 #relationship covariance
	sab <- C - (sccs*(n-1))/(n*(n-2)) - scc/(n*(n-2)) #actor-partner covariance
	saa <- A - (scc*(n-1))/(n*(n-2)) - sccs/(n*(n-2)) #actor variance
	sbb <- B - (scc*(n-1))/(n*(n-2)) - sccs/(n*(n-2)) #partner variance
	
	saa2 <- ifelse(saa>=0, saa, NaN)
	sbb2 <- ifelse(sbb>=0, sbb, NaN)
	scc2 <- ifelse(scc>=0, scc, NaN)
	
	raa <- saa2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized actor variance
	rbb <- sbb2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized partner variance
	rcc <- scc2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized relationship variance
	rab <- ifelse(saa>0 & sbb>0,sab/sqrt(saa*sbb),NaN) #actor-partner correlation
	rccs <- sccs/scc2 #relationship correlation
	
	
	w = (n^2 - 3*n + 6) * (n^2 - 3*n + 4)

	sesaa <- sqrt((2*saa^2) / (n+1) + (2*(n^6 - 7*n^5 + 28*n^4 - 66*n^3 + 102*n^2 - 84*n + 32)* scc^2)/ (w*(n+1)*n^2*(n-2)^2)
	        + (2*(n^3-n^2-2*n+16)*(n^2-2*n+2)*sccs^2) / (w*(n+1)*n^2*(n-2)^2)
	        + (4*saa*((n-1) * scc + sccs)) / ((n+1)*n*(n-2))
	        + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2))

	sesbb <- sqrt((2*sbb^2) / (n+1) + (2*(n^6 - 7*n^5 + 28*n^4 - 66*n^3 + 102*n^2 - 84*n + 32)* scc^2)/ (w*(n+1)*n^2*(n-2)^2)
	        + (2*(n^3-n^2-2*n+16)*(n^2-2*n+2)*sccs^2) / (w*(n+1)*n^2*(n-2)^2)
	        + (4*sbb*((n-1) * scc + sccs)) / ((n+1)*n*(n-2))
	        + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2))

	sescc <- sqrt((2*(n^2-3*n+5)* ((scc^2 + sccs^2))) / w + (4*scc*sccs)/ w)

	sesab <- sqrt(((n-3)*sab^2) / ((n+1)*(n-2)) + ((n^6 - 5*n^5 + 19*n^4 - 45*n^3 + 90*n^2 - 96*n + 64)* scc^2)/ (w*(n+1)*n^2*(n-2)^2)
	        + ((n^6 - 7*n^5 + 31*n^4 - 83*n^3 + 150*n^2 - 144*n + 64)* sccs^2)/ (w*(n+1)*n^2*(n-2)^2)
	        + ((n-1)*saa*sbb)/((n+1)*(n-2))
	        + ((n-1)*(saa+sbb)*((n-1)*scc+sccs))/((n+1)*n*(n-2)^2)
	        + (2*(n-3)*sab*(scc+(n-1)*sccs))/((n+1)*n*(n-2)^2)
	        + (4*(n^5-5*n^4+20*n^3-42*n^2+60*n-32)*scc*sccs) / (w*(n+1)*n^2*(n-2)^2))

	sesccs <- sqrt((2*(n^2-3*n+5)*((scc^2+sccs^2)))/w + (4*scc*sccs)/w)


	# error variance is NA if only one group is present
	estimate <- c(saa,sbb,scc,NA,sab,sccs)
	standardized <- clamp(raa,rbb,rcc,NA,rab,rccs)
	
	se <- c(ifelse(estimate[1:3]>=0,c(sesaa,sesbb,sescc),NaN),NA,sesab,sesccs)
	t.value <- estimate/se
	p.value <- 1-pt(abs(t.value), n-1)
	# Kovarianzen werden zweiseitig getestet:
	p.value[4:5] <- p.value[4:5]*2
	
	
	# calculate reliability for actor and partner effects
	rel.a <- saa / (saa + scc*(n-1)/(n*(n-2)) + sccs/(n*(n-2)))
	if (saa < 0) rel.a <- NaN
	rel.b <- sbb / (sbb + scc*(n-1)/(n*(n-2)) + sccs/(n*(n-2)))
	if (sbb < 0) rel.b <- NaN
	
	attr(eff$eff[,2], "reliability") <- rel.a
	attr(eff$eff[,3], "reliability") <- rel.b
	
	# join everything in one dataframe
	univariate <- data.frame(type=unilabels_b, estimate, standardized, se, t.value, p.value)
	
	# if one variance component is below zero: erase covariances
	# erase indices for negative variances
	univariate[1:3,][univariate$estimate[1:3]<0,3:6] <- NaN
	if (saa <= 0 | sbb <= 0) {univariate[5,3:6] <- NaN}
	if (scc <= 0) {univariate[6,3:6] <- NaN}

	res <- list(effects = eff$eff, effectsRel = eff$effRel, effects.gm = eff$eff.gm, varComp = univariate, relMat.av=e, relMat.diff=d, group.size=n, latent=FALSE, anal.type="Univariate analysis of one round robin variable", n.NA = n.NA)
	class(res) <- "RRuni"
	attr(res, "group.size") <- n
	attr(res, "varname") <- attr(RRMatrix, "varname")
	attr(res, "self") <- attr(eff, "self")
	
	# if self ratings are present: add to results object
	#print(attr(RRMatrix, "group.id"))
	dummy <- capture.output(self <- invisible(selfCor(res)))
	if (!is.null(self)) {res[["selfCor"]] <- self}
	
	return(res)
}



# combines two RR-matrices, depending on parameter 'latent'
# latent = TRUE: both matrices are treated as two measures for one underlying construct
# latent = FALSE: both matrices are treated as independent variables
# noCorrection = TRUE: even if univariate estimates are negative, bivariate covariances are NOT set to NA (this is necessary, when the manifest bivariat results are transferred into the bivariate latent analysis, see TAG1)

RR.bivariate <- function(RRMatrix1, RRMatrix2, analysis="manifest", na.rm=FALSE, verbose=TRUE, noCorrection=FALSE, index="", varname=NA) {
	
	if (!(analysis %in% c("latent", "manifest"))) stop("Parameter 'analysis' must either be 'latent' or 'manifest'. Calculations aborted.")
	
	
	dif1 <- union(setdiff(rownames(RRMatrix1), rownames(RRMatrix2)), setdiff(rownames(RRMatrix2), rownames(RRMatrix1)))
	if (length(dif1)>0 & verbose==TRUE) {
		warning(paste(length(dif1),"participant(s) have been excluded from the bivariate/latent analysis due to missing data in one of both variables"), call.=FALSE)
	}
	
	
	# Clean up data: only participants are allowed that are in BOTH matrices
	allparticipants <- intersect(rownames(RRMatrix1), rownames(RRMatrix2))
	
	a1 <- attributes(RRMatrix1)
	a1$self.ratings <- a1$self.ratings[allparticipants]
	a1$dimnames <- NULL
	a2 <- attributes(RRMatrix2)
	a2$self.ratings <- a2$self.ratings[allparticipants]
	a2$dimnames <- NULL
	a1$dim <- rep(length(allparticipants), 2)
	a2$dim <- rep(length(allparticipants), 2)

	RRMatrix1 <- RRMatrix1[allparticipants,allparticipants]
	RRMatrix2 <- RRMatrix2[allparticipants,allparticipants]
	dimn <- rownames(RRMatrix1)
	attributes(RRMatrix1) <- a1
	attributes(RRMatrix2) <- a2
	rownames(RRMatrix1) <- rownames(RRMatrix2) <- colnames(RRMatrix1) <- colnames(RRMatrix2) <- dimn
	
	
	RR.1 <- RR.univariate(RRMatrix1, na.rm, verbose, index=index, varname=varname)
	RR.2 <- RR.univariate(RRMatrix2, na.rm, verbose, index=index, varname=varname)	
	varComp.1 <- RR.1$varComp$estimate
	varComp.2 <- RR.2$varComp$estimate
	n <- nrow(RRMatrix1)
	
	

	#Bivariate Relationships
	A <- sum(RR.1$effects[,2]*RR.2$effects[,2])/(n-1)
	B <- sum(RR.1$effects[,2]*RR.2$effects[,3])/(n-1)
	C <- sum(RR.1$effects[,3]*RR.2$effects[,2])/(n-1)
	D <- sum(RR.1$effects[,3]*RR.2$effects[,3])/(n-1)
	E <- sum(RR.1$relMat.av*RR.2$relMat.av,na.rm=TRUE)/(((n-1)*(n-2)/2)-1)
	F <- (sum(RR.1$relMat.diff*RR.2$relMat.diff,na.rm=TRUE)/2)/((n-1)*(n-2))
	sch <- (E+F)/2 #intrapersonal relationship covariance
	schs <- (E-F)/2 #interpersonal relationship covariance
	saf <- A - (sch*(n-1))/(n*(n-2)) - schs/(n*(n-2)) #actor-actor covariance
	sag <- B - (schs*(n-1))/(n*(n-2)) - sch/(n*(n-2)) #actor-partner covariance
	sbf <- C - (schs*(n-1))/(n*(n-2)) - sch/(n*(n-2)) #partner-actor covariance
	sbg <- D - (sch*(n-1))/(n*(n-2)) - schs/(n*(n-2)) #partner-partner covariance
	
	
	# standardized covariances (=correlations), bivariate case
	#standardized <- clamp(raf,rbg,rag,rbf,rch,rchs)
	w <- getOption("warn")
	options(warn=-1)
		raf <-  saf/(sqrt(varComp.1[1])*sqrt(varComp.2[1])) # bivariate correlations
		rbg <-  sbg/(sqrt(varComp.1[2])*sqrt(varComp.2[2]))
		rag <-  sag/(sqrt(varComp.1[1])*sqrt(varComp.2[2]))
		rbf <-  sbf/(sqrt(varComp.1[2])*sqrt(varComp.2[1]))
		rch <-  sch/(sqrt(varComp.1[3])*sqrt(varComp.2[3]))
		rchs <- schs/(sqrt(varComp.1[3])*sqrt(varComp.2[3]))
	options(warn=w)
	
	
	if (analysis=="latent") {
		stabpervar1 <- saf
		stabtarvar1 <- sbg
		stabrelvar1 <- sch
		
    	stabapcov1 <- (sag + sbf)/2		
    	stabdycov1 <- schs
		unstabper1 <- (varComp.1[1] + varComp.2[1])/2 - saf
		unstabtar1 <- (varComp.1[2] + varComp.2[2])/2 - sbg
		unstabrel1 <- (varComp.1[3] + varComp.2[3]) / 2 - sch
		
		
		saf2 <- max(saf, 0)
		sbg2 <- max(sbg, 0)
		sch2 <- max(sch, 0)
		
		stable1 <- saf2 + sbg2 + sch2
		unstable1 <- max(unstabper1, 0) + max(unstabtar1, 0) + max (unstabrel1, 0) 
		unstable.raw <- sum(unstabper1, unstabtar1, unstabrel1)
		stabler1 <- stable1 / (stable1 + unstable1)
		unstabler1 <- unstable1 / (stable1 + unstable1)
		stabperr1 <- saf2/(stable1 + unstable1)
		stabtarr1 <- sbg2/(stable1+unstable1)
		stabrelr1 <- sch2/(stable1+unstable1)
		stabdycor1 <- ifelse(sch2>0, stabdycov1/sch2, NaN)
		stabapcor1 <- ifelse(saf>0 & sbg>0, stabapcov1 / sqrt(saf*sbg), NaN)
	}
	

	# Standard errors (se) und t-values (t) of bivariate srm-parameters

	coef1 <- n^11-14*n^10+89*n^9-342*n^8+872*n^7-1505*n^6+1698*n^5-1063*n^4+116*n^3+292*n^2-224*n+64
	coef2 <- n^10-12*n^9+66*n^8-227*n^7+534*n^6-857*n^5+883*n^4-416*n^3-148*n^2+224*n-64
	coef3 <- n^10-11*n^9+48*n^8-93*n^7-2*n^6+388*n^5-763*n^4+572*n^3+4*n^2-224*n+64
	coef4 <- n^11-16*n^10+117*n^9-520*n^8+1540*n^7-3083*n^6+3970*n^5-2689*n^4-4*n^3+1212*n^2-544*n-64
	coef5 <- n^10-11*n^9+42*n^8-33*n^7-258*n^6+976*n^5-1453*n^4+788*n^3+348*n^2-544*n-64
	coef6 <- 2*(n^10-14*n^9+92*n^8-383*n^7+1074*n^6-1963*n^5+2101*n^4-752*n^3-780*n^2+544*n+64)
	coef7 <- (n-3)*n*(n^6-9*n^5+35*n^4-75*n^3+76*n^2-12*n-48)

	suppressWarnings(
		sesaf <- sqrt(((n-1)*varComp.1[1]*varComp.2[1]+(n-3)*saf^2)/((n-2)*(n+1))
		         + ((n-1)^2*(varComp.1[1]*varComp.2[3]+varComp.1[3]*varComp.2[1])+(n-1)*(varComp.1[1]*varComp.2[5]+varComp.1[5]*varComp.2[1])+2*(n-3)*saf*((n-1)*sch+schs))/((n-2)^2*n*(n+1))
		         + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+coef3*varComp.1[5]*varComp.2[5]+coef4*sch^2+coef5*schs^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7))
		)

	suppressWarnings(
		sesbg <- sqrt(((n-1)*varComp.1[2]*varComp.2[2]+(n-3)*sbg^2)/((n-2)*(n+1))
		         + ((n-1)^2*(varComp.1[2]*varComp.2[3]+varComp.1[3]*varComp.2[2])+(n-1)*(varComp.1[2]*varComp.2[5]+varComp.1[5]*varComp.2[2])+2*(n-3)*sbg*((n-1)*sch+schs))/((n-2)^2*n*(n+1))
		         + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+coef3*varComp.1[5]*varComp.2[5]+coef4*sch^2+coef5*schs^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7))
		)
		
		
	suppressWarnings(
	sesch <- sqrt(((n^6-9*n^5+32*n^4-57*n^3+43*n^2+6*n-8)*(varComp.1[3]*varComp.2[3]+varComp.1[5]*varComp.2[5])+2*(n^4-6*n^3+3*n^2+18*n-8)*sch*schs)/coef7
	         + ((n^4-6*n^3+11*n^2-6*n+8)*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+(n^6-9*n^5+28*n^4-33*n^3-9*n^2+54*n+8)*(sch^2 + schs^2))/coef7)
			)
			
	suppressWarnings(
	sesag <- sqrt(((n-1)*varComp.1[1]*varComp.2[2]+(n-3)*sag^2)/((n-2)*(n+1))
	         + ((n-1)^2*(varComp.1[1]*varComp.2[3]+varComp.1[3]*varComp.2[2])+(n-1)*(varComp.1[1]*varComp.2[5]+varComp.1[5]*varComp.2[2])+2*(n-3)*sag*((n-1)*schs+sch))/((n-2)^2*n*(n+1))
	         + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+coef3*varComp.1[5]*varComp.2[5]+coef4*schs^2+coef5*sch^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7))
	)

	suppressWarnings(
	sesbf <- sqrt(((n-1)*varComp.1[2]*varComp.2[1]+(n-3)*sbf^2)/((n-2)*(n+1))
	         + ((n-1)^2*(varComp.1[2]*varComp.2[3]+varComp.1[3]*varComp.2[1])+(n-1)*(varComp.1[2]*varComp.2[5]+varComp.1[5]*varComp.2[1])+2*(n-3)*sbf*((n-1)*schs+sch))/((n-2)^2*n*(n+1))
	         + (coef1*varComp.1[3]*varComp.2[3]+coef2*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+coef3*varComp.1[5]*varComp.2[5]+coef4*schs^2+coef5*sch^2+coef6*sch*schs)/((n-2)^3*n^2*(n+1)*coef7))
	)


	suppressWarnings(
	seschs <- sqrt(((n^6-9*n^5+32*n^4-57*n^3+43*n^2+6*n-8)*(varComp.1[3]*varComp.2[3]+varComp.1[5]*varComp.2[5])+2*(n^4-6*n^3+3*n^2+18*n-8)*sch*schs)/coef7
	         + ((n^4-6*n^3+11*n^2-6*n+8)*(varComp.1[3]*varComp.2[5]+varComp.1[5]*varComp.2[3])+(n^6-9*n^5+28*n^4-33*n^3-9*n^2+54*n+8)*(sch^2 + schs^2))/coef7)
	)


	taf <- saf/sesaf
	tbg <- sbg/sesbg
	tch <- sch/sesch
	tag <- sag/sesag
	tbf <- sbf/sesbf
	tchs <- schs/seschs
	
	if (analysis=="latent") {
		sestabpervar1 <- sesaf
		sestabtarvar1 <- sesbg
		sestabrelvar1 <- sesch

		tstabpervar1 <- ifelse(saf>=0, saf/sesaf, NaN)
		tstabtarvar1 <- ifelse(sbg>=0, sbg/sesbg, NaN)
		tstabrelvar1 <- ifelse(sch>=0, sch/sesch, NaN)
	}

	#########################################Result Matrix


# Result Matrix for two independent variables

if (analysis=="manifest") {

	univariate <- list()
	univariate[[1]] <- RR.1
	univariate[[2]] <- RR.2
	
	estimate <- c(saf,sbg,sag,sbf,sch,schs)
	standardized <- clamp(raf,rbg,rag,rbf,rch,rchs)
	
	se <- c(sesaf,sesbg,sesag,sesbf,sesch,seschs)
	t.value <- c(taf,tbg,tag,tbf,tch,tchs)
	pvalues <- (1-pt(abs(t.value), n-1))*2 	# alles Kovarianzen, daher zweiseitig testen!
	bivariate <- data.frame(type=bilabels_bb, estimate, standardized, se, t.value, p.value=pvalues)
	
	if (noCorrection==FALSE) {
		# erase covariances if one variance component is < 0
		if (RR.1$varComp[1,2] <= 0) bivariate[c(1,3),3:6] <- NaN
		if (RR.1$varComp[2,2] <= 0) bivariate[c(2,4),3:6] <- NaN
		if (RR.2$varComp[1,2] <= 0) bivariate[c(1,4),3:6] <- NaN
		if (RR.2$varComp[2,2] <= 0) bivariate[c(2,3),3:6] <- NaN
		if (RR.1$varComp[3,2] <= 0) bivariate[c(5,6),3:6] <- NaN
		if (RR.2$varComp[3,2] <= 0) bivariate[c(5,6),3:6] <- NaN
	}	
	
	res <- list(univariate = univariate, bivariate = bivariate, latent=FALSE, anal.type="Bivariate analysis of two variables, each measured by one round robin variable")
	attr(res, "group.size") <- n
	class(res) <- "RRbi"
	
} else 
{
	# Result matrix for latent analysis 

	unstand <- c(stabpervar1,stabtarvar1,stabrelvar1,unstable1,stabapcov1,stabdycov1)
	stand <- clamp(stabperr1, stabtarr1, stabrelr1, unstabler1,stabapcor1,stabdycor1)
	stand[is.infinite(stand)] <- NaN
	
	# se, t, und p werden aus dem bivariaten Fall uebernommen
	se <- c(sestabpervar1,sestabtarvar1,sestabrelvar1,NA,(sesag+sesbf)/2,seschs)
	tvalues<-c(tstabpervar1,tstabtarvar1,tstabrelvar1,NA,stabapcov1/((sesag+sesbf)/2),tchs)
	pvalues <- (1-pt(abs(tvalues), n-1))
	pvalues[4:5] <- pvalues[4:5]*2
	
	results <- data.frame(type=unilabels_b, estimate=unstand,standardized=stand,se=se,t.value=tvalues, p.value=pvalues)
	
	# erase indices for negative variances
	results[1:3,][results$estimate[1:3]<0,3:6] <- NaN
	
	#-------------------------------
	# calculate reliability for actor and partner effects
	
	unstabdycov1 <- (varComp.1[5] + varComp.2[5]) / 2 - schs
	
	r <- 2 # r = number of replications - in our case, it's always 2
	rel.a <- stabpervar1 / ((stabpervar1 + (unstabper1/r)) + (stabrelvar1+(unstabrel1/r))*(n-1)/(n*(n-2)) + (stabdycov1+(unstabdycov1/r))/(n*(n-2)))
	if (stabpervar1 < 0) rel.a <- NaN
	
	rel.p <- stabtarvar1 / ((stabtarvar1 + (unstabtar1/r)) + (stabrelvar1+(unstabrel1/r))*(n-1)/(n*(n-2)) + (stabdycov1+(unstabdycov1/r))/(n*(n-2)))
	if (stabtarvar1 < 0) rel.p <- NaN
	
	rel.r <- stabrelvar1 / (stabrelvar1+(unstabrel1/r))
	if (stabrelvar1 < 0) rel.r <- NaN
	
	#-------------------------------

	## average effects of both latent indicators
	eff2 <- (as.matrix(RR.1$effects[,-1]) + as.matrix(RR.2$effects[,-1])) / 2
	eff3 <- data.frame(RR.1$effects$id, eff2)
	colnames(eff3) <- colnames(RR.1$effects)
	eff <- eff3
	
	attr(eff[,2], "type") <- "actor"
	attr(eff[,3], "type") <- "partner"
	if (ncol(eff)>=4) attr(eff[,4], "type") <- "self"
	
	eff2 <- (as.matrix(RR.1$effects.gm[,-1]) + as.matrix(RR.2$effects.gm[,-1])) / 2
	eff3 <- data.frame(RR.1$effects.gm$id, eff2)
	colnames(eff3) <- colnames(RR.1$effects.gm)
	eff.gm <- eff3
	attr(eff.gm[,2], "type") <- "actor"
	attr(eff.gm[,3], "type") <- "partner"
	
	
	effRel2 <- merge(RR.1$effectsRel, RR.2$effectsRel, by=c("actor.id", "partner.id", "dyad"))
	effRel2$relationship <- apply(effRel2[,c("relationship.x", "relationship.y")], 1, mean, na.rm=TRUE)
	effRel <- effRel2[,c("actor.id", "partner.id", "dyad", "relationship")]
	effRel <- effRel[order(effRel$dyad),]
	
	
	attr(eff[,grep(localOptions$suffixes[1], colnames(eff), fixed=TRUE)], "reliability") <- rel.a
	attr(eff[,grep(localOptions$suffixes[2], colnames(eff), fixed=TRUE)], "reliability") <- rel.p
	attr(effRel$relationship, "reliability") <- rel.r
	
	res <- list(effects = eff, effects.gm=eff.gm, effectsRel=effRel, varComp=results, unstabdycov1=unstabdycov1, unstabper1=unstabper1, unstabtar1=unstabtar1, unstabrel1=unstabrel1, unstable.raw=unstable.raw, latent=TRUE, anal.type="Latent construct analysis of one construct measured by two round robin variables")
	attr(res, "group.size") <- n
	attr(res, "varname") <- paste(attr(RR.1, "varname"), attr(RR.2, "varname"), sep="/")
	if ((attr(RR.1, "self") == TRUE) & (attr(RR.2, "self") == TRUE)) {attr(res, "self") <- TRUE} else {attr(res, "self") <- FALSE}
	class(res) <- "RRuni"
}

	
	return(res)
}




# helper function
ifg <- function(g) {
	if (is.null(g)) {
		return("");
	} else {
		return(paste("|",g));
	}
}


checkVar <- function(x, minVar=0) {
	if (is.null(minVar)) return(FALSE)
	if (is.na(minVar)) return(FALSE)
	if (is.null(x)) return(TRUE)
	if (is.nan(x)) return(TRUE)
	if (is.na(x)) return(TRUE)
	if (x < minVar) return(TRUE)
	return(FALSE)
}


# Wrapper function: depending on parameters, different results are calculated:
#' @export
#' @importFrom reshape2 melt
#' @importFrom reshape2 acast
#' @importFrom plyr ldply
#' @importFrom plyr laply
RR <- function(formula, data, na.rm=FALSE, minData = 1, verbose=TRUE, g.id=NULL, index="", exclude.ids="", varname=NA, minVar=localOptions$minVar, ...) {

	extra <- list(...)

	# set default
	analysis <- "manifest"
	RRMatrix2 <- RRMatrix3 <- RRMatrix4 <- NULL

	# transform long format (formula) to quadratic matrices
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");

	#remove spaces from formula
	f1 <- formula
	lhs <- strsplit(gsub(" ","",as.character(f1)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(f1)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]

	actor.id <- rhs[1]
	partner.id <- rhs[2]
	if (length(rhs)>=3) {group.id <- rhs[3]} else {group.id=NULL}


	# if a grouping factor is provided: forward to RR.multi
	if (!is.null(group.id)) {
		res <- RR.multi(f1, data=data, na.rm=na.rm, verbose=verbose, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, ...)
		
		if (!is.null(res$univariate)) {
			
			# bivariate case
			
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				if (checkVar(res$univariate[[1]]$varComp[1, 3], minVar)) {
					res$univariate[[1]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[1]]$varComp[2, 3], minVar)) {
					res$univariate[[1]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
				if (checkVar(res$univariate[[2]]$varComp[1, 3], minVar)) {
					res$univariate[[2]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[2]]$varComp[2, 3], minVar)) {
					res$univariate[[2]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
			}
			
			
		} else {
				# if variance < minVar: set effects to NA
				if (!is.na(minVar)) {
					if (checkVar(res$varComp[1, 3], minVar)) {
						res$effects[,3][1:nrow(res$effects)] <- NA
						res$effects.gm[,3][1:nrow(res$effects)] <- NA
					}
					if (checkVar(res$varComp[2, 3], minVar)) {
						res$effects[,4][1:nrow(res$effects)] <- NA
						res$effects.gm[,4][1:nrow(res$effects)] <- NA

					}
				}
		}
		
		res$minVar <- minVar
		for (g in 1:length(res$groups)) res$groups[[g]]$minVar <- minVar
		
		return(res)
		
	}



	# univariater Fall:
	if (length(lhs)==1) {
		lhs1 <- strsplit(lhs, "/")[[1]]
	
		# manifester vs. latenter Fall
		if (length(lhs1)==1) {
			RRMatrix1 <- long2matrix(formula(paste(lhs1,"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			analysis <- "manifest"
			
		} else if (length(lhs1)==2) {
			
			# if (!is.null(extra[["bistyle"]])) {v2 <- FALSE} else {v2=verbose}
			
			ex1 <- attr(long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2 <- attr(long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex3 <- Reduce(union, list(exclude.ids, ex1, ex2))
			
			
			RRMatrix1 <- long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, bistyle=TRUE)[[1]]
			RRMatrix2 <- long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, ...)[[1]]
			analysis <- "latent"
		}
	
	} else 
	# bivariater Fall
	if (length(lhs)==2) {
	
		lhs1 <- strsplit(lhs[1], "/")[[1]]
	
		# manifester vs. latenter Fall
		if (length(lhs1)==1) {
			
			RRMatrix1 <- long2matrix(formula(paste(lhs[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			
			RRMatrix2 <- long2matrix(formula(paste(lhs[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			analysis <- "manifest"
		} else if (length(lhs1)==2) {
			lhs2 <- strsplit(lhs[2], "/")[[1]]
			
			
			# exclude participants
			
			ex1a <- attr(long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex1b <- attr(long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2a <- attr(long2matrix(formula(paste(lhs2[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2b <- attr(long2matrix(formula(paste(lhs2[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex3 <- Reduce(union, list(exclude.ids, ex1a, ex1b, ex2a, ex2b))
			
			RRMatrix1 <- long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, bistyle=TRUE)[[1]]
			RRMatrix2 <- long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			RRMatrix3 <- long2matrix(formula(paste(lhs2[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			RRMatrix4 <- long2matrix(formula(paste(lhs2[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			analysis <- "latent"
		}
	} else {stop("Error: Unknown term in formula.")}
	



## if all RRMatrices are NULL: stop
if (is.null(RRMatrix1) & is.null(RRMatrix2) & is.null(RRMatrix3) & is.null(RRMatrix4)) {
	return(NULL);
}
	
# depending on given parameters different results are calculated

#-----------------------------
#- One group

	if (is.null(RRMatrix2) & is.null(RRMatrix3) & is.null(RRMatrix4)) {
		if (analysis=="latent") {
			return(NULL);
			# warning("Warning: analysis='latent' only is valid, when two different RRMatrices for one latent construct are given")
		}
		
		res <- RR.univariate(RRMatrix1, na.rm, verbose, index=index, varname=varname)
		
		# if variance < minVar: set effects to NA
		if (!is.na(minVar)) {
			if (checkVar(res$varComp[1, 3], minVar)) {
				res$effects[,2][1:nrow(res$effects)] <- NA
				res$effects.gm[,2][1:nrow(res$effects)] <- NA
			}
			if (checkVar(res$varComp[2, 3], minVar)) {
				res$effects[,3][1:nrow(res$effects)] <- NA
				res$effects.gm[,3][1:nrow(res$effects)] <- NA
				
			}
		}
		
		res$minVar <- minVar
		
		return(res)
	}
	
#-----------------------------
#- Two groups, independent or latent constructs

	if (is.null(RRMatrix3) & is.null(RRMatrix4)) {
		
		if (is.null(RRMatrix1) | is.null(RRMatrix2)) {
			# if (verbose) {warning("Error: One of both round robin matrices has to few participants!", call.=FALSE)}
			return();
		}
		res <- RR.bivariate(RRMatrix1, RRMatrix2, analysis=analysis, na.rm=na.rm, verbose=verbose, index=index, varname=varname)
		
		if (!is.null(res$univariate)) {
			
			# bivariate case
			
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				if (checkVar(res$univariate[[1]]$varComp[1, 3], minVar)) {
					res$univariate[[1]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[1]]$varComp[2, 3], minVar)) {
					res$univariate[[1]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
				if (checkVar(res$univariate[[2]]$varComp[1, 3], minVar)) {
					res$univariate[[2]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[2]]$varComp[2, 3], minVar)) {
					res$univariate[[2]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
			}
			
			
		} else {
				# if variance < minVar: set effects to NA
				if (!is.na(minVar)) {
					if (checkVar(res$varComp[1, 3], minVar)) {
						res$effects[,3][1:nrow(res$effects)] <- NA
						res$effects.gm[,3][1:nrow(res$effects)] <- NA
					}
					if (checkVar(res$varComp[2, 3], minVar)) {
						res$effects[,4][1:nrow(res$effects)] <- NA
						res$effects.gm[,4][1:nrow(res$effects)] <- NA

					}
				}
		}
		
		res$minVar <- minVar
		return(res);
	}
	
	
#-----------------------------
#- four groups: two constructs measured with each two variables

	if (!is.null(RRMatrix1) & !is.null(RRMatrix2) & !is.null(RRMatrix3) & !is.null(RRMatrix4)) {
		
		
		# calculate latent effects for both constructs
		lat1.full <- RR.bivariate(RRMatrix1, RRMatrix2, analysis="latent", na.rm=na.rm, verbose=FALSE, index=index, varname=varname)
		lat2.full <- RR.bivariate(RRMatrix3, RRMatrix4, analysis="latent", na.rm=na.rm, verbose=FALSE, index=index, varname=varname)
		lat1 <- lat1.full$varComp$estimate
		lat2 <- lat2.full$varComp$estimate
				
		
		# calculate new raw data: mean of both indicators; a manifest bivariate analysis is conducted on the mean variable
		# --> all calculations are correct, except the standardization --> this has to be done on the latent data
		RR12 <- (RRMatrix1+RRMatrix2)/2
		RR34 <- (RRMatrix3+RRMatrix4)/2
		
		#TAG1 <-- this is a bookmark, do not remove
		bivariate <- RR.bivariate(RR12, RR34, na.rm=na.rm, verbose=FALSE, noCorrection=TRUE)$bivariate
		
		#Estimation of bivariate relations on construct level
		
		w <- getOption("warn")
		options(warn=-1)
		denom <- c(
			sqrt(lat1[1]*lat2[1]),
			sqrt(lat1[2]*lat2[2]),
			sqrt(lat1[1]*lat2[2]),
			sqrt(lat1[2]*lat2[1]),
			sqrt (lat1[3]*lat2[3]),
			sqrt (lat1[3]*lat2[3])
		)
		options(warn=w)
		
		bivariate$standardized <- clamp(bivariate$estimate / denom)

		# erase covariances if one variance component is < 0
		if (lat1[1] <= 0) bivariate[c(1,3),3:6] <- NaN
		if (lat1[2] <= 0) bivariate[c(2,4),3:6] <- NaN
		if (lat2[1] <= 0) bivariate[c(1,4),3:6] <- NaN
		if (lat2[2] <= 0) bivariate[c(2,3),3:6] <- NaN
		if (lat1[3] <= 0) bivariate[c(5,6),3:6] <- NaN
		if (lat2[3] <= 0) bivariate[c(5,6),3:6] <- NaN
		
		
		univariate <- list()
		univariate[[1]] <- lat1.full
		univariate[[2]] <- lat2.full
		
		grandres <- list(univariate = univariate, bivariate = bivariate)
		class(grandres) <- "RR"
		grandres$anal.type <- "Bivariate analysis of two constructs, each measured by two round robin variables"
		attr(grandres, "group.size") <- nrow(RRMatrix2)
		
		
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				if (checkVar(grandres$univariate[[1]]$varComp[1, 3], minVar)) {
					grandres$univariate[[1]]$effects[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[1]]$effects.gm[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(grandres$univariate[[1]]$varComp[2, 3], minVar)) {
					grandres$univariate[[1]]$effects[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[1]]$effects.gm[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA

				}
				if (checkVar(grandres$univariate[[2]]$varComp[1, 3], minVar)) {
					grandres$univariate[[2]]$effects[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[2]]$effects.gm[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(grandres$univariate[[2]]$varComp[2, 3], minVar)) {
					grandres$univariate[[2]]$effects[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[2]]$effects.gm[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA

				}
			}

		grandres$minVar <- minVar
		return(grandres)		
	} else {
		# warning("Error: One of the round robin matrices has to few participants!", call.=FALSE)
	}
	
	return(NULL);
}


 


# weighted variance, inspired by a function from Gavin Simpson on R-Help
var.wt <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    return((sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2)))
}




weighted.t.test <- function(x, w, mu, conf.level = 0.95, alternative="two.sided", na.rm=TRUE) {
	
	if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
	
	if (na.rm) { 
		w <- w[i <- !is.na(x)] 
		x <- x[i] 
	}
	
	# to achieve consistent behavior in loops, return NA-structure in case of complete missings
	if (sum(is.na(x)) == length(x)) return(list(estimate=NA, se=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
	
	# if only one value is present: this is the best estimate, no significance test provided
	if (sum(!is.na(x)) == 1) {
		warning("Warning weighted.t.test: only one value provided; this value is returned without test of significance!", call.=FALSE)
		return(list(estimate=x[which(!is.na(x))], se=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
	}
	
	x.w <- weighted.mean(x,w, na.rm=na.rm)
	var.w <- var.wt(x,w, na.rm=na.rm)
	df <- length(x)-1
	t.value <- sqrt(length(x))*((x.w-mu)/sqrt(var.w))
	se <- sqrt(var.w)/sqrt(length(x))
	
	if (alternative == "less") {
		pval <- pt(t.value, df)
		cint <- c(-Inf, x.w + se*qt(conf.level, df) )
    }
    else if (alternative == "greater") {
		pval <- pt(t.value, df, lower.tail = FALSE)
		cint <- c(x.w - se * qt(conf.level, df), Inf)
    }
    else {
		pval <- 2 * pt(-abs(t.value), df)
		alpha <- 1 - conf.level
        cint <- x.w + se*qt(1 - alpha/2, df)*c(-1,1)
    }
	
	names(t.value) <- "t"
	return(list(estimate=x.w, se=se, conf.int=cint, statistic=t.value, df=df, p.value=pval))
}



# helper-functions
posOrNA <- function(x) {
	return(ifelse(x>=0, x, NA))
}



# uni1, uni2: univariate Analysen der beiden Konstrukte (Daten werden zum Standardisieren der bivariaten Koeffizienten gebraucht)

getWTest <- function(RR0, res1, typ="univariate", uni1=NA, uni2=NA, unstable=NA) {
	
	if (is.null(RR0)) return();
	
	if (typ=="univariate") {
		if (length(RR0$univariate)==2) {varComp <- RR0$univariate[[1]]$varComp} else {varComp <- RR0$varComp}
		varComp$p.value <- NA
		varComp$estimate <- NA

		for (v in names(table(res1$type))) {
			w.t <- weighted.t.test(res1$estimate[res1$type == v], res1$group.size[res1$type == v]-1, mu=0)
	
			varComp[varComp$type==v,]$estimate <- w.t$estimate
			varComp[varComp$type==v,]$se <- w.t$se
			varComp[varComp$type==v,]$t.value <- w.t$statistic
			varComp[varComp$type==v,]$p.value <- w.t$p.value
		}
		varComp$p.value[1:4] <- varComp$p.value[1:4] / 2
		# Varianzen nur einseitig testen (Voreinstellung bei weighted.t.test ist zweiseitig)

		
		
		# unstable variance im latent bivariaten Fall wird von aussen in die Funktion gegeben
		if (!is.na(unstable)) {
			varComp[4,2] <- unstable
		}


		#standardized coefficients need special treatment ...
		varComp[1,3] <- posOrNA(varComp[1,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[2,3] <- posOrNA(varComp[2,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[3,3] <- posOrNA(varComp[3,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[4,3] <- posOrNA(varComp[4,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		w <- getOption("warn")
		options(warn=-1)
			varComp[5,3] <- varComp[5,2]/ sqrt(varComp[1,2]*varComp[2,2])
			varComp[6,3] <- varComp[6,2]/ varComp[3,2]
		options(warn=w)
		
		varComp[,3] <- clamp(varComp[,3])
		
		# variance below zero: erase all other indices
		bz <- which(varComp[1:3, 2]<0)
		if (length(bz)>0) {
			varComp[bz, 3:6] <- NaN
			if (varComp[1,2]<0 | varComp[2,2]<0) varComp[5,3:6] <- NaN
		}
		
		# error variance: set se, t, p to NA (instead of NaN)
		varComp[4,4:6] <- NA
	
		return(varComp)
	}
	
	if (typ=="bivariate") {
		
		
		bivariate <- RR0$bivariate
		bivariate$p.value <- NA

		for (v in names(table(res1$type))) {
			w.t <- weighted.t.test(res1$estimate[res1$type == v], res1$group.size[res1$type == v]-1, mu=0)
	
			bivariate[bivariate$type==v,]$estimate <- w.t$estimate
			bivariate[bivariate$type==v,]$se <- w.t$se
			bivariate[bivariate$type==v,]$t.value <- w.t$statistic
			bivariate[bivariate$type==v,]$p.value <- w.t$p.value
		}

		#standardized coefficients need special treatment ...
		w <- getOption("warn")
		options(warn=-1)
			bivariate[1,3] <- bivariate[1,2]/ sqrt(uni1[1,2]*uni2[1,2])
			bivariate[2,3] <- bivariate[2,2]/ sqrt(uni1[2,2]*uni2[2,2])
			bivariate[3,3] <- bivariate[3,2]/ sqrt(uni1[1,2]*uni2[2,2])
			bivariate[4,3] <- bivariate[4,2]/ sqrt(uni1[2,2]*uni2[1,2])
			bivariate[5,3] <- bivariate[5,2]/ sqrt(uni1[3,2]*uni2[3,2])
			bivariate[6,3] <- bivariate[6,2]/ sqrt(uni1[3,2]*uni2[3,2])
		options(warn=w)
		
		
		# erase covariances if one variance component is < 0
		if (uni1[1,2] <= 0) bivariate[c(1,3),3:6] <- NaN
		if (uni2[1,2] <= 0) bivariate[c(1,4),3:6] <- NaN
		if (uni1[2,2] <= 0) bivariate[c(2,4),3:6] <- NaN
		if (uni2[2,2] <= 0) bivariate[c(2,3),3:6] <- NaN

		bivariate[,3] <- clamp(bivariate[,3])
	
		return(bivariate)
	}
}





RR.multi.uni <- function(formule, data, na.rm=FALSE, verbose=TRUE, index="", minData=1, exclude.ids="", varname=NA, ...) {

	# this function needs data in long format ...
	extra <- list(...)
		
	# parse formula
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");
	
	# f1 = formula without grouping factor
	fstring <- paste(as.character(formule)[c(2,1,3)], collapse=" ")
	f0 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"\\|", perl=TRUE)[[1]]
	f1 <- formula(f0[1])
	f3 <- strsplit(strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][1], "+", fixed=TRUE)[[1]]
	group.id <- f0[2]

	mode <- ifelse(length(f3)==2,"bi","uni")
		
	res <- data.frame()
	res.bi <- data.frame()
	g.uni <- list()
	saa <- sbb <- scc <- sccs <- n.m <- c()
	undc1 <- unp1 <- unt1 <- unr1 <- un.raw  <- c()
	
	self <- FALSE	# are self ratings present?
	
	for (g in names(table(data[,group.id]))) {
		
		#print(g)
		RR0 <- RR(f1, data=data[data[,group.id] == g,], verbose=verbose, na.rm=na.rm, g.id=group.id, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, minVar=NA, ...)
		
		#print(str(RR0))

		if (is.null(RR0)) {next;} else {RR1 <- RR0}
		if (attr(RR0, "self") == TRUE) {self <- TRUE}
		g.id <- g
		
		RR0$group.id <- g.id

		# if (RR0$latent==FALSE) {
# 			RR0$effects.gm$group.id <- g.id
# 		} else {
# 			eff.gm <- list(relationship=NA)
# 		}
		
		g.uni[[g]] <- RR0
		
		if (RR0$latent==FALSE) {
			saa <- c(saa, RR0$varComp[1,2])
			sbb <- c(sbb, RR0$varComp[2,2])
			scc <- c(scc, RR0$varComp[3,2])
			sccs <- c(sccs, RR0$varComp[6,2])
		} else {
			undc1 <- c(undc1, RR0$unstabdycov1)
			unp1 <- c(unp1, RR0$unstabper1)
			unt1 <- c(unt1, RR0$unstabtar1)
			unr1 <- c(unr1, RR0$unstabrel1)
		}
		
		n.m <- c(n.m, attr(RR0, "group.size"))
		
		u1 <- RR0$varComp
		u1$variable <- 1

		u1$group.size <-  attr(RR0, "group.size")
		u1$group.id <- g.id
		
		res <- rbind(res, u1)
		
	}

	# aus der liste die Effekte extrahieren und zusammenfuegen
	effect <- ldply(g.uni, function(x) {return(x$effects)})
	effect[,1:2] <- effect[,2:1]
	colnames(effect)[1:2] <- c("id", "group.id")
	effect[,1] <- factor(effect[,1])
	effect[,2] <- factor(effect[,2])
	
	type <- c("actor", "partner", "self")
	for (ty in 3:ncol(effect)) {
		attr(effect[,ty], "type") <- type[ty-2]
	}

	eff.gm <- ldply(g.uni, function(x) {return(x$effects.gm)})
	eff.gm[,1:2] <- eff.gm[,2:1]
	colnames(eff.gm)[1:2] <- c("id", "group.id")
	eff.gm[,1] <- factor(eff.gm[,1])
	eff.gm[,2] <- factor(eff.gm[,2])
	
	effectRel <- ldply(g.uni, function(x) {return(x$effectsRel)})
	colnames(effectRel)[1:3] <- c("group.id", all.vars(f1)[2:3])
	
	effectRel[,1] <- factor(effectRel[,1])
	effectRel[,2] <- factor(effectRel[,2])
	effectRel[,3] <- factor(effectRel[,3])
	effectRel[,4] <- factor(effectRel[,4])

	# im latenten Fall: die Error variance erst am Ende berechnen (d.h., alle error componenten ueber alle Gruppen mitteln, die unter NUll auf Null setzen, dann addieren)
	
	unstable.raw.m <- NA
	if (RR1$latent==TRUE) {
		unstable.raw.m <- max(weighted.mean(unp1, n.m), 0) + max(weighted.mean(unt1, n.m), 0) + max(weighted.mean(unr1, n.m), 0)
	}	

	if (length(effect) == 0) {
		effect <- data.frame(actor=NA, partner=NA, relationship=NA)
	}
	# get weighted variance components
	varComp <- getWTest(RR1, res, unstable=ifelse(is.null(unstable.raw.m), NULL, unstable.raw.m))


	# calculate reliability for actor and partner effects, and variance of group means
	group.var <- NA
	
	if (!is.null(n.m)) {

		n <- mean(n.m)
		
		if (RR1$latent==FALSE) {
			saa.m <- weighted.mean(saa, n.m-1)
			sbb.m <- weighted.mean(sbb, n.m-1)
			scc.m <- weighted.mean(scc, n.m-1)
			sccs.m <- weighted.mean(sccs, n.m-1)
	
			rel.a <- saa.m / (saa.m + scc.m*(n-1)/(n*(n-2)) + sccs.m/(n*(n-2)))
			if (saa.m < 0) rel.a <- NaN
			rel.p <- sbb.m / (sbb.m + scc.m*(n-1)/(n*(n-2)) + sccs.m/(n*(n-2)))
			if (sbb.m < 0) rel.p <- NaN
		} else {
			
			
			unp1.m <- weighted.mean(unp1, n.m-1, na.rm=TRUE)
			unt1.m <- weighted.mean(unt1, n.m-1, na.rm=TRUE)
			unr1.m <- weighted.mean(unr1, n.m-1, na.rm=TRUE)
			undc1.m <- weighted.mean(undc1, n.m-1, na.rm=TRUE)
			
			
			r <- 2 # r = number of replications - in our case, it's always 2
			rel.a <- varComp$estimate[1] / ((varComp$estimate[1] + (unp1.m/r)) + (varComp$estimate[3]+(unr1.m/r))*(n-1)/(n*(n-2)) + (varComp$estimate[6]+(undc1.m/r))/(n*(n-2)))
			if (varComp$estimate[1] < 0) rel.a <- NaN

			
			rel.p <- varComp$estimate[2] / ((varComp$estimate[2] + (unt1.m/r)) + (varComp$estimate[3]+(unr1.m/r))*(n-1)/(n*(n-2)) + (varComp$estimate[6]+(undc1.m/r))/(n*(n-2)))
			if (varComp$estimate[2] < 0) rel.p <- NaN

			rel.r <- varComp$estimate[3] / (varComp$estimate[3]+(unr1.m/r))
			if (varComp$estimate[3] < 0) rel.r <- NaN
			
			attr(effectRel$relationship, "reliability") <- clamp(rel.r)
		}

		attr(effect[,grep(localOptions$suffixes[1], colnames(effect), fixed=TRUE)], "reliability") <- clamp(rel.a)
		attr(effect[,grep(localOptions$suffixes[2], colnames(effect), fixed=TRUE)], "reliability") <- clamp(rel.p)
		
		# group mean & group variance		
		group.means <- tapply(data[,all.vars(f1)[1]], data[,group.id], mean, na.rm=TRUE)
		group.var <- var(group.means) - ((varComp$estimate[1] + varComp$estimate[2] + 2*varComp$estimate[5])/n + (varComp$estimate[3]+varComp$estimate[6])/(n*(n-1)))
		
		
	}	
	
	#------------------------------------
	#--  Variance Components: calculate weighted mean and weighted between groups t-test
	#------------------------------------
	
	anal.type <- paste(RR1$anal.type, " in multiple groups",sep="")
	
	if (!is.null(varComp)) {
	
		res2 <- list(effects = effect, effectsRel = effectRel, effects.gm = eff.gm, varComp = varComp, groups = g.uni, varComp.groups=res, group.var=group.var, anal.type=anal.type)
		class(res2) <- "RRmulti"
		attr(res2, "varname") <- attr(g.uni[[1]], "varname")
		attr(res2, "self") <- self
		
		# # noch rausfinden, welche Teilnehmer ausgeschlossen wurden
		# l1 <- long2matrix(formule, data, verbose=FALSE)
		# attr(res2, "excluded.participants") <- attr(l1, "excluded.participants")
		# attr(res2, "excluded.groups") <- attr(l1, "excluded.groups")
		
		return(res2)
	} else {return();}
	
}





RR.multi <- function(formule, data, na.rm=FALSE, verbose=TRUE, index="", minData=1, exclude.ids="", varname=NA, ...) {

	# this function needs data in long format ...
	extra <- list(...)
	
	# parse formula
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");
	
	# f1 = formula without grouping factor
	fstring <- paste(as.character(formule)[c(2,1,3)], collapse=" ")
	f0 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"\\|", perl=TRUE)[[1]]
	f1 <- formula(f0[1])
	group.id <- f0[2]
	
	f3 <- strsplit(strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][1], "+", fixed=TRUE)[[1]]
	f4 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][2]
	
	if (sum(grepl("/", f3, fixed=TRUE))>1) {analysis <- "latent"} else {analysis <- "manifest"}

	# Vorlage fuer die Varianzkomponente erstellen
	df <- data[data[,group.id]==data[1,group.id],]
	mode <- ifelse(length(f3)==2,"bi","uni")

	if (mode=="uni") return(RR.multi.uni(formule, data, na.rm, verbose, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, ...))

	# ... ansonsten bi-mode durchfuehren
	
	# find out, which participants are excluded and exclude them from all variables
	if (analysis=="manifest") {
		ex1 <- attr(long2matrix(formula(paste(f3[1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2 <- attr(long2matrix(formula(paste(f3[2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex3 <- Reduce(union, list(ex1, ex2, exclude.ids))
	} else {
		ex1a <- attr(long2matrix(formula(paste(strsplit(f3[1], "/", fixed=TRUE)[[1]][1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex1b <- attr(long2matrix(formula(paste(strsplit(f3[1], "/", fixed=TRUE)[[1]][2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2a <- attr(long2matrix(formula(paste(strsplit(f3[2], "/", fixed=TRUE)[[1]][1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2b <- attr(long2matrix(formula(paste(strsplit(f3[2], "/", fixed=TRUE)[[1]][2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex3 <- Reduce(union, list(ex1a, ex1b, ex2a, ex2b, exclude.ids))
	}
	
	V1 <- RR.multi.uni(formula(paste(f3[1], "~", f4)), data, na.rm, verbose=verbose, index=index, minData=minData, exclude.ids=ex3, bistyle=TRUE)
	V2 <- RR.multi.uni(formula(paste(f3[2], "~", f4)), data, na.rm, verbose=FALSE, index=index, minData=minData, exclude.ids=ex3)
	V2$varComp.groups$variable <- 2
		
	res <- data.frame()
	res.bi <- data.frame()
	bi.groups <- list()
	
	for (g in names(table(data[,group.id]))) {
		
			RR0 <- RR(f1, data=data[data[,group.id] == g,], verbose=FALSE, na.rm=na.rm, minData=minData, exclude.ids=ex3, minVar=NA)
			
			if (is.null(RR0)) {next;} else
			{RR1 <- bi.groups[[g]]  <- RR0}
			
			if (!is.null(RR1$bivariate)) {
				res.bi <- rbind(res.bi, data.frame(RR1$bivariate, group.size=attr(RR1, "group.size"), group=g))
			}
		
	}

	
	anal.type <- paste(RR1$anal.type, "in multiple groups")
	
	bivariate <- getWTest(RR1, res.bi, typ="bivariate", V1$varComp, V2$varComp)
	
	res <- list(univariate = list(V1, V2), bivariate = bivariate, anal.type=anal.type, groups=bi.groups)
	class(res) <- "RRmulti"

	return(res)
}





# Workaround: When @export does not recognize that this is a S3-method, you need the extra @method statement
#' @method print RRmulti
#' @export 
print.RRmulti <- function(x, ...) {
	print.RR(x, ...)
}

#' @method print RRuni
#' @export
print.RRuni <- function(x, ...) {
	print.RR(x, ...)
}

#' @method print RRbi
#' @export
print.RRbi <- function(x, ...) {
	print.RR(x, ...)
}

# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=3, prepoint=0) {
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", ".", sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", ".", sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}


# x muss hier direkt auf das univariate-Objekt verweisen
print.uni <- function(x, ..., measure=NA, digits=3, r.names=NULL, minVar=0) {
	
	# print descriptivers for multi group
	if (length(x$groups) > 1) {
		groupsizes <- laply(x$groups, function(y) return(attr(y, "group.size")))
		av.groupsize <- round(mean(groupsizes), 2)
		
		 print(paste("Group descriptives: n = ",length(x$groups),"; average group size = ",av.groupsize, "; range: ", min(groupsizes), "-", max(groupsizes)))
	}
	
	uni <- round(x$varComp[,2:ncol(x$varComp)], digits)	
	
	if (checkVar(uni[1, 2], minVar)) {uni[5, 2:5] <- NA}
	if (checkVar(uni[2, 2], minVar)) {uni[5, 2:5] <- NA}
	if (checkVar(uni[3, 2], minVar)) {uni[6, 2:5] <- NA}
		
	if (is.na(measure)) {
		measure <- localOptions$style
	} else {
		measure <- match.arg(measure, c("behavior", "perception", "metaperception"))
	}
	
	if (!is.null(r.names)) {rownames(uni) <- r.names} else {
		if (measure == "behavior") rownames(uni) <- unilabels_b
		if (measure == "perception") rownames(uni) <- unilabels_p
		if (measure == "metaperception") {
			warning("Warning: the current RR-object only consists of a single variable. Labels for metaperception are only provided when two variables are specified.", call.=FALSE)
			rownames(uni) <- unilabels_b
		}
	}
	
	print(uni)
	
	
	# Actor effect reliability
	if (!is.null(x$effects[,grep(localOptions$suffixes[1], colnames(x$effects), fixed=TRUE)])) print(paste(role[[measure]][1], "effect reliability:",f2(attr(x$effects[,grep(localOptions$suffixes[1], colnames(x$effects), fixed=TRUE)], "reliability"), 3)))
	
	# Partner effect reliability
	if (!is.null(x$effects[,grep(localOptions$suffixes[2], colnames(x$effects), fixed=TRUE)])) print(paste(role[[measure]][2], "effect reliability:",f2(attr(x$effects[,grep(localOptions$suffixes[2], colnames(x$effects), fixed=TRUE)], "reliability"), 3)))
	
	# Relationship effect reliability
	if (!is.null(attr(x$effectsRel$relationship, "reliability"))) print(paste(role[[measure]][3], "effect reliability:",f2(attr(x$effectsRel$relationship, "reliability"), 3)))
	
	selfCor(x, measure=measure)
}



# Here the default print method for RR-objects gets overwritten, so that 
# the information in the RR-class is displayed in a convenient way
print.RR <- function(x, ..., measure1=NA, measure2=NA, digits=3, measure=NULL) {
	
	if (is.na(measure1)) {
		measure1 <- localOptions$style
	} else {
		measure1 <- match.arg(measure1, c("behavior", "perception", "metaperception"))
	}
	if (is.na(measure2)) {
		measure2 <- measure1
	} else {
		measure2 <- match.arg(measure2, c("behavior", "perception", "metaperception"))
	}
	
	print("Round-Robin object ('RR'), calculated by TripleR")
	print(x$anal.type)
	
	
	if (!is.null(measure)) {measure1 <- measure}
	
	# bivariate case
	if (length(x$univariate) == 2) {
		
		uni <- lapply(x$univariate, function(x) return(x))
		bi <- round(x$bivariate[,2:ncol(x$bivariate)], digits)
				
		# Erase bivariate correlations for variance components < minVar
		if (checkVar(uni[[1]]$varComp[1, 2], x$minVar)) {bi[c(1,3), 2:5] <- NA}
		if (checkVar(uni[[1]]$varComp[2, 2], x$minVar)) {bi[c(2,4), 2:5] <- NA}
		if (checkVar(uni[[1]]$varComp[3, 2], x$minVar)) {bi[c(5,6), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[1, 2], x$minVar)) {bi[c(1,4), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[2, 2], x$minVar)) {bi[c(2,3), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[3, 2], x$minVar)) {bi[c(5,6), 2:5] <- NA}
		                                     
		r.names1 <- r.names2 <- NULL
		if (measure1 == "behavior" & measure2 == "behavior") {
			rownames(bi) <- bilabels_bb
		} else
   		if (measure1 == "behavior" & measure2 == "perception") {
				rownames(bi) <- bilabels_bp
		} else
		if (measure1 == "perception" & measure2 == "behavior") {
				rownames(bi) <- bilabels_pb
		} else
		if (measure1 == "perception" & measure2 == "perception") {
				rownames(bi) <- bilabels_pp
		} else
		if (measure1 == "perception" & measure2 == "metaperception") {
			r.names1 <- unilabels_b_meta1
			r.names2 <- unilabels_b_meta2
			rownames(bi) <- bilabels_meta
		} else {
			stop("This combination of measurement labels does not fit.")
		}
		print(paste("Univariate analyses for:", attr(uni[[1]], "varname")))
		print.uni(uni[[1]], measure=measure1, r.names=r.names1, minVar=x$minVar)
		cat("\n")
		print(paste("Univariate analyses for:", attr(uni[[2]], "varname")))
		print.uni(uni[[2]], measure=measure2, r.names=r.names2, minVar=x$minVar)
		cat("\n")
		print("Bivariate analyses:")
		
		print(bi)
		
		if (length(uni[[1]]$groups) != length(uni[[2]]$groups)) {
			warning(paste("Note: Univariate analyses of both variables are based on different numbers of groups. Bivariate analyses therefore are based on the common groups of both variables (n=",min(length(uni[[1]]$groups), length(uni[[2]]$groups)),")", sep=""), call.=FALSE)
		}
		
	} else
	
	# univariate case
	{
		print(paste("Univariate analyses for:", attr(x, "varname")))
		print.uni(x, measure=measure1, minVar=x$minVar)
	}
	
	
}

