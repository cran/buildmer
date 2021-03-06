get2LL <- function (m) as.numeric(-2*stats::logLik(m))
getdf  <- function (m) attr(stats::logLik(m),'df') #note: for GAMs this automatically uses the sum of edf2
getdev <- function (m) {
	if (all(c('deviance','null.deviance') %in% names(m))) return(1-m$deviance/m$null.deviance)
	if (!is.null(summary(m)$r.squared)) return(1-summary(m)$r.squared)
	ff <- fitted(m)
	if (attr(terms(formula(m)),'intercept')) ff <- ff - mean(ff)
	rr <- resid(m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}

crit.AIC <- function (p,ref,alt) if (is.null(ref)) stats::AIC(alt) else stats::AIC(alt) - stats::AIC(ref)
crit.BIC <- function (p,ref,alt) if (is.null(ref)) stats::BIC(alt) else stats::BIC(alt) - stats::BIC(ref)
crit.F <- function (p,ref,alt) {
	if (!inherits(alt,'gam')) {
		stop('crit.F currently only works with gam/bam models')
	}
	r2_alt  <- summary(alt,re.test=FALSE)$r.sq
	ddf_alt <- alt$df.residual
	ndf_alt <- sum(alt$edf2)
	if (is.null(ref)) {
		r2_ref  <- 0
		ddf_ref <- nobs(alt)
		ndf_ref <- 0
	} else {
		if (!inherits(ref,'gam')) {
			stop('crit.F currently only works with gam/bam models')
		}
		r2_ref  <- summary(ref,re.test=FALSE)$r.sq
		ddf_ref <- ref$df.residual
		ndf_ref <- sum(ref$edf2)
	}
	if (is.null(r2_alt) || is.null(r2_ref)) {
		stop('r^2 not available for this family, cannot compute the F criterion!')
	}
	Fval <- (r2_alt - r2_ref) / ((1 - r2_alt) / (ddf_alt))
	ndf  <- ndf_alt - ndf_ref
	ddf  <- ddf_alt
	if (is.na(Fval)) {
		return(log(1))
	}
	if (Fval <= 0 || ndf <= 0) {
		return(log1p(abs(Fval))) #gives the order step some idea of which model is the least unhelpful
	}
	if (alt$scale.estimated) {
		pf(Fval,ndf,ddf,lower.tail=FALSE,log.p=TRUE)
	} else {
		pchisq(ndf*Fval,ndf,lower.tail=FALSE,log.p=TRUE)
	}
}
crit.LRT <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- getdf(alt)
		f1   <- ~0
	} else {
		chLL <- get2LL(ref) - get2LL(alt)
		chdf <- getdf(alt)  - getdf(ref)
		f1   <- formula(ref)
	}
	if (chdf <= 0) {
		return(0)
	}

	# If the two models differ in lme4-style random effects, we need to correct the p-value
	# We cannot use the formula stored in p here because that may still go through re2mgcv
	f2 <- formula(alt)
	tab1 <- tabulate.formula(f1)
	tab2 <- tabulate.formula(f2)
	fe.same <- isTRUE(all.equal(tab1[ is.na(tab1$grouping),],tab2[ is.na(tab2$grouping),]))
	re.same <- isTRUE(all.equal(tab1[!is.na(tab1$grouping),],tab2[!is.na(tab2$grouping),]))
	if (fe.same && !re.same) {
		# Stram & Lee (1994): mixture of chisq(chdf) and chisq(chdf-1)
		p1 <- stats::pchisq(chLL,chdf  ,lower.tail=FALSE,log.p=TRUE) - log(2)
		p2 <- stats::pchisq(chLL,chdf-1,lower.tail=FALSE,log.p=TRUE) - log(2)
		log(exp(p1) + exp(p2))
	} else {
		stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
	}
}
crit.2LL <- function (p,ref,alt) if (is.null(ref)) get2LL(alt) else get2LL(alt) - get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (p,ref,alt) if (is.null(ref)) getdev(alt) else getdev(alt) - getdev(ref)
crit.deviance <- crit.devexp

elim.AIC <- function (diff) diff > -.001
elim.BIC <- elim.AIC
elim.F   <- function (logp) exp(logp) >= .05
elim.LRT <- elim.F
elim.2LL <- elim.AIC
elim.LL  <- elim.AIC
elim.devexp <- elim.AIC
elim.deviance <- elim.AIC
