#' Adds terms to a formula
#' @param formula The formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @examples
#' library(buildmer)
#' form <- Reaction ~ Days + (1|Subject)
#' add.terms(form,'Days|Subject')
#' add.terms(form,'(0+Days|Subject)')
#' add.terms(form,c('many','more|terms','to|terms','(be|added)','to|test'))
#' @export
add.terms <- function(formula,add) {
	dep       <- if (length(formula) < 3) '' else as.character(formula[2])
	terms     <- terms(formula,keep.order=TRUE)
	intercept <- attr(terms,'intercept')
	offset    <- attr(terms,'offset')
	vars      <- attr(terms,'variables')
	terms     <- attr(terms,'term.labels')
	terms     <- c(intercept,terms)
	if (!is.null(offset)) {
		offset.term <- as.character(list(vars[[1+offset]]))
		terms <- c(offset.term,terms)
	}

	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	if (length(random.terms)) {
		random.terms <- sapply(random.terms,function(x) if (mkTerm(x)[[1]] != '(') paste0('(',x,')') else x)
	}

	for (term in add) {
		if (is.random.term(term)) {
			bar <- mkTerm(term)
			if (bar[[1]] == '(') {
				# independent term: tack it on at the end
				random.terms <- c(random.terms,term)
				next
			}
			bar.grouping <- as.character(bar[3])
			bar.terms <- bar[[2]]
			# Find suitable terms for 'intruding', i.e.: can we add the term requested to a pre-existing term group?
			suitable <- if (length(random.terms)) which(sapply(random.terms,function(term) mkTerm(term)[[2]][[3]] == bar.grouping)) else NULL
			if (length(suitable)) {
				random.terms[[suitable[1]]] <- sapply(random.terms[[suitable[1]]],function(x) {
					bar <- mkTerm(x)[[2]]
					grouping <- as.character(bar[3])
					terms <- as.character(bar[2])
					form <- mkForm(terms)
					terms <- terms(form,keep.order=TRUE)
					intercept <- attr(terms,'intercept')
					terms <- attr(terms,'term.labels')
					terms <- c(terms,bar.terms)
					if (intercept) terms <- c('1',terms)
					else terms[[1]] <- paste0('0 + ',terms[[1]])
					terms <- paste(terms,collapse=' + ')
					paste0('(',terms,' | ',grouping,')')
				})
			} else {
				#still have to tack it on at the end in the end...
				term <- paste(bar.terms,collapse=' + ')
				if (!any(bar.terms == '1')) term <- paste0('0 + ',term) #for some reason, "!'1' %in% bar.terms" complains about requiring vector arguments... well duh?
				random.terms <- c(random.terms,paste0('(',term,' | ',bar.grouping,')'))
			}
		} else fixed.terms <- c(fixed.terms,term)
	}
	terms <- c(fixed.terms,random.terms)
	stats::as.formula(paste0(dep,'~',paste0(terms,collapse='+')),environment(formula))
}

#' Converts a buildmer term list into a proper model formula
#' @param dep The dependent variable.
#' @param terms The term list.
#' @param env The environment of the formula to return.
#' @return A formula.
#' @examples
#' library(buildmer)
#' form1 <- Reaction ~ Days + (Days|Subject)
#' terms <- tabulate.formula(form1)
#' form2 <- build.formula(dep='Reaction',terms)
#' 
#' # check that the two formulas give the same results
#' library(lme4)
#' check <- function(f) resid(lmer(f,sleepstudy))
#' all.equal(check(form1),check(form2))
#' 
#' # can also do double bars now
#' form1 <- Reaction ~ Days + (Days||Subject)
#' terms <- tabulate.formula(form1)
#' form2 <- build.formula(dep='Reaction',terms)
#' all.equal(check(form1),check(form2))
#' @export
build.formula <- function(dep,terms,env=parent.frame()) {
	fixed.intercept <- is.na(terms$grouping) & terms$term == '1'
	if (any(fixed.intercept)) {
		form <- stats::as.formula(paste(dep,'~1'))
		terms <- terms[!fixed.intercept,]
	} else  form <- stats::as.formula(paste(dep,'~0'))
	while (nrow(terms)) {
		# we can't use a simple for loop: the data frame will mutate in-place when we encounter grouping factors
		if (is.na(terms[1,'index'])) {
			cur <- terms[1,'term']
			terms <- terms[-1,]
		} else {
			ix <- terms[1,'index']
			ii <- !is.na(terms$index) & terms$index == ix
			cur <- terms[ii,]
			termlist <- cur$term
			if (!'1' %in% termlist) termlist <- c('0',termlist)
			termlist <- paste0(termlist,collapse='+')
			cur <- paste0('(',paste0(termlist,'|',unique(cur$grouping)),')')
			terms <- terms[!ii,]
		}
		form <- add.terms(form,cur)
	}
	environment(form) <- env
	form
}

#' Tests a model for convergence
#' @param model The model object to test.
#' @param singular.ok A logical indicating whether singular fits are accepted as `converged' or not. Relevant only for lme4 models.
#' @param grad.tol The tolerance to use for checking the gradient. This is currently only used by mgcv, glmmTMB, and clm(m) models.
#' @param hess.tol The tolerance to use for checking the Hessian for negative eigenvalues. This is currently only used by mgcv, glmmTMB, and cl(m)m models.
#' @return Logical indicating whether the model converged.
#' @examples
#' library(buildmer)
#' library(lme4)
#' good1 <- lm(Reaction ~ Days,sleepstudy)
#' good2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
#' bad <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,control=lmerControl(
#'             optimizer='bobyqa',optCtrl=list(maxfun=1)))
#' sapply(list(good1,good2,bad),converged)
#' @export
converged <- function(model,singular.ok=FALSE,grad.tol=if (inherits(model,'bam')) 10 else .1,hess.tol=if (inherits(model,'bam')) 10 else .01) {
	grad.tol <- eval(eval(grad.tol)) #the second eval evals `if`,
	hess.tol <- eval(eval(hess.tol)) #the first eval evals `quote` from buildmerControl
	setattr <- function(x,msg,err=NULL) {
		if (!is.null(err)) msg <- paste0(msg,' (',err,')')
		attr(x,'reason') <- msg
		x
	}
	failure <- function(msg,err=NULL) setattr(FALSE,msg,err)
	success <- function(msg,err=NULL) setattr(TRUE,msg,err)
	#eigen sometimes returns complex vectors; https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2013-November/000719.html
	eigval  <- function(...) suppressWarnings(as.numeric(eigen(...)$values))
	if (inherits(model,'try-error')) return(failure(model))
	if (inherits(model,'buildmer')) return(converged(model@model))
	if (inherits(model,'gam')) {
		if (!is.null(model$mgcv.conv)) {
			if (is.logical(model$mgcv.conv)) {
				if (!model$mgcv.conv) return(failure('mgcv reports nonconvergence (mgcv.conv)',model$mgcv.conv))
			} else {
				if (!is.null(model$mgcv.conv$fully.converged) && !model$mgcv.conv$fully.converged) return(failure('mgcv reports nonconvergence (mgcv.conv$fully.converged)',model$mgcv.conv$fully.converged))
				if (model$mgcv.conv$message != 'full convergence') return(failure('mgcv reports nonconvergence (mgcv.conv$message)',model$mgcv.conv$fully.converged))
			}
		}
		if (!is.null(model$converged) && !model$converged) return(failure('mgcv reports nonconvergence (converged)',model$converged))
		if (is.null(model$outer.info)) {
			if (is.null(model$mgcv.conv)) return(success('No smoothing-parameter selection necessary'))
			if ((err <- model$mgcv.conv$rms.grad) > grad.tol) return(failure(paste0('mgcv reports absolute gradient containing values >',grad.tol),err))
			if (!model$mgcv.conv$hess.pos.def) return(failure('mgcv reports non-positive-definite Hessian'))
			return(success('All buildmer checks passed (gam with PQL)'))
		} else {
			if (!is.null(model$outer.info$conv) && (err <- model$outer.info$conv) != 'full convergence') return(failure('mgcv reports nonconvergence (outer.info$conv)',err))
			if (!is.null(model$outer.info$grad) && (err <- max(abs(model$outer.info$grad))) > grad.tol)  return(failure(paste0('Absolute gradient contains values >',grad.tol),err))
			if (!is.null(model$outer.info$hess)) {
				err <- try(min(eigval(model$outer.info$hess)),silent=TRUE)
				if (inherits(err,'try-error')) return(failure('Eigenvalue decomposition of Hessian failed',err))
				if (err < -hess.tol) return(failure(paste0('Hessian contains negative eigenvalues <',-hess.tol),err))
			}
			return(success('All buildmer checks passed (gam with outer iteration)'))
		}
	}
	if (inherits(model,'merMod')) {
		if ((err <- model@optinfo$conv$opt) != 0) return(failure('Optimizer reports not having finished',err))
		if (lme4::isSingular(model) && !singular.ok) {
			return(failure('Singular fit'))
		}
		if (length(model@optinfo$conv$lme4)) {
			if (is.null(model@optinfo$conv$lme4$code) && !singular.ok) return(failure('Singular fit (due to absence of optinfo convergence code)'))
			if (any((err <- model@optinfo$conv$lme4$code) != 0)) return(failure('lme4 reports not having converged',err))
		}
		if (is.null(model@optinfo$derivs)) return(success('No derivative information available -- succeeding by default (dangerous!)'))
		grad <- model@optinfo$derivs$gradient
		hess <- model@optinfo$derivs$Hessian
		scaled.grad <- try(solve(hess,grad),silent=TRUE)
		if (inherits(scaled.grad,'try-error')) return(failure('Manual gradient checks were unable to compute scaled gradient',ev))
		if ((err <- max(pmin(abs(scaled.grad),abs(grad)))) > grad.tol) return(failure(paste0('Absolute gradient contains values >',grad.tol),err))
		err <- try(min(eigval(hess)),silent=TRUE)
		if (inherits(err,'try-error')) return(failure('Eigenvalue decomposition of Hessian failed',err))
		if (err < -hess.tol) return(failure(paste0('Hessian contains negative eigenvalues <',-hess.tol),err))
		return(success('All buildmer checks passed (merMod)'))
	}
	if (inherits(model,'glmmTMB')) {
		if (!is.null(model$fit$convergence) && (err <- model$fit$convergence) != 0) return(failure('glmmTMB reports nonconvergence',err))
		if (!is.null(model$sdr$pdHess)) {
			if (!model$sdr$pdHess) return(failure('glmmTMB reports non-positive-definite Hessian'))
			if (sum(dim(model$sdr$cov.fixed))) {
				if ((err <- max(abs(model$sdr$gradient.fixed))) > grad.tol) return(failure(paste0('Absolute gradient contains values >',grad.tol),err))
				ev <- try(1/eigval(model$sdr$cov.fixed),silent=TRUE)
				if (inherits(ev,'try-error')) return(failure('Eigenvalue decomposition of Hessian failed',ev))
				if ((err <- min(ev)) < -hess.tol) return(failure(paste0('Hessian contains negative eigenvalues <',-hess.tol),err))
			}
		}
		return(success('All buildmer checks passed (glmmTMB)'))
	}
	if (inherits(model,'nnet'))   return(if (model$convergence == 0) success('All buildmer checks passed (nnet)'  ) else failure('nnet reports nonconvergence',model$convergence))
	if (inherits(model,'MixMod')) return(if (model$converged) success('All buildmer checks passed (MixMod)') else failure('GLMMadaptive reports nonconvergence'))
	if (inherits(model,'clm') || inherits(model,'clmm')) {
		if (inherits(model,'clm')) {
			if (model$convergence$code != 0) return(failure('clm reports nonconvergence',model$messages))
			if (model$maxGradient > grad.tol) return(failure(paste0('Absolute gradient contains values >',grad.tol),model$maxGradient))
		}
		if (inherits(model,'clmm')) {
			if (model$optRes$convergence != 0) return(failure('clmm reports nonconvergence',model$optRes$message))
			if ((err <- max(abs(model$gradient))) > grad.tol) return(failure(paste0('Absolute gradient contains values >',grad.tol),err))
		}
		ev <- try(eigval(model$Hessian),silent=TRUE)
		if (inherits(ev,'try-error')) return(failure('Eigenvalue decomposition of Hessian failed',ev))
		if ((err <- min(ev)) < -hess.tol) return(failure(paste0('Hessian contains negative eigenvalues <',-hess.tol),err))
		return(success('All buildmer checks passed (clm/clmm)'))
	}
	if (inherits(model,'lme')) {
		if (!singular.ok && inherits(model$apVar,'character')) {
			if (model$apVar == 'Approximate variance-covariance matrix not available') {
				stop('model$apVar is not available; refit model with apVar control option set to TRUE to enable convergence checking')
			}
			return(failure('lme reports effectively singular convergence',model$apVar))
		}
		return(success('All buildmer checks passed (lme)'))
	}
	if (inherits(model,'glmertree')) return(converged(model$glmer))
	if (inherits(model,'lmertree'))  return(converged(model$lmer))
	success('No checks needed or known for this model type',class(model))
}

#' Sets up deviation contrasts
#' 
#' This function is an alias to \code{\link{contr.deviation}} intended to allow the \code{C} command to work for deviation-coding the same way as it does for other, base-R, contrasts.
#' 
#' @param n See \code{\link{contr.sum}}.
#' @param contrasts See \code{\link{contr.sum}}.
#' @param sparse See \code{\link{contr.sum}}.
#' @importFrom stats contr.sum
#' @return A matrix of deviation-coded contrasts.
#' @seealso \code{\link{contr.deviation}}, \code{\link{contr.sum}}
#' @examples
#' iris$Species <- C(iris$Species,deviation)
#' @export
deviation <- function(n,contrasts=TRUE,sparse=FALSE) {
	contr.deviation(n,contrasts,sparse)
}

#' Sets up deviation contrasts
#' 
#' Deviation-coding is a contrast-coding scheme that compares each level to the average of the other levels. It is often confused with sum-coding, which compares each level to the average of \emph{all} levels. In the binary case, deviation coding works out to \( (+0.5,-0.5) \); in the ternary case and beyond, the coding is a little bit more involved.
#' 
#' \code{contr.deviation} starts by calling \code{\link{contr.sum}} and hence accepts the same arguments as it does.
#' @param n See \code{\link{contr.sum}}.
#' @param contrasts See \code{\link{contr.sum}}.
#' @param sparse See \code{\link{contr.sum}}.
#' @importFrom stats contr.sum
#' @return A matrix of deviation-coded contrasts.
#' @seealso \code{\link{contr.sum}}
#' @examples
#' iris$Species <- contr.deviation(iris$Species)
#' @export
contr.deviation <- function(n,contrasts=TRUE,sparse=FALSE) {
	c  <- contr.sum(n,contrasts,sparse)
	nr <- nrow(c)
	c[c == 0] <- -1
	c[c == 1] <- nr - 1
	c / nr
}

#' Converts lme4 random-effect terms to mgcv 're' smooths
#' @param formula The lme4 formula.
#' @param data The data.
#' @param drop Logical indicating whether constant, non-intercept columns should be dropped. Default \code{TRUE}. A warning is issued if a column needed to be dropped. Note that repeated intercept columns are silently merged without a warning.
#' @return A list containing a new formula (has a dependent variable but loses track of what terms belong together in a single block), a buildmer term list (lacks the dependent variable but represents what terms belong together), and a new data set.
#' @details
#' To ensure that the various buildmer functions (e.g. buildgam, buildbam) correctly apply marginality constraints between parametric terms and random-effect terms, the parametric terms are renamed using the same name-mangling scheme that is automatically applied to the random effects. To ensure that the mangled names will be grouped together, instead of passing the returned \code{formula}, use the returned \code{termlist} with a \code{dep} argument. See the description of the \code{dep} argument in \code{\link{buildmerControl}}.
#' @examples
#' library(buildmer)
#' re <- re2mgcv(temp ~ angle + (1|replicate) + (1|recipe),lme4::cake)
#' model <- buildgam(re$formula,re$data)
#' # note: the below does NOT work, as the dependent variable is looked up in the data by name!
#' \dontshow{if (FALSE)}
#' re <- re2mgcv(log(Reaction) ~ Days + (Days|Subject),lme4::sleepstudy)
#' @export
re2mgcv <- function(formula,data,drop=TRUE) {
	re2mgcv.internal(formula,data,drop,FALSE)
}

#' Converts lme4 random-effect terms to uncorrelated lme4 random-effect terms, including factor random slopes
#' 
#' lme4 and similar packages (e.g. glmmTMB) have an issue where correlations between \emph{factor} random slopes are not dropped, even if using the double-bar syntax or \code{\link{diag}} to diagonalize the random-effects formula. \code{re2uncorr} works around this issue by explicitly converting factor random slopes to numerics, which do not suffer from this problem.
#' @param formula The lme4 formula.
#' @param data The data.
#' @param drop Logical indicating whether constant, non-intercept columns should be dropped. Default \code{TRUE}. A warning is issued if a column needed to be dropped. Note that repeated intercept columns are silently merged without a warning.
#' @return A list containing a new formula (has a dependent variable but loses track of what terms belong together in a single block), a buildmer term list (lacks the dependent variable but represents what terms belong together), and a new data set.
#' @details
#' To ensure that the various buildmer functions (e.g. buildmer, buildglmmTMB) correctly apply marginality constraints between fixed effects and random-effect terms, the fixed effects are renamed using the same name-mangling scheme that is automatically applied to the random effects. To ensure that the mangled names will be grouped together, instead of passing the returned \code{formula}, use the returned \code{termlist} with a \code{dep} argument. See the description of the \code{dep} argument in \code{\link{buildmerControl}}.
#' @examples
#' library(buildmer)
#' re <- re2uncorr(f1 ~ vowel*timepoint*following +
#' 	(vowel*timepoint*following|participant) + (timepoint|word),vowels)
#' \dontshow{if (FALSE)}
#' model <- buildmer(re$formula,re$data)
#' @export
re2uncorr <- function(formula,data,drop=TRUE) {
	re2mgcv.internal(formula,data,drop,TRUE)
}

#' Removes terms from a formula
#' @param formula The formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `(term|group)' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param check A logical indicating whether effects should be checked for marginality. If \code{TRUE} (default), effects will not be removed if doing so would violate marginality. Setting \code{check} to \code{FALSE} will remove terms unconditionally.
#' @examples
#' library(buildmer)
#' remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)')
#' # illustration of the marginality checking mechanism:
#' # this refuses to remove the term:
#' remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)')
#' # so does this, because marginality is checked before removal:
#' remove.terms(Reaction ~ Days + (Days|Subject),c('(Days|Subject)','(1|Subject)'))
#' # but it works with check=FALSE
#' remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)',check=FALSE)
#' @export
remove.terms <- function(formula,remove,check=TRUE) {
	marginality.ok <- function(remove,have) {
		forbidden <- if (!all(have == '1')) '1' else NULL
		for (x in have) {
			x.star <- if (has.smooth.terms(stats::as.formula(paste0('~',x)))) paste(unpack.smooth.terms(x),collapse='*') else gsub(':','*',x) #replace any interaction by the star operator, which will cause as.formula to pull in all lower-order terms necessary without any more work from us!
			partterms <- attr(terms(stats::as.formula(paste0('~',x.star))),'term.labels')
			forbidden <- c(forbidden,partterms[partterms != x])
		}
		!remove %in% forbidden
	}
	unwrap.terms <- function(terms,inner=FALSE,intercept=FALSE) {
		form <- stats::as.formula(paste0('~',terms))
		terms <- terms(form,keep.order=TRUE)
		if (intercept) intercept <- attr(terms,'intercept')
		if (inner) return(terms[[2]])
		terms <- attr(terms,'term.labels')
		if (intercept) terms <- c('1',terms)
		terms
	}

	dep <- as.character(formula[2])
	terms <- terms(formula,keep.order=TRUE)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)

	# Build lists to check which terms are currently present.
	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- decompose.random.terms(random.terms)

	# Protect the intercept: do not remove the fixed intercept if we have a random intercept.
	# This is handled here because the intercept is handled by a separate 'intercept' variable, rather than being in the 'remove' vector.
	if ('1' %in% remove) {
	        if (!check || !'1' %in% unlist(random.terms)) {
			intercept <- 0
		}
	}
	remove <- remove[remove != '1']

	# Prepare the final list of terms for which removal was requested
	if (!length(remove)) {
		remove <- '0'
	}
	remove <- stats::as.formula(paste0('~',paste(remove,collapse='+')))
	remove.random <- get.random.list(remove)
	remove.fixed <- attr(terms(remove),'term.labels')
	remove.fixed <- Filter(Negate(is.random.term),remove.fixed)

	if (check) {
		# Do not remove fixed effects if they have corresponding random effects
		remove.fixed <- remove.fixed[!remove.fixed %in% unlist(random.terms)]

		# Do not remove effects participating in higher-order interactions
		remove.fixed <- remove.fixed[marginality.ok(remove.fixed,fixed.terms)]
		if (length(remove.random)) for (i in 1:length(remove.random)) {
			nm <- names(remove.random)[i]
			have <- unlist(random.terms[names(random.terms) == nm])
			remove.random[[i]] <- remove.random[[i]][marginality.ok(remove.random[[i]],have)]
		}
	}

	# Perform actual removal
	terms <- lapply(terms,function(term) {
		if (is.random.term(term)) {
			terms <- decompose.random.terms(term)
			terms <- lapply(1:length(terms),function(i) {
				g <- names(terms)[i]
				terms <- terms[[i]]
				terms <- terms[!terms %in% remove.random[[g]]]
				if (!length(terms)) return(NULL)
				if (!'1' %in% terms) terms <- c('0',terms)
				paste0('(',paste0(terms,collapse=' + '),'|',g,')')
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) return(NULL)
			paste0(terms,collapse=' + ')
		} else {
			if (term %in% remove.fixed) return(NULL)
			term
		}
	})
	terms <- Filter(Negate(is.null),terms)

	# Wrap up
	if (length(terms)) {
		if (!intercept) {
			terms <- c('0',terms)
		}
		return(stats::as.formula(paste0(dep,'~',paste(terms,collapse='+')),environment(formula)))
	}
	stats::as.formula(if (intercept) '~1' else '~0',environment(formula))
}

#' Parses a formula into a buildmer term list
#' @param formula A formula.
#' @param group A character vector of regular expressions. Terms matching the same regular expression are assigned the same block, and will be evaluated together in buildmer functions.
#' @return A buildmer term list, which is just a normal data frame.
#' @examples
#' form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following +
#'              ((vowel1+vowel2+vowel3+vowel4)*timepoint*following|participant) + (timepoint|word))
#' tabulate.formula(form)
#' tabulate.formula(form,group='vowel[1-4]')
#' @seealso buildmer-package
#' @export
tabulate.formula <- function(formula,group=NULL) {
	mkGroups <- function(t) {
		for (x in group) {
			t <- gsub(x,x,t,perl=TRUE)
		}
		t
	}

	dep       <- as.character(formula[2])
	terms     <- terms(formula,keep.order=TRUE)
	intercept <- attr(terms,'intercept')
	offset    <- attr(terms,'offset')
	vars      <- attr(terms,'variables')
	terms     <- attr(terms,'term.labels')
	if (intercept) {
		terms <- c('1',terms)
	}
	if (!is.null(offset)) {
		offset.term <- as.character(list(vars[[1+offset]]))
		terms <- c(offset.term,terms)
	}
	if (!length(terms)) {
		return(data.frame(index=character(),grouping=character(),term=character(),code=character(),block=character(),stringsAsFactors=FALSE))
	}

	# Build lists to check which terms are currently present.
	fixed.terms  <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- decompose.random.terms(random.terms)

	terms <- lapply(1:length(terms),function(i) {
		term <- terms[[i]]
		if (is.random.term(term)) {
			terms <- decompose.random.terms(term)
			terms <- lapply(1:length(terms),function(j) {
				g <- names(terms)[j]
				terms <- terms[[j]]
				if (!length(terms)) {
					return(NULL)
				}
				ix <- paste(i,j)
				data.frame(index=ix,grouping=g,term=terms,code=paste(ix,g,terms),block=paste(NA,g,mkGroups(terms)),stringsAsFactors=FALSE)
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) {
				return(NULL)
			}
			do.call(rbind,terms)
		} else {
			data.frame(index=NA,grouping=NA,term=term,code=term,block=paste(NA,NA,mkGroups(term)),stringsAsFactors=FALSE)
		}
	})
	terms <- Filter(Negate(is.null),terms)
	tab <- do.call(rbind,terms)
	environment(tab) <- environment(formula)
	tab
}

#' Generates an LRT elimination function with custom alpha level
#' 
#' The \code{elim} argument in \code{\link{buildmerControl}} can take any user-specified elimination function. \code{LRTalpha} generates such a function that uses the likelihood-ratio test, based on a user-specified alpha level. (For the default alpha of .05, one can also simply specify the string \code{'LRT'} or the function \code{buildmer:::elim.LRT}).
#' @param alpha The alpha level for the likelihood-ratio test.
#' @seealso \code{\link{buildmerControl}}
#' @export
LRTalpha <- function(alpha) {
	if (alpha <= 0 || alpha >= 1) {
		stop("'alpha' should be in (0,1)")
	}
	ret <- function(logp) exp(logp) >= alpha
	attr(ret,'elim.name') <- 'LRT'
	ret
}
