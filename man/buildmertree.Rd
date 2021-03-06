% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildmer.R
\name{buildmertree}
\alias{buildmertree}
\title{Use \code{buildmer} to perform stepwise elimination for \code{lmertree()} and \code{glmertree()} models from package \code{glmertree}}
\usage{
buildmertree(
  formula,
  data = NULL,
  family = gaussian(),
  buildmerControl = buildmerControl(crit = "AIC"),
  ...
)
}
\arguments{
\item{formula}{Either a \code{glmertree} formula, looking like \code{dep ~ left | middle | right} where the \code{middle} part is an \code{lme4}-style random-effects specification, or an ordinary formula (or buildmer term list thereof) specifying only the dependent variable and the fixed and random effects for the regression part. In the latter case, the additional argument \code{partitioning} must be specified as a one-sided formula containing the partitioning part of the model.}

\item{data}{See the general documentation under \code{\link{buildmer-package}}}

\item{family}{See the general documentation under \code{\link{buildmer-package}}}

\item{buildmerControl}{Control arguments for buildmer --- see the general documentation under \code{\link{buildmerControl}}}

\item{...}{Additional options to be passed to \code{lmertree} or \code{glmertree}. (They will also be passed to \code{(g)lmtree} in so far as they're applicable. The single exception is the \code{control} argument, which is assumed to be meant only for \code{(g)lmertree} and not for \code{(g)lmtree}, and will \emph{not} be passed on to \code{(g)lmtree}). For backward-compatibility reasons, will also accept buildmer control parameters, although those specified in \code{buildmerControl} will take precedence.}
}
\description{
Use \code{buildmer} to perform stepwise elimination for \code{lmertree()} and \code{glmertree()} models from package \code{glmertree}
}
\details{
Note that the likelihood-ratio test is not available for \code{glmertree} models, as it cannot be assured that the models being compared are nested. The default is thus to use AIC.
In the generalized case or when testing many partitioning variables, it is recommended to pass \code{joint=FALSE}, as this results in a dramatic speed gain and reduces the odds of the final \code{glmer} model failing to converge or converging singularly.
}
\examples{
if (requireNamespace('glmertree')) {
	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
		buildmerControl=buildmerControl(crit='LL',direction='order'),
        data=lme4::sleepstudy)
\donttest{
	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
		buildmerControl=buildmerControl(crit='LL',direction='order'),
	        data=lme4::sleepstudy,family=Gamma(link=identity),joint=FALSE)
}
}
}
\seealso{
\code{\link{buildmer-package}}
}
