library(buildmer)
library(testthat)
test_that('converged',{
	skip_on_cran()
	library(lme4)
	good1 <- lm(Reaction ~ Days,sleepstudy)
	good2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
	bad <- suppressWarnings(lmer(Reaction ~ Days + (Days|Subject),sleepstudy,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=1))))
	expect_equal(sapply(list(good1,good2,bad),converged),c(TRUE,TRUE,FALSE))
})

test_that('converged should reject singular fits',{
	skip_on_cran()
	library(lme4)
	re <- re2uncorr(following ~ vowel*timepoint + (vowel*timepoint||participant),vowels)
	model <- buildmer(re$formula,family=binomial,data=re$data,buildmerControl=list(direction='order',args=list(verbose=2,control=lmerControl(optimizer='bobyqa'))))
	tab <- tabulate.formula(formula(model@model))
	buildmer:::testthat.compare.df(tab,'converged')
})
