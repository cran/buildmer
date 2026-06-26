library(buildmer)
library(testthat)
test_that('re2uncorr',{
	re <- re2uncorr(f1 ~ vowel*timepoint*following + (vowel*timepoint*following|participant) + (timepoint|word),vowels)
	buildmer:::testthat.compare.df(re$termlist,'re2uncorr')
})
