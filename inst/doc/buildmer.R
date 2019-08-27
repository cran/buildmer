## ------------------------------------------------------------------------
library(buildmer)
head(vowels)

## ------------------------------------------------------------------------
f <- f1 ~ vowel*timepoint*following * neighborhood*information*stress + 
	 (vowel*timepoint*following * neighborhood+information+stress | participant) +
	 (timepoint | word)

## ------------------------------------------------------------------------
f <- f1 ~ vowel*timepoint*following +
	 (vowel*timepoint*following | participant) +
	 (timepoint | word)

## ----eval=F--------------------------------------------------------------
#  library(lme4)
#  m <- buildmer(f,data=vowels,direction='order',control=lmerControl(optimizer='bobyqa'))

## ----echo=F--------------------------------------------------------------
cat('Determining predictor order
Currently evaluating LRT for: vowel, timepoint, following
Fitting as (g)lm: f1 ~ vowel
Fitting as (g)lm: f1 ~ timepoint
Fitting as (g)lm: f1 ~ following
Updating formula: f1 ~ following
Currently evaluating LRT for: vowel, timepoint
Fitting as (g)lm: f1 ~ following + vowel
Fitting as (g)lm: f1 ~ following + timepoint
Updating formula: f1 ~ following + vowel
Currently evaluating LRT for: timepoint, vowel:following
[...]
Currently evaluating LRT for: vowel | participant
Fitting via lme4, with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following + vowel | participant)
boundary (singular) fit: see ?isSingular
None of the models converged - giving up ordering attempt.')

## ----include=F-----------------------------------------------------------
library(lme4)
#hack for consistency with actual output without actually fitting the model every time I change something in the vignette
m <- buildmer:::mkBuildmer(model=list(formula=(function () as.formula('f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following + timepoint:following | participant) + (1 + timepoint | word)',.GlobalEnv))()))

## ------------------------------------------------------------------------
(f <- formula(m@model))

## ----eval=F--------------------------------------------------------------
#  m <- buildmer(f,data=vowels,direction='backward',control=lmerControl(optimizer='bobyqa'))

## ----echo=F--------------------------------------------------------------
cat('Fitting ML and REML reference models
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Testing terms
Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following | participant)
      grouping                      term
1         <NA>                         1
2         <NA>                 following
3         <NA>                     vowel
4         <NA>                 timepoint
5         <NA>           vowel:timepoint
6         <NA>       following:timepoint
7         <NA>           following:vowel
8         <NA> following:vowel:timepoint
9  participant                         1
10 participant                 timepoint
11 participant                 following
12 participant       timepoint:following
13        word                         1
14        word                 timepoint
                                block           LRT Iteration
 1                            NA NA 1            NA         1
 2                    NA NA following            NA         1
 3                        NA NA vowel            NA         1
 4                    NA NA timepoint            NA         1
 5              NA NA vowel:timepoint            NA         1
 6          NA NA following:timepoint            NA         1
 7              NA NA following:vowel            NA         1
 8    NA NA following:vowel:timepoint  4.967287e-01         1
 9                   NA participant 1            NA         1
10           NA participant timepoint            NA         1
11           NA participant following            NA         1
12 NA participant timepoint:following  8.319280e-11         1
13                          NA word 1            NA         1
14                  NA word timepoint 2.134598e-153         1
Updating formula: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting ML and REML reference models
[...]
All terms are significant
Finalizing by converting the model to lmerTest')

## ----echo=F,message=F----------------------------------------------------
f2 <- as.formula('f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)',.GlobalEnv)
m <- buildmer(f2,vowels,direction=NULL,control=lmerControl(optimizer='bobyqa'))

## ------------------------------------------------------------------------
summary(m)

## ------------------------------------------------------------------------
tabulate.formula(f)

## ------------------------------------------------------------------------
vowels <- cbind(vowels,model.matrix(~vowel,vowels))

## ------------------------------------------------------------------------
form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
	     ((vowel1+vowel2+vowel3+vowel4)*timepoint*following | participant) +
	     (timepoint | word))
terms <- tabulate.formula(form,group='vowel[^:]')

## ----eval=F--------------------------------------------------------------
#  m <- buildmer(terms,data=vowels,dep='f1',control=lmerControl(optimizer='bobyqa'))

