## ----echo=FALSE---------------------------------------------------------------------------------------------
options(width=110)

## -----------------------------------------------------------------------------------------------------------
library(buildmer)
head(vowels)

## -----------------------------------------------------------------------------------------------------------
f <- f1 ~ vowel*timepoint*following * neighborhood*information*stress + 
	 (vowel*timepoint*following * neighborhood+information+stress | participant) +
	 (timepoint | word)

## -----------------------------------------------------------------------------------------------------------
f <- f1 ~ vowel*timepoint*following +
	 (vowel*timepoint*following | participant) +
	 (timepoint | word)

## ----eval=FALSE---------------------------------------------------------------------------------------------
# library(lme4)
# m <- buildmer(f,data=vowels,buildmerControl=buildmerControl(direction='order',
# 	      args=list(control=lmerControl(optimizer='bobyqa'))))

## ----echo=FALSE---------------------------------------------------------------------------------------------
cat('Determining predictor order
Fitting via lm: f1 ~ 1
Currently evaluating LRT for: following, timepoint, vowel
Fitting via lm: f1 ~ 1 + following
Fitting via lm: f1 ~ 1 + timepoint
Fitting via lm: f1 ~ 1 + vowel
Updating formula: f1 ~ 1 + vowel
Currently evaluating LRT for: following, timepoint
Fitting via lm: f1 ~ 1 + vowel + following
Fitting via lm: f1 ~ 1 + vowel + timepoint
Updating formula: f1 ~ 1 + vowel + timepoint
Currently evaluating LRT for: following, vowel:timepoint
Fitting via lm: f1 ~ 1 + vowel + timepoint + following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint
Currently evaluating LRT for: following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following
Currently evaluating LRT for: timepoint:following, vowel:following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + vowel:following
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following
Currently evaluating LRT for: vowel:following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following
Currently evaluating LRT for: vowel:timepoint:following
Fitting via lm: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following
Fitting via gam, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following
Currently evaluating LRT for: 1 | participant, 1 | word
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 | participant)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 | word)
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 | participant)
Currently evaluating LRT for: following | participant, timepoint | participant, vowel |
    participant, 1 | word
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + following |
    participant)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint |
    participant)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + vowel | participant)
boundary (singular) fit: see ?isSingular
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 | participant) + (1 |
    word)
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 | participant) + (1 | word)
Currently evaluating LRT for: following | participant, timepoint | participant, vowel |
    participant, timepoint | word
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + following |
    participant) + (1 | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint |
    participant) + (1 | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + vowel | participant)
    + (1 | word)
boundary (singular) fit: see ?isSingular
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 | participant) + (1 +
    timepoint | word)
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 | participant) + (1 + timepoint | word)
Currently evaluating LRT for: following | participant, timepoint | participant, vowel |
    participant
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + following |
    participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint |
    participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + vowel | participant)
    + (1 + timepoint | word)
boundary (singular) fit: see ?isSingular
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 + timepoint | participant) + (1 + timepoint |
    word)
Currently evaluating LRT for: following | participant, vowel | participant
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint + following
    | participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint + vowel |
    participant) + (1 + timepoint | word)
boundary (singular) fit: see ?isSingular
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 + timepoint + following | participant) + (1 +
    timepoint | word)
Currently evaluating LRT for: timepoint:following | participant, vowel | participant
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint + following
    + timepoint:following | participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint + following
    + vowel | participant) + (1 + timepoint | word)
boundary (singular) fit: see ?isSingular
Updating formula: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following + timepoint:following +
    vowel:following + vowel:timepoint:following + (1 + timepoint + following + timepoint:following
    | participant) + (1 + timepoint | word)
Currently evaluating LRT for: vowel | participant
Fitting via lmer, with REML: f1 ~ 1 + vowel + timepoint + vowel:timepoint + following +
    timepoint:following + vowel:following + vowel:timepoint:following + (1 + timepoint + following
    + timepoint:following + vowel | participant) + (1 + timepoint | word)
boundary (singular) fit: see ?isSingular
Ending the ordering procedure due to having reached the maximal feasible model - all higher models
    failed to converge. The types of convergence failure are: Singular fit
Finalizing by converting the model to lmerTest')

## ----include=F----------------------------------------------------------------------------------------------
library(lme4)
#hack for consistency with actual output without actually fitting the model every time I change something in the vignette
m <- buildmer:::mkBuildmer(model=list(formula=(function () as.formula('f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following + timepoint:following | participant) + (1 + timepoint | word)',.GlobalEnv))()))

## -----------------------------------------------------------------------------------------------------------
(f <- formula(m@model))

## ----eval=FALSE---------------------------------------------------------------------------------------------
# m <- buildmer(f,data=vowels,buildmerControl=list(direction='backward',
# 	      args=list(control=lmerControl(optimizer='bobyqa'))))

## ----echo=FALSE---------------------------------------------------------------------------------------------
cat('Fitting ML and REML reference models
Fitting via lmer, with REML: f1 ~ following + vowel + timepoint + vowel:timepoint +
    following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following
    + timepoint:following | participant) + (1 + timepoint | word)

Fitting via lmer, with REML: f1 ~ following + vowel + timepoint + vowel:timepoint +
    following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following
    + timepoint:following | participant) + (1 + timepoint | word)
Testing terms
Fitting via lmer, with ML: f1 ~ 1 + following + vowel + timepoint + vowel:timepoint +
    following:timepoint + following:vowel + (1 + timepoint + following + timepoint:following |
    participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + following + vowel + timepoint + vowel:timepoint +
    following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following
    | participant) + (1 + timepoint | word)
Fitting via lmer, with REML: f1 ~ 1 + following + vowel + timepoint + vowel:timepoint +
    following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint + following
    + timepoint:following | participant) + (1 | word)
      grouping                      term                              block Iteration           LRT
1         <NA>                         1                            NA NA 1         1            NA
2         <NA>                 following                    NA NA following         1            NA
3         <NA>                     vowel                        NA NA vowel         1            NA
4         <NA>                 timepoint                    NA NA timepoint         1            NA
5         <NA>           vowel:timepoint              NA NA vowel:timepoint         1            NA
6         <NA>       following:timepoint          NA NA following:timepoint         1            NA
7         <NA>           following:vowel              NA NA following:vowel         1            NA
8         <NA> following:vowel:timepoint    NA NA following:vowel:timepoint         1  3.609316e-30
9  participant                         1                   NA participant 1         1            NA
10 participant                 timepoint           NA participant timepoint         1            NA
11 participant                 following           NA participant following         1            NA
12 participant       timepoint:following NA participant timepoint:following         1  1.013211e-10
13        word                         1                          NA word 1         1            NA
14        word                 timepoint                  NA word timepoint         1 2.198802e-153
All terms are significant
Finalizing by converting the model to lmerTest')

## ----echo=FALSE,message=FALSE-------------------------------------------------------------------------------
f2 <- as.formula('f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)',.GlobalEnv)
m <- buildmer(f2,vowels,buildmerControl=list(direction=NULL,args=list(control=lmerControl(optimizer='bobyqa'))))

## -----------------------------------------------------------------------------------------------------------
summary(m)

## -----------------------------------------------------------------------------------------------------------
tabulate.formula(f)

## -----------------------------------------------------------------------------------------------------------
vowels <- cbind(vowels,model.matrix(~vowel,vowels))

## -----------------------------------------------------------------------------------------------------------
form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
	     ((vowel1+vowel2+vowel3+vowel4)*timepoint*following | participant) +
	     (timepoint | word))
terms <- tabulate.formula(form,group='vowel[^:]')

## ----eval=FALSE---------------------------------------------------------------------------------------------
# m <- buildmer(terms,data=vowels,buildmerControl=buildmerControl(dep='f1',
# 	      args=list(control=lmerControl(optimizer='bobyqa'))))
# ## (output not shown)

