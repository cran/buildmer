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

## ----include=F-----------------------------------------------------------
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
m <- buildmer:::mkBuildmer(model=list(formula=(function ()
f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
)()))

## ------------------------------------------------------------------------
(f <- formula(m@model))

## ----eval=F--------------------------------------------------------------
#  m <- buildmer(f,data=vowels,direction='backward',control=lmerControl(optimizer='bobyqa'))

## ----include=F-----------------------------------------------------------
cat('Fitting ML and REML reference models
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Testing terms
Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + following:vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following | participant)
      grouping                      term block           LRT Iteration
1         <NA>                         1     1            NA         1
2         <NA>                 following     2            NA         1
3         <NA>                     vowel     3            NA         1
4         <NA>                 timepoint     4            NA         1
5         <NA>           vowel:timepoint     5            NA         1
6         <NA>       following:timepoint     6            NA         1
7         <NA>           following:vowel     7            NA         1
8         <NA> following:vowel:timepoint     8  4.967286e-01         1
9         word                         1     9            NA         1
10        word                 timepoint    10 2.134598e-153         1
11 participant                         1    11            NA         1
12 participant                 timepoint    12            NA         1
13 participant                 following    13            NA         1
14 participant       timepoint:following    14  8.319280e-11         1
Updating formula: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
Fitting ML and REML reference models
[...]
All terms are significant')
#Updating formula: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting ML and REML reference models
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Testing terms
#Fitting with ML: f1 ~ following + vowel + timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + following:vowel + (1 + timepoint | word) + (1 + timepoint + following | participant)
#   index    grouping                term                                 code block           LRT Iteration
#1   <NA>        <NA>                   1                              NA NA 1     1            NA         1
#2   <NA>        <NA>           following                      NA NA following     2            NA         1
#3   <NA>        <NA>               vowel                          NA NA vowel     3            NA         1
#4   <NA>        <NA>           timepoint                      NA NA timepoint     4            NA         1
#5   <NA>        <NA>     vowel:timepoint                NA NA vowel:timepoint     5  2.820922e-10         1
#6   <NA>        <NA> following:timepoint            NA NA following:timepoint     6  6.701786e-04         1
#7   <NA>        <NA>     following:vowel                NA NA following:vowel     7  1.104266e-01         1
#9    9 1        word                   1                           9 1 word 1     9            NA         1
#10   9 1        word           timepoint                   9 1 word timepoint    10 1.218638e-157         1
#11  10 1 participant                   1                   10 1 participant 1    11            NA         1
#12  10 1 participant           timepoint           10 1 participant timepoint    12            NA         1
#13  10 1 participant           following           10 1 participant following    13            NA         1
#14  10 1 participant timepoint:following 10 1 participant timepoint:following    14  8.298328e-11         1
#Updating formula: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting ML and REML reference models
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Testing terms
#Fitting with ML: f1 ~ following + vowel + timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with ML: f1 ~ following + vowel + timepoint + vowel:timepoint + (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 | word) + (1 + timepoint + following + timepoint:following | participant)
#Fitting with REML: f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + (1 + timepoint | word) + (1 + timepoint + following | participant)
#   index    grouping                term                                 code block           LRT Iteration
#1   <NA>        <NA>                   1                              NA NA 1     1            NA         2
#2   <NA>        <NA>           following                      NA NA following     2            NA         2
#3   <NA>        <NA>               vowel                          NA NA vowel     3            NA         2
#4   <NA>        <NA>           timepoint                      NA NA timepoint     4            NA         2
#5   <NA>        <NA>     vowel:timepoint                NA NA vowel:timepoint     5  2.827268e-10         2
#6   <NA>        <NA> following:timepoint            NA NA following:timepoint     6  6.714423e-04         2
#9    9 1        word                   1                           9 1 word 1     9            NA         2
#10   9 1        word           timepoint                   9 1 word timepoint    10 6.501098e-158         2
#11  10 1 participant                   1                   10 1 participant 1    11            NA         2
#12  10 1 participant           timepoint           10 1 participant timepoint    12            NA         2
#13  10 1 participant           following           10 1 participant following    13            NA         2
#14  10 1 participant timepoint:following 10 1 participant timepoint:following    14  8.450196e-11         2
#All terms are significant
f2 <- f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + 
    (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | 
    participant)
m <- buildmer(f2,vowels,direction=NULL)

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
(terms <- tabulate.formula(form))

## ------------------------------------------------------------------------
terms[ 2: 5,'block'] <- 'same1'
terms[ 7:10,'block'] <- 'same2'
terms[12:15,'block'] <- 'same3'
terms[17:20,'block'] <- 'same4'
terms[22:25,'block'] <- 'same5'
terms[27:30,'block'] <- 'same6'
terms[32:35,'block'] <- 'same7'
terms[37:40,'block'] <- 'same8'

## ----eval=F--------------------------------------------------------------
#  m <- buildmer(terms,data=vowels,dep='f1',control=lmerControl(optimizer='bobyqa'))

