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
f2 <- f1 ~ following + vowel + timepoint + vowel:timepoint + following:timepoint + 
    (1 + timepoint | word) + (1 + timepoint + following + timepoint:following | 
    participant)
m <- buildmer(f2,vowels,direction=NULL)

## ------------------------------------------------------------------------
(f <- formula(m@model))

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

