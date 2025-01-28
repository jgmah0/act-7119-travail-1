###
### Travail 1, ACT-7119
### Exemple synthèse des modèles de base
###
###
###
###
##

source("modele-comonotonicite-tot/calculer_ES.R")
source("modele-comonotonicite-tot/calculer_frep_XcondN.R")
source("modele-comonotonicite-tot/calculer_mesures_risque_S.R")
source("modele-comonotonicite-tot/echantillonner_model_3.R")

source("modele-antimonotone-tot/calculer_ES.R")
source("modele-antimonotone-tot/calculer_frep_XcondN.R")
source("modele-antimonotone-tot/calculer_mesures_risque_S.R")
source("modele-antimonotone-tot/echantillonner_model_4.R")
lam <- 2
bet <- 1 / 100

kappa <- 0.95

rho <- 0.001


### Modèle classique (indépendance)

## TVaR
1 / (1 - kappa) * sum(dpois(1:15, lam) * (1/bet) * (1:15) *
                          pgamma(qgamma(kappa, (1:15), bet),  (1:15) + 1, bet,
                                 lower.tail = FALSE))
## Mesure entropique

1 / rho * log(exp(2 * (bet / (bet - rho) - 1)))



### Modèle avec indépendance et comonotonicité
sum(dpois(1:15, lam) * (1:15) * (qexp(kappa, bet) + 1 / bet))

 sum(dpois(1:15, lam) * (qexp(kappa, bet / k) + 1 / (bet / k)))

### Modèle avec composantes comonotones
cbind(0:15, ppois(0:15, lam), 1- ppois(0:15, lam))


calculer_ES_crm_comonotonicite(function(x) ppois(x, lam),
                               function(x) qexp(x, bet),
                               kmax=15)

# validation
FN <- ppois(0:15, lam)

i <- 1:12

antideriv <- function(u) ((u - 1) * log(1 - u) - u)

-100 * sum(i * (antideriv(ppois(1:12, lam)) -
                    antideriv(ppois(0:11, lam))))  #ok comme integrate

## EX|N
k <- 2
tvarX <- function(k, bet) qexp(k, bet) + 1/bet

1 / (ppois(kk, lam) - ppois(kk - 1, lam)) * ((1 - ppois(k-1, lam)) *
                                                 tvarX(ppois(k-1, lam), bet) -
                                                 (1 - ppois(k, lam)) *
                                                 tvarX(ppois(k, lam), bet))

integrate(function(y) qexp(y, bet), ppois(1, lam), ppois(2, lam))$value / (ppois(kk, lam) - ppois(kk - 1, lam))


## TVARS

calculer_TVaRS_crm_comonotonicite(0.95, function(x) qpois(x, lam),
                                  function(x) qexp(x, bet))

calculer_TVaRS_crm_comonotonicite(0.95, function(x) qpois(x, lam),
                                 function(x) qexp(x, bet))

## Mesure entropique
calculer_entropique_crm_comonotonicite(0.001, function(x) ppois(x, lam),
                                       function(x) qpois(x, lam),
                                       function(x) qexp(x, bet), 15)

### Modèle avec composantes antimonotones

calculer_ES_crm_antimonotones(function(x) ppois(x, lam),
                              function(x) qexp(x, bet),
                              15)

1 / (ppois(2, lam) - ppois(1, lam)) * (ppois(2, lam) * tvarX(1 - ppois(2, lam), bet) -
                                           ppois(1, lam) * tvarX(1 - ppois(1, lam), bet))

calculer_TVaRS_crm_antimonotones(0.95, function(x) qpois(x, lam),
                                 function(x) qexp(x, bet))

calculer_entropique_crm_antimonotones(0.001, function(x) ppois(x, lam),
                                      function(x) qpois(x, lam),
                                      function(x) qexp(x, bet), 15)
