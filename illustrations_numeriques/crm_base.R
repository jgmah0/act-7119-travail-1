###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec modèles de dépendance de base
###
##

source("../modele-comonotonicite-tot/calculer_ES.R")
source("../modele-comonotonicite-tot/calculer_frep_XcondN.R")
source("../modele-comonotonicite-tot/calculer_mesures_risque_S.R")
source("../modele-comonotonicite-tot/echantillonner_model_3.R")

source("../modele-antimonotone-tot/calculer_ES.R")
source("../modele-antimonotone-tot/calculer_frep_XcondN.R")
source("../modele-antimonotone-tot/calculer_mesures_risque_S.R")
source("../modele-antimonotone-tot/echantillonner_model_4.R")

source("../modele-independance/echantillonner_modele_2.R")
source("../modele-independance-comonotone/echantillonner_modele_1.R")
source("../modele-independance/echantillonner_modele_2.R")

source("../utilitaires/approximer_kmax.R")


lambdaPois <- 2
alphaGa <- 3
betaGa <- 0.002

seuil <- 10^(-7)

kmax <- approximer_kmax_version_2(function(x) dpois(x, lambdaPois), 0, seuil)

# Code de l'exemple synthèse de base
FS <- function(x) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * pgamma(x, (1:kmax) * alphaGa, betaGa))

FS(500)
FS(1000)

kap <- c(0.99, 0.995)

VaR <- function(kap) optimize(function(a) abs(FS(a) - kap), c(0, 20000))$minimum

VaR(0.99)
VaR(0.995)

# Avec la formule de [Cossette et Marceau, 2022]
TVaR <- function(kap) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * (1 / (1 - kap)) * (((1:kmax) * alphaGa) / betaGa) * pgamma(VaR(kap), ((1:kmax) * alphaGa) + 1, betaGa, lower.tail = FALSE) )



1 / (1 - kappa) * sum(dpois(1:20, lam) * (1:20)/bet *
                          pgamma(v,  (1:20) + 1, bet,
                                 lower.tail = FALSE))
## Mesure entropique
rho <- c(0.0001, 0.0002)

1 / rho * log(exp(lambdaPois * (((betaGa / (betaGa - rho))^alphaGa) - 1)))
