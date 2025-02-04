###
### Travail 1, ACT-7119
###
###
### Exemples numériques avec le modèle 1
###
###
###
source("utilitaires/utilitaires.R")
source("modele-independance-comonotone/echantillonner_modele_1.R")
source("modele-independance-comonotone/calculer_ES.R")
source("modele-independance-comonotone/calculer_fS.R")
#### Validation des résultats simulation versus numérique
### echantillon
n <- 1000000
lam <- 3
lam2 <- 100
set.seed(20250202)
sim <- echantillonner_modele_1.R(n,
                                 qDistN = function(x) qpois(x, lam),
                                 qDistX = function(x) qpois(x, lam2))

simul <- structurer_echantillon(sim$realisations, sim$max_N_i)


S <- rowSums(simul[, -1])

## ES
mean(S)
calculer_ES_indep(dDistX = function(x) dpois(x, lam2),
                  qDistX = function(x) qpois(x, lam2),
                  dDistN = function(x) dpois(x, lam),
                  qDistN = function(x) qpois(x, lam))

nfft <- 2^20
fs <- calculer_fS_indep(nfft = nfft,
                        dDistX = function(x) dpois(x, lam2),
                        fgpN = function(phi) exp(lam * (phi - 1)))
sum(fs)
sum((seq(nfft) - 1) * fs)

## VaR_k(S)
kap = 0.9

VaR_kapp_discr(kapp = kap,
               h = 1,
               fx = fs)
sort(S)[n*kap]

## TVaR_k(S)
kap = 0.9

mean(sort(S)[(n*kap + 1):n])
TVaR_kapp_discr(kapp = kap,
                h = 1,
                fx = fs) # reste à ajuster pour considérer k réellement

# Mesure entropique
rho = 0.00001 # prendre petite valeur parce que ça peut rapidement exploser

1/rho * log(mean(exp(rho*S)))
Mesure_entropique_discr(rho = rho,
                        h = 1,
                        fx = fs)





