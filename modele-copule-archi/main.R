###
### Travail 1, ACT-7119
###
###
### Exemples numériques avec les copules archimédiennes (facteur simple)
###
###
###
source("utilitaires/utilitaires.R")
source("modele-copule-archi/calculer_ES.R")
source("modele-copule-archi/calculer_EXw_condN.R")
source("modele-copule-archi/calculer_fX_condN.R")
source("modele-copule-archi/calculer_varX_condN.R")
source("modele-copule-archi/echantillonner_crm_cop_archi.R")
source("modele-copule-archi/geom_shifted.R")
source("modele-copule-archi/calculer_fS.R")
#### Validation des résultats simulation versus numérique
### echantillon
n = 1000000
alpha <- 0.4
lam <- 3
lam2 <- 100
set.seed(20250202)
sim <- echantillonner_crm_cop_archi(n,
                             qDistTheta = function(x) qgeo_shifted(x, alpha),
                             tlsTheta = function(x) (1-alpha) * exp(-x) / (1 - alpha * exp(-x)),
                             qDistN = function(x) qpois(x, lam),
                             qDistX = function(x) qpois(x, lam2))

simul <- structurer_echantillon(sim$realisations, sim$max_N_i)


S <- rowSums(simul[, -1])


## ES
mean(S)
calculer_ES_crm_archi_simple(alpha,
                             tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                             dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                             qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                             pDistX = function(x) ppois(x, lam2),
                             qDistX = function(x) qpois(x, lam2),
                             dDistN = function(x) dpois(x, lam),
                             pDistN = function(x) ppois(x, lam),
                             qDistN = function(x) qpois(x, lam))

Exw <- numeric(16)
for (k in 1:16)
{
    Exw[k] <- calculer_EXw_condN_crm_archi_simple(k, 1, alpha,
                                        tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                        dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                        qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                        pDistX = function(x) ppois(x, lam2),
                                        qDistX = function(x) qpois(x, lam2),
                                        dDistN = function(x) dpois(x, lam),
                                        pDistN = function(x) ppois(x, lam))
}

sum((1:16) * dpois(1:16, lam) * Exw) # donne comme le calcul de ES

nfft = 2^15

fs <- calculer_fS_archi_simple(nfft = nfft, alpha = alpha,
                               tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                               dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                               qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                               pDistX = function(x) ppois(x, lam2),
                               qDistN = function(x) qpois(x, lam),
                               pDistN = function(x) ppois(x, lam))
sum(fs)
sum(0:(nfft-1) * fs)
## EX|N , Frep cond X | N
k <- 4
X_N <- simul[simul[, 1] == k, 2:(k+1)]

mean(X_N)
Exw[k]


fxn <- numeric(200)
for (i in 1:200)
{
    fxn[i] <- calculer_fx_condN_crm_archi_simple(i, k, alpha,
                                                 tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                 dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                                 qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                                 pDistX = function(x) ppois(x, lam2),
                                                 dDistN = function(x) dpois(x, lam),
                                                 pDistN = function(x) ppois(x, lam))
}

sum(fxn)
xx <- 110
mean(X_N <= xx)
cumsum(fxn)[xx]

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
rho = 0.001 # prendre petite valeur parce que ça peut rapidement exploser

1/rho * log(mean(exp(rho*S)))
Mesure_entropique_discr(rho = rho,
                        fx = fs)






