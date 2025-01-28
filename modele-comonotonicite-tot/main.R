options(scipen = 1000)

source("echantillonner_model_3.R")
source("calculer_FS.R")
source("calculer_ES.R")
source("calculer_VarS.R")
source("calculer_mesures_risque_S.R")
source("calculer_frep_XcondN.R")
source("../utilitaires/approximer_kmax.R")


# Test - model 3 de [Cossette et al., 2019]
n <- 500000
ech_l <- echantillonner_modele_3(n, function(x) qpois(x, 10), function(x) qexp(x, 0.001))
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)
ech[1,]

S <- apply(ech[, -1], 1, sum)

mean(S)
var(S)

10 * (1 / 0.001)
10 * (1 / (0.001^2)) + ((1 / 0.001)^2) * 10


# Validation de FS
ktest <- approximer_kmax(function(x) ppois(x, 10), 0.0001)

calculer_FS_crm_comonotonicite(0,
                               function(x) ppois(x, 10),
                               function(x) qexp(x, 0.001),
                               ktest) == ppois(0, 10)
calculer_FS_crm_comonotonicite(1000, function(x) ppois(x, 10), function(x) qexp(x, 0.001), ktest)
calculer_FS_crm_comonotonicite(100000000, function(x) ppois(x, 10), function(x) qexp(x, 0.001), ktest)
# Semble tendre vers 1

c(empirique = mean(S), theorique = calculer_ES_crm_comonotonicite(function(x) ppois(x, 10), function(x) qexp(x, 0.001), ktest))
c(empirique = var(S), theorique = calculer_VarS_crm_comonotonicite(function(x) ppois(x, 10), function(x) qexp(x, 0.001), ktest))

calculer_VaRS_crm_comonotonicite(0.99, function(x) qpois(x, 10), function(x) qexp(x, 0.001))
calculer_TVaRS_crm_comonotonicite(0.99, function(x) qpois(x, 10), function(x) qexp(x, 0.001))


calculer_frep_XcondN_crm_comonotonicite(seq(100, 500, by = 100),
                                        8,
                                        function(x) ppois(x, 10),
                                        function(x) pexp(x, 0.001))

plot(1:20,
     calculer_frep_XcondN_crm_comonotonicite(500,
                                             1:20,
                                             function(x) ppois(x, 10),
                                             function(x) pexp(x, 0.001)))
