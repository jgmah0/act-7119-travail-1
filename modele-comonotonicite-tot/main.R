options(scipen = 1000)

source("echantillonner_model_3.R")
source("calculer_FS.R")


# Test - model 3 de [Cossette et al., 2019]
ech_l <- echantillonner_modele_3(1000, function(x) qpois(x, 10), function(x) qexp(x, 0.001))
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)
ech[1,]

S <- apply(ech[, -1], 1, sum)
mean(S)
var(S)

10 * (1 / 0.001)
10 * (1 / (0.001^2)) + ((1 / 0.001)^2) * 10


# Test 2
ktest <- approximer_kmax(function(x) ppois(x, 4), 0.0001)

calculer_FS_crm_comonotonicite(1000, function(x) ppois(x, 4), function(x) qexp(x, 0.001), ktest)
