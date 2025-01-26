###
### Travail 1
### Code source principal pour appeler les fonctions développées
###
##

options(scipen = 1000)

source("echantillonner_model_3.R")
source("echantillonner_crm_cop_archi_hiera.R")
source("utilitaires.R")


# Test - model 3 de [Cossette et al., 2019]
ech_l <- echantillonner_modele_3(1000)
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)
ech[1,]

S <- apply(ech[, -1], 1, sum)
mean(S)
var(S)

10 * (1 / 0.001)
10 * (1 / (0.001^2)) + ((1 / 0.001)^2) * 10


# Test - crm avec dépendance selon des copules Archimédiennes hiérarchiques de
# [Cossette et al., 2019]
# Exemple 14 de [Cossette et al., 2019]
# ...


