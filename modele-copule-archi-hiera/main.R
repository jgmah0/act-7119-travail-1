###
### Travail 1, ACT-7119
### Script principal
### CRM avec dépendance selon une copule Archimédienne hiérarchique
###
##


source("calculer_ES.R")
source("../utilitaires/approximer_kmax.R")
source("../utilitaires/dist_binomiale_negative_essais.R")
source("../modele-copule-archi/geom_shifted.R")


# Exemple test
# N ~ Pois(3)
# X ~ Pois(10)
# Theta0 ~ Geom(0.2)
# Theta1 ~ Geom(0.8)

kmax <- approximer_kmax_version_2(function(x) dpois(x, 3), threshold = 0.00001)
imax <- approximer_kmax_version_2(function(x) dpois(x, 10), threshold = 0.00001)
h <- 1
q0 <- 0.2
q1 <- 0.8

calculer_ES_archi_hiera(kmax, imax, h,
                        function(x) tls_inv_geom_essais(x, q0),
                        function(x) ppois(x, 3),
                        function(x) dpois(x, 3),
                        function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
                        function(x) ppois(x, 10),
                        function(x) dgeo_shifted(x, 1 - q0),
                        function(x, th0) dnbinom_essais(x, th0, q1),
                        0.00001)
