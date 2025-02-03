###
### Travail 1, ACT-7119
### Script principal
### CRM avec dépendance selon une copule Archimédienne hiérarchique
###
##


source("calculer_ES.R")
source("../utilitaires/approximer_kmax.R")
source("../utilitaires/dist_binomiale_negative_essais.R")
source("../utilitaires/utilitaires.R")
source("../modele-copule-archi/geom_shifted.R")
source("../modele-copule-archi/geom_shifted.R")


# Exemple test
# N ~ Pois(3)
# X ~ Pois(10)
# Theta0 ~ Geom(0.2)
# Theta1 ~ Geom(0.8)

kmax <- approximer_kmax_version_2(function(x) dpois(x, 1.2),
                                  borne_inf_support = 0,
                                  threshold = 0.0001)
h <- 1000
imax <- ceiling(qgamma(1 - 0.00001, 3, 0.001) / h)
q0 <- 0.2
q1 <- 0.2

calculer_ES_archi_hiera(kmax, imax, h,
                        function(x) tls_inv_geom_essais(x, q0),
                        function(x) ppois(x, 1.2),
                        function(x) dpois(x, 1.2),
                        function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
                        function(x) pgamma(x, 3, 0.001),
                        function(x) dgeo_shifted(x, 1 - q0),
                        function(x, th0) dnbinom_essais(x, th0, q1),
                        0.00001)
