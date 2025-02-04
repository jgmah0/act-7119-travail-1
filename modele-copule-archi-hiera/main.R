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
source("echantillonner_crm_cop_archi_hiera.R")
source("../modele-comonotonicite-tot/echantillonner_model_3.R")
source("calculer_fS.R")
source("calculer_VarS.R")

# Exemple test
# N ~ Pois(3)
# X ~ Pois(10)
# Theta0 ~ Geom(0.2)
# Theta1 ~ Geom(0.8)

kmax <- approximer_kmax_version_2(function(x) dpois(x, 1.2),
                                  borne_inf_support = 0,
                                  threshold = 0.0001)
h <- 1
imax <- ceiling(qgamma(1 - 0.00001, 3, 0.001) / h)
q0 <- 0.5
q1 <- 0.5

calculer_ES_archi_hiera(kmax, imax, h,
                        function(x) tls_inv_geom_essais(x, q0),
                        function(x) ppois(x, 1.2),
                        function(x) dpois(x, 1.2),
                        function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
                        function(x) pgamma(x, 3, 0.001),
                        function(x) dgeo_shifted(x, 1 - q0),
                        function(x, th0) dnbinom_essais(x, th0, q1),
                        0.00001) # 3 898.88

# À tester
calculer_VarS_archi_hiera(kmax, imax, h,
                          function(x) tls_inv_geom_essais(x, q0),
                          function(x) ppois(x, 1.2),
                          function(x) dpois(x, 1.2),
                          function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
                          function(x) pgamma(x, 3, 0.001),
                          function(x) dgeo_shifted(x, 1 - q0),
                          function(x, th0) dnbinom_essais(x, th0, q1),
                          0.00001)


n <- 1000000
ech_l <- echantillonner_crm_cop_archi_hiera(n,
                                            function(x) qgeo_shifted(x, 1 - q0),
                                            function(x) tls_geom_essais(x, q0),
                                            function(x) qgeo_shifted(x, 1 - q1),
                                            function(x) fgp_geom_essais(tls_geom_essais(x, q1), q0),
                                            function(x) qpois(x, 1.2),
                                            function(x) qgamma(x, 3, 0.001))
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)


mean(rowSums(ech[, -1])) # E[S] empirique
var(rowSums(ech[, -1])) # Var(S) empirique


fS <- calculer_fS_archi_hiera(2^14, kmax, h,
                              function(x) tls_inv_geom_essais(x, q0),
                              function(x) ppois(x, 1.2),
                              function(x) dpois(x, 1.2),
                              function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
                              function(x) pgamma(x, 3, 0.001),
                              function(x) dgeo_shifted(x, 1 - q0),
                              function(x, th0) dnbinom_essais(x, th0, q1),
                              0.0000001)

sum(fS)
# Somme à 1
