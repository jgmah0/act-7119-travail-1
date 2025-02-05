###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec dépendance selon une copule Archimédienne hiérarchique AMH-AMH
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

seuil <- 0.0000001

lambdaPois <- 2

h <- 1
alphaGa <- 3
betaGa <- 0.002
alpha0 <- c(0.3, 0.8)
alpha1 <- c(0.3, 0.8)
q0 <- 1 - alpha0
q1 <- 1 - alpha1

kmax <- approximer_kmax_version_2(function(x) dpois(x, lambdaPois),
                                  borne_inf_support = 0,
                                  threshold = seuil)
imax <- ceiling(qgamma(1 - seuil, alphaGa, betaGa) / h)

# À tester
# varS <- calculer_VarS_archi_hiera(kmax, imax, h,
#                                   function(x) tls_inv_geom_essais(x, q0),
#                                   function(x) ppois(x, 1.2),
#                                   function(x) dpois(x, 1.2),
#                                   function(x) tls_inv_geom_essais_comp_geom_essais(x, q0, q1),
#                                   function(x) pgamma(x, 3, 0.001),
#                                   function(x) dgeo_shifted(x, 1 - q0),
#                                   function(x, th0) dnbinom_essais(x, th0, q1[1]),
#                                   seuil)
#
#
n <- 1000000
ech_l <- echantillonner_crm_cop_archi_hiera(n,
                                            function(x) qgeo_shifted(x, 1 - q0[1]),
                                            function(x) tls_geom_essais(x, q0[1]),
                                            function(x) qgeo_shifted(x, 1 - q1[1]),
                                            function(x) fgp_geom_essais(tls_geom_essais(x, q1[1]), q0[1]),
                                            function(x) qpois(x, lambdaPois),
                                            function(x) qgamma(x, alphaGa, betaGa))
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)


mean(rowSums(ech[, -1])) # E[S] empirique
var(rowSums(ech[, -1])) # Var(S) empirique


nfft <- 2^16
fS_cas11 <- calculer_fS_archi_hiera(nfft, kmax, h,
                                    function(x) tls_inv_geom_essais(x, q0[1]),
                                    function(x) ppois(x, lambdaPois),
                                    function(x) dpois(x, lambdaPois),
                                    function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[1], q1[1]),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) dgeo_shifted(x, 1 - q0[1]),
                                    function(x, th0) dnbinom_essais(x, th0, q1[1]),
                                    seuil)

sum(fS_cas11)
sum(fS_cas11) + seuil

# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas11)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas11) - (sum(((0:(nfft - 1)) * h) * fS_cas11)^2)

sum(fS_cas11[seq((500 / h) + 1)])
sum(fS_cas11[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas11)
VaR_kapp_discr(0.995, h, fS_cas11)

TVaR_kapp_discr_v2(0.99, h, fS_cas11)
TVaR_kapp_discr_v2(0.995, h, fS_cas11)

# Mesure_entropique_discr_v2(0.001, h, fS_cas11)

Mesure_entropique_discr_v2(0.001, h, pmax(fS_cas11, 0)) # valeurs négatives à la fin du vecteur fS.

# Mesure_entropique_discr_v2(0.0001, h, fS_cas11)

Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas11, 0)) # valeurs négatives à la fin du vecteur fS.

fS_cas12 <- calculer_fS_archi_hiera(nfft, kmax, h,
                                    function(x) tls_inv_geom_essais(x, q0[1]),
                                    function(x) ppois(x, lambdaPois),
                                    function(x) dpois(x, lambdaPois),
                                    function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[1], q1[2]),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) dgeo_shifted(x, 1 - q0[1]),
                                    function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                    seuil)

sum(fS_cas12[seq((500 / h) + 1)])
sum(fS_cas12[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas12)
VaR_kapp_discr(0.995, h, fS_cas12)

TVaR_kapp_discr_v2(0.99, h, fS_cas12)
TVaR_kapp_discr_v2(0.995, h, fS_cas12)

Mesure_entropique_discr_v2(0.001, h, pmax(fS_cas12, 0)) # valeurs négatives à la fin du vecteur fS.
Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas12, 0)) # valeurs négatives à la fin du vecteur fS.

fS_cas21 <- calculer_fS_archi_hiera(nfft, kmax, h,
                                    function(x) tls_inv_geom_essais(x, q0[2]),
                                    function(x) ppois(x, lambdaPois),
                                    function(x) dpois(x, lambdaPois),
                                    function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[1]),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) dgeo_shifted(x, 1 - q0[2]),
                                    function(x, th0) dnbinom_essais(x, th0, q1[1]),
                                    seuil)

sum(fS_cas21[seq((500 / h) + 1)])
sum(fS_cas21[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas21)
VaR_kapp_discr(0.995, h, fS_cas21)

TVaR_kapp_discr_v2(0.99, h, fS_cas21)
TVaR_kapp_discr_v2(0.995, h, fS_cas21)

Mesure_entropique_discr_v2(0.001, h, pmax(fS_cas21, 0)) # valeurs négatives à la fin du vecteur fS.
Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas21, 0)) # valeurs négatives à la fin du vecteur fS.

fS_cas22 <- calculer_fS_archi_hiera(nfft, kmax, h,
                                    function(x) tls_inv_geom_essais(x, q0[2]),
                                    function(x) ppois(x, lambdaPois),
                                    function(x) dpois(x, lambdaPois),
                                    function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) dgeo_shifted(x, 1 - q0[2]),
                                    function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                    seuil)

sum(fS_cas22[seq((500 / h) + 1)])
sum(fS_cas22[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas22)
VaR_kapp_discr(0.995, h, fS_cas22)

TVaR_kapp_discr_v2(0.99, h, fS_cas22)
TVaR_kapp_discr_v2(0.995, h, fS_cas22)

Mesure_entropique_discr_v2(0.001, h, pmax(fS_cas22, 0)) # valeurs négatives à la fin du vecteur fS.
Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas22, 0)) # valeurs négatives à la fin du vecteur fS.

# valider fS en calculant espS à la mitaine et vérifier avec les valeurs
# trouvées ci-dessus.


# sum(fS)
# Somme à 1