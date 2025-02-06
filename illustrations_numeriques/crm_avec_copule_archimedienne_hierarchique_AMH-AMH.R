###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec dépendance selon une copule Archimédienne hiérarchique AMH-AMH
###
##

source("../modele-copule-archi-hiera/calculer_fmp_conjointe_composantes.R")
source("../modele-copule-archi-hiera/calculer_ES.R")
source("../utilitaires/approximer_kmax.R")
source("../utilitaires/dist_binomiale_negative_essais.R")
source("../utilitaires/utilitaires.R")
source("../modele-copule-archi/geom_shifted.R")
source("../modele-copule-archi/geom_shifted.R")
source("../modele-copule-archi-hiera/echantillonner_crm_cop_archi_hiera.R")
source("../modele-copule-archi-hiera/calculer_fS.R")
source("../modele-copule-archi-hiera/calculer_VarS.R")

seuil <- 0.0000001

lambdaPois <- 2

h <- 1
alphaGa <- 3
betaGa <- 0.002
alpha0 <- c(0.3, 0.55)
alpha1 <- c(0.3, 0.55)
q0 <- 1 - alpha0
q1 <- 1 - alpha1

kmax <- approximer_kmax_version_2(function(x) dpois(x, lambdaPois),
                                  borne_inf_support = 0,
                                  threshold = seuil)
imax <- ceiling(qgamma(1 - seuil, alphaGa, betaGa) / h)

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
sum(((0:(nfft - 1)) * h) * fS_cas11) # validé par simulation (approximation ok)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas11) - (sum(((0:(nfft - 1)) * h) * fS_cas11)^2)
# Valeur obtenue avec 1 000 000 de réalisations : 7 528 584


sum(fS_cas11[seq((500 / h) + 1)])
sum(fS_cas11[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas11)
VaR_kapp_discr(0.995, h, fS_cas11)

TVaR_kapp_discr_v2(0.99, h, fS_cas11)
TVaR_kapp_discr_v2(0.995, h, fS_cas11)

# Mesure_entropique_discr_v2(0.001, h, fS_cas11)

Mesure_entropique_discr_v2(0.0002, h, pmax(fS_cas11, 0)) # valeurs négatives à la fin du vecteur fS.

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

# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas12)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas12) - (sum(((0:(nfft - 1)) * h) * fS_cas12)^2)

sum(fS_cas12[seq((500 / h) + 1)])
sum(fS_cas12[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas12)
VaR_kapp_discr(0.995, h, fS_cas12)

TVaR_kapp_discr_v2(0.99, h, fS_cas12)
TVaR_kapp_discr_v2(0.995, h, fS_cas12)

Mesure_entropique_discr_v2(0.0002, h, pmax(fS_cas12, 0)) # valeurs négatives à la fin du vecteur fS.
Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas12, 0)) # valeurs négatives à la fin du vecteur fS.

saveRDS(fS_cas11, "fS_cas11.rds")
saveRDS(fS_cas12, "fS_cas12.rds")

fS_cas21 <- calculer_fS_archi_hiera(nfft, kmax, h,
                                    function(x) tls_inv_geom_essais(x, q0[2]),
                                    function(x) ppois(x, lambdaPois),
                                    function(x) dpois(x, lambdaPois),
                                    function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[1]),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) dgeo_shifted(x, 1 - q0[2]),
                                    function(x, th0) dnbinom_essais(x, th0, q1[1]),
                                    seuil)

saveRDS(fS_cas21, "fS_cas21.rds")

# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas21)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas21) - (sum(((0:(nfft - 1)) * h) * fS_cas21)^2)

sum(fS_cas21[seq((500 / h) + 1)])
sum(fS_cas21[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas21)
VaR_kapp_discr(0.995, h, fS_cas21)

TVaR_kapp_discr_v2(0.99, h, fS_cas21)
TVaR_kapp_discr_v2(0.995, h, fS_cas21)

Mesure_entropique_discr_v2(0.0002, h, pmax(fS_cas21, 0)) # valeurs négatives à la fin du vecteur fS.
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

saveRDS(fS_cas22, "fS_cas22.rds")


# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas22)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas22) - (sum(((0:(nfft - 1)) * h) * fS_cas22)^2)

sum(fS_cas22[seq((500 / h) + 1)])
sum(fS_cas22[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas22)
VaR_kapp_discr(0.995, h, fS_cas22)

TVaR_kapp_discr_v2(0.99, h, fS_cas22)
TVaR_kapp_discr_v2(0.995, h, fS_cas22)

Mesure_entropique_discr_v2(0.0002, h, pmax(fS_cas22, 0)) # valeurs négatives à la fin du vecteur fS.
Mesure_entropique_discr_v2(0.0001, h, pmax(fS_cas22, 0)) # valeurs négatives à la fin du vecteur fS.

fS_cas11 <- readRDS("fS_cas11.rds")
fS_cas12 <- readRDS("fS_cas12.rds")
fS_cas21 <- readRDS("fS_cas21.rds")
fS_cas22 <- readRDS("fS_cas22.rds")

test_comparaison <- matrix(c(cumsum(fS_cas11[1:8000]), cumsum(fS_cas12[1:8000])), ncol = 2)

test_comparaison[7900:8000,]

plot((0:15000) * h, cumsum(fS_cas11[1:15001]), lwd = 2, type = "l", col = "green", xlab = "ih", ylab = "F_S (ih)")
title("Fonctions de répartition de la v.a. S pour diverses valeurs de alpha_0\net alpha_1 dans le contexte d'un CRM avec dépendance\nselon une copule Archimédienne hiérarchique AMH-AMH")
lines((0:15000) * h, cumsum(fS_cas12[1:15001]), lwd = 1, col = "blue")
lines((0:15000) * h, cumsum(fS_cas21[1:15001]), lwd = 2, col = "purple")
lines((0:15000) * h, cumsum(fS_cas22[1:15001]), lwd = 1, col = "orange")
legend(8000, 0.6, c("alpha0 = 0.3, alpha1 = 0.3",
               "alpha0 = 0.3, alpha1 = 0.55",
               "alpha0 = 0.55, alpha1 = 0.3",
               "alpha0 = 0.55, alpha1 = 0.55"),
       col = c("green", "blue", "purple", "orange"),
       lwd = rep(2.5, 4))


sup_fx_cond <- seq(0, 6000, 20)
fX_cond_N_k1 <- sapply(sup_fx_cond,
                       function(y) calculer_fX_cond_N_archi_hiera(1, y, h,
                                                                  function(x) tls_inv_geom_essais(x, q0[1]),
                                                                  function(x) ppois(x, lambdaPois),
                                                                  function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                                                  function(x) pgamma(x, alphaGa, betaGa),
                                                                  function(x) dgeo_shifted(x, 1 - q0[1]),
                                                                  function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                                                  seuil))

fX_cond_N_k2 <- sapply(sup_fx_cond,
                       function(y) calculer_fX_cond_N_archi_hiera(2, y, h,
                                                                  function(x) tls_inv_geom_essais(x, q0[1]),
                                                                  function(x) ppois(x, lambdaPois),
                                                                  function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                                                  function(x) pgamma(x, alphaGa, betaGa),
                                                                  function(x) dgeo_shifted(x, 1 - q0[1]),
                                                                  function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                                                  seuil))

fX_cond_N_k3 <- sapply(sup_fx_cond,
                       function(y) calculer_fX_cond_N_archi_hiera(3, y, h,
                                                                  function(x) tls_inv_geom_essais(x, q0[1]),
                                                                  function(x) ppois(x, lambdaPois),
                                                                  function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                                                  function(x) pgamma(x, alphaGa, betaGa),
                                                                  function(x) dgeo_shifted(x, 1 - q0[1]),
                                                                  function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                                                  seuil))

fX_cond_N_k4 <- sapply(sup_fx_cond,
                       function(y) calculer_fX_cond_N_archi_hiera(4, y, h,
                                                                  function(x) tls_inv_geom_essais(x, q0[1]),
                                                                  function(x) ppois(x, lambdaPois),
                                                                  function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                                                  function(x) pgamma(x, alphaGa, betaGa),
                                                                  function(x) dgeo_shifted(x, 1 - q0[1]),
                                                                  function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                                                  seuil))

fX_cond_N_k5 <- sapply(sup_fx_cond,
                       function(y) calculer_fX_cond_N_archi_hiera(5, y, h,
                                                                  function(x) tls_inv_geom_essais(x, q0[1]),
                                                                  function(x) ppois(x, lambdaPois),
                                                                  function(x) tls_inv_geom_essais_comp_geom_essais(x, q0[2], q1[2]),
                                                                  function(x) pgamma(x, alphaGa, betaGa),
                                                                  function(x) dgeo_shifted(x, 1 - q0[1]),
                                                                  function(x, th0) dnbinom_essais(x, th0, q1[2]),
                                                                  seuil))

plot(sup_fx_cond, fX_cond_N_k1, type = "l", lwd = 2, col = "blue", xlab = "ih", ylab = "F.m.p. de X sachant N = k")
title("F.m.p. de X conditionnelle à N = k pour quelques valeurs de k\npour un CRM avec copule Archimédienne hiérarchique AMH-AMH\nde paramètres alpha0 = 0.3 et alpha1 = 0.8")
lines(sup_fx_cond, fX_cond_N_k2, lwd = 2, col = "green")
lines(sup_fx_cond, fX_cond_N_k3, lwd = 2, col = "orange")
lines(sup_fx_cond, fX_cond_N_k4, lwd = 2, col = "purple")
lines(sup_fx_cond, fX_cond_N_k5, lwd = 2, col = "red")
legend(4000, 0.00065, c("f_{X | N = 1} (ih)",
               "f_{X | N = 2} (ih)",
               "f_{X | N = 3} (ih)",
               "f_{X | N = 4} (ih)",
               "f_{X | N = 5} (ih)"),
       col = c("blue", "green", "orange", "purple", "red"),
       lwd = rep(2, 5))


# Échantillonnage
n <- 1000000
ech_l <- echantillonner_crm_cop_archi_hiera(n,
                                            function(x) qgeo_shifted(x, 1 - q0[1]),
                                            function(x) tls_geom_essais(x, q0[1]),
                                            function(x) qgeo_shifted(x, 1 - q1[2]),
                                            function(x) fgp_geom_essais(tls_geom_essais(x, q1[2]), q0[1]),
                                            function(x) qpois(x, lambdaPois),
                                            function(x) qgamma(x, alphaGa, betaGa))
ech <- structurer_echantillon(ech_l$realisations, ech_l$max_N_i)


mean(rowSums(ech[, -1])) # E[S] empirique
var(rowSums(ech[, -1])) # Var(S) empirique

ech_S <- rowSums(ech[, -1])
sorted_ech_S <- sort(ech_S)
sorted_ech_S[0.99 * n]
sorted_ech_S[0.995 * n]


# Échantillonnage - dépendance négative N, vecteur X
n <- 1000000
ech_l_rc <- echantillonner_crm_cop_archi_hiera_rc(n,
                                                  function(x) qgeo_shifted(x, 1 - q0[1]),
                                                  function(x) tls_geom_essais(x, q0[1]),
                                                  function(x) qgeo_shifted(x, 1 - q1[2]),
                                                  function(x) fgp_geom_essais(tls_geom_essais(x, q1[2]), q0[1]),
                                                  function(x) qpois(x, lambdaPois),
                                                  function(x) qgamma(x, alphaGa, betaGa))
ech_rc <- structurer_echantillon(ech_l_rc$realisations, ech_l_rc$max_N_i)


mean(rowSums(ech_rc[, -1])) # E[S] empirique
var(rowSums(ech_rc[, -1])) # Var(S) empirique

ech_S_rc <- rowSums(ech_rc[, -1])
sorted_ech_S <- sort(ech_S_rc)
sorted_ech_S[0.99 * n]
sorted_ech_S[0.995 * n]


# Espérance de la v.a. S, modèle classique
lambdaPois * (alphaGa / betaGa)
