###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec dépendance selon une copule Archimédienne AMH (1 facteur)
###
##

# On utilise le même code que celui utilisé pour le CRM avec la copule de
# Frank.

source("../modele-copule-archi/geom_shifted.R")
source("../utilitaires/utilitaires.R")
source("../modele-copule-archi/calculer_fS.R")
source("../modele-copule-archi/calculer_fX_condN.R")
source("../modele-copule-archi/echantillonner_crm_cop_archi.R")


lambdaPois <- 2

alphaGa <- 3
betaGa <- 0.002
al_theta <- c(0.1, 0.3, 0.5, 0.7)
nfft <- 2^15


fS_cas1 <- calculer_fS_archi_simple(nfft,
                                    al_theta[1],
                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                    function(x, al) dgeo_shifted(x, al),
                                    function(x, al) qgeo_shifted(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas1)
fS_cas2 <- calculer_fS_archi_simple(nfft,
                                    al_theta[2],
                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                    function(x, al) dgeo_shifted(x, al),
                                    function(x, al) qgeo_shifted(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas2)
fS_cas3 <- calculer_fS_archi_simple(nfft,
                                    al_theta[3],
                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                    function(x, al) dgeo_shifted(x, al),
                                    function(x, al) qgeo_shifted(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas3)
fS_cas4 <- calculer_fS_archi_simple(nfft,
                                    al_theta[4],
                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                    function(x, al) dgeo_shifted(x, al),
                                    function(x, al) qgeo_shifted(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas4)
# Un pas de discrétisation de 1 semble être utilisé dans
# "calculer_fS_archi_simple".
h <- 1


# **** CAS 1
# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas1)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas1) - (sum(((0:(nfft - 1)) * h) * fS_cas1)^2)


sum(fS_cas1[seq((500 / h) + 1)])
sum(fS_cas1[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas1)
VaR_kapp_discr(0.995, h, fS_cas1)

TVaR_kapp_discr_v2(0.99, h, fS_cas1)
TVaR_kapp_discr_v2(0.995, h, fS_cas1)

Mesure_entropique_discr_v2(0.0001, h, fS_cas1)
Mesure_entropique_discr_v2(0.0002, h, fS_cas1)



# **** CAS 2
# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas2)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas2) - (sum(((0:(nfft - 1)) * h) * fS_cas2)^2)


sum(fS_cas2[seq((500 / h) + 1)])
sum(fS_cas2[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas2)
VaR_kapp_discr(0.995, h, fS_cas2)

TVaR_kapp_discr_v2(0.99, h, fS_cas2)
TVaR_kapp_discr_v2(0.995, h, fS_cas2)

Mesure_entropique_discr_v2(0.0001, h, fS_cas2)
Mesure_entropique_discr_v2(0.0002, h, fS_cas2)


# **** CAS 3
# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas3)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas3) - (sum(((0:(nfft - 1)) * h) * fS_cas3)^2)


sum(fS_cas3[seq((500 / h) + 1)])
sum(fS_cas3[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas3)
VaR_kapp_discr(0.995, h, fS_cas3)

TVaR_kapp_discr_v2(0.99, h, fS_cas3)
TVaR_kapp_discr_v2(0.995, h, fS_cas3)

Mesure_entropique_discr_v2(0.0001, h, fS_cas3)
Mesure_entropique_discr_v2(0.0002, h, fS_cas3)


# **** CAS 4
# E[S]
sum(((0:(nfft - 1)) * h) * fS_cas4)

# Var(S)
sum((((0:(nfft - 1)) * h)^2) * fS_cas4) - (sum(((0:(nfft - 1)) * h) * fS_cas4)^2)


sum(fS_cas4[seq((500 / h) + 1)])
sum(fS_cas4[seq((1000 / h) + 1)])

VaR_kapp_discr(0.99, h, fS_cas4)
VaR_kapp_discr(0.995, h, fS_cas4)

TVaR_kapp_discr_v2(0.99, h, fS_cas4)
TVaR_kapp_discr_v2(0.995, h, fS_cas4)

Mesure_entropique_discr_v2(0.0001, h, fS_cas4)
Mesure_entropique_discr_v2(0.0002, h, fS_cas4)


# Fonction de répartition de la v.a. S
plot((0:15000) * h, cumsum(fS_cas1[1:15001]), lwd = 1.4, type = "l", col = "green", xlab = "ih", ylab = "F_S (ih)")
title("Fonctions de répartition de la v.a. S d'un CRM\navec copule AMH de paramètre alpha")
lines((0:15000) * h, cumsum(fS_cas2[1:15001]), lwd = 1.4, col = "blue")
lines((0:15000) * h, cumsum(fS_cas3[1:15001]), lwd = 1.4, col = "purple")
lines((0:15000) * h, cumsum(fS_cas4[1:15001]), lwd = 1.4, col = "orange")
legend(10000, 0.6, c("alpha = 0.1",
                     "alpha = 0.3",
                     "alpha = 0.5",
                     "alpha = 0.7"),
       col = c("green", "blue", "purple", "orange"),
       lwd = rep(1.4, 4))


# f.m.p. conditionnelle à N = k de X
fx_cond_N1 <- sapply(0:6000,
                     function(y) calculer_fx_condN_crm_archi_simple(y, 1, al_theta[2],
                                                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                                                    function(x, al) dgeo_shifted(x, al),
                                                                    function(x, al) qgeo_shifted(x, al),
                                                                    function(x) pgamma(x, alphaGa, betaGa),
                                                                    function(x) dpois(x, lambdaPois),
                                                                    function(x) ppois(x, lambdaPois)))
fx_cond_N2 <- sapply(0:6000,
                     function(y) calculer_fx_condN_crm_archi_simple(y, 2, al_theta[2],
                                                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                                                    function(x, al) dgeo_shifted(x, al),
                                                                    function(x, al) qgeo_shifted(x, al),
                                                                    function(x) pgamma(x, alphaGa, betaGa),
                                                                    function(x) dpois(x, lambdaPois),
                                                                    function(x) ppois(x, lambdaPois)))
fx_cond_N3 <- sapply(0:6000,
                     function(y) calculer_fx_condN_crm_archi_simple(y, 3, al_theta[2],
                                                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                                                    function(x, al) dgeo_shifted(x, al),
                                                                    function(x, al) qgeo_shifted(x, al),
                                                                    function(x) pgamma(x, alphaGa, betaGa),
                                                                    function(x) dpois(x, lambdaPois),
                                                                    function(x) ppois(x, lambdaPois)))

fx_cond_N4 <- sapply(0:6000,
                     function(y) calculer_fx_condN_crm_archi_simple(y, 4, al_theta[2],
                                                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                                                    function(x, al) dgeo_shifted(x, al),
                                                                    function(x, al) qgeo_shifted(x, al),
                                                                    function(x) pgamma(x, alphaGa, betaGa),
                                                                    function(x) dpois(x, lambdaPois),
                                                                    function(x) ppois(x, lambdaPois)))
fx_cond_N5 <- sapply(0:6000,
                     function(y) calculer_fx_condN_crm_archi_simple(y, 5, al_theta[2],
                                                                    function(x, al) tls_inv_geom_essais_alpha(x, al),
                                                                    function(x, al) dgeo_shifted(x, al),
                                                                    function(x, al) qgeo_shifted(x, al),
                                                                    function(x) pgamma(x, alphaGa, betaGa),
                                                                    function(x) dpois(x, lambdaPois),
                                                                    function(x) ppois(x, lambdaPois)))

sup_fx_cond <- 0:6000

plot(sup_fx_cond, fx_cond_N1, type = "l", lwd = 2, col = "blue", xlab = "ih", ylab = "F.m.p. de X sachant N = k")
title("F.m.p. de X conditionnelle à N = k (pour k = 1, ..., 5)\npour un CRM avec copule AMH de paramètre 0.3")
lines(sup_fx_cond, fx_cond_N2, lwd = 2, col = "green")
lines(sup_fx_cond, fx_cond_N3, lwd = 2, col = "orange")
lines(sup_fx_cond, fx_cond_N4, lwd = 2, col = "purple")
lines(sup_fx_cond, fx_cond_N5, lwd = 2, col = "red")
legend(4000, 0.00055, c("f_{X | N = 1} (ih)",
                        "f_{X | N = 2} (ih)",
                        "f_{X | N = 3} (ih)",
                        "f_{X | N = 4} (ih)",
                        "f_{X | N = 5} (ih)"),
       col = c("blue", "green", "orange", "purple", "red"),
       lwd = rep(2, 5))



# Échantillonnage
n <- 1000000
ech_non_structuree <- echantillonner_crm_cop_archi(n,
                                                   function(x) qgeo_shifted(x, al_theta[2]),
                                                   function(x) tls_geom_essais_alpha(x, al_theta[2]),
                                                   function(x) qpois(x, lambdaPois),
                                                   function(x) qgamma(x, alphaGa, betaGa),
                                                   q_theta_not_vectorized = TRUE)

ech_structuree <- structurer_echantillon(ech_non_structuree$realisations,
                                         ech_non_structuree$max_N_i)


mean(rowSums(ech_structuree[, -1])) # E[S] empirique
var(rowSums(ech_structuree[, -1])) # Var(S) empirique

ech_S <- rowSums(ech_structuree[, -1])
sorted_ech_S <- sort(ech_S)
sorted_ech_S[0.99 * n]
sorted_ech_S[0.995 * n]
