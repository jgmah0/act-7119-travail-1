###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec dépendance selon une copule Archimédienne de Frank (1 facteur)
###
##

source("../modele-copule-archi/geom_shifted.R")
source("../utilitaires/utilitaires.R")
source("../modele-copule-archi/calculer_fS.R")


lambdaPois <- 2

alphaGa <- 3
betaGa <- 0.002
al_theta <- c(0.1, 0.2, 0.5, 0.8)
nfft <- 2^15


fS_cas1 <- calculer_fS_archi_simple(nfft,
                                    al_theta[1],
                                    function(x, al) tls_inv_logarithmic(x, al),
                                    function(x, al) dlogarithmique_alpha(x, al),
                                    function(x, al) qlogarithmique_alpha(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas1)
fS_cas2 <- calculer_fS_archi_simple(nfft,
                                    al_theta[2],
                                    function(x, al) tls_inv_logarithmic(x, al),
                                    function(x, al) dlogarithmique_alpha(x, al),
                                    function(x, al) qlogarithmique_alpha(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas2)
fS_cas3 <- calculer_fS_archi_simple(nfft,
                                    al_theta[3],
                                    function(x, al) tls_inv_logarithmic(x, al),
                                    function(x, al) dlogarithmique_alpha(x, al),
                                    function(x, al) qlogarithmique_alpha(x, al),
                                    function(x) pgamma(x, alphaGa, betaGa),
                                    function(x) qpois(x, lambdaPois),
                                    function(x) ppois(x, lambdaPois))
sum(fS_cas3)
fS_cas4 <- calculer_fS_archi_simple(nfft,
                                    al_theta[4],
                                    function(x, al) tls_inv_logarithmic(x, al),
                                    function(x, al) dlogarithmique_alpha(x, al),
                                    function(x, al) qlogarithmique_alpha(x, al),
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

Mesure_entropique_discr_v2(0.001, h, fS_cas1)
Mesure_entropique_discr_v2(0.0001, h, fS_cas1)



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

Mesure_entropique_discr_v2(0.001, h, fS_cas2)
Mesure_entropique_discr_v2(0.0001, h, fS_cas2)


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

Mesure_entropique_discr_v2(0.001, h, fS_cas3)
Mesure_entropique_discr_v2(0.0001, h, fS_cas3)


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

Mesure_entropique_discr_v2(0.001, h, fS_cas4)
Mesure_entropique_discr_v2(0.0001, h, fS_cas4)


