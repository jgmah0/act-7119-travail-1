###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de f_{N, X1, ..., Xk} (n, i1*h, ..., ik*h)
### CRM avec copule Archimédienne hiérarchique
###
##


source("../utilitaires/approximer_kmax.R")


#
# Fonction pour calculer la valeur de f_{N | Theta0 = theta0} (k)
#
# Exemple :
# calculer_fN_cond_th0(3,
#                      2
#                      function(x) tlsInv(x, param1),
#                      function(x) ppois(x, 2))
#
calculer_fN_cond_th0_archi_hiera <- function(k, th0, tlsInvTh0, pDistN)
{
    exp(-th0 * tlsInvTh0(pDistN(k))) -
        (k != 0) * exp(-th0 * tlsInvTh0(pDistN(k - 1)))
}


#
# Fonction pour calculer la valeur de
#   f_{X | Theta0 = theta0, Theta01 = theta01} (ih)
#
# Exemple :
# calculer_fX_cond_th0_th01(1,
#                           0.1,
#                           3,
#                           function(x) tlsInv(x, param1),
#                           function(x) ppois(x, 2))
#
calculer_fX_cond_th0_th01_archi_hiera <- function(i, h, th01, tlsInvTh1,
                                                  pDistX)
{
    exp(-th01 * tlsInvTh1(pDistN(i * h))) -
        (i != 0) * exp(-th01 * tlsInvTh1(pDistN((i - 1) * h)))
}


#
# Fonction pour calculer la valeur de
#   f_{N, X1, ..., Xk} (n, i1*h, ..., ik*h)
#
# Exemple :
#
#
calculer_fNX1Xk_archi_hiera <- function(k, vec_i, h,
                                        tlsInvTh0, pDistN,
                                        tlsInvTh1, pDistX,
                                        dDistTh0, dDistTh01,
                                        seuil)
{
    dimension_va <- length(vec_i) + 1
    kmax_th0 <- approximer_kmax_version_2(dDistTh0, seuil)
    kmax_th01 <- approximer_kmax_version_2(dDistTh01, seuil)

    f_conj_comp_partielle <- 0
    for (th01 in seq(kmax_th01))
    {
        f_conj_comp_partielle <- f_conj_comp_partielle +
            dDistTh01(th01) *
            prod(calculer_fX_cond_th0_th01_archi_hiera(vec_i,
                                                       h,
                                                       th01,
                                                       tlsInvTh1,
                                                       pDistX))
    }

    f_conj_comp <- 0
    for (th0 in seq(kmax_th0))
    {
        f_conj_comp <- f_conj_comp +
            dDistTh0(th0) *
            calculer_fN_cond_th0_archi_hiera(k,
                                             th0,
                                             tlsInvTh0,
                                             pDistN)
    }

    f_conj_comp * f_conj_comp_partielle
}

