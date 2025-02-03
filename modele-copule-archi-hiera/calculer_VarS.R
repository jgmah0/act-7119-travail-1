###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de Var(S)
### CRM avec copule Archimédienne hiérarchique
###
### Exemple :
### **À faire
###
###
##

source("calculer_fmp_conjointe_composantes.R")
source("calculer_ES.R")

calculer_VarX_cond_N_archi_hiera <- function(k, h, imax,
                                             tlsInvTh0, pDistN, dDistN,
                                             tlsInvTh1, pDistX,
                                             dDistTh0, dDistTh01,
                                             seuil)
{
    calculer_EX2_cond_N_archi_hiera(k, h, imax,
                                    tlsInvTh0, pDistN, dDistN,
                                    tlsInvTh1, pDistX,
                                    dDistTh0, dDistTh01,
                                    seuil) -
        calculer_EX_cond_N_archi_hiera(k, h, imax,
                                       tlsInvTh0, pDistN, dDistN,
                                       tlsInvTh1, pDistX,
                                       dDistTh0, dDistTh01,
                                       seuil)
}

calculer_CovX1X2_cond_N_archi_hiera <- function(k, h, imax, tlsInvTh0, pDistN,
                                                dDistN, tlsInvTh1, pDistX,
                                                dDistTh0, dDistTh01, seuil)
{

    EX1X2_cond_N <- 0

    for (i in seq(imax))
    {
        for (j in seq(imax))
        {
            EX1X2_cond_N <- EX1X2_cond_N +
                (i * h) * (j * h) * calculer_fNX1Xk_archi_hiera(k, c(i, j), h,
                                                                tlsInvTh0,
                                                                pDistN,
                                                                tlsInvTh1,
                                                                pDistX,
                                                                dDistTh0,
                                                                dDistTh01,
                                                                seuil)
        }
    }

    (1 / dDistN(k)) * EX1X2_cond_N - (calculer_EX_cond_N_archi_hiera(k, h,
                                                                     imax,
                                                                     tlsInvTh0,
                                                                     pDistN,
                                                                     dDistN,
                                                                     tlsInvTh1,
                                                                     pDistX,
                                                                     dDistTh0,
                                                                     dDistTh01,
                                                                     seuil))^2
}


calculer_VarS_archi_hiera <- function(kmax, imax, h, tlsInvTh0, pDistN, dDistN,
                                      tlsInvTh1, pDistX, dDistTh0, dDistTh01,
                                      seuil)
{
    VarS <- 0
    var_N_EX_cond_N_terme1 <- 0

    for (k in 0:kmax)
    {
        VarS <- VarS +
            dDistN(k) *
            (k *
                 calculer_VarX_cond_N_archi_hiera(k, h, imax,
                                                  tlsInvTh0, pDistN, dDistN,
                                                  tlsInvTh1, pDistX,
                                                  dDistTh0, dDistTh01,
                                                  seuil) +
                 k *
                 (k - 1) *
                 calculer_CovX1X2_cond_N_archi_hiera(k, h, imax, tlsInvTh0,
                                                     pDistN, dDistN, tlsInvTh1,
                                                     pDistX, dDistTh0,
                                                     dDistTh01, seuil))

        var_N_EX_cond_N_terme1 <- var_N_EX_cond_N_terme1 +
            dDistN(k) *
            (k^2) * calculer_EX_cond_N_archi_hiera(k, h, imax,
                                                   tlsInvTh0, pDistN, dDistN,
                                                   tlsInvTh1, pDistX,
                                                   dDistTh0, dDistTh01,
                                                   seuil)

    }

    VarS +
        var_N_EX_cond_N_terme1 -
        (calculer_ES_archi_hiera(kmax, imax, h, tlsInvTh0, pDistN, dDistN,
                                 tlsInvTh1, pDistX, dDistTh0, dDistTh01,
                                 seuil)^2)
}
