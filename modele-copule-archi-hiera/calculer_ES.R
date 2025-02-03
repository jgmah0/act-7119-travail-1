###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de E[S]
### CRM avec copule Archimédienne hiérarchique
###
### Exemple :
### **À faire
###
###
##

source("calculer_fmp_conjointe_composantes.R")

calculer_EX_cond_N_archi_hiera <- function(k, h, imax,
                                           tlsInvTh0, pDistN, dDistN,
                                           tlsInvTh1, pDistX,
                                           dDistTh0, dDistTh01,
                                           seuil)
{
    EX_cond_N <- 0

    for (i in 0:imax)
    {
        print(i)
        EX_cond_N <- EX_cond_N + i * h * calculer_fNX1Xk_archi_hiera(k, i, h,
                                                                     tlsInvTh0,
                                                                     pDistN,
                                                                     tlsInvTh1,
                                                                     pDistX,
                                                                     dDistTh0,
                                                                     dDistTh01,
                                                                     seuil)
    }

    (1 / dDistN(k)) * EX_cond_N
}

calculer_EX2_cond_N_archi_hiera <- function(k, h, imax,
                                            tlsInvTh0, pDistN, dDistN,
                                            tlsInvTh1, pDistX,
                                            dDistTh0, dDistTh01,
                                            seuil)
{
    EX_cond_N <- 0

    for (i in 0:imax)
    {
        EX_cond_N <- EX_cond_N + ((i * h)^2) *
            calculer_fNX1Xk_archi_hiera(k, i, h,
                                        tlsInvTh0,
                                        pDistN,
                                        tlsInvTh1,
                                        pDistX,
                                        dDistTh0,
                                        dDistTh01,
                                        seuil)
    }

    (1 / dDistN(k)) * EX_cond_N
}


calculer_ES_archi_hiera <- function(kmax, imax, h, tlsInvTh0, pDistN, dDistN,
                                    tlsInvTh1, pDistX, dDistTh0, dDistTh01,
                                    seuil)
{
    ES <- 0

    for (k in 0:kmax)
    {
        print(paste0("k = ", k))
        ES <- ES + k * dDistN(k) * calculer_EX_cond_N_archi_hiera(k, h, imax,
                                                                  tlsInvTh0,
                                                                  pDistN,
                                                                  dDistN,
                                                                  tlsInvTh1,
                                                                  pDistX,
                                                                  dDistTh0,
                                                                  dDistTh01,
                                                                  seuil)
    }

    ES
}
