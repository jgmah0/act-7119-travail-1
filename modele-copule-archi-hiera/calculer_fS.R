###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de f_S (ih)
### CRM avec copule Archimédienne hiérarchique
###
### Exemple :
### **À faire
###
###
##

# source("calculer_fmp_conjointe_composantes.R")

# Fmp de S au complet selon taille nfft.
calculer_fS_archi_hiera <- function(nfft, kmax, h, tlsInvTh0, pDistN,
                                    dDistN, tlsInvTh1, pDistX, dDistTh0,
                                    dDistTh01, seuil)
{
    fS <- numeric(nfft)
    support_i <- 0:(nfft - 1)

    kmax_th0 <- approximer_kmax_version_2(dDistTh0, borne_inf_support = 1, threshold = seuil)

    for (k in 0:kmax)
    {
        print(k)
        for (th0 in seq(kmax_th0))
        {
            prN_k_sachant_th0 <- calculer_fN_cond_th0_archi_hiera(k,
                                                                  th0,
                                                                  tlsInvTh0,
                                                                  pDistN)
            pr_th0 <- dDistTh0(th0)

            # Est-ce que cela est toujours vrai que le support commence à th0 pour Theta01 toujours?
            # Je crois que oui, car B est strictement positif et "theta0" termes
            # sont sommés.
            kmax_th01 <- approximer_kmax_version_2(function(x) dDistTh01(x, th0), borne_inf_support = th0, threshold = seuil)
            for (th01 in seq(kmax_th01))
            {
                pr_th01 <- dDistTh01(th01, th0)

                fx_sachant_th0_th01 <- calculer_fX_cond_th0_th01_archi_hiera(support_i, h, th01, tlsInvTh1, pDistX)
                fx_sachant_th0_th01_convolution_k <- Re(fft(fft(fx_sachant_th0_th01)^k, inverse = TRUE)) / nfft

                fS <- fS +
                    fx_sachant_th0_th01_convolution_k *
                    pr_th0 *
                    pr_th01 *
                    prN_k_sachant_th0
            }
        }
    }

    fS
}

# calculer_fN_cond_th0_archi_hiera(k, th0, tlsInvTh0, pDistN)
# calculer_fX_cond_th0_th01_archi_hiera(i, h, th01, tlsInvTh1, pDistX)
