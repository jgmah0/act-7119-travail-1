###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de F_{X | N = k} (x)
### CRM avec composantes comonotones
###
### Exemple :
### calculer_frep_XcondN_crm_comonotonicite(1000,
###                                         3,
###                                         function(x) ppois(x, 2),
###                                         function(x) pexp(x, 0.01))
###
##

calculer_frep_XcondN_crm_comonotonicite <- function(x, k, pDistN, pDistX)
{
    (pDistX(x) < pDistN(k - 1)) * 0 +
        ( (pDistN(k - 1) <= pDistX(x)) &
              (pDistX(x) <= pDistN(k)) ) *
        ( (pDistX(x) - pDistN(k - 1)) / (pDistN(k) - pDistN(k - 1)) ) +
        (pDistN(k) < pDistX(x)) * 1
}
