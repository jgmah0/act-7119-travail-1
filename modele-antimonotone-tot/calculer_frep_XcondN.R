###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de F_{X | N = k} (x)
### CRM avec composantes antimonotones
###
### Exemple :
### calculer_frep_XcondN_crm_antimonotones(1000,
###                                         3,
###                                         function(x) ppois(x, 2),
###                                         function(x) pexp(x, 0.01))
###
##

calculer_frep_XcondN_crm_antimonotones <- function(x, k, pDistN, pDistX)
{
    (1 - pDistX(x) <= pDistN(k - 1)) * 1 +
        ( (pDistN(k - 1) < 1 - pDistX(x)) &
              (1 - pDistX(x) <= pDistN(k)) ) *
        ( (pDistN(k) + pDistX(x) - 1) / (pDistN(k) - pDistN(k - 1)) ) +
        (pDistN(k) < 1 - pDistX(x)) * 0
}
