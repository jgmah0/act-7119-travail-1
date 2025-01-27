###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de F_S (x)
### CRM avec composantes comonotones
###
### Exemple :
### calculer_FS_crm_comonotonicite(x,
###                                function(x) ppois(x, 2),
###                                function(x) qexp(x, 0.01),
###                                8)
###
##

calculer_FS_crm_comonotonicite <- function(x, pDistN, qDistX, kmax)
{
    FS <- pDistN(0)

    for (i in seq(kmax))
    {
        FS <- FS + integrate(function(y) as.numeric((i * qDistX(y)) <= x), pDistN(i - 1), pDistN(i))$value
    }

    FS
}
