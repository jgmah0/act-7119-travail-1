###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de Var(S)
### CRM avec composantes comonotones
###
### Exemple :
### calculer_VarS_crm_comonotonicite(function(x) ppois(x, 2),
###                                  function(x) qexp(x, 0.01),
###                                  8)
###
##

calculer_VarS_crm_comonotonicite <- function(pDistN, qDistX, kmax)
{
    ES <- 0
    ES2 <- 0

    for (i in seq(kmax))
    {
        ES <- ES + integrate(function(y) i * qDistX(y), pDistN(i - 1), pDistN(i))$value
        ES2 <- ES2 + integrate(function(y) (i * qDistX(y))^2, pDistN(i - 1), pDistN(i))$value
    }

    ES2 - (ES^2)
}

