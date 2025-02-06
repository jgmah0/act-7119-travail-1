###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de E[S]
### CRM avec composantes antimonotones
###
### Exemple :
### calculer_ES_crm_antimonotones(function(x) ppois(x, 2),
###                                function(x) qexp(x, 0.01),
###                                8)
###
##

calculer_ES_crm_antimonotones <- function(pDistN, qDistX, kmax)
{
    ES <- 0

    for (i in seq(kmax))
    {
        ES <- ES + integrate(function(y) i * qDistX(1 - y), pDistN(i - 1), pDistN(i))$value
    }

    ES
}

calculer_EStronq_crm_antimonotones <- function(d, pDistN, qDistX, kmax)
{
    ES <- 0

    for (i in seq(kmax))
    {
        ES <- ES + integrate(function(y) (i * qDistX(1 - y)) * ((i * qDistX(1 - y)) > d ), pDistN(i - 1), pDistN(i))$value
    }

    ES
}

