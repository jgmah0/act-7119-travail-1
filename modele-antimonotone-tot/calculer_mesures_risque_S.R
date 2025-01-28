###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de mesures de risques
### CRM avec composantes antimonotones
###
### Exemples :
### calculer_VaRS_crm_antimonotones(0.99,
###                                  function(x) qpois(x, 2),
###                                  function(x) qexp(x, 0.01))
###
### calculer_TVaRS_crm_antimonotones(0.99,
###                                   function(x) qpois(x, 2),
###                                   function(x) qexp(x, 0.01))
###
##

calculer_VaRS_crm_antimonotones <- function(k, qDistN, qDistX)
{
    qDistN(k) * qDistX(1 - k)
}


calculer_TVaRS_crm_antimonotones <- function(k, qDistN, qDistX)
{
    (1 / (1 - k)) *
        integrate(function(x) calculer_VaRS_crm_antimonotones(x,
                                                               qDistN,
                                                               qDistX),
                  k,
                  1)$value
}

calculer_entropique_crm_antimonotones <- function(rho, pDistN, qDistN, qDistX, kmax)
{
    MS <- pDistN(0)

    for (i in seq(kmax))
    {
        MS <- MS + integrate(function(y) exp(rho * qDistN(y) * qDistX(1 - y)),
                             pDistN(i - 1), pDistN(i))$value
    }

    1 / rho * log(MS)
}
