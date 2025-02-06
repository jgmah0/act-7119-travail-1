###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de mesures de risques
### CRM avec composantes antimonotones
###
### Exemples :
### calculer_VaRS_crm_antimonotones(0.99,
###                                  function(x) ppois(x, 2),
###                                  function(x) qexp(x, 0.01), 15)
###
### calculer_TVaRS_crm_antimonotones(0.99,
###                                   function(x) ppois(x, 2),
###                                   function(x) qexp(x, 0.01))
###
##
source("modele-antimonotone-tot/calculer_FS.R")
source("modele-antimonotone-tot/calculer_ES.R")
calculer_VaRS_crm_antimonotones <- function(k, pDistN, qDistX, kmax)
{
    optimize(function(a) abs(calculer_FS_crm_antimonotones(a, pDistN, qDistX, kmax) - k),
             c(0, 4000))$minimum
}


calculer_TVaRS_crm_antimonotones <- function(k, pDistN, qDistX)
{
    vv <- optimize(function(a) abs(calculer_FS_crm_antimonotones(a,
                                                                 pDistN,
                                                                 qDistX,
                                            15) - k), c(0, 4000))$minimum
    (1 / (1 - k)) *
        calculer_EStronq_crm_antimonotones(vv, pDistN,
                                           qDistX,
                                           15)
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
