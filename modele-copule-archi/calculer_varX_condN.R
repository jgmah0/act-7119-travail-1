###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de Var(X | N = k)
### CRM avec copules archimédiennes (modèle à facteur simple)
###
### Exemple :
### calculer_varX_condN_crm_archi_simple(k = 5,
###                                      alpha = 0.2,
###                                      tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
###                                      dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
###                                      qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
###                                      pDistX = function(x) ppois(x, 10),
###                                      qDistX = function(x) qpois(x, 10),
###                                      dDistN = function(x) dpois(x, 2),
###                                      pDistN = function(x) ppois(x, 2))
###

source("modele-copule-archi/calculer_EXw_condN.R")
calculer_varX_condN_crm_archi_simple <- function(k, alpha, tlsinvDistTheta,
                                                 dDistTheta, qDistTheta,
                                                 pDistX, qDistX, dDistN, pDistN)
{
  calculer_EXw_condN_crm_archi_simple(k,
                                      w = 2,
                                      alpha,
                                      tlsinvDistTheta,
                                      dDistTheta,
                                      qDistTheta,
                                      pDistX,
                                      qDistX,
                                      dDistN,
                                      pDistN) -
    calculer_EXw_condN_crm_archi_simple(k,
                                        w = 1,
                                        alpha,
                                        tlsinvDistTheta,
                                        dDistTheta,
                                        qDistTheta,
                                        pDistX,
                                        qDistX,
                                        dDistN,
                                        pDistN)^2
}

# Test
calculer_varX_condN_crm_archi_simple(k = 5,
                                     alpha = 0.2,
                                     tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                     dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                     qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                     pDistX = function(x) ppois(x, 10),
                                     qDistX = function(x) qpois(x, 10),
                                     dDistN = function(x) dpois(x, 2),
                                     pDistN = function(x) ppois(x, 2))



