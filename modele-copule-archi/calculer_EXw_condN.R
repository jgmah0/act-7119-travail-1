###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de E[X^w|N=k]
### CRM avec copules archimédiennes (modèle à facteur simple)
###
### Exemple :
### calculer_EXw_condN_crm_archi_simple(k = 10,
###                                     w = 2,
###                                     alpha = 0.2,
###                                     tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
###                                     dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
###                                     qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
###                                     pDistX = function(x) ppois(x, 10),
###                                     qDistX = function(x) qpois(x, 10),
###                                     dDistN = function(x) dpois(x, 2),
###                                     pDistN = function(x) ppois(x, 2))
###                              

# E[X^w|N=k]
calculer_EXw_condN_crm_archi_simple <- function(k, w, alpha, tlsinvDistTheta, dDistTheta, qDistTheta, pDistX, qDistX, dDistN, pDistN)
{
  epsilon = 1/1000000000 # Erreur acceptable
  
  # On veut sommer suffisamment de termes
  jmax = qDistTheta(1 - epsilon, alpha)
  imax = qDistX(1 - epsilon)
  
  EXw_condN = 0 # Initialise E[X^w|N=k]
  
  for (j in seq(jmax))
  {
    for (i in seq(imax))
    {
      EXw_condN = EXw_condN + dDistTheta(j, alpha) * i^w * (exp(-j * tlsinvDistTheta(pDistN(k), alpha)) - 
                                                   exp(-j * tlsinvDistTheta(pDistN(k - 1), alpha))) * (exp(-j * tlsinvDistTheta(pDistX(i), alpha)) - 
                                                                                                  exp(-j * tlsinvDistTheta(pDistX(i - 1), alpha)))
    }
  }
  EXw_condN / dDistN(k)
}

# Test
source("modele-copule-archi/geom_shifted.R")
calculer_EXw_condN_crm_archi_simple(k = 10,
                                    w = 2,
                                    alpha = 0.2,
                                    tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                    dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                    qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                    pDistX = function(x) ppois(x, 10),
                                    qDistX = function(x) qpois(x, 10),
                                    dDistN = function(x) dpois(x, 2),
                                    pDistN = function(x) ppois(x, 2))




