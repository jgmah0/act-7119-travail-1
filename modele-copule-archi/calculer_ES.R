###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de E[S]
### CRM avec copules archimédiennes (modèle à facteur simple)
###
### Exemple :
### calculer_ES_crm_archi_simple(alpha = 0.2, # 0 <= alpha < 1
###                              tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha), # TLS inverse d'une loi géométrique
###                              dDistTheta = function(x, alpha) dgeo_shifted(x, alpha), # Copule AMH
###                              qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
###                              pDistX = function(x) ppois(x, 10), # Sévérité discrète
###                              qDistX = function(x) qpois(x, 10),
###                              dDistN = function(x) dpois(x, 2), # Fréquence poisson
###                              pDistN = function(x) ppois(x, 2),
###                              qDistN = function(x) qpois(x, 2))
###

# E[S]
calculer_ES_crm_archi_simple <- function(alpha, tlsinvDistTheta, dDistTheta,
                                         qDistTheta, pDistX, qDistX, dDistN,
                                         pDistN, qDistN)
{
  epsilon = 1/1000000000 # Erreur acceptable

  # On veut sommer suffisamment de termes
  kmax = qDistN(1 - epsilon)
  jmax = qDistTheta(1 - epsilon, alpha)
  imax = qDistX(1 - epsilon)

  ES = 0 # Initialise E[S]

  for (k in seq(kmax))
  {
    for (j in seq(jmax))
    {
      for (i in seq(imax))
      {
        ES = ES + k *
            dDistTheta(j, alpha) *
            (i * (exp(-j *
                          tlsinvDistTheta(pDistN(k), alpha)) -

                      exp(-j * tlsinvDistTheta(pDistN(k - 1), alpha))) *
                 (exp(-j * tlsinvDistTheta(pDistX(i), alpha)) -

                      exp(-j * tlsinvDistTheta(pDistX(i - 1), alpha))))
      }
    }
  }
  ES
}

# Test
source("modele-copule-archi/geom_shifted.R")
calculer_ES_crm_archi_simple(alpha = 0.2,
                             tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                             dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                             qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                             pDistX = function(x) ppois(x, 10),
                             qDistX = function(x) qpois(x, 10),
                             dDistN = function(x) dpois(x, 2),
                             pDistN = function(x) ppois(x, 2),
                             qDistN = function(x) qpois(x, 2))




