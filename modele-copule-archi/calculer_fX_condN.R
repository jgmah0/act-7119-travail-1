###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de f_{X | N = k} (ih)
### CRM avec copules archimédiennes (modèle à facteur simple)
###
### Exemple :
### calculer_fx_condN_crm_archi_simple(i = 1,
###                                    k = 3,
###                                    alpha, # 0 <= alpha < 1
###                                    tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha), # TLS inverse d'une loi géométrique
###                                    dDistTheta = function(x) dgeo_shifted(x, alpha), # Copule AMH
###                                    qDistTheta = function(x) qgeo_shifted(x, alpha),
###                                    pDistX = function(x) ppois(x, 10), # Sévérité discrète
###                                    dDistN = function(x) dpois(x, 2), # Fréquence poisson
###                                    pDistN = function(x) ppois(x, 2))
###

# f_{X | N = k} (ih)
calculer_fx_condN_crm_archi_simple <- function(i, k, alpha, tlsinvDistTheta, dDistTheta, qDistTheta, pDistX, dDistN, pDistN)
{

  epsilon = 1/1000000000 # Erreur acceptable
  jmax = qDistTheta(1 - epsilon, alpha)

  fx_condN = 0 # Initialise f_{X | N = k} (i)

  for (j in seq(jmax))
  {
    fx_condN = fx_condN +
      dDistTheta(j, alpha) * (exp(-j * tlsinvDistTheta(pDistN(k), alpha)) -
                         exp(-j * tlsinvDistTheta(pDistN(k - 1), alpha))) * (exp(-j * tlsinvDistTheta(pDistX(i), alpha)) -
                                                                        exp(-j * tlsinvDistTheta(pDistX(i - 1), alpha)))
  }

  fx_condN = fx_condN / dDistN(k)

  fx_condN

}

# Vérification qu'on somme à 1
# source("modele-copule-archi/geom_shifted.R")
xx = seq(10000) - 1
sum(sapply(xx, function(y) calculer_fx_condN_crm_archi_simple(i = y,
                                                         k = 3,
                                                         alpha = 0.2,
                                                         tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                         dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                                         qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                                         pDistX = function(x) ppois(x, 10),
                                                         dDistN = function(x) dpois(x, 2),
                                                         pDistN = function(x) ppois(x, 2))))
# Somme à 1

### Autres vérifications ----
calculer_fx_condTheta_crm_archi_simple <- function(k, theta, alpha, tlsinvDistTheta, pDistX)
{

  (k == 0)*exp(-theta * tlsinvDistTheta(pDistX(0), alpha)) +
    (k != 0)*(exp(-theta * tlsinvDistTheta(pDistX(k), alpha)) - exp(-theta * tlsinvDistTheta(pDistX(k - 1), alpha)))

}

xx = seq(10000) - 1
sum(sapply(xx, function(y) calculer_fx_condTheta_crm_archi_simple(k = y,
                                                                  theta = 1,
                                                                  alpha = 0.2,
                                                                  tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                                  pDistX = function(x) ppois(x, 10))))
# Somme à 1
calculer_fN_condTheta_crm_archi_simple <- function(k, theta, alpha, tlsinvDistTheta, pDistN)
{

  (k == 0)*exp(-theta * tlsinvDistTheta(pDistN(0), alpha)) +
    (k != 0)*(exp(-theta * tlsinvDistTheta(pDistN(k), alpha)) - exp(-theta * tlsinvDistTheta(pDistN(k - 1), alpha)))

}

xx = seq(10000) - 1
sum(sapply(xx, function(y) calculer_fN_condTheta_crm_archi_simple(k = y,
                                                                  theta = 1,
                                                                  alpha = 0.2,
                                                                  tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                                  pDistN = function(x) ppois(x, 10))))
# Somme à 1

# Version alternative
calculer_fx_condN_crm_archi_simple_alt <- function(i, k, alpha, tlsinvDistTheta, dDistTheta, qDistTheta, pDistX, dDistN, pDistN)
{

  epsilon = 1/1000000000 # Erreur acceptable

  jmax = qDistTheta(1 - epsilon, alpha)

  fx_condN = 0 # Initialise f_{X | N = k} (i)

  for (j in seq(jmax))
  {
    fx_condN = fx_condN + dDistTheta(j, alpha)*calculer_fx_condTheta_crm_archi_simple(i,
                                                                 theta = j,
                                                                 alpha,
                                                                 tlsinvDistTheta,
                                                                 pDistX)*calculer_fN_condTheta_crm_archi_simple(k,
                                                                                                                theta = j,
                                                                                                                alpha,
                                                                                                                tlsinvDistTheta,
                                                                                                                pDistN)
  }

  fx_condN = fx_condN / dDistN(k)

  fx_condN

}

xx = seq(10000) - 1
sum(sapply(xx, function(y) calculer_fx_condN_crm_archi_simple_alt(i = y,
                                                                  k = 3,
                                                                  alpha = 0.2,
                                                                  tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                                  dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                                                  qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                                                  pDistX = function(x) ppois(x, 10),
                                                                  dDistN = function(x) dpois(x, 2),
                                                                  pDistN = function(x) ppois(x, 2))))
# Somme à 1



