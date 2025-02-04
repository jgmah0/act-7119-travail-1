###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de E[S]
### CRM indépendance
###
### Exemple :
### calculer_ES_indep(dDistX = function(x) dpois(x, 3),
###                   qDistX = function(x) qpois(x, 3),
###                   dDistN = function(x) dpois(x, 4),
###                   qDistN = function(x) qpois(x, 4))
###

# E[S]
calculer_ES_indep <- function(dDistX, qDistX, dDistN, qDistN)
{
  epsilon = 1/1000000000 # Erreur acceptable

  # On veut sommer suffisamment de termes
  kmax = qDistN(1 - epsilon)
  imax = qDistX(1 - epsilon)
  
  EN = sum(seq(kmax) * sapply(seq(kmax), dDistN))
  EX = sum(seq(imax) * sapply(seq(imax), dDistX))
  
  ES = EN * EX
  
  ES
}

# Test
calculer_ES_indep(dDistX = function(x) dpois(x, 3),
                  qDistX = function(x) qpois(x, 3),
                  dDistN = function(x) dpois(x, 4),
                  qDistN = function(x) qpois(x, 4))




