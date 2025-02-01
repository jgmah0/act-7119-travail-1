###
### Travail 1, ACT-7119
### Algorithme 13 de [Cossette et al., 2019] version non hierarchique
###
###


# À tester.
echantillonner_crm_cop_archi <- function(n, qDistTheta, tlsTheta, qDistB, qDistN, qDistX)
{
  realisations <- list()
  realisation_i <- vector(mode = "numeric")
  
  theta <- qDistTheta(runif(n)) # Échantillonner theta
  R <- rexp(n) # Échantillonner R
  U <- tlsTheta(R / theta) # Calculer U
  
  realisations_N <- qDistN(U) # Simuler les réalisations de N
  
  # Initialiser les valeurs
  R_X <- 0
  
  max_N_i <- 0
  
  for (i in seq(n))
  {
    # Stocker la réalisation i de N
    realisation_i[1] <- realisations_N[i]
    
    if (realisation_i[1] != 0)
    {
      
      # Pour i = 1,...,N
      for (j in seq(realisations_N[i]))
      {
        
        R_X <- rexp(1) # Échantilonner R_i
        X_X = qDistX(U) # Calculer X_I avec le même U que pour la simulation de N
        
        realisation_i <- c(realisation_i, X_X) # Créer le vecteur c(N, X_1,...,X_N)
      }
    }
    
    if (realisation_i[1] > max_N_i)
    {
      max_N_i <- realisation_i[1]
    }
    
    realisations[[i]] <- realisation_i
  }
  
  list(realisations = realisations,
       max_N_i = max_N_i)
}
