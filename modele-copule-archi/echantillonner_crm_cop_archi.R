###
### Travail 1, ACT-7119
### Algorithme 13 de [Cossette et al., 2019] version non hierarchique
###
###


echantillonner_crm_cop_archi <- function(n, qDistTheta, tlsTheta, qDistN, qDistX, q_theta_not_vectorized = FALSE)
{
  realisations <- list()
  realisation_i <- vector(mode = "numeric")
  theta <- numeric(0)

  if (q_theta_not_vectorized)
  {
      theta <- sapply(runif(n), function(y) qDistTheta(y)) # Échantillonner theta
      print(summary(theta))
  }
  else
  {
      theta <- qDistTheta(runif(n)) # Échantillonner theta
  }

  R <- rexp(n) # Échantillonner R
  U <- tlsTheta(R / theta) # Calculer U

  realisations_N <- qDistN(U) # Simuler les réalisations de N

  # Initialiser les valeurs

  max_N_i <- 0

  for (i in seq(n))
  {
    print(i)
    # Stocker la réalisation i de N
    realisation_i <- vector(mode = "numeric")
    realisation_i[1] <- realisations_N[i]

    if (realisation_i[1] != 0)
    {

      # Pour i = 1,...,N
        X_X = qDistX(tlsTheta(rexp(realisations_N[i])/theta[i])) # Calculer X_I avec le même Theta que pour la simulation de N

        realisation_i <- c(realisation_i, X_X) # Créer le vecteur c(N, X_1,...,X_N)

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


structurer_echantillon <- function(ech_list, max_N_i)
{
    n <- length(ech_list)

    realisations <- matrix(rep(0, n * (max_N_i + 1)), nrow = n, ncol = max_N_i + 1)

    for (i in seq(n))
    {
        realisations[i, 1:(ech_list[[i]][1] + 1)] <- ech_list[[i]]
    }

    realisations
}
