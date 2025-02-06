###
### Travail 1, ACT-7119
### Algorithme 13 de [Cossette et al., 2019]
###
##


echantillonner_crm_cop_archi_hiera <- function(n, qDistTheta0, tlsTheta0, qDistB, tlsTheta1, qDistN, qDistX)
{
  realisations <- list()
  realisation_i <- vector(mode = "numeric")

  theta0 <- qDistTheta0(runif(n)) # Échantillonner theta0
  R <- rexp(n) # Échantillonner R
  U <- tlsTheta0(R / theta0) # Calculer U

  realisations_N <- qDistN(U) # Simuler les réalisations de N

  # Initialiser les valeurs
  R_X <- 0
  U_X <- 0
  Theta01 <- 0

  max_N_i <- 0

  for (i in seq(n))
  {
    realisation_i <- numeric(0)
    print(i)
    # Stocker la réalisation i de N
    realisation_i[1] <- realisations_N[i]

    if (realisation_i[1] != 0)
    {

      Theta01 <- sum(qDistB(runif(theta0[i]))) # Échantilonner Theta01
      R_X <- rexp(realisations_N[i]) # Échantilonner R_i
      U_X <- tlsTheta1(R_X / Theta01) # Calculer U_i
      X_X = qDistX(U_X) # Calculer X_i

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


echantillonner_crm_cop_archi_hiera_rc <- function(n, qDistTheta0, tlsTheta0, qDistB, tlsTheta1, qDistN, qDistX)
{
    realisations <- list()
    realisation_i <- vector(mode = "numeric")

    theta0 <- qDistTheta0(runif(n)) # Échantillonner theta0
    R <- rexp(n) # Échantillonner R
    U <- tlsTheta0(R / theta0) # Calculer U

    realisations_N <- qDistN(1 - U) # Simuler les réalisations de N

    # **La différence avec la fonction ci-dessus est dans la simulation des
    # réalisations de N. "1 - U" à la place de "U".


    # Initialiser les valeurs
    R_X <- 0
    U_X <- 0
    Theta01 <- 0

    max_N_i <- 0

    for (i in seq(n))
    {
        realisation_i <- numeric(0)
        print(i)
        # Stocker la réalisation i de N
        realisation_i[1] <- realisations_N[i]

        if (realisation_i[1] != 0)
        {

            Theta01 <- sum(qDistB(runif(theta0[i]))) # Échantilonner Theta01
            R_X <- rexp(realisations_N[i]) # Échantilonner R_i
            U_X <- tlsTheta1(R_X / Theta01) # Calculer U_i
            X_X = qDistX(U_X) # Calculer X_i

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
