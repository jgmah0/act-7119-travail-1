###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec dépendance selon une copule Archimédienne de Clayton (1 facteur)
###
##

source("../modele-copule-archi/geom_shifted.R")
source("../utilitaires/utilitaires.R")
source("../modele-copule-archi/calculer_fS.R")
source("../modele-copule-archi/calculer_fX_condN.R")
source("../modele-copule-archi/echantillonner_crm_cop_archi.R")



# Fonctions pour un exemple numérique avec simulation
esp_empirique <- function(realisations)
{
    mean(realisations)
}

variance_empirique <- function(realisations)
{
    mean(realisations^2) - (mean(realisations)^2)
}

fdensite_empirique <- function(x, realisations)
{
    sum((x - (10^-1) <= realisations) & (x + (10^-1) >= realisations)) / length(realisations)
}

frep_empirique <- function(x, realisations)
{
    sum(realisations <= x) / length(realisations)
}

VaR_empirique <- function(kap, realisations)
{
    realisations_ordonnees <- sort(realisations)

    realisations_ordonnees[ceiling(kap*length(realisations))]
}

TVaR_empirique <- function(kap, realisations)
{
    VaR_kappa <- VaR_empirique(kap, realisations)

    VaR_kappa + (1 / (1 - kap)) * mean(pmax(realisations - VaR_kappa, 0))
}

mesure_entropique_empirique <- function(rho, realisations)
{
    (1 / rho) * log(mean(exp(rho * realisations)))
}



fx_cond_n_preparation_donnees <- function(k, n, realisations_non_structurees)
{
    realisations_n_eq_k <- numeric(0)
    for (i in seq(n))
    {
        if (unlist(realisations_non_structurees[[1]][i])[1] == k)
        {
            realisations_n_eq_k <- c(realisations_n_eq_k, unlist(realisations_non_structurees[[1]][i])[-1])
        }
    }

    realisations_n_eq_k
}# utiliser ensuite f_empirique


fx_cond_n_preparation_donnees <- function(k, n, realisations_non_structurees)
{
    realisations_n_eq_k <- numeric(0)
    for (i in seq(n))
    {
        if (unlist(realisations_non_structurees[[1]][i])[1] == k)
        {
            realisations_n_eq_k <- c(realisations_n_eq_k, unlist(realisations_non_structurees[[1]][i])[-1])
        }
    }

    realisations_n_eq_k
}



# Paramètres
lambdaPois <- 2
alphaGa <- 3
betaGa <- 0.002
al_theta <- c(0.2, 0.5, 1, 3)
n <- 1000000


# Échantillonnage 1
ech_non_structuree_cas1 <- echantillonner_crm_cop_archi(n,
                                                        function(x) qgamma(x, 1 / al_theta[1], 1),
                                                        function(x) tls_gamma(x, 1 / al_theta[1], 1),
                                                        function(x) qpois(x, lambdaPois),
                                                        function(x) qgamma(x, alphaGa, betaGa),
                                                        q_theta_not_vectorized = TRUE)

ech_structuree_cas1 <- structurer_echantillon(ech_non_structuree_cas1$realisations,
                                              ech_non_structuree_cas1$max_N_i)

ech_S_cas1 <- rowSums(ech_structuree_cas1[, -1])

esp_empirique(ech_S_cas1)
variance_empirique(ech_S_cas1)
frep_empirique(500, ech_S_cas1)
frep_empirique(1000, ech_S_cas1)
VaR_empirique(0.99, ech_S_cas1)
VaR_empirique(0.995, ech_S_cas1)
TVaR_empirique(0.99, ech_S_cas1)
TVaR_empirique(0.995, ech_S_cas1)
mesure_entropique_empirique(0.0001, ech_S_cas1)
mesure_entropique_empirique(0.0002, ech_S_cas1)


# Échantillonnage 2
ech_non_structuree_cas2 <- echantillonner_crm_cop_archi(n,
                                                        function(x) qgamma(x, 1 / al_theta[2], 1),
                                                        function(x) tls_gamma(x, 1 / al_theta[2], 1),
                                                        function(x) qpois(x, lambdaPois),
                                                        function(x) qgamma(x, alphaGa, betaGa),
                                                        q_theta_not_vectorized = TRUE)

ech_structuree_cas2 <- structurer_echantillon(ech_non_structuree_cas2$realisations,
                                              ech_non_structuree_cas2$max_N_i)

ech_S_cas2 <- rowSums(ech_structuree_cas2[, -1])

esp_empirique(ech_S_cas2)
variance_empirique(ech_S_cas2)
frep_empirique(500, ech_S_cas2)
frep_empirique(1000, ech_S_cas2)
VaR_empirique(0.99, ech_S_cas2)
VaR_empirique(0.995, ech_S_cas2)
TVaR_empirique(0.99, ech_S_cas2)
TVaR_empirique(0.995, ech_S_cas2)
mesure_entropique_empirique(0.0001, ech_S_cas2)
mesure_entropique_empirique(0.0002, ech_S_cas2)


# Échantillonnage 3
ech_non_structuree_cas3 <- echantillonner_crm_cop_archi(n,
                                                        function(x) qgamma(x, 1 / al_theta[3], 1),
                                                        function(x) tls_gamma(x, 1 / al_theta[3], 1),
                                                        function(x) qpois(x, lambdaPois),
                                                        function(x) qgamma(x, alphaGa, betaGa),
                                                        q_theta_not_vectorized = TRUE)

ech_structuree_cas3 <- structurer_echantillon(ech_non_structuree_cas3$realisations,
                                              ech_non_structuree_cas3$max_N_i)

ech_S_cas3 <- rowSums(ech_structuree_cas3[, -1])



esp_empirique(ech_S_cas3)
variance_empirique(ech_S_cas3)
frep_empirique(500, ech_S_cas3)
frep_empirique(1000, ech_S_cas3)
VaR_empirique(0.99, ech_S_cas3)
VaR_empirique(0.995, ech_S_cas3)
TVaR_empirique(0.99, ech_S_cas3)
TVaR_empirique(0.995, ech_S_cas3)
mesure_entropique_empirique(0.0001, ech_S_cas3)
mesure_entropique_empirique(0.0002, ech_S_cas3)



# Échantillonnage 4
ech_non_structuree_cas4 <- echantillonner_crm_cop_archi(n,
                                                        function(x) qgamma(x, 1 / al_theta[4], 1),
                                                        function(x) tls_gamma(x, 1 / al_theta[4], 1),
                                                        function(x) qpois(x, lambdaPois),
                                                        function(x) qgamma(x, alphaGa, betaGa),
                                                        q_theta_not_vectorized = TRUE)

ech_structuree_cas4 <- structurer_echantillon(ech_non_structuree_cas4$realisations,
                                              ech_non_structuree_cas4$max_N_i)

ech_S_cas4 <- rowSums(ech_structuree_cas4[, -1])


esp_empirique(ech_S_cas4)
variance_empirique(ech_S_cas4)
frep_empirique(500, ech_S_cas4)
frep_empirique(1000, ech_S_cas4)
VaR_empirique(0.99, ech_S_cas4)
VaR_empirique(0.995, ech_S_cas4)
TVaR_empirique(0.99, ech_S_cas4)
TVaR_empirique(0.995, ech_S_cas4)
mesure_entropique_empirique(0.0001, ech_S_cas4)
mesure_entropique_empirique(0.0002, ech_S_cas4)


# Fonction de répartion de la v.a. S
supS <- seq(0, 15000, by = 50)
plot(supS, sapply(supS, function(y) frep_empirique(y, ech_S_cas1)), lwd = 1.4, type = "l", col = "green", xlab = "x", ylab = "F_S (x)")
title("Fonctions de répartition de la v.a. S d'un CRM avec dépendance\nselon une copule de Clayton pour quatre valeur de paramètres de dépendance")
lines(supS, sapply(supS, function(y) frep_empirique(y, ech_S_cas2)), lwd = 1.4, col = "blue")
lines(supS, sapply(supS, function(y) frep_empirique(y, ech_S_cas3)), lwd = 1.4, col = "purple")
lines(supS, sapply(supS, function(y) frep_empirique(y, ech_S_cas4)), lwd = 1.4, col = "orange")
legend(12000, 0.6, c("alpha = 0.2",
                     "alpha = 0.5",
                     "alpha = 1",
                     "alpha = 3"),
       col = c("green", "blue", "purple", "orange"),
       lwd = rep(1.4, 4))



# fX sachant N = k
donnees_X_sachant_N1 <- fx_cond_n_preparation_donnees(1, n, ech_non_structuree_cas2)
donnees_X_sachant_N2 <- fx_cond_n_preparation_donnees(2, n, ech_non_structuree_cas2)
donnees_X_sachant_N3 <- fx_cond_n_preparation_donnees(3, n, ech_non_structuree_cas2)
donnees_X_sachant_N4 <- fx_cond_n_preparation_donnees(4, n, ech_non_structuree_cas2)
donnees_X_sachant_N5 <- fx_cond_n_preparation_donnees(5, n, ech_non_structuree_cas2)

plot(density(donnees_X_sachant_N1, bw = 50, kernel = "gaussian"), type = "l", lwd = 2, col = "blue", xlab = "x", xlim = c(0, 6000),
     ylab = "F.m.p. de X sachant N = k",
     main = "Fonctions de densité (estimées par la méthode des noyaux) de X conditionnelles\nà N = k (pour k = 1, ..., 5) pour un CRM avec dépendance modélisée selon\nune copule de Clayton de paramètre alpha = 0.5")
lines(density(donnees_X_sachant_N2, bw = 50, kernel = "gaussian"), lwd = 2, col = "green")
lines(density(donnees_X_sachant_N3, bw = 50, kernel = "gaussian"), lwd = 2, col = "orange")
lines(density(donnees_X_sachant_N4, bw = 50, kernel = "gaussian"), lwd = 2, col = "purple")
lines(density(donnees_X_sachant_N5, bw = 50, kernel = "gaussian"), lwd = 2, col = "red")
legend(4000, 0.00055, c("f_{X | N = 1} (ih)",
                        "f_{X | N = 2} (ih)",
                        "f_{X | N = 3} (ih)",
                        "f_{X | N = 4} (ih)",
                        "f_{X | N = 5} (ih)"),
       col = c("blue", "green", "orange", "purple", "red"),
       lwd = rep(2, 5))
