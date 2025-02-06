###
### Travail 1, ACT-7119
### Illustration numérique
### CRM avec modèles de dépendance de base
###
##

source("modele-comonotonicite-tot/calculer_ES.R")
source("modele-comonotonicite-tot/calculer_VarS.R")
source("modele-comonotonicite-tot/calculer_FS.R")
source("modele-comonotonicite-tot/calculer_frep_XcondN.R")
source("modele-comonotonicite-tot/calculer_mesures_risque_S.R")
source("modele-comonotonicite-tot/echantillonner_model_3.R")

source("modele-antimonotone-tot/calculer_ES.R")
source("modele-antimonotone-tot/calculer_VarS.R")
source("modele-antimonotone-tot/calculer_FS.R")
source("modele-antimonotone-tot/calculer_frep_XcondN.R")
source("modele-antimonotone-tot/calculer_mesures_risque_S.R")
source("modele-antimonotone-tot/echantillonner_model_4.R")

source("modele-independance/echantillonner_modele_2.R")
source("modele-independance-comonotone/echantillonner_modele_1.R")
source("modele-independance/echantillonner_modele_2.R")

source("utilitaires/approximer_kmax.R")


lambdaPois <- 2
alphaGa <- 3
betaGa <- 0.002

seuil <- 10^(-7)

kmax <- approximer_kmax_version_2(function(x) dpois(x, lambdaPois), 0, seuil)
### CRM Classique
# Espérance
lambdaPois * alphaGa / betaGa

# Variance
lambdaPois * alphaGa / betaGa^2 + (alphaGa / betaGa)^2 * lambdaPois

# Répart
FS <- function(x) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * pgamma(x, (1:kmax) * alphaGa, betaGa))

FS(500)
FS(1000)

# VaR
kap <- c(0.99, 0.995)

VaR <- function(kap) optimize(function(a) abs(FS(a) - kap), c(0, 20000))$minimum

VaR(0.99)
VaR(0.995)

FS(10356.77) #ok
FS(11433.53) #ok

# TVaR Avec la formule de [Cossette et Marceau, 2022]
TVaR <- function(kap)
{
    dpois(0, lambdaPois) +
        sum(dpois(1:kmax, lambdaPois) * (1 / (1 - kap)) *
                (((1:kmax) * alphaGa) / betaGa) *
                pgamma(VaR(kap), ((1:kmax) * alphaGa) + 1,
                       betaGa, lower.tail = FALSE) )

}

sapply(kap, TVaR)
## Mesure entropique
rho <- c(0.0001)

1 / rho * log(exp(lambdaPois * (((betaGa / (betaGa - rho))^alphaGa) - 1)))

## validation par simulation
set.seed(743)
nsim <- 1000000
simul <- echantillonner_modele_2(nsim, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

sim <- structurer_echantillon(simul$realisations, simul$max_N_i)

S <- rowSums(sim[, -1])

mean(S) # E[S]
mean(S^2) - mean(S)^2

mean(S <= 500)
mean(S <= 1000)

So <- sort(S)

v <- c(So[nsim * 0.99], So[nsim * 0.995])
v[1] + 1 / (1 - 0.99) * mean(pmax(S - v[1], 0))
v[2] + 1 / (1 - 0.995) * mean(pmax(S - v[2], 0))

rho <- c(0.0001)
log(mean(exp(rho[1] * S)))/rho[1]

### Modèle avec indépendance et comonotonicité

# Espérance comme dans le modèle classique

# Variance
(lambdaPois + lambdaPois^2) * alphaGa / betaGa^2 + (alphaGa / betaGa)^2 * lambdaPois

# Fonction de répartition
FS <- function(x) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * pgamma(x / 1:kmax, alphaGa, betaGa))

FS(500)
FS(1000)

# VaR
(v0.99 <- optimize(function(a) abs(FS(a) - 0.99), c(0, 20000))$minimum)
(v0.995 <- optimize(function(a) abs(FS(a) - 0.995), c(0, 20000))$minimum)

estronq_d <- function(d) alphaGa / betaGa * pgamma(d, alphaGa + 1, betaGa,
                                                   lower.tail = FALSE)

sum(dpois(1:kmax, lambdaPois) * (1:kmax) *
        estronq_d(v0.99/(1:kmax))) / (1 - 0.99)

sum(dpois(1:kmax, lambdaPois) * (1:kmax) *
        estronq_d(v0.995/(1:kmax))) / (1 - 0.995)

ms <- function(t) sum(dpois(0:kmax, lambdaPois) *
                          (betaGa / (betaGa - t * (0:kmax)))^{alphaGa})

log(ms(0.0001)) /  0.0001

## validation par simulation
set.seed(743)
nsim <- 1000000
simul <- echantillonner_modele_1(nsim, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

sim <- structurer_echantillon(simul$realisations, simul$max_N_i)

S <- rowSums(sim[, -1])

mean(S) # E[S]
mean(S^2) - mean(S)^2

mean(S <= 500)
mean(S <= 1000)

So <- sort(S)

v <- c(So[nsim * 0.99], So[nsim * 0.995])
v[1] + 1 / (1 - 0.99) * mean(pmax(S - v[1], 0))
v[2] + 1 / (1 - 0.995) * mean(pmax(S - v[2], 0))

rho <- c(0.0001)
log(mean(exp(rho[1] * S)))/rho[1]


### Modèle avec composantes comonotones
calculer_ES_crm_comonotonicite(function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_VarS_crm_comonotonicite(function(x) ppois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa),
                                 kmax=kmax)

calculer_FS_crm_comonotonicite(500, function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_FS_crm_comonotonicite(1000, function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_VaRS_crm_comonotonicite(0.99, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))
calculer_VaRS_crm_comonotonicite(0.995, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

calculer_TVaRS_crm_comonotonicite(0.99, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))
calculer_TVaRS_crm_comonotonicite(0.995, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

calculer_entropique_crm_comonotonicite(0.0001, function(x) ppois(x, lambdaPois),
                                       function(x) qpois(x, lambdaPois),
                                       function(x) qgamma(x, alphaGa, betaGa),
                                       kmax)
## validation par simulation
set.seed(743)
nsim <- 1000000
simul <- echantillonner_modele_3(nsim, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

sim <- structurer_echantillon(simul$realisations, simul$max_N_i)

S <- rowSums(sim[, -1])

mean(S) # E[S]
mean(S^2) - mean(S)^2

mean(S <= 500)
mean(S <= 1000)

So <- sort(S)

v <- c(So[nsim * 0.99], So[nsim * 0.995])
v[1] + 1 / (1 - 0.99) * mean(pmax(S - v[1], 0))
v[2] + 1 / (1 - 0.995) * mean(pmax(S - v[2], 0))

rho <- c(0.0001)
log(mean(exp(rho[1] * S)))/rho[1]


### Modèle avec composantes antimonotones
calculer_ES_crm_antimonotones(function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_VarS_crm_antimonotones(function(x) ppois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa),
                                 kmax=kmax)

calculer_FS_crm_antimonotones(500, function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_FS_crm_antimonotones(1000, function(x) ppois(x, lambdaPois),
                               function(x) qgamma(x, alphaGa, betaGa),
                               kmax=kmax)

calculer_VaRS_crm_antimonotones(0.99, function(x) ppois(x, lambdaPois),
                                function(x) qgamma(x, alphaGa, betaGa),
                                kmax=kmax)

calculer_VaRS_crm_antimonotones(0.995, function(x) ppois(x, lambdaPois),
                                function(x) qgamma(x, alphaGa, betaGa),
                                kmax=kmax)

calculer_TVaRS_crm_antimonotones(0.99, function(x) ppois(x, lambdaPois),
                                  function(x) qgamma(x, alphaGa, betaGa))

calculer_TVaRS_crm_antimonotones(0.995, function(x) ppois(x, lambdaPois),
                                  function(x) qgamma(x, alphaGa, betaGa))

calculer_entropique_crm_antimonotones(0.0001, function(x) ppois(x, lambdaPois),
                                       function(x) qpois(x, lambdaPois),
                                       function(x) qgamma(x, alphaGa, betaGa),
                                       kmax)

## validation par simulation
set.seed(743)
nsim <- 1000000
simul <- echantillonner_modele_4(nsim, function(x) qpois(x, lambdaPois),
                                 function(x) qgamma(x, alphaGa, betaGa))

sim <- structurer_echantillon(simul$realisations, simul$max_N_i)

S <- rowSums(sim[, -1])

mean(S) # E[S]
mean(S^2) - mean(S)^2

mean(S <= 500)
mean(S <= 1000)

So <- sort(S)

v <- c(So[nsim * 0.99], So[nsim * 0.995])
v[1] + 1 / (1 - 0.99) * mean(pmax(S - v[1], 0))
v[2] + 1 / (1 - 0.995) * mean(pmax(S - v[2], 0))

rho <- c(0.0001)
log(mean(exp(rho[1] * S)))/rho[1]


##### Densités conditionnelles #####
### Avec composantes comonotones
fxsachantn <- function(x, k, lam, alp, bet)
{
    dgamma(x, alp, bet) / (ppois(k, lam) - ppois(k - 1, lam)) *
        (qgamma(ppois(k - 1, lam), alp, bet) < x) *
        (qgamma(ppois(k, lam), alp, bet) > x)
}

# Couleurs
colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "turquoise",
            "darkgreen", "black")

# Première
curve(fxsachantn(x, 1, lambdaPois, alphaGa, betaGa), xlim = c(0, 8000), n = 1e3, col = colors[1],
      ylab = "Densité",
      xlab = "Montant d'un sinistre x",
      main = "Densité du montant de sinistre X sachant N = k\n(Comonotonicité entre la fréquence et la sévérité)")

# Autres
for (i in 2:9) {
    curve(fxsachantn(x, i, lambdaPois, alphaGa, betaGa), xlim = c(0, 8000), n = 1e3, col = colors[i], add = TRUE)
}

curve(dgamma(x, alphaGa, betaGa), xlim = c(0, 8000), col = "black", add = TRUE)
# Légende
legend_labels <- sapply(1:10, function(i) bquote(f[X|N == .(i)]))

# Légende
legend("topright", legend = c(paste0("k = ", 1:9), "indep"), col = colors, lty = 1, cex = 0.8)

### Avec composantes comonotones
Fxsachantn <- function(x, k, lam, alp, bet)
{

    if (x > qgamma(ppois(k, lam), alp,bet))
        return(1)

    (pgamma(x, alp, bet) - ppois(k-1, lam)) / (ppois(k, lam) - ppois(k - 1, lam)) *
        (qgamma(ppois(k - 1, lam), alp, bet) < x)
}

# Couleurs
colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "turquoise",
            "darkgreen", "black")

# Première
curve(sapply(x, function(x) Fxsachantn(x, 1, lambdaPois, alphaGa, betaGa)), xlim = c(0, 10000), ylim = c(0, 1), n = 1e3, col = colors[1],
      ylab = "Fonction de répartition",
      xlab = "Montant d'un sinistre x",
      main = "Répartition montant de sinistre X sachant N = k\n(Comonotonicité entre la fréquence et la sévérité)")

# Autres
for (i in 2:9) {
    curve(sapply(x, function(x) Fxsachantn(x, i, lambdaPois, alphaGa, betaGa)), xlim = c(0, 10000), n = 1e3, col = colors[i], add = TRUE)
}
curve(pgamma(x, alphaGa, betaGa), xlim = c(0, 10000), col = "black", add = TRUE)
# Légende
legend_labels <- sapply(1:9, function(i) bquote(f[X|N == .(i)]))

# Légende
legend("topright", legend = c(paste0("k = ", 1:9), "indep"), col = colors, lty = 1, cex = 0.8)

### Avec composantes antimonotones
fxsachantn <- function(x, k, lam, alp, bet)
{
    dgamma(x, alp, bet) / (ppois(k, lam) - ppois(k - 1, lam)) *
        (qgamma(1 - ppois(k, lam), alp, bet) <= x) *
        (qgamma(1 - ppois(k - 1, lam), alp, bet) > x)
}

# Couleurs
colors <- c("red", "blue", "darkgreen", "purple", "black")

# Première
curve(fxsachantn(x, 4, lambdaPois, alphaGa, betaGa), xlim = c(0, 8000), n = 1e4, col = colors[4],
      ylab = "Densité",
      xlab = "Montant d'un sinistre x",
      main = "Densité du montant de sinistre X sachant N = k\n(Antimonotonicité entre la fréquence et la sévérité)",
      lwd = 2)


# Autres
for (i in 1:3) {
    curve(fxsachantn(x, i, lambdaPois, alphaGa, betaGa), xlim = c(0, 8000), n = 1e4, col = colors[i],
          add = TRUE, lwd = 2)
}

curve(dgamma(x, alphaGa, betaGa), xlim = c(0, 8000), col = "black", add = TRUE, lwd = 2)

# Légende
legend_labels <- sapply(1:4, function(i) bquote(f[X|N == .(i)]))

# Légende
legend("topright", legend = c(paste0("k = ", 1:4), "indep"), col = colors, lty = 1, cex = 0.8,
       lwd =2)
### Répartition conditionnelle
Fxsachantn <- function(x, k, lam, alp, bet)
{
    if (x > qgamma(1 - ppois(k - 1, lam), alp, bet))
        return(1)

    (pgamma(x, alp, bet) + ppois(k, lam) - 1) / (ppois(k, lam) - ppois(k - 1, lam)) *
        (qgamma(1 - ppois(k, lam), alp, bet) < x)
}

# Couleurs
colors <- c("red", "blue", "darkgreen", "purple", "black")

# Première
curve(sapply(x, function(x) Fxsachantn(x, 4, lambdaPois, alphaGa, betaGa)), xlim = c(0, 8000), n = 1e4, col = colors[4],
      ylab = "Densité",
      xlab = "Montant d'un sinistre x",
      main = "Répartition du montant de sinistre X sachant N = k\n(Antimonotonicité entre la fréquence et la sévérité)",
      lwd = 2)

# Autres
for (i in 1:3) {
    curve(sapply(x, function(x) Fxsachantn(x, i, lambdaPois, alphaGa, betaGa)), xlim = c(0, 8000), n = 1e4, col = colors[i],
          add = TRUE, lwd = 2)
}
curve(pgamma(x, alphaGa, betaGa), xlim = c(0, 8000), col = "black", add = TRUE, lwd = 2)

# Légende
legend_labels <- sapply(1:4, function(i) bquote(f[X|N == .(i)]))

# Légende
legend("topright", legend = c(paste0("k = ", 1:4), "indep"), col = colors, lty = 1, cex = 0.8,
       lwd =2)

###
FS_classique <- function(x) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * pgamma(x, (1:kmax) * alphaGa, betaGa))
FS_classique_comono <- function(x) dpois(0, lambdaPois) + sum(dpois(1:kmax, lambdaPois) * pgamma(x / 1:kmax, alphaGa, betaGa))

curve(sapply(x, function(x) FS_classique(x)), xlim = c(0, 30000), col = "black", ylim= c(0, 1),
      main = "Fonctions de répartition de S", ylab = "Probabilité", xlab = "Perte totale s", lwd = 2)

curve(sapply(x, function(x) FS_classique_comono(x)), xlim = c(0, 30000), add = TRUE, col = "red", lwd = 2)

curve(sapply(x, function(x) calculer_FS_crm_comonotonicite(x, function(a) ppois(a, lambdaPois),
                                                           function(a) qgamma(a, alphaGa, betaGa),
                                                           kmax=kmax)), xlim = c(0, 30000), add = TRUE, col = "purple", lwd = 2)

curve(sapply(x, function(x) calculer_FS_crm_antimonotones(x, function(a) ppois(a, lambdaPois),
                                                           function(a) qgamma(a, alphaGa, betaGa),
                                                           kmax=kmax)), xlim = c(0, 30000), add = TRUE, col = "green", lwd = 2)

abline(h = 0.99, lty = 4)

abline(v = 3032, lty = 4, col = "green")
abline(v = 10357, lty = 4, col = "black")
abline(v = 13866, lty = 4, col = "red")
abline(v = 25218, lty = 4, col = "purple")
# Légende
legend("bottomright", legend = c("Classique", "X+", "S+", "S-", "VaR_0.99"), col = c("black", "red", "purple", "green"), lty = c(1, 1, 1, 1, 4), cex = 0.8,
       lwd =c(2, 2, 2, 2, 1))
