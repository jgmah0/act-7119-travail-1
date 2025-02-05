###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de f_S (ih)
### CRM avec copules archimédiennes (modèle à facteur simple)
###
### Exemple :
### calculer_fS_archi_simple(nfft = 2^9,
###                                    alpha, # 0 <= alpha < 1
###                                    tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha), # TLS inverse d'une loi géométrique
###                                    dDistTheta = function(x) dgeo_shifted(x, alpha), # Copule AMH
###                                    qDistTheta = function(x) qgeo_shifted(x, alpha),
###                                    pDistX = function(x) ppois(x, 10), # Sévérité discrète
###                                    dDistN = function(x) dpois(x, 2), # Fréquence poisson
###                                    pDistN = function(x) ppois(x, 2))
###

# f_{S} (ih)
calculer_fS_archi_simple <- function(nfft, alpha, tlsinvDistTheta, dDistTheta, qDistTheta, pDistX, qDistN, pDistN)
{

  epsilon = 1/1000000000 # Erreur acceptable
  jmax = qDistTheta(1 - epsilon, alpha)
  kmax <- qDistN(1 - epsilon)
  fS <- numeric(nfft)
  ii <- 1:(nfft - 1)

  for (j in seq(jmax))
  {
      # Pr[X | theta]
      fxj <- numeric(nfft)
      fxj[1] <- exp(-j * tlsinvDistTheta(pDistX(0), alpha))
      fxj[2:(nfft)] <- (exp(-j * tlsinvDistTheta(pDistX(ii), alpha)) -
                  exp(-j * tlsinvDistTheta(pDistX(ii - 1), alpha)))
      ffxj <- fft(fxj)

      # Pr[N | theta]
      fnj <- numeric(kmax + 1)
      fnj[1] <- exp(-j * tlsinvDistTheta(pDistN(0), alpha))
      fnj[2:(kmax+1)] <- (exp(-j * tlsinvDistTheta(pDistN(1:kmax), alpha)) -
                  exp(-j * tlsinvDistTheta(pDistN(1:kmax - 1), alpha)))

      fxk <- matrix(numeric((kmax + 1) * nfft), ncol = kmax + 1)
      fxk[1, 1] <- fnj[1]  # masse à 0

      for (k in 1:kmax)
      {
        fxk[, k + 1] <- Re(fft(ffxj^k, inverse = TRUE)) / nfft * fnj[k + 1]
      }
      fsj <- rowSums(fxk)

      fS <-  fS + dDistTheta(j, alpha) * fsj
  }

    fS

}


# f_{S} (ih)
calculer_fSfft_archi_simple <- function(nfft, alpha, tlsinvDistTheta, dDistTheta, qDistTheta, pDistX, qDistN, pDistN)
{

    epsilon = 1/1000000000 # Erreur acceptable
    jmax = qDistTheta(1 - epsilon, alpha)
    kmax <- qDistN(1 - epsilon)
    fS <- numeric(nfft)
    ii <- 1:(nfft - 1)

    for (j in seq(jmax))
    {
        # Pr[X | theta]
        fxj <- numeric(nfft)
        fxj[1] <- exp(-j * tlsinvDistTheta(pDistX(0), alpha))
        fxj[2:(nfft)] <- (exp(-j * tlsinvDistTheta(pDistX(ii), alpha)) -
                              exp(-j * tlsinvDistTheta(pDistX(ii - 1), alpha)))
        ffxj <- fft(fxj)

        # Pr[N | theta]
        fnj <- numeric(kmax + 1)
        fnj[1] <- exp(-j * tlsinvDistTheta(pDistN(0), alpha))
        fnj[2:(kmax+1)] <- (exp(-j * tlsinvDistTheta(pDistN(1:kmax), alpha)) -
                                exp(-j * tlsinvDistTheta(pDistN(1:kmax - 1), alpha)))

        fxk <- matrix(numeric((kmax + 1) * nfft), ncol = kmax + 1)
        fxk[1, 1] <- fnj[1]  # masse à 0
        fnj_nfft <- numeric(nfft)
        fnj_nfft[1:(kmax + 1)] <- fnj
        fgpn <- sapply(ffxj, function(a) sum(a^(0:(nfft-1)) * fnj_nfft))

        fsj <- Re(fft(fgpn, inverse = TRUE)) / nfft
        print(round(sum(fsj), 6) == 1)
        fS <-  fS + dDistTheta(j, alpha) * fsj
    }

    fS

}



# Vérification qu'on somme à 1
source("geom_shifted.R")
fs <- calculer_fS_archi_simple(nfft = 2^9, alpha = 0.2,
                                                         tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                                                         dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                                                         qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                                                         pDistX = function(x) ppois(x, 10),
                                                         qDistN = function(x) qpois(x, 2),
                                                         pDistN = function(x) ppois(x, 2))

sum(fs)


fsfft <- calculer_fSfft_archi_simple(nfft = 2^9, alpha = 0.2,
                               tlsinvDistTheta = function(x, alpha) log((1-alpha)/x + alpha),
                               dDistTheta = function(x, alpha) dgeo_shifted(x, alpha),
                               qDistTheta = function(x, alpha) qgeo_shifted(x, alpha),
                               pDistX = function(x) ppois(x, 10),
                               qDistN = function(x) qpois(x, 2),
                               pDistN = function(x) ppois(x, 2))

sum(fs)

