# Modèle indépendance
# calculer_fS_modele_classique <- function(nfft, fgp_N, h, pDistX, method = "lower")
# {
#     fmpX <- discr_v2(nfft, h, pDistX, method = method)
#     # if (length(fmpX) < nfft)
#     #     stop("Erreur : length(fmpX) < nfft")
#
#     fmpX <- fmpX[1:nfft]
#
#     fft_fgp_N_fX <- fgp_N(fft(fmpX))
#
#
#     Re(fft(fft_fgp_N_fX, inverse = TRUE)) / nfft
# }
#
# fgp_N_poisson <- function(t, lam)
# {
#     exp(lam * (t - 1))
# }
# nfft <- 2^15
# fS_classique <- calculer_fS_modele_classique(nfft,
#                                              function(x) fgp_N_poisson(x, lambdaPois),
#                                              h,
#                                              function(x) pgamma(x, alphaGa, betaGa))
#
# sum(fS_classique)
# lambdaPois * (alphaGa / betaGa) # validation
# sum((0:(nfft - 1)) * fS_classique)