###
### Travail 1, ACT-7119
### Fonction pour calculer la valeur de f_S (ih)
### CRM indépendance
###
### Exemple :
###

# f_{S} (ih)
calculer_fS_indep <- function(nfft, dDistX, fgpN)
{
  fx = sapply(seq(nfft) - 1, dDistX)
  phix = fft(fx)
  phis = fgpN(phix)
  fs = Re(fft(phis, inverse = T))/nfft
  
  fs
}

# Vérification qu'on somme à 1
fs <- calculer_fS_indep(nfft = 2^20,
                        dDistX = function(x) dpois(x, 3),
                        fgpN = function(phi) exp(4 * (phi - 1)))

sum(fs)
