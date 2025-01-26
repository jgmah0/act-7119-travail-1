### ACT-7119
## TP1

## Reproduction exemple 7 de Cossette et al. (2018)
# pour s'échauffer

qi <- 0.1 * 1:4
nfft <- 2^{5}
n <- 10

fx1 <- dbinom(0:(nfft-1), n, qi[1])
fx2 <- dbinom(0:(nfft-1), n, qi[2])
fx3 <- dbinom(0:(nfft-1), n, qi[3])
fx4 <- dbinom(0:(nfft-1), n, qi[4])

sum(fx1); sum(fx2); sum(fx3); sum(fx4)

fxi_theta <- function(theta, fxi, alpha)
{
    # copule Frank
    FXi <- cumsum(fxi)
    inverselst <- function(u) -log((1-exp(-alpha * u))/(1 -exp(-alpha)))
    FXI_theta <- exp(-theta * inverselst(FXi))

    head(c(FXI_theta, 0) - c(0, FXI_theta), -1)
}

alpha <- 1

ftheta_frank <- (1-exp(-alpha))^(1:50)/((1:50) * alpha)

fx1theta <- matrix(numeric(50*nfft), nrow=nfft)

for (i in 1:50)
{
    fx1theta[, i] <- fxi_theta(i, fx1, alpha)
}

fx2theta <- matrix(numeric(50*nfft), nrow=nfft)

for (i in 1:50)
{
    fx2theta[, i] <- fxi_theta(i, fx2, alpha)
}


fx3theta <- matrix(numeric(50*nfft), nrow=nfft)

for (i in 1:50)
{
    fx3theta[, i] <- fxi_theta(i, fx3, alpha)
}


fx4theta <- matrix(numeric(50*nfft), nrow=nfft)

for (i in 1:50)
{
    fx4theta[, i] <- fxi_theta(i, fx4, alpha)
}

colSums(fx1theta)
colSums(fx2theta)
colSums(fx3theta)
colSums(fx4theta)

ffx1 <- apply(fx1theta, 2, fft)
ffx2 <- apply(fx2theta, 2, fft)
ffx3 <- apply(fx3theta, 2, fft)
ffx4 <- apply(fx4theta, 2, fft)

ffs <- matrix(numeric(50*nfft), nrow=nfft)

for (i in 1:50)
{
    ffs[, i] <- ffx1[, i] * ffx2[, i] * ffx3[, i] * ffx4[, i]
}

fstheta <- apply(ffs, 2, function(a) Re(fft(a, inverse=TRUE))/nfft)


colSums(fstheta)


fs <- numeric(32)

for (i in 1:32)
{
    fs[i] <- sum(fstheta[i, ] * ftheta_frank)
}

# Tableau 2 de Cossette (2018)
fs[1]
fs[6]
fs[11]
