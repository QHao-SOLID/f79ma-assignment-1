# R-code to generate data for F79MA Assignment 1. This code MUST be 
# placed at the start of your own R-script. You must edit
# the argument to the set.seed( ) function to fit your own
# registration number
#
# Please also ensure you have created all required variables
# (Quantity1, Quantity2, Quantity3, Quantity4, Quantity5, Quantity6)
# and include code for all of your plots
#

RNGkind(sample.kind = "Rejection")
lambda0 <- 2
lambda1 <- 3
lambda2 <- 4
#
simulate_YA_YB <- function(n){ # n is the sample size
  set.seed(7830) ## Change this number
  X0 <- rpois(n,lambda0)
  X1 <- rpois(n,lambda1)
  X2 <- rpois(n,lambda2)
  YA <- X0 + X1 # sample of tremor numbers in zone A
  YB <- X0 + X2 # sample of tremor numbers in zone B
  return(list(YA = YA,YB = YB))
}
#####################################################################
# Please insert your R code after this line
#####################################################################

#####################################################################
# Task 3 (theory, carried out in R)
#
# For one day, Y^A = X0 + X1 and Y^B = X0 + X2, where X0, X1, X2 are
# independent Poisson variables with rates lambda0, lambda1, lambda2.
# Both zone counts contain X0, so Y^A and Y^B are not independent.
#
# Condition on X0 = k.  Then X1 = a - k and X2 = b - k, so k can only
# take the values 0, 1, ..., min(a, b); X0 cannot exceed either observed
# count.  Multiplying the three Poisson probabilities and summing over k
# gives the joint distribution of (Y^A, Y^B):
#
#   P(a, b) = sum_{k=0}^{min(a,b)} l0^k l1^(a-k) l2^(b-k)
#                                 / ( k! (a-k)! (b-k)! )                 (1)
#
# P(a, b) is the unnormalised weight of the pair (a, b): the joint
# probability is exp(-(l0 + l1 + l2)) * P(a, b).  That exponential factor
# is identical for every (a, b), which is why it appears once, separately,
# in the log-likelihood below.
#
# For an independent sample of n days:
#
#   l(l0) = -n (l0 + l1 + l2) + sum_i log P(Y_i^A, Y_i^B)                (2)
#
#   U(l0) = dl/dl0 = -n + sum_i P(Y_i^A - 1, Y_i^B - 1)
#                              / P(Y_i^A, Y_i^B)                         (3)
#
# Equation (3) uses  d/dl0 P(a, b) = P(a - 1, b - 1),  which holds because
#   d/dl0 [ l0^k / k! ] = k l0^(k-1) / k! = l0^(k-1) / (k-1)! .
# The k = 0 term of P(a, b) does not survive the derivative (it has no l0
# in it), and P(a, b) = 0 whenever a < 0 or b < 0, so days with an
# observed count of zero contribute nothing to the sum in (3).
#
# U is strictly decreasing in l0 (the log-likelihood is strictly concave),
# so U(l0) = 0 has exactly one root, and that root is the MLE.
#####################################################################

# (1) The joint distribution P(a, b) for one day, a and b being the two
#     observed counts.  Written to follow the report line by line.
P <- function(a, b, l0, l1 = lambda1, l2 = lambda2) {
  total <- 0
  for (k in 0:min(a, b)) {
    total <- total + l0^k * l1^(a - k) * l2^(b - k) /
      (factorial(k) * factorial(a - k) * factorial(b - k))
  }
  total
}

# (1) applied to every day of the sample at once
P_sample <- function(yA, yB, l0, l1 = lambda1, l2 = lambda2) {
  mapply(P, a = yA, b = yB, MoreArgs = list(l0 = l0, l1 = l1, l2 = l2))
}

# (2) log-likelihood of the sample
loglik <- function(l0, yA, yB, l1 = lambda1, l2 = lambda2) {
  -length(yA) * (l0 + l1 + l2) + sum(log(P_sample(yA, yB, l0, l1, l2)))
}

# (3) score function of the sample
score <- function(l0, yA, yB, l1 = lambda1, l2 = lambda2) {
  ratio <- numeric(length(yA))
  for (i in seq_along(yA)) {
    if (yA[i] >= 1 && yB[i] >= 1) {
      ratio[i] <- P(yA[i] - 1, yB[i] - 1, l0, l1, l2) /
        P(yA[i], yB[i], l0, l1, l2)
    }
  }
  -length(yA) + sum(ratio)
}

#####################################################################
# Task 4 - a single sample of size n = 30, and the MLE of lambda0
#####################################################################

n <- 30
sample30 <- simulate_YA_YB(n)
YA <- sample30$YA
YB <- sample30$YB

cat("Task 4:  single simulated sample of size n =", n, "\n")
cat("  Y^A :", YA, "\n")
cat("  Y^B :", YB, "\n")
cat("  mean(Y^A) =", mean(YA), " (E[Y^A] = lambda0 + lambda1 =",
    lambda0 + lambda1, ")\n")
cat("  mean(Y^B) =", mean(YB), " (E[Y^B] = lambda0 + lambda2 =",
    lambda0 + lambda2, ")\n\n")

# values of lambda0 at which the score and the log-likelihood are evaluated
l0_grid <- seq(0.1, 6, by = 0.005)
score_grid <- sapply(l0_grid, score, yA = YA, yB = YB)
loglik_grid <- sapply(l0_grid, loglik, yA = YA, yB = YB)

# the MLE is the unique root of the observed score function
Quantity1 <- uniroot(score, interval = c(0.001, 100),
                     yA = YA, yB = YB, tol = 1e-12)$root

# plot: log-likelihood on the left, score function on the right
png("task4_plots.png", width = 1700, height = 850, res = 145)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

plot(l0_grid, loglik_grid, type = "l", lwd = 2,
     xlab = expression(lambda[0]),
     ylab = expression("log-likelihood" ~ italic(l)(lambda[0])),
     main = "Observed log-likelihood")
abline(v = Quantity1, col = "red", lty = 3, lwd = 2)
points(Quantity1, loglik(Quantity1, YA, YB), col = "red", pch = 19, cex = 1.2)
legend("topright", bty = "n", lty = 3, lwd = 2, col = "red",
       legend = bquote("MLE =" ~ .(round(Quantity1, 4))))

plot(l0_grid, score_grid, type = "l", lwd = 2,
     xlab = expression(lambda[0]), ylab = expression(U(lambda[0])),
     main = "Observed score function")
abline(h = 0, col = "grey40", lty = 2)
abline(v = Quantity1, col = "red", lty = 3, lwd = 2)
points(Quantity1, 0, col = "red", pch = 19, cex = 1.2)
legend("topright", bty = "n", lty = c(1, 2, 3), lwd = c(2, 1, 2),
       col = c("black", "grey40", "red"),
       legend = c(expression(U(lambda[0])), "U = 0",
                  bquote("MLE =" ~ .(round(Quantity1, 4)))))
dev.off()

cat("  Quantity1 = MLE of lambda0 =", Quantity1, "\n")
cat("  check: U(Quantity1)        =", score(Quantity1, YA, YB), "\n")
cat("  check: l(Quantity1)        =", loglik(Quantity1, YA, YB), "\n")
cat("  check: l at the true value =", loglik(lambda0, YA, YB),
    " (must be no larger)\n")

#####################################################################
# Task 5 (theory, used again in Tasks 7 to 10)
#
#   T = ( mean(Y^A) + mean(Y^B) - lambda1 - lambda2 ) / 2
#
# Bias: the days are identically distributed, so
#   E[T] = ( E[Y^A] + E[Y^B] - l1 - l2 ) / 2
#        = ( (l0 + l1) + (l0 + l2) - l1 - l2 ) / 2
#        = l0 ,
# so T is an unbiased estimator of lambda0.
#
# Variance: the n days are independent, so the variance of each mean is
# the variance of one day divided by n, and
#   Var(T) = (1/4) [ Var(Ybar^A) + Var(Ybar^B) + 2 Cov(Ybar^A, Ybar^B) ]
#          = (1/4) [ (l0 + l1)/n + (l0 + l2)/n + 2 l0/n ]
#          = (4 l0 + l1 + l2) / (4 n) .
# With l0 = 2, l1 = 3, l2 = 4 and n = 30 this is 15/120 = 0.125.
#####################################################################

T_estimator <- function(yA, yB, l1 = lambda1, l2 = lambda2) {
  (mean(yA) + mean(yB) - l1 - l2) / 2
}

#####################################################################
# Task 6 - T for the same sample of size n = 30, stored as Quantity2
#####################################################################

Quantity2 <- T_estimator(YA, YB)

cat("\nTask 6:  alternative estimator T, same sample of size", n, "\n")
cat("  mean(Y^A) =", mean(YA), " and mean(Y^B) =", mean(YB), "\n")
cat("  Quantity2 = T   =", Quantity2, "\n")
cat("  Quantity1 = MLE =", Quantity1, "\n")
cat("  T - MLE         =", Quantity2 - Quantity1, "\n")
cat("  theoretical Var(T) = (4*l0 + l1 + l2) / (4n) =",
    (4 * lambda0 + lambda1 + lambda2) / (4 * n),
    ", so sd(T) =",
    sqrt((4 * lambda0 + lambda1 + lambda2) / (4 * n)), "\n")
cat("  both estimates fall below the true lambda0 =", lambda0,
    "because both sample means are below their expectations\n")

#####################################################################
# Task 7 - histogram of T from repeated sampling, sample size n = 30
#
# simulate_YA_YB() calls set.seed() internally, so calling it b times in a
# loop would return the SAME sample b times, not b independent samples.
# The function does return n i.i.d. days, so the way to obtain repeated
# independent samples is to simulate n * b days in a single call and cut
# the result into b consecutive blocks of n days.  The days are i.i.d.,
# so those blocks are independent samples of size n.
#####################################################################

n_samples <- 10000                          # number of independent samples
days7 <- simulate_YA_YB(n * n_samples)      # n * b i.i.d. days, one call

YA7 <- matrix(days7$YA, nrow = n)           # column b holds sample b
YB7 <- matrix(days7$YB, nrow = n)

T7 <- (colMeans(YA7) + colMeans(YB7) - lambda1 - lambda2) / 2

# theoretical values of T from Task 5
varT <- (4 * lambda0 + lambda1 + lambda2) / (4 * n)
sdT <- sqrt(varT)
# skewness of T: it is a linear combination of three independent means,
# so its third central moment is  l0/n^2 + (l1 + l2)/(8 n^2)
skewT <- (lambda0 + (lambda1 + lambda2) / 8) / n^2 / sdT^3

cat("\nTask 7:", n_samples, "independent samples of size n =", n, "\n")
cat("  mean(T)     =", mean(T7), " (theory lambda0 =", lambda0, ")\n")
cat("  median(T)   =", median(T7), "\n")
cat("  sd(T)       =", sd(T7), " (theory", sdT, ")\n")
cat("  var(T)      =", var(T7), " (theory", varT, ")\n")
cat("  skewness(T) =", mean((T7 - mean(T7))^3) / sd(T7)^3,
    " (theory", skewT, ")\n")
cat("  range(T)    =", range(T7), "\n")
cat("  quantiles   =",
    quantile(T7, c(0.025, 0.25, 0.5, 0.75, 0.975)), "\n")
cat("  P(|T - 2| < 1 sd) =", mean(abs(T7 - lambda0) < sdT),
    " (normal gives 0.683)\n")
cat("  P(|T - 2| < 2 sd) =", mean(abs(T7 - lambda0) < 2 * sdT),
    " (normal gives 0.954)\n")

# T is lattice valued with spacing 1/60 (T = (sum(YA) + sum(YB) - 210)/60),
# so bin edges are placed on that lattice with a width of 4 lattice steps
# to avoid a comb artefact in the histogram.
breaks7 <- seq(0.4, 3.6, by = 1 / 15)

png("task7_hist.png", width = 1500, height = 950, res = 145)
hist(T7, breaks = breaks7, probability = TRUE,
     col = "grey85", border = "white",
     xlab = expression(T), ylab = "Density",
     main = paste0("Histogram of T from ", n_samples,
                   " independent samples of size n = ", n))
curve(dnorm(x, mean = lambda0, sd = sdT), add = TRUE, col = "red", lwd = 2)
abline(v = lambda0, col = "blue", lty = 2, lwd = 2)
abline(v = mean(T7), col = "darkgreen", lty = 3, lwd = 2)
legend("topright", bty = "n", lwd = 2, lty = c(1, 2, 3),
       col = c("red", "blue", "darkgreen"),
       legend = c(paste0("Normal(2, ", round(varT, 3), ")"),
                  "true mean of T = 2",
                  paste0("sample mean = ", round(mean(T7), 4))))
dev.off()

#####################################################################
# Task 8 - approximate 95% confidence interval for T
#
# Task 7 gives the sampling distribution of T, but it is centred on the
# lambda0 = 2 that was used to simulate it.  We want an interval centred
# on our own estimate Quantity2, so we keep only the SHAPE of that
# distribution: subtract its centre and add Quantity2 instead.
#####################################################################

centred_T7 <- T7 - mean(T7)                 # shape of the distribution of T
q8 <- quantile(centred_T7, c(0.025, 0.975))

Quantity3 <- Quantity2 + q8[1]              # lower 95% limit
Quantity4 <- Quantity2 + q8[2]              # upper 95% limit

cat("\nTask 8:  approximate 95% confidence interval for lambda0\n")
cat("  Quantity2 (the centre)  =", Quantity2, "\n")
cat("  2.5% and 97.5% points of T - mean(T) from Task 7:",
    q8[1], "and", q8[2], "\n")
cat("  Quantity3 (lower limit) =", Quantity3, "\n")
cat("  Quantity4 (upper limit) =", Quantity4, "\n")
cat("  interval width          =", Quantity4 - Quantity3, "\n")
cat("  covers the true lambda0 = 2 ?",
    Quantity3 <= lambda0 && lambda0 <= Quantity4, "\n")
cat("  covers the MLE Quantity1     ?",
    Quantity3 <= Quantity1 && Quantity1 <= Quantity4, "\n")
cat("  (the symmetric normal interval +-1.96*sd(T) would be",
    Quantity2 - 1.96 * sd(T7), "to", Quantity2 + 1.96 * sd(T7), ")\n")

#####################################################################
# Task 9 - Fisher information and the Cramer-Rao lower bound
#
# For one day the score is  U1(l0) = -1 + P(a-1,b-1)/P(a,b),  and when l0
# is the true value E[U1(l0)] = 0, so the information contributed by one
# observation is
#     i(l0) = Var[U1(l0)] = E[U1(l0)^2].
# For n independent days the score is U(l0) = sum_i U1(l0), which is the
# sample score function of Task 4, and
#     I(l0) = n * i(l0) = Var[U(l0)] = E[U(l0)^2].
# Estimation: simulate many samples, evaluate U(2) on each, average U(2)^2.
#####################################################################

B9 <- 200000
days9 <- simulate_YA_YB(n * B9)
YA9 <- matrix(days9$YA, nrow = n)
YB9 <- matrix(days9$YB, nrow = n)

U9 <- vapply(seq_len(B9), function(b) score(lambda0, YA9[, b], YB9[, b]), 0)

Quantity5 <- mean(U9^2)          # estimated Fisher information I(2), sample size n
Quantity6 <- 1 / Quantity5       # Cramer-Rao lower bound for an unbiased estimator

cat("\nTask 9:  Fisher information and Cramer-Rao lower bound, l0 = 2, n =",
    n, "\n")
cat("  simulated mean of U(2) =", mean(U9), " (should be 0)\n")
cat("  Quantity5 = I(2)       =", Quantity5, "\n")
cat("  Quantity6 = 1 / I(2)   =", Quantity6, "\n")
cat("  information per day i(2) = I(2)/n =", Quantity5 / n, "\n")
cat("  Monte Carlo standard error of Quantity5 =",
    sd(U9^2) / sqrt(B9), "\n")

# independent check: average U1(l0)^2 over the EXACT joint distribution of
# one day instead of over simulated days
grid9 <- expand.grid(a = 0:40, b = 0:40)
grid9$p <- exp(-(lambda0 + lambda1 + lambda2)) *
  mapply(P, grid9$a, grid9$b, MoreArgs = list(l0 = lambda0))
pos9 <- grid9$a >= 1 & grid9$b >= 1
grid9$u1 <- -1
grid9$u1[pos9] <- mapply(P, grid9$a[pos9] - 1, grid9$b[pos9] - 1,
                         MoreArgs = list(l0 = lambda0)) /
  mapply(P, grid9$a[pos9], grid9$b[pos9], MoreArgs = list(l0 = lambda0)) - 1
i9_exact <- sum(grid9$p * grid9$u1^2)

cat("  exact i(2) by summing over the joint distribution =", i9_exact, "\n")
cat("  exact I(2) = n * i(2)                             =", n * i9_exact, "\n")
cat("  exact mean of U1(2), which must be 0              =",
    sum(grid9$p * grid9$u1), "\n")
cat("  (sum of that joint distribution =", sum(grid9$p), ")\n")

#####################################################################
# Task 10 - variance of T against the CRLB as n runs from 10 to 50
#
# From Task 5:   Var(T) = (4 l0 + l1 + l2) / (4 n)
# From Task 9:   I_n(l0) = n * i(l0), so CRLB(n) = 1 / (n * i(l0)),
#                with i(2) = Quantity5 / 30 the information per day.
#####################################################################

i2 <- Quantity5 / n
ns <- 10:50

varT_n <- (4 * lambda0 + lambda1 + lambda2) / (4 * ns)   # theory
crlb_n <- 1 / (ns * i2)                                  # CRLB

B10 <- 5000                                              # simulation check
varT_sim <- sapply(ns, function(nn) {
  d <- simulate_YA_YB(nn * B10)
  tt <- (colMeans(matrix(d$YA, nrow = nn)) +
           colMeans(matrix(d$YB, nrow = nn)) - lambda1 - lambda2) / 2
  var(tt)
})

cat("\nTask 10:  Var(T) and the CRLB for n = 10 to 50\n")
cat("  ratio Var(T) / CRLB is constant at", mean(varT_n / crlb_n),
    "for every n from 10 to 50\n")
cat("  i.e. T is about", round(100 * (mean(varT_n / crlb_n) - 1), 1),
    "% less efficient than the bound, at every n\n")
cat("  n = 10 : Var(T) =", varT_n[1], " CRLB =", crlb_n[1],
    " simulated =", varT_sim[1], "\n")
cat("  n = 30 : Var(T) =", varT_n[21], " CRLB =", crlb_n[21],
    " simulated =", varT_sim[21], "\n")
cat("  n = 50 : Var(T) =", varT_n[41], " CRLB =", crlb_n[41],
    " simulated =", varT_sim[41], "\n")

png("task10_crlb.png", width = 1700, height = 850, res = 145)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

plot(ns, varT_n, type = "l", lwd = 2, ylim = c(0, max(varT_n)),
     xlab = "n", ylab = "variance",
     main = "Variance of T and the CRLB")
lines(ns, crlb_n, col = "red", lwd = 2, lty = 2)
points(ns, varT_sim, pch = 20, col = "grey35")
legend("topright", bty = "n", lwd = c(2, 2, 1), lty = c(1, 2, NA),
       pch = c(NA, NA, 20), col = c("black", "red", "grey35"),
       legend = c(paste0("Var(T) = ", 4 * lambda0 + lambda1 + lambda2, " / (4n)"),
                  "CRLB = 1 / (n i(2))", "simulated Var(T)"))

plot(ns, varT_n / crlb_n, type = "l", lwd = 2, ylim = c(0.9, 1.2),
     xlab = "n", ylab = "Var(T) / CRLB",
     main = "Efficiency of T does not improve with n")
abline(h = 1, col = "red", lty = 2, lwd = 2)
dev.off()
