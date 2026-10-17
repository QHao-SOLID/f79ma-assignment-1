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

P <- function(a, b, l0, l1 = lambda1, l2 = lambda2) {
  total <- 0
  for (k in 0:min(a, b)) {
    total <- total + l0^k * l1^(a - k) * l2^(b - k) /
      (factorial(k) * factorial(a - k) * factorial(b - k))
  }
  total
}

P_sample <- function(yA, yB, l0, l1 = lambda1, l2 = lambda2) {
  mapply(P, a = yA, b = yB, MoreArgs = list(l0 = l0, l1 = l1, l2 = l2))
}

loglik <- function(l0, yA, yB, l1 = lambda1, l2 = lambda2) {
  -length(yA) * (l0 + l1 + l2) + sum(log(P_sample(yA, yB, l0, l1, l2)))
}

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

T_estimator <- function(yA, yB, l1 = lambda1, l2 = lambda2) {
  (mean(yA) + mean(yB) - l1 - l2) / 2
}

# Task 4

n <- 30
sample30 <- simulate_YA_YB(n)
YA <- sample30$YA
YB <- sample30$YB

cat("Y^A:", YA, "\n")
cat("Y^B:", YB, "\n")
cat("mean(Y^A) =", mean(YA), "\n")
cat("mean(Y^B) =", mean(YB), "\n")

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

plot(
  l0_grid,
  loglik_grid,
  type = "l",
  lwd = 2,
  xlab = expression(lambda[0]),
  ylab = expression("log-likelihood" ~ italic(l)(lambda[0])),
  main = "Observed log-likelihood"
)
abline(
  v = Quantity1,
  col = "red",
  lty = 3,
  lwd = 2
)
points(
  Quantity1,
  loglik(Quantity1, YA, YB),
  col = "red",
  pch = 19,
  cex = 1.2
)
legend(
  "topright",
  bty = "n",
  lty = 3,
  lwd = 2,
  col = "red",
  legend = bquote("MLE =" ~ .(round(Quantity1, 4)))
)

plot(
  l0_grid,
  score_grid,
  type = "l",
  lwd = 2,
  xlab = expression(lambda[0]), ylab = expression(U(lambda[0])),
  main = "Observed score function"
)
abline(
  h = 0,
  col = "grey40",
  lty = 2
)
abline(
  v = Quantity1,
  col = "red",
  lty = 3,
  lwd = 2
)
points(
  Quantity1,
  0,
  col = "red",
  pch = 19,
  cex = 1.2
)
legend(
  "topright",
  bty = "n",
  lty = c(1, 2, 3),
  lwd = c(2, 1, 2),
  col = c("black", "grey40", "red"),
  legend = c(
    expression(U(lambda[0])),
    "U = 0",
    bquote("MLE =" ~ .(round(Quantity1, 4)))
  )
)
dev.off()

cat("lambda0 =", Quantity1, "\n")

# Task 6
Quantity2 <- (mean(YA) + mean(YB) - lambda1 - lambda2) / 2
cat("Quantity2 =", Quantity2, "\n")
cat("Var(T) = (4*l0 + l1 + l2) / (4n) =", (4 * lambda0 + lambda1 + lambda2) / (4 * n), "\n")
cat("SD(T) = sqrt(Var(T)) =", sqrt((4 * lambda0 + lambda1 + lambda2) / (4 * n)), "\n")

# Task 7

n_samples <- 10000
days7 <- simulate_YA_YB(n * n_samples)

YA7 <- matrix(days7$YA, nrow = n)
YB7 <- matrix(days7$YB, nrow = n)

T7 <- (colMeans(YA7) + colMeans(YB7) - lambda1 - lambda2) / 2

cat("mean(T) =", mean(T7), "\n")
cat("median(T) =", median(T7), "\n")
cat("sd(T) =", sd(T7), " (theory", sqrt((4 * lambda0 + lambda1 + lambda2) / (4 * n)), ")\n")
cat("var(T) =", var(T7), " (theory", (4 * lambda0 + lambda1 + lambda2) / (4 * n), ")\n")
cat("range(T) =", range(T7), "\n")
cat("quantiles =", quantile(T7, c(0.025, 0.25, 0.5, 0.75, 0.975)), "\n")

breaks7 <- seq(0.4, 3.6, by = 1 / 15)

png("task7_hist.png", width = 1500, height = 950, res = 145)
hist(
  T7,
  breaks = breaks7,
  probability = TRUE,
  col = "grey85",
  border = "white",
  xlab = expression(T),
  ylab = "Density",
  main = paste0(
    "Histogram of T from ", n_samples, " independent samples of size n = ", n
  )
)
curve(
  dnorm(
    x,
    mean = lambda0,
    sd = sqrt((4 * lambda0 + lambda1 + lambda2) / (4 * n))
  ),
  add = TRUE,
  col = "red",
  lwd = 2
)
abline(
  v = lambda0,
  col = "blue",
  lty = 2,
  lwd = 2
)
abline(
  v = mean(T7),
  col = "darkgreen",
  lty = 3,
  lwd = 2
)
legend(
  "topright",
  bty = "n",
  lwd = 2,
  lty = c(1, 2, 3),
  col = c("red", "blue", "darkgreen"),
  legend = c(
    paste0("Normal(2, ", round((4 * lambda0 + lambda1 + lambda2) / (4 * n), 3), ")"),
    "true mean of T = 2",
    paste0("sample mean = ", round(mean(T7), 4))
  )
)
dev.off()

# Task 8

centred_T7 <- T7 - mean(T7)
q8 <- quantile(centred_T7, c(0.025, 0.975))

Quantity3 <- Quantity2 + q8[1]
Quantity4 <- Quantity2 + q8[2]

cat("Quantity3 (lower limit) =", Quantity3, "\n")
cat("Quantity4 (upper limit) =", Quantity4, "\n")
cat("interval width =", Quantity4 - Quantity3, "\n")

# Task 9

B9 <- 200000
days9 <- simulate_YA_YB(n * B9)
YA9 <- matrix(days9$YA, nrow = n)
YB9 <- matrix(days9$YB, nrow = n)

U9 <- vapply(
  seq_len(B9),
  function(b) score(lambda0, YA9[, b], YB9[, b]),
  0
)

Quantity5 <- mean(U9^2)
Quantity6 <- 1 / Quantity5

cat("simulated mean of U(2) =", mean(U9), " (should be 0)\n")
cat("Quantity5 = I(2) =", Quantity5, "\n")
cat("Quantity6 = 1 / I(2) =", Quantity6, "\n")

# Task 10

i2 <- Quantity5 / n
ns <- 10:50

varT_n <- (4 * lambda0 + lambda1 + lambda2) / (4 * ns)
crlb_n <- 1 / (ns * i2)

B10 <- 5000
varT_sim <- sapply(ns, function(nn) {
  d <- simulate_YA_YB(nn * B10)
  tt <- (colMeans(matrix(d$YA, nrow = nn)) +
           colMeans(matrix(d$YB, nrow = nn)) - lambda1 - lambda2) / 2
  var(tt)
})

cat("\nTask 10:  Var(T) and the CRLB for n = 10 to 50\n")
cat("ratio Var(T) / CRLB is constant at", mean(varT_n / crlb_n), "for every n from 10 to 50\n")
cat("n = 10 : Var(T) =", varT_n[1], " CRLB =", crlb_n[1], " simulated =", varT_sim[1], "\n")
cat("n = 30 : Var(T) =", varT_n[21], " CRLB =", crlb_n[21], " simulated =", varT_sim[21], "\n")
cat("n = 50 : Var(T) =", varT_n[41], " CRLB =", crlb_n[41], " simulated =", varT_sim[41], "\n")

png("task10_crlb.png", width = 1700, height = 850, res = 145)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

plot(
  ns,
  varT_n,
  type = "l",
  lwd = 2,
  ylim = c(0, max(varT_n)),
  xlab = "n", ylab = "variance",
  main = "Variance of T and the CRLB"
)
lines(
  ns,
  crlb_n,
  col = "red",
  lwd = 2,
  lty = 2
)
points(
  ns,
  varT_sim,
  pch = 20,
  col = "grey35"
)
legend(
  "topright",
  bty = "n",
  lwd = c(2, 2, 1), lty = c(1, 2, NA),
  pch = c(NA, NA, 20), col = c("black", "red", "grey35"),
  legend = c(
    paste0("Var(T) = ", 4 * lambda0 + lambda1 + lambda2, " / (4n)"),
    "CRLB = 1 / (n i(2))", "simulated Var(T)"
  )
)

plot(
  ns,
  varT_n / crlb_n,
  type = "l",
  lwd = 2,
  ylim = c(0.9, 1.2),
  xlab = "n",
  ylab = "Var(T) / CRLB",
  main = "Efficiency of T does not improve with n"
)
abline(
  h = 1,
  col = "red",
  lty = 2,
  lwd = 2
)
dev.off()
