# R-code to generate data for F79MA Assignment 1. This code MUST be 
# placed at the start of your own R-script. You must edit
# the argument to the set.seed( ) function to fit your own
# registration number
#
# Please also ensure you have created all required variables
# (Quantity1, Quantity2, Quantity3, Quantity4, Quantity5, Quantity6)
# and include code for all of your plots
#

set.seed(7830)

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
# Task 3 (theory, implemented numerically) - joint distribution of
# (Y^A, Y^B), log-likelihood and score function.
#
# For one day, Y^A = X0 + X1 and Y^B = X0 + X2, where X0, X1, X2 are
# independent Poisson(lambda0), Poisson(lambda1), Poisson(lambda2).
# The pair is NOT independent because both counts contain X0.  Writing
# a and b for the observed counts, k for the (unknown) value of X0, and
# noting X0 <= min(a, b) since X1 = a - k >= 0 and X2 = b - k >= 0:
#
#   p(a, b) = exp(-(l0 + l1 + l2)) * S(a, b)                        (1)
#   S(a, b) = sum_{k=0}^{min(a,b)} l0^k l1^(a-k) l2^(b-k) /
#                                  ( k! (a-k)! (b-k)! )
#
# which is the bivariate Poisson distribution (Holgate, 1964).
# For an i.i.d. sample of n days:
#
#   l(l0)  = -n (l0 + l1 + l2) + sum_i log S(yA_i, yB_i)            (2)
#   U(l0)  = dl/dl0 = -n + sum_i S(yA_i - 1, yB_i - 1)/S(yA_i,yB_i)  (3)
#          = -n + (1/l0) * sum_i E[ X0_i | yA_i, yB_i ]
#
# (3) is used because U(l0) = 0 is exactly the EM fixed point
# l0 = (1/n) sum_i E[X0_i | data].  Since d2l/dl0^2 < 0 the likelihood
# is strictly concave and U is strictly decreasing, so the root of U is
# the unique MLE.
#####################################################################

# log S(a, b), vectorised over a and b, evaluated with a log-sum-exp so
# that the factorials and powers never overflow or underflow.
logS <- function(a, b, l0, l1, l2) {
  n <- max(length(a), length(b))
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  acc <- rep(-Inf, n)
  for (k in 0:max(pmin(a, b))) {
    ok <- k <= pmin(a, b)              # X0 = k is only possible if k <= min(a,b)
    if (!any(ok)) next
    lt <- rep(-Inf, n)
    lt[ok] <- k * log(l0) + (a[ok] - k) * log(l1) + (b[ok] - k) * log(l2) -
      lgamma(k + 1) - lgamma(a[ok] - k + 1) - lgamma(b[ok] - k + 1)
    m <- pmax(acc, lt)
    acc <- ifelse(is.finite(m), m + log(exp(acc - m) + exp(lt - m)), -Inf)
  }
  acc                                  # returns log S(a, b)
}

# log-likelihood (2) for a sample
loglik_l0 <- function(l0, yA, yB, l1, l2)
  -length(yA) * (l0 + l1 + l2) + sum(logS(yA, yB, l0, l1, l2))

# score function (3) for a sample
score_l0 <- function(l0, yA, yB, l1, l2) {
  pos <- yA >= 1 & yB >= 1             # terms with yA = 0 or yB = 0 contribute 0
  -length(yA) + sum(ifelse(pos, exp(logS(yA - 1, yB - 1, l0, l1, l2) -
                                      logS(yA, yB, l0, l1, l2)), 0))
}

#####################################################################
# Task 4 - plot the observed score function for one sample of n = 30
# and find the MLE of lambda0.  Store the MLE as Quantity1.
#####################################################################

n4 <- 30
sample4 <- simulate_YA_YB(n4)
YA4 <- sample4$YA
YB4 <- sample4$YB

cat("Task 4 - simulated sample of size 30\n")
cat("  yA:", YA4, "\n")
cat("  yB:", YB4, "\n")
cat("  mean(yA) =", mean(YA4), " (E[YA] =", lambda0 + lambda1, ")\n")
cat("  mean(yB) =", mean(YB4), " (E[YB] =", lambda0 + lambda2, ")\n\n")

# observed score function U(l0) evaluated over a grid of l0
lambda_grid <- seq(0.1, 6, by = 0.002)
score_grid <- vapply(lambda_grid, score_l0, 0,
                     yA = YA4, yB = YB4, l1 = lambda1, l2 = lambda2)

# MLE: the unique root of U(l0) = 0 (U is strictly decreasing)
Quantity1 <- uniroot(score_l0, interval = c(0.001, 100),
                     yA = YA4, yB = YB4, l1 = lambda1, l2 = lambda2,
                     tol = 1e-12)$root

# independent check: maximise the log-likelihood directly
mle_check <- optim(1, function(l) -loglik_l0(l, YA4, YB4, lambda1, lambda2),
                   method = "Brent", lower = 0.001, upper = 100)$par

# plot of the observed score function
png("task4_score.png", width = 1500, height = 950, res = 150)
plot(lambda_grid, score_grid, type = "l", lwd = 2,
     xlab = expression(lambda[0]),
     ylab = expression(U(lambda[0])),
     main = "Observed score function for the simulated sample (n = 30)")
abline(h = 0, col = "grey40", lty = 2)
abline(v = Quantity1, col = "red", lty = 3, lwd = 2)
points(Quantity1, 0, col = "red", pch = 19, cex = 1.3)
legend("topright", bty = "n", lty = c(1, 2, 3), lwd = c(2, 1, 2),
       col = c("black", "grey40", "red"),
       legend = c(expression(U(lambda[0])), "U = 0",
                  bquote("MLE =" ~ .(round(Quantity1, 4)))))
dev.off()

cat("  MLE of lambda0 (Quantity1)  =", Quantity1, "\n")
cat("  check by direct maximisation =", mle_check, "\n")
cat("  U(Quantity1)                 =",
    score_l0(Quantity1, YA4, YB4, lambda1, lambda2), "\n")
cat("  quantity T of Task 5/6       =",
    (mean(YA4) + mean(YB4) - lambda1 - lambda2) / 2, "\n")