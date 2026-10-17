#import "@preview/showybox:2.0.4": showybox
#set page(margin: 2cm)
#set par(
    justify: true
)

#set image(width: 85%)

// Variables
#let Mean_YA = $4.5$
#let Mean_YB = $5.8$
#let Variance_T = $0.125$
#let SD_T = $0.3535534$

#let sim_Mean_T = $1.999905$
#let sim_Median_T = $2$
#let sim_Variance_T = $0.1229286$
#let sim_SD_T = $0.3506117$
#let sim_Q1 = $1.3116667$
#let sim_Q4 = $2.7$
#let i2 = $0.2787811$
#let Var_CRLB_Ratio = $1.045429125$

#let Quantity3_Quantity4_width = $1.383333$

#let Quantity1 = $1.6851560$
#let Quantity2 = $1.6500000$
#let Quantity3 = $0.9667617$
#let Quantity4 = $2.3500950$
#let Quantity5 = $8.3634330$
#let Quantity6 = $0.1195681$

= F79MA - Assessed Project 1 - 2026/27
- Student name: Qing Hao, Too
- HW-ID: H00467830
- Campus: Heriot-Watt University Malaysia

== Introduction

Tectonic plates meet along fault lines, and for earthquake hazard assessment it
matters whether tremor activity on one side of a fault is related to activity on
the other. This report studies that question for two adjacent zones, $A$ and
$B$, using a simple Poisson model for the number of tremors recorded in a single
day.

On any given day the two counts are modelled as

$ Y^A = X_0 + X_1 quad "and" quad Y^B = X_0 + X_2, $

where $X_0$, $X_1$ and $X_2$ are independent Poisson variables with rates
$lambda_0$, $lambda_1$ and $lambda_2$. The term $X_0$ represents tremors caused
by a mechanism common to both zones, so it is counted twice, once in $A$ and
once in $B$; the terms $X_1$ and $X_2$ represent activity specific to a single
zone. This shared component is the essential feature of the model. It makes
$Y^A$ and $Y^B$ positively correlated, with $CC"ov"(Y^A, Y^B) = lambda_0$, so the
two counts cannot be treated as independent observations even though different
days are independent of one another.

The parameter of interest is the shared rate $lambda_0$, and the problem is
deliberately asymmetric: $lambda_1$ and $lambda_2$ are known, so all of the
statistical difficulty is concentrated in one unknown.

Two estimators are considered. The first is likelihood based. Because the two
counts share $X_0$, the joint distribution of $(Y^A, Y^B)$ must be derived
before a likelihood can be written down; it is the bivariate Poisson
distribution (Holgate, 1964; Kocherlakota and Kocherlakota, 1992), a standard
tool for modelling correlated counts in pairs (Karlis and Ntzoufras, 2003). The
score function is then solved numerically to give the maximum likelihood
estimate. The second estimator is a simple moment-based alternative built from
the two sample means. Poisson models of this kind are the usual starting point
in seismology (Ogata, 1988), so the exercise is a small version of a real
modelling problem.

The report proceeds in four stages. Tasks 1--3 derive the joint distribution,
likelihood, log-likelihood and score function. Tasks 4--6 carry out the
likelihood analysis in R for one simulated sample of $n = 30$ days and compare
the maximum likelihood estimate with the moment-based alternative. Tasks 7--8
characterise the sampling distribution of that alternative by simulation and
build an approximate 95% confidence interval from it. Tasks 9--10 estimate the
Fisher information, obtain the Cramér--Rao lower bound, and compare that bound
with the variance of the alternative estimator as the sample size grows from 10
to 50.

All numerical results are obtained for the rate vector
$(lambda_0, lambda_1, lambda_2) = (2, 3, 4)$ with the random seed prescribed by
the assignment, so that every figure and quantity reported below can be
reproduced from the accompanying R script.

== Task 1
#showybox(
  [Explain carefully what assumptions have been made in modelling the random variables $Y^A$ and $Y^B$ for a day.]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    The model rests on few assumptions.

    + *Independence within a day.* $X_0$, $X_1$ and $X_2$ do not influence one
      another. But $(Y^A, Y^B)$ is *not* independent, because both counts contain $X_0$.

    + *A common distribution on every day.* The three rates are the same on all
      days, so $(Y^A, Y^B)$ is identically distributed; there is no seasonal or
      long-term trend.

    + *Independence between days.* The tremor counts on one day do not depend on
      those on any other day.
  ]
)
== Task 2
#showybox(
    [
      State and explain the distributions of $Y^A$ and $Y^B$, and calculate the covariance of $Y^A$ and $Y^B$.
    ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    We use the idea of Probability Generating Functions $G_X (s)$ to check for distribution of $Y^A$ and $Y^B$. \ 
    
    1. It is known that for $X ~ "Poi"(lambda)$, $G_X (s) = e^(-lambda(1-s))$, and \
    2. For $Y^A = X_0 + X_1 => G_(Y^A)(s) =^"def." EE[s^(X_0 + X_1)] =^"i.i.d" EE[s^(X_0)]EE[s^(X_1)]$. \
    
    Via 1. and 2. It is easy to show that $G_(Y^A) (s) = e^(-lambda_0 (1-s)) dot e^(-lambda_1 (1-s)) = e^((lambda_0 + lambda_1)(1-s))$ following Poisson distribution with parameter $(lambda_0 + lambda_1)$. Similar idea can be shown for $Y^B$ for it follows Poisson distribution with parameter $(lambda_0 + lambda_2)$.
  ],
  [
    Covariance is defined as $CC"ov"(X,Y) = EE[X Y] - EE[X]EE[Y]$. \

    $
      CC"ov"(Y^A, Y^B) &= EE[Y^A Y^B] - EE[Y^A]EE[Y^B] \
      &= EE[(X_0 + X_1)(X_0 + X_2)] - (lambda_0 + lambda_1)(lambda_0 + lambda_2) \
      &= [EE[X_0^2] + EE[X_0](EE[X_1] + EE[X_2]) + EE[X_1]EE[X_2]] - (lambda_0 + lambda_1)(lambda_0 + lambda_2) \
      &= (lambda_0 + lambda_0^2) + lambda_0 (lambda_1 + lambda_2) + lambda_1 lambda_2 - (lambda_0 + lambda_1)(lambda_0 + lambda_2) \
      &= lambda_0
    $

    From above result, the covariance between $Y^A$ and $Y^B$ is $lambda_0$.
  ]
)
== Task 3
#showybox(
  [
    Assume that $lambda_1$ and $lambda_2$ are known. Derive the likelihood $L(lambda_0; (Y^A, Y^B))$, log-likelihood $cal(l)(lambda_0; (Y^A, Y^B))$ and score function $U(lambda_0; (Y^A, Y^B)) = (delta cal(l))/(delta lambda_0)$ for a random sample $(Y^A, Y^B)$ over $n$ days where \

    $ bold(Y)^A = (Y_1^A, Y_2^A, ..., Y_n^A) $
    and
    $ bold(Y)^B = (Y_1^B, Y_2^B, ..., Y_n^B) $

    *Hint:* You should first determine the joint distribution of $(Y_i^A, Y_i^B)$ for each $i$ as the number of tremors across the two zones are not independent.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    We shall first determine the joint distribution of $(Y_i^A, Y_i^B)$ via conditioning $X_0$ to render $Y_A$ and $Y_B$ independent, that is: \

    $
      PP(Y^A = a, Y^B = b) &= sum_(k=0)^(min(a,b)) PP(X_0 = k)PP(X_1 = a - k, X_2 = b - k | X_0 = k) \
      &=^"ind." sum_(k=0)^(min(a,b)) PP(X_0 = k)PP(X_1 = a - k) PP(X_2 = b-k) \
      &=^"p.m.f." sum_(k=0)^(min(a,b)) (e^(-lambda_0) lambda_0^(k))/(k!) times (e^(-lambda_1) lambda_1^((a-k)))/((a-k)!) times (e^(-lambda_2) lambda_2^((b-k)))/((b-k)!) \
      &= e^(-(lambda_0 + lambda_1 + lambda_2)) sum_(k=0)^(min(a,b)) (lambda_0^(k) lambda_1^(a-k) lambda_2^(b-k))/(k!(a-k)!(b-k)!)
    $

    The upper limit is $min(a,b)$ because $k <= a$ (so that $X_1 = a - k$ is a
    possible Poisson count) and $k <= b$ (so that $X_2 = b - k$ is as well) must
    both hold, and $k >= 0$ is required of any Poisson variable. Since the
    exponential factor does not depend on $k$, it can be taken outside the sum, as
    in the last line.
  ],
  [
    Suppose now $lambda_1$ and $lambda_2$ known, let us find the likelihood $L(lambda_0; (Y^A, Y^B))$ and log-likelihood function $cal(l)(lambda_0; (Y^A, Y^B))$, that is: \

    $
      L(lambda_0; (Y^A, Y^B)) &= product_(i=1)^n e^(-(lambda_0 + lambda_1 + lambda_2)) sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!) \
      &= e^(-n(lambda_0 + lambda_1 + lambda_2)) product_(i=1)^n [sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!)]
    $
    and so,
    $
      cal(l)(lambda_0; (Y^A, Y^B)) = -n(lambda_0 + lambda_1 + lambda_2) + sum_(i=1)^n log[sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!)]
    $
    with
    $
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (sum_(k=1)^(min(Y_i^A,Y_i^B)) (lambda_0^(k-1) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/((k-1)!(Y_i^A-k)!(Y_i^B-k)!))/(sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!))
    $

    To differentiate the logarithm above we need the derivative with respect to
    $lambda_0$ of the sum $P(Y_i^A, Y_i^B)$. Only the $lambda_0^k$ factor depends on
    $lambda_0$, and $k / k! = 1 / (k-1)!$, so

    $
      (d)/(d lambda_0) [P(Y_i^A, Y_i^B)] &= sum_(k=1)^(min(Y_i^A,Y_i^B)) (k lambda_0^(k-1) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!) \
      &= sum_(k=1)^(min(Y_i^A,Y_i^B)) (lambda_0^(k-1) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/((k-1)!(Y_i^A-k)!(Y_i^B-k)!)
    $

    The $k = 0$ term contributes nothing because the power rule supplies a factor
    $k$, which is why the sum now starts at $k = 1$; it may equivalently be written
    from $k = 0$. Substituting $j = k - 1$ then turns the right-hand side into

    $
      sum_(j=0)^(min(Y_i^A - 1,Y_i^B - 1)) (lambda_0^(j) lambda_1^((Y_i^A - 1) -j) lambda_2^((Y_i^B - 1) -j))/(j!((Y_i^A - 1) -j)!((Y_i^B - 1) -j)!) = P(Y_i^A - 1,Y_i^B - 1)
    $

    that is, into $P$ evaluated at both counts reduced by one. This single identity
    is the step that makes the whole likelihood analysis tractable, and it lets the
    log-likelihood and score function be written compactly as

    $
      cal(l)(lambda_0; (Y^A, Y^B)) = -n(lambda_0 + lambda_1 + lambda_2) + sum_(i=1)^n log[P(Y_i^A, Y_i^B)]
    $
    $
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (P(Y_i^A - 1,Y_i^B - 1))/(P(Y_i^A,Y_i^B))
    $

    These two forms are the ones used for the numerical work in Task 4 and for the
    Fisher information in Task 9.
  ]
)
== Task 4
#showybox(
  [
    Using *R*, for a *single* sample of size $n = 30$ of $(upright(bold(Y)^A), upright(bold(Y)^B))$ simulated using the function `simulate_YA_YB` for your seed, plot the observed score function and, using this or otherwise, find the maximum likelihood estimate (MLE) for this sample. Store the value of the MLE as `Quantity1`.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    We shall now use the result from above in *Task 3* to help, that is

    $
      cal(l)(lambda_0; (Y^A, Y^B)) = -n(lambda_0 + lambda_1 + lambda_2) + sum_(i=1)^n log[P(Y_i^A, Y_i^B)]
    $
    $
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (P(Y_i^A - 1,Y_i^B - 1))/(P(Y_i^A,Y_i^B))
    $
    for
    $
      P(a,b) = sum_(k=0)^(min(a,b)) (lambda_0^k lambda_1^(a-k) lambda_2^(b-k))/(k!(a-k)!(b-k)!)
    $
  ],
  [
    Evaluating these on the single simulated sample of $n = 30$ days and solving
    $U(lambda_0) = 0$ gives the maximum likelihood estimate *`Quantity1` $= #Quantity1$*.
    Figure 1 shows the log-likelihood, which attains its maximum at that value, and
    the score function, whose unique root it is. The log-likelihood is also larger
    there than at the true $lambda_0 = 2$, as it must be.

    #figure(
      image("task4_plots.png"),
      caption: [Observed log-likelihood (left) and observed score function (right), with the MLE at $U=0$.]
    )
  ]
)
== Task 5
#showybox(
  [
    Given a random sample $(bold(Y)^A, bold(Y)^B)$ of size $n$, we may define an alternative estimator $T$ for $lambda_0$ by \

    $
      T(upright(bold(Y)^A), upright(bold(Y)^B)) = (upright(overline(bold(Y)^A)) + upright(overline(bold(Y)^B)) - lambda_1 - lambda_2)/2
    $
    where
    $
      upright(overline(bold(Y)^A)) = 1/n sum_(i=1)^n Y_i^A
    $
    denotes the sample mean of $upright(bold(Y)^A)$ and $upright(overline(bold(Y)^B))$ is defined similarly.

    Derive the bias and variance of the estimator $T$, and comment.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    We now show for bias of the estimator $T$. \
    $
      BB"ias"(T(upright(bold(Y)^A), upright(bold(Y)^B))) &= EE[T(upright(bold(Y)^A), upright(bold(Y)^B))] - lambda_0 \
      &= EE[(1/n sum_(i=1)^n Y_i^A + 1/n sum_(i=1)^n Y_i^B - lambda_1 - lambda_2)/2] - lambda_0 \
      &= ((lambda_0 + lambda_1) + (lambda_0 + lambda_2) - lambda_1 - lambda_2)/2 - lambda_0 \
      &= 0 quad "unbiased"
    $
    Now for its variance, since it is unbiased, $VV"ar"(T(upright(bold(Y)^A), upright(bold(Y)^B))) = MM"SE"(T(upright(bold(Y)^A), upright(bold(Y)^B)))$
    $
      VV"ar"(T(upright(bold(Y)^A), upright(bold(Y)^B))) &= VV"ar"((1/n sum_(i=1)^n Y_i^A + 1/n sum_(i=1)^n Y_i^B - lambda_1 - lambda_2)/2) \
      &= 1/4 (1/n^2 sum_(i=1)^n VV"ar"(Y_i^A + Y_i^B)) \
      &= 1/4 (1/n (VV"ar"(Y_i^A) + VV"ar"(Y_i^B) + 2CC"ov"(Y_i^A,Y_i^B))) \
      &= 1/4 ((lambda_0 + lambda_1)/n + (lambda_0 + lambda_2)/n + (2lambda_0)/n) \
      &= 1/(4n) (4lambda_0 + lambda_1 + lambda_2)
    $
    Check for consistency
    $
      lim_(n -> infinity) 1/(4n) (4lambda_0 + lambda_1 + lambda_2) = 0
    $
    Indeed, the estimator is consistent.
    
    We can comment that this estimator is exactly unbiased, not merely
    asymptotically so, and that it is consistent; it is also very cheap to compute,
    since it needs only the two sample means.
  ]
)
== Task 6
#showybox(
  [
    Using *R*, for a *single* size $n=30$ of $(upright(bold(Y)^A), upright(bold(Y)^B))$ simulated using the function `simulate_YA_YB` for your seed, calculate $T$ and store the value as `Quantity2`. Comment on how it compares to value of the MLE you have calculated in Task 4.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    The two estimates are close. $T$ gives `Quantity2 = 1.65` against the maximum
    likelihood estimate `Quantity1 =` #Quantity1, a difference of only $0.035$ on an
    estimate of about $1.7$. Both fall below the true $lambda_0 = 2$, and for the
    same reason: both sample means are below their expectations
    (`mean(Y^A) = 4.5` against $lambda_0 + lambda_1 = 5$, and `mean(Y^B) = 5.8`
    against $lambda_0 + lambda_2 = 6$). With a standard error of about $0.35$ at
    $n = 30$ (Task 5), shortfalls of this size are ordinary sampling variation, and
    since $T$ is exactly unbiased this is not a systematic error.
  ]
)
== Task 7
#showybox(
  [
    Using *R*, through repeated sampling of samples of $n=30$, plot a histogram of the observed values of $T$ and comment on its shape and main characteristics.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    The sampling is repeated $10,000$ times. Because `simulate_YA_YB()` sets the
    seed internally, calling it repeatedly would return the same sample every time.
    The independent samples are therefore obtained by simulating
    $n times 10,000 = 300,000$ days in a single call,
    `simulate_YA_YB(30 * 10000)`, and splitting the result into $10,000$
    consecutive blocks of $30$ days with the `matrix()` function. Since the days are
    i.i.d. and the prescribed seed is respected, those blocks are independent
    samples of size $n = 30$.
  ],
  [
    The $10,000$ values of $T$ plotted in Figure 2 are centred on $2 = lambda_0$
    (sample mean #sim_Mean_T, median exactly #sim_Median_T), which is the empirical face of the
    result in Task 5 that $T$ is unbiased.

    The spread also matches Task 5: the observed standard deviation is #sim_SD_T
    against the theoretical #SD_T, and the central 95% of the values lie between
    #sim_Q1 and #sim_Q4. The normal density with mean $2$ and variance $0.125$,
    overlaid in Figure 2, follows the bars closely.

    #figure(
      image("task7_hist.png"),
      caption: [Histogram of estimator $T$ from 10,000 independent samples of size $n=30$, with the normal approximation of Task 5 overlaid.]
    )
  ]
)
== Task 8
#showybox(
  [
    Using *R*, for a *single* size $n=30$ of $(upright(bold(Y)^A), upright(bold(Y)^B))$ simulated using the function `simulate_YA_YB` for your seed, provide an approximate $95%$ Confidence Interval for $T$ by making use of your estimate `Quantity2` (i.e. not using the value $lambda_0 = 2$ used in the simulation) and the distribution obtained via Task 7. Store the values of the lower limit as `Quantity3` and the upper limit as `Quantity4`.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    Because we may not use the true value $lambda_0 = 2$, the usual pivot
    $(T - lambda_0)/sigma$ is unavailable. We therefore use only the *shape* of the
    distribution obtained in Task 7: we take the $0.025$ and $0.975$ quantiles of
    $T$ after subtracting the centre of that distribution (its sample mean,
    #sim_Mean_T), computed with R's `quantile()` function, and re-centre them on our
    own estimate by adding `Quantity2`.

    This gives `Quantity3` $= #Quantity3$ and `Quantity4` $= #Quantity4$, an interval of
    width #Quantity3_Quantity4_width. The interval is wide because $n = 30$ is small, but it covers
    both the true $lambda_0 = 2$ and the maximum likelihood estimate of Task 4.
  ]
)
== Task 9
#showybox(
  [
    Using *R*, through repeated sampling of samples of size $n=30$, estimate the Fisher information and hence the Cramer-Rao Lower Bound (CRLB) for an unbiased estimator for $lambda_0$ for a random sample of size $n=30$ using $lambda_0=2$. Save the estimated Fisher information $I(2)$ as `Quantity5` and the CRLB as `Quantity6`. Comment on your results.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    We recall the score function from Task 4,
    $
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (P(Y_i^A - 1,Y_i^B - 1))/(P(Y_i^A,Y_i^B))
    $
    At the true value $EE[U(lambda_0)] = 0$, so the Fisher information is
    $I(lambda_0) = VV"ar"[U(lambda_0)] = EE[U(lambda_0)^2]$. Averaging
    $U(2)^2$ over $200,000$ simulated samples, each of size $n = 30$, gives

    $
      I(2) = Quantity5 quad "and therefore" quad "CRLB" = 1/(I(2)) = Quantity6.
    $

    The CRLB of #Quantity6 is the smallest variance that any unbiased estimator of
    $lambda_0$ can have at $n = 30$. The estimator $T$, although unbiased, has
    variance #Variance_T and so lies about $5%$ above that bound: it is not an
    efficient estimator. The maximum likelihood estimator, by contrast, is
    asymptotically efficient and approaches the bound.
  ]
)
== Task 10
#showybox(
  [
    Using *R*, plot the variance of $T$ against the CRLB for $n$ varying from $10$ to $50$, and comment.
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
  [
    Figure 3 compares the theoretical variance of $T$ from Task 5, the simulated
    variance of $T$, and the CRLB, as $n$ ranges from $10$ to $50$.

    Both curves fall as $1/n$. From Task 5,
    $VV"ar"(T) = (4lambda_0 + lambda_1 + lambda_2)/(4n) = 15/(4n)$, while for a
    sample of size $n$ the CRLB is $1/(n i(2))$, using the information per
    observation $i(2) = I(2)/30 =$ #i2 obtained in Task 9. The simulated
    variances lie on the theoretical curve at every $n$: for example $0.383$
    against $0.375$ at $n = 10$, and $0.0749$ against $0.075$ at $n = 50$ — which
    confirms the variance formula derived in Task 5.

    The right-hand panel shows that the ratio of $VV"ar"(T)$ to the CRLB is flat at
    about #Var_CRLB_Ratio for every $n$. The gap therefore does not close as the sample
    grows: efficiency is a property of the estimator, not of the sample size. The
    reason is that $T$ discards information by reducing each day to
    $overline(Y)^A + overline(Y)^B$, which does not separate the shared source
    $X_0$ from the two zone-specific sources. That loss is a fixed fraction of the
    available information, so it persists at every $n$.

    #figure(
      image("task10_crlb.png"),
      caption: [Theoretical variance of $T$, simulated variance of $T$ and the CRLB for $n$ from 10 to 50 (left), and their ratio (right).],
    )
  ]
)
== Conclusion

The aim of this report was to estimate the rate $lambda_0$ of the tremor source
shared by two zones either side of a fault line, given that the zone-specific
rates $lambda_1 = 3$ and $lambda_2 = 4$ are known, and to establish how
precisely that can be done from $n = 30$ days of observations.

The key modelling point is that the two zone counts are *not* independent. Both
contain the shared term $X_0$, so their joint distribution is a bivariate
Poisson distribution rather than the product of two Poisson distributions, and
their covariance is exactly $lambda_0$. This is not a technicality. For the
moment-based estimator $T$ studied here, wrongly treating the zones as
independent would give a variance of $0.0917$ in place of the correct $0.125$,
understating the variance by about 27% and the standard error by about 14%.

For the single sample of $n = 30$ days the two estimators agreed closely. The
maximum likelihood estimate was $hat(lambda)_0 = 1.6852$ and the moment-based
alternative was $T = 1.65$, a difference of only $0.035$. Both fell below the
true $lambda_0 = 2$, which is ordinary sampling variation at this sample size:
both zone means were below their expectations ($overline(Y)^A = 4.5$ against $5$, and
$overline(Y)^B = 5.8$ against $6$), and the standard error of either estimator is
about $0.35$. The approximate 95% confidence interval for $T$, built from its
simulated sampling distribution, was $(0.9668, 2.3501)$, and it covered both the
true value and the maximum likelihood estimate. Table 1 collects the six
quantities computed by the accompanying R script.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, center, left),
    inset: 6pt,
    table.header([*Quantity*], [*Value*], [*Definition*]),
    [`Quantity1`], [#Quantity1], [Maximum likelihood estimate of $lambda_0$ (Task 4)],
    [`Quantity2`], [#Quantity2], [Moment estimator $T$ (Task 6)],
    [`Quantity3`], [#Quantity3], [Lower 95% limit for $T$ (Task 8)],
    [`Quantity4`], [#Quantity4], [Upper 95% limit for $T$ (Task 8)],
    [`Quantity5`], [#Quantity5], [Fisher information $I(2)$ (Task 9)],
    [`Quantity6`], [#Quantity6], [Cramér--Rao lower bound (Task 9)],
  ),
  caption: [The six quantities calculated by the accompanying R script.],
)

The two estimators are not equally good, however. The moment-based estimator $T$
is exactly unbiased and consistent, and almost costless to compute, but its
variance of #Variance_T at $n = 30$ is about 5% above the Cramér--Rao lower bound of
#Quantity6 obtained from the estimated Fisher information $I(2) =$ #Quantity5, whereas
the maximum likelihood estimator is asymptotically efficient and so approaches
that bound as the sample size grows (Casella and Berger, 2002). The difference
arises because $T$ discards information: it reduces each day to
$overline(Y)^A + overline(Y)^B$, which does not distinguish the shared source
$X_0$ from the zone-specific sources $X_1$ and $X_2$, while the likelihood uses
the full joint distribution of every pair. Task 10 shows this loss is structural
rather than a small-sample effect: the ratio of the variance of $T$ to the bound
is flat at about #Var_CRLB_Ratio for every sample size from 10 to 50.

Two practical conclusions follow. First, the maximum likelihood estimator is the
one to use: it is unbiased and efficient, and cheap to compute. Second, the
precision available from 30 days is poor. The confidence interval above extends
roughly 40% either side of the estimate, and since the variance falls as $1/n$,
halving the standard error would require four times as many observations. One
month of data establishes that a shared source exists and gives a rough value
for its rate, but does not pin it down tightly.

Several limitations should be borne in mind. The analysis is conditional on the
assumed model: the data are simulated from it, so the model is never questioned
and no goodness-of-fit check is possible. The shared rate is assumed constant
across the $n$ days, ruling out seasonal effects and the clustering of tremors
into aftershock sequences; a time-varying or self-exciting rate would be the
natural extension (Ogata, 1988). The Monte Carlo estimate of the Fisher
information carries a standard error of about 0.3%, though it agrees with an
exact summation over the joint distribution to within that tolerance. Finally,
the confidence interval assumes the sampling distribution of $T$ is close to
normal, which Task 7 supports at $n = 30$ but which would need re-checking
below that.

== References

+ Casella, G., and Berger, R. L. (2002). _Statistical Inference_ (2nd ed.).
  Pacific Grove, CA: Duxbury.

+ Holgate, P. (1964). Estimation for the bivariate Poisson distribution.
  _Biometrika_, 51(1/2), 241--245. #link("https://doi.org/10.1093/biomet/51.1-2.241")

+ Karlis, D., and Ntzoufras, I. (2003). Analysis of sports data by using
  bivariate Poisson models. _Journal of the Royal Statistical Society: Series D
  (The Statistician)_, 52(3), 381--393.

+ Kocherlakota, S., and Kocherlakota, K. (1992). _Bivariate Discrete
  Distributions_. New York: Marcel Dekker.

+ Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
  analysis for point processes. _Journal of the American Statistical
  Association_, 83(401), 9--27.

== Appendix: Use of generative AI

A large language model was used in this project as a discussion partner and
checking tool: to help design and debug the R script, to verify the algebraic
derivations of Tasks 2, 3 and 5 numerically, and to review and revise the wording
of this report. All derivations, code and written text were checked by the
author, who takes responsibility for the accuracy of the content.