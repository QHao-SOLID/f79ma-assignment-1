#import "@preview/showybox:2.0.4": showybox
#set par(
    justify: true
)

= F79MA - Assessed Project 1 - 2026/27
- Student name: Qing Hao, Too
- HW-ID: H00467830
- Campus: Heriot-Watt University Malaysia

== Introduction
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
    For the following model, we've made a few assumptions, that is:

    1. $X_0, X_1$ and $X_2$ are independent, that means, what happens in the future does not depend on what is happening currently.
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
    2. For $Y^A = X_0 + X_1 => G_(Y^A)(s) =^"def." EE[S^(X_0 + X_1)] =^"i.i.d" EE[S^(X_0)]EE[S^(X_1)]$. \
    
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
      PP(Y^A = a, Y^B = b) &= sum_(k=0)^(min(a,b)) PP(X_0 = k)PP(X_1 = a - X_0, X_2 = b - X_0 | X_0 = k) \
      &=^"ind." sum_(k=0)^(min(a,b)) PP(X_0 = k)PP(X_1 = a - k) PP(X_2 = b-k) \
      &=^"p.m.f." sum_(k=0)^(min(a,b)) (e^(-lambda_0) lambda_0^(k))/(k!) times (e^(-lambda_1) lambda_1^((a-k)))/((a-k)!) times (e^(-lambda_2) lambda_2^((b-k)))/((b-k)!) \
      &= e^(-(lambda_0 + lambda_1 + lambda_2)) sum_(k=0)^(min(a,b)) (lambda_0^(k) lambda_1^(a-k) lambda_2^(b-k))/(k!(a-k)!(b-k)!)
    $

    Observe that with use $min(a,b)$ as the conditioning limit, that is due to satisfying the condition $a-k>=0, b-k>=0$ and $k>=0$ for Poisson only accepts non-negative inputs.
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

    #showybox(
      title: "Optional info.",
      frame: (
        border-color: blue,
        title-color: blue.lighten(30%),
        body-color: blue.lighten(95%)
      ),
      [
        Some interesting notes from $U(lambda_0; (Y^A, Y^B))$ is the summation below:

        $
          P(Y_i^A, Y_i^B) = sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!)
        $
        When we differentiate the summation (w.r.t $lambda_0$), we get
        $
          (d)/(d lambda_0) [P(Y_i^A, Y_i^B)] = sum_(k=1)^(min(Y_i^A,Y_i^B)) (lambda_0^(k-1) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/((k-1)!(Y_i^A-k)!(Y_i^B-k)!)
        $
        Let $j = k-1$, we now have
        $
          sum_(j=0)^(min(Y_i^A - 1,Y_i^B - 1)) (lambda_0^(j) lambda_1^((Y_i^A - 1) -j) lambda_2^((Y_i^B - 1) -j))/(j!((Y_i^A - 1) -j)!((Y_i^B - 1) -j)!) = P(Y_i^A - 1,Y_i^B - 1)
        $
        So
        $
          cal(l)(lambda_0; (Y^A, Y^B)) = -n(lambda_0 + lambda_1 + lambda_2) + sum_(i=1)^n log[P(Y_i^A, Y_i^B)]
        $
        $
          U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (P(Y_i^A - 1,Y_i^B - 1))/(P(Y_i^A,Y_i^B))
        $

        This will be very useful for the tasks ahead!
      ]
    )
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
      P(a,b) = sum_(k=0)^(min(a,b)) (lambda_0^k lambda_1^(a-k) lambda_2^(b-k))/(k!(a-k)!(b-k!))
    $
  ],
  [
    Below is the plot for the log-likelihood function $cal(l)(lambda_0; (Y^A, Y^B))$ as well as the score function $U(lambda_0; (Y^A, Y^B))$.

    #figure(
      image("task4_plots.png"),
      caption: [Observed log-likelihood (left) and Observed score function (right) with MLE at $U=0$.]
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
      &= EE[(1/n sum_(i=1)^n Y_i^A + 1/n sum_(i=1)^n Y_i^A - lambda_1 - lambda_2)/2] - lambda_0 \
      &= ((lambda_0 + lambda_1) + (lambda_0 + lambda_2) - lambda_1 - lambda_2)/2 - lambda_0 \
      &= 0 quad "unbiased"
    $
    Now for its variance, since it is unbiased, $VV"ar"(T(upright(bold(Y)^A), upright(bold(Y)^B))) = MM"SE"(T(upright(bold(Y)^A), upright(bold(Y)^B)))$
    $
      VV"ar"(T(upright(bold(Y)^A), upright(bold(Y)^B))) &= VV"ar"((1/n sum_(i=1)^n Y_i^A + 1/n sum_(i=1)^n Y_i^A - lambda_1 - lambda_2)/2) \
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
    
    We can comment that the following estimator is asymptotically unbiased, as its bias is zero and is consistent.
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
    Comparing $T$ (`Quantity2 = 1.65`) and MLE (`Quantity1 = 1.6852`), differ only by $0.0352$, but both falls below the true $lambda_0 = 2$ because both sample shows mean below the true mean (`mean(Y^A) = 4.5` vs. $lambda_0 + lambda_1 = 5$ and `mean(Y^B) = 5.8` vs. $lambda_0 + lambda_2 = 6$).
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
    The sampling process is repeated for total of $10,000$ times under the same seed.
    
    The way such sample is obtained is done through `simulate_YA_YB(n * n_samples) = simulate_YA_YB(30 * 10000)` as so to comply with fixed seed rule. Then we split to $10,000$ group via the `matrix()` function, each entry have exactly $30$ samples.
  ],
  [
    Below shows a plot of a histogram of the observed values of $T$. With comparison to Normal distribution of parameter $mu = 2, sigma^2 = 0.125$ ($sigma^2$ value obtained through Task 5's $VV"ar"(T(upright(bold(Y)^A), upright(bold(Y)^B)))$).

    #figure(
      image("task7_hist.png"),
      caption: [Histogram plot of estimator $T$ from 10,000 independent samples of size $n=30$ with comparison to Normal distribution.]
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
    As not permitted to use true parameter $lambda_0 = 2$, we cannot use the standard pivot variable $(T - mu)/(sigma)$, so we use utilise R's `quantile()` function to obtain the $0.025$ and $0.975$ percentile. But before that, we want an interval centred on our own estimate $T$, so we keep only the shape of that distribution: subtract Task 7's sample data mean then add `Quantity2`.

    We will then obtain the interval $(0.966762, 2.350095)$, the width of $1.383333$ is fair for we are only given $n=30$.
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
    We recall the score function in Task 4,
    $
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n (P(Y_i^A - 1,Y_i^B - 1))/(P(Y_i^A,Y_i^B))
    $
    The Fisher information is given by $I(theta) = EE[U^2]$. We will estimate it via computing the average score squared of $200,000$ samples with each having $n = 30$.

    The result estimates $I(2) = 8.363433$ and so $"CRLB" = 1/(I(2)) = 0.119568$. The estimated CRLB tells us that $0.119568$ is the smallest variance for any unbiased estimators of $lambda_0$ at $n=30$. This tells us that the estimator $T$ despite being unbiased, it is not the *most* efficient estimator as it's variance sits right above CRLB at $0.125$.
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
    Below illustrates the plot of theoretical variance of $T$, CRLB as well as the simulated variance of $T$.

    #figure(
      image("task10_crlb.png"),
      caption: "Comparison of Theoretical Variance, Simulated Variance and CRLB for n = 10 to n = 50 (left), Efficiency as n increases (right)",
    )

    Both curve fall as $1/n$, as $VV"ar"(T) = 15/(4n)$ and $"CRLB" = 1/(n dot I(2))$. The plot on the right shows that the eficiency is flat at $approx 1.049$ using the exact $I(2)$, this means that the gap does not close, because $T$ does not include information of the relationship between the two, that is the shared $X_0$. The loss of that information is a fixed fraction of the available information, so it persists at every $n$.
  ]
)
== Conclusion