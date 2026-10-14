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

    1. $X_0, X_1$ and $X_2$ are independent
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
    We use the idea of Probability Generating Functions $G_x (s)$ to check for distribution of $Y^A$ and $Y^B$. \ 
    
    1. It is known that for $X ~ "Poi"(lambda)$, $G_x (s) = e^(-lambda(1-s))$, and \
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
      U(lambda_0; (Y^A, Y^B)) = -n + sum_(i=1)^n ([sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k-1) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/((k-1)!(Y_i^A-k)!(Y_i^B-k)!)])/(sum_(k=0)^(min(Y_i^A,Y_i^B)) (lambda_0^(k) lambda_1^(Y_i^A-k) lambda_2^(Y_i^B-k))/(k!(Y_i^A-k)!(Y_i^B-k)!))
    $
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
      &= 1/4 (1/n^2 sum_(i=1)^n VV"ar"(Y_i^A) + 1/n^2 sum_(i=1)^n VV"ar"(Y_i^B)) \
      &= 1/4 ((lambda_0 + lambda_1)/n + (lambda_0 + lambda_2)/n) \
      &= 1/(4n) (2lambda_0 + lambda_1 + lambda_2)
    $
    Check for consistency
    $
      lim_(n -> infinity) 1/(4n) (2lambda_0 + lambda_1 + lambda_2) = 0
    $
    Indeed, the estimator is consistent.
    
    We can comment that the following estimator is asymptotically unbiased, as its bias is zero and is consistent.
  ]
)
== Task 6
#showybox(
  [
    
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
)
== Task 7
#showybox(
  [
    
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
)
== Task 8
#showybox(
  [
    
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
)
== Task 9
#showybox(
  [
    
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
)
== Task 10
#showybox(
  [
    
  ]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
  breakable: true,
)
== Conclusion