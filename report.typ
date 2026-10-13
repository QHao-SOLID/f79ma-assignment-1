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
== Task 2
#showybox(
    [State and explain the distributions of $Y^A$ and $Y^B$, and calculate the covariance of $Y^A$ and $Y^B$.]
)

#showybox(
  shadow: (
    offset: 3pt,
  ),
  title: "Answer",
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
== Task 4
== Task 5
== Task 6
== Task 7
== Task 8
== Task 9
== Task 10
== Conclusion