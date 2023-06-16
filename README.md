# torustest

**PCIdep** is an $\texttt{R}$ package implementing the approaches introduced in [1] to perform selective inference after hierarchical or $k$-means clustering when the observations are drawn from a general matrix normal model:

$$
\mathbf{X}\sim\mathcal{MN}_{n\times p}(\boldsymbol\mu,\mathbf{U},\mathbf{\Sigma}),
$$

where $\boldsymbol\mu\in\mathcal{M}_{n\times p}(\mathbb{R})$ and $\mathbf{U}\in\mathcal{M}_{n\times n}(\mathbb{R})$ and $\mathbf{\Sigma}\in\mathcal{M}_{p\times p}(\mathbb{R})$ are positive definite matrices encoding the depence structure between observations and features respectively.

**PCIdep** is the natural extension to the general matrix normal model of the work in [2] ([clusterpval](https://github.com/lucylgao/clusterpval)) and [3] ([KMeansInference](https://github.com/yiqunchen/KmeansInference)) where the framework for selective inference after hierarchical clustering and $k$-means respectively is presented when $\mathbf{U}=\mathbf{I}_n$ and $\mathbf{\Sigma}=\sigma\mathbf{I}_p$.

### Installing PCIdep

**PCIdep** can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")
```

For any inquires, please file an [issue](https://github.com/gonzalez-delgado/PCIdep/issues) or [contact us](mailto:javier.gonzalez-delgado@math.univ-toulouse.fr).

### Two-sample goodness-of-fit tests

Given two measures $P,Q\in\mathcal{P}(\mathbb{T}^2)$, we consider the null hypothesis $H_0:P=Q$. The package **torustest** allows the assessment of $H_0$ through two different procedures, detailed below.

#### 1. Projection into one-dimensional closed geodesics

The first approach tests the equality of the one-dimensional projected distributions of $P$ and $Q$ into $N_g$ closed geodesics of $\mathbb{T}^2$, which are isomorphic to the circle $S^1$. Given a pair of samples in the periodic $[0,1]\times[0,1]$, the function [twosample.geodesic.torus.test](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.geodesic.torus.test.R) returns a $p$-value for $H_0$ after specifying the number of projections $N_g$, which may be randomly chosen.

#### Example

```
n <- 100 # Sample size
 
# Simulate the null distribution of the circle test statistic
sim_free_null <- torustest::sim.null.stat(500, NC = 2)

# Bivariate von Mises distributions
samp_1 <- BAMBI::rvmcos(n) / (2 * pi) 
samp_2 <- BAMBI::rvmcos(n) / (2 * pi)

# 4 geodesics are chosen randomly
torustest::twosample.geodesic.torus.test(samp_1, samp_2, n_geodesics = 3, NC_geodesic = 2, sim_null = sim_free_null) 

# 4 geodesics are chosen a priori
glist <- list(c(1, 0), c(0, 1), c(1, 1), c(2, 3))
torustest::twosample.geodesic.torus.test(samp_1, samp_2, geodesic_list = glist, NC_geodesic = 2, sim_null = sim_free_null) 

```

#### Two-sample goodness-of-fit for measures on $S^1$

For each pair of projected measures, a Wasserstein two-sample goodness-of-fit test for measures supported on the circle is performed. Such test is also available for any pair of input samples in the periodic $[0,1]$, via the function [twosample.test.s1.R](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.test.s1.R).

#### Example

```
n <- 50 # Sample size
 
# Simulate the statistic null distribution
NR <- 100
sim_free_null <- torustest::sim.null.stat(500, NC = 2)

x <- runif(n, 0, 1)
y <- runif(n, 0, 1)
torustest::twosample.test.s1(x, y, sim_free_null) 

x <- as.numeric(circular::rvonmises(n, pi, 1)/(2*pi))
y <- as.numeric(circular::rvonmises(n, pi, 0)/(2*pi))
torustest::twosample.test.s1(x, y, sim_free_null) 
```

### 2. $p$-value upper bound

The second approach computes an upper bound for the $p$-value $\mathbb{P}_{H_0}(\mathcal{T}_2(P_n,Q_m) \geq t)$, where $P_n,Q_m$ are the empirical probability measures of $P,Q$ respectively, the statistic $\mathcal{T}_2(P_n,Q_m)$ denotes their squared $2$-Wasserstein distance and $t$ the statistic realization for the pair of samples considered. $H_0$ is rejected at level $\alpha>0$ if the upper bound is less than $\alpha$. This test is asymptotically consistent at level $\alpha$, for any $\alpha>0$, and has reasonable power for large sample sizes. The function [twosample.ubound.torus.test function](https://github.com/gonzalez-delgado/torustest/blob/master/R/twosample.ubound.torus.test.R) returns the upper-bound for a pair of samples in the periodic $[0,1]\times[0,1]$.

#### Example

```
n <- 2000 # Sample size
 
# Bivariate von Mises distribution
samp_1 <- BAMBI::rvmcos(n, kappa1 = 1, kappa2 = 1, mu1 = 0, mu2 = 0)/(2*pi)
samp_2 <- BAMBI::rvmcos(n, kappa1 = 1, kappa2 = 1, mu1 = 0, mu2 = 0)/(2*pi)
torustest::twosample.ubound.torus.test(samp_1, samp_2) 

samp_1 <- BAMBI::rvmcos(n ,kappa1 = 0, kappa2 = 0, mu1 = 0.5, mu2 = 0.5)/(2*pi)
samp_2 <- BAMBI::rvmcos(n, kappa1 = 1, kappa2 = 1, mu1 = 0.5, mu2 = 0.5)/(2*pi)
torustest::twosample.ubound.torus.test(samp_1, samp_2) 
```

### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. arXiv:2108.00165. [[arxiv]](https://arxiv.org/abs/2108.00165)[[HAL]](https://hal.archives-ouvertes.fr/hal-03369795v2).
