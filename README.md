# PCIdep

**PCIdep** is an $\texttt{R}$ package implementing the approaches introduced in [1] to perform selective inference after hierarchical or $k$-means clustering when the observations are drawn from a general matrix normal model:

$$
\mathbf{X}\sim\mathcal{MN}_{n\times p}(\boldsymbol\mu,\mathbf{U},\mathbf{\Sigma}),
$$

where  $`\boldsymbol\mu\in\mathcal{M}_{n\times p}(\mathbb{R})`$ and $`\mathbf{U}\in\mathcal{M}_{n\times n}(\mathbb{R})`$, $`\mathbf{\Sigma}\in\mathcal{M}_{p\times p}(\mathbb{R})`$ are positive definite matrices encoding the depence structure between observations and features respectively.

**PCIdep** is the natural extension to the general matrix normal model of the work in [clusterpval](https://github.com/lucylgao/clusterpval) [2] and [KMeansInference](https://github.com/yiqunchen/KmeansInference) [3] where the framework for selective inference after hierarchical clustering and $k$-means respectively is presented for $\mathbf{U}=\mathbf{I}_n$ and $\mathbf{\Sigma}=\sigma\mathbf{I}_p$.

### Installing PCIdep

**PCIdep** can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")
```

For any inquires, please file an [issue](https://github.com/gonzalez-delgado/PCIdep/issues) or [contact us](mailto:javier.gonzalez-delgado@math.univ-toulouse.fr).

### Selective inference after hierarchical clustering

The function [test.clusters.hc](https://github.com/gonzalez-delgado/PCIdep/blob/main/R/test.clusters.hc.R) adapts the framework presented in [clusterpval](https://github.com/lucylgao/clusterpval) to the general matrix normal model. It allows selective inference after hierarchical clustering (HC) with multiple types of linkage. Over-estimation of $`\mathbf{\Sigma}`$ for known $`\mathbf{U}`$ is possible while asymptotically respecting the selective type I error (under some conditions on $`\mathbf{U}`$).

#### Example

```
# Model parameters
n <- 50
p <- 20
M <- Matrix::Matrix(0, nrow = n , ncol = p) # Mean matrix
Sigma <- stats::toeplitz(seq(1, 0.1, length = p)) # Sigma: dependence between features
U <- matrixNormal::I(n) # U: dependence between observations

# Simulate matrix normal samples
X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma) # i.i.d. copy of X

# HC with average linkage under a global null hypothesis
test.clusters.hc(X, U, Sigma, NC = 3, clusters = sample(1:3, 2), plot = TRUE, linkage = "average")

# HC with complete linkage under the global null hypothesis and over-estimation of Sigma
test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, clusters = sample(1:3, 2), plot = TRUE, linkage = "complete")
```

### Selective inference after $k$-means clustering

The function [test.clusters.km](https://github.com/gonzalez-delgado/PCIdep/blob/main/R/test.clusters.km.R) adapts the framework presented in [KMeansInference](https://github.com/yiqunchen/KmeansInference) to the general matrix normal model. It allows selective inference after $k$-means clustering. Over-estimation of $`\mathbf{\Sigma}`$ for known $`\mathbf{U}`$ is possible while asymptotically respecting the selective type I error (under some conditions on $`\mathbf{U}`$). 

#### Example

```
# Model parameters
n <- 50
p <- 20
M <- Matrix::Matrix(0, nrow = n , ncol = p) # Mean matrix
Sigma <- stats::toeplitz(seq(1, 0.1, length = p)) # Sigma: dependence between features
U <- matrixNormal::I(n) # U: dependence between observations

# Simulate matrix normal samples
X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma) # i.i.d. copy of X

# k-means under the global null hypothesis
test.clusters.km(X, U, Sigma, NC = 3, clusters = sample(1:3, 2))

# k-means under the global null hypothesis and over-estimation of Sigma
test.clusters.km(X, U, Sigma = NULL, Y = Y, NC = 3, clusters = sample(1:3, 2))
```

### References

[2] L. L. Gao, J. Bien, and D. Witten. Selective inference for hierarchical clustering. Journal of the American Statistical Association, 0(0):1–11, 2022. 

[3] Y. T. Chen and D. M. Witten. Selective inference for k-means clustering, 2022. arXiv:2203.15267.
