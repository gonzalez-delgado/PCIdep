# PCIdep: Post-Clustering Inference under DEPendence

$\texttt{PCIdep}$ is an $\texttt{R}$ package implementing the approaches introduced in [González-Delgado _et al._ (2023)](http://arxiv.org/abs/2310.11822) to perform selective inference after clustering when the observations are drawn from a general matrix normal model:

$$
\mathbf{X}\sim\mathcal{MN}_{n\times p}(\boldsymbol\mu,\mathbf{U},\mathbf{\Sigma}),
$$

where  $`\boldsymbol\mu\in\mathcal{M}_{n\times p}(\mathbb{R})`$ and $`\mathbf{U}\in\mathcal{M}_{n\times n}(\mathbb{R})`$, $`\mathbf{\Sigma}\in\mathcal{M}_{p\times p}(\mathbb{R})`$ are positive definite matrices encoding the dependence structure between observations and features respectively.

$\texttt{PCIdep}$ is the natural extension to the general matrix normal model of the work in [clusterpval](https://github.com/lucylgao/clusterpval) and [KMeansInference](https://github.com/yiqunchen/KmeansInference), where the framework for selective inference after hierarchical clustering and $k$-means respectively is presented for $\mathbf{U}=\mathbf{I}_n$ and $\mathbf{\Sigma}=\sigma\mathbf{I}_p$, allowing for the estimation of the unknown parameter $\sigma$. Inference after any user-specified clustering algorithm can be performed via Monte Carlo approximation. 

### User guide

For a comprehensive introduction to the package and its functionalities, please refer to the [user guide](https://gonzalez-delgado.github.io/PCIdep/index.html). 

### Installing PCIdep

$\texttt{PCIdep}$ can be installed using

```
devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")
```

For any inquires, please file an [issue](https://github.com/gonzalez-delgado/PCIdep/issues) or [contact us](mailto:javier.gonzalez-delgado@ensai.fr).
