# Assess whether a pair of clusters remains invariant across two partitions.

Assess whether a pair of clusters remains invariant across two
partitions.

## Usage

``` r
preserve.cl(cl, cl_phi, clusters)
```

## Arguments

- cl:

  A vector of size \\n\\ with integers in \\1,\ldots,n\\ defining the
  classes of the first partition.

- cl_phi:

  A vector of size \\n\\ with integers in \\1,\ldots,n\\ defining the
  classes of the second partition.

- clusters:

  A vector of two integers chosen among the coordinates of `cl`,
  denoting the pair of clusters in `cl` to be analyzed.

## Value

Whether the cluster `clusters`\[1\] and the cluster `clusters`\[2\] in
`cl` remain invariant in `cl_phi`, independently of label modifications.
