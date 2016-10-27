# bbac: Bregman block average co-clustering

This is an R implementation of the [Bregman Block Average Co-clustering (BBAC) algorithm (Banerjee et al., 2007a)](http://www.jmlr.org/papers/volume8/banerjee07a/banerjee07a.pdf). The [BBAC code in Matlab](http://www.ideal.ece.utexas.edu/software.html) written by Srujana Merugu and Arindam Banerjee was used as a reference.

## Example

```{r}
# Load bbac package
source("bbac.R")
set.seed(1)

# Generate synthetic data
Z <- matrix(rep(1:4, 25), 10, 10)

# Run Cheng-Church and Information-Theoretic co-clustering algorithms
cheng_church  <- bbac(Z, k = 2, l = 2, distance = "e", scheme = 6)
inf_theoretic <- bbac(Z, k = 2, l = 2, distance = "d", scheme = 5)

# Show co-clusters
par(mfrow=c(1, 2))
plot_coclusters(Z, cheng_church$R, cheng_church$C)
title(paste("Cheng-Church algorithm", cheng_church$status))
plot_coclusters(Z, inf_theoretic$R, inf_theoretic$C)
title(paste("Information-Theoretic algorithm", inf_theoretic$status))
```