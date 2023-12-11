# LS3MU: Laplacian Spike and Slab Selector for Matrix with Uncertainty

Spike and slab regression methods have been widely used in sparse model selection tasks. Considering high-dimensional sparse regression with measurement error in variables, 
we assume the errors are Normally distributed with mean 0 and errors for each variable have a specific variance. $LS3MU$ package function is developed to conduct variable 
selection tasks for such problem setting. The main function $lsum$ works under EM framework, where the unknown true design matrix and spike-slab prior indicators of regression parameters are 
treated as latent variables. The final output coefficients are obtained following a decreasing sequence of spike parameters and the paths can be displayed by calling functions in 
$LS3MU$ package.

## Installation
install_github("ShuyuZoeyGuo/LS3MU")

## Manual
[Manual](https://github.com/ShuyuZoeyGuo/LS3MU/blob/main/doc/LS3MU-manual.pdf)

## Vignettes
[Vignettes](https://github.com/ShuyuZoeyGuo/LS3MU/blob/main/doc/vignettes.pdf)
