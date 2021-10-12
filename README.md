# sprintr

The `sprintr` package contains implementation of a computationally efficient method to fit large-scale interaction models based on a reluctant interaction selection principle.
The details of the method can be found in 
[Yu, Bien, and Tibshirani (2021) *Reluctant Interaction Modeling*](https://arxiv.org/abs/1907.08414), v2.

To install `sprintr` from [github](http://github.com), type in R console
```R
devtools::install_github("hugogogo/sprintr", build_vignettes = TRUE)
```
Note that the installation above requires using R package [devtools](https://CRAN.R-project.org/package=devtools)
(which can be installed using `install.packages("devtools")`).

A [vignette](https://github.com/hugogogo/sprintr/blob/master/vignettes/using_sprintr.pdf) contains examples of how to use the package.

**Note that the sprinter method is changed in v2 of the manuscript. This package reflects this change.** 
