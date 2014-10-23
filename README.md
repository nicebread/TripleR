# TripleR

Social Relation Model (SRM) analyses for single or multiple round-robin groups are performed. These analyses are either based on one manifest variable, one latent construct measured by two manifest variables, two manifest variables and their bivariate relations, or two latent constructs each measured by two manifest variables. Within-group t-tests for variance components and covariances are provided for single groups, between-groups t-tests for multiple groups. Handling for missing values is provided.

## Installation

The stable version can be installed from [CRAN](http://cran.r-project.org/web/packages/TripleR/index.html):

    install.packages("TripleR")

The current development version can be installed from this repository:

    install.packages(c("devtools", "plyr", "reshape2", "ggplot2"), dependencies=TRUE)
    library(devtools)
    install_github("nicebread/TripleR", subdir="package")	

