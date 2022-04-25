<!-- README.md is generated from README.Rmd. Please edit that file 
library(badger)
rmarkdown::render("README.rmd")
-->

[![Build
Status](https://travis-ci.org/nicebread/TripleR.svg)](https://travis-ci.org/nicebread/TripleR)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--04--25-yellowgreen.svg)](/commits/master)
[![](https://www.r-pkg.org/badges/version/TripleR?color=orange)](https://cran.r-project.org/package=TripleR)
[![packageversion](https://img.shields.io/badge/Package%20version-1.5.4-orange.svg?style=flat-square)](commits/master)
![](http://cranlogs.r-pkg.org/badges/TripleR)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:e90f8d849adb0ba98566e7670dd953f90125e259/)](https://archive.softwareheritage.org/swh:1:dir:e90f8d849adb0ba98566e7670dd953f90125e259;origin=https://github.com/nicebread/TripleR;visit=swh:1:snp:5252dbbb315b873bf5f438914c6802a736934196;anchor=swh:1:rev:7e506c6d81dac91f239af4c8db8dea5dfd1ebd1f)

# TripleR: Social Relation Model (SRM) Analyses for Single or Multiple Groups

Social Relation Model (SRM) analyses for single or multiple round-robin
groups are performed. These analyses are either based on one manifest
variable, one latent construct measured by two manifest variables, two
manifest variables and their bivariate relations, or two latent
constructs each measured by two manifest variables. Within-group t-tests
for variance components and covariances are provided for single groups.
For multiple groups two types of significance tests are provided:
between-groups t-tests (as in SOREMO) and enhanced standard errors based
on [Lashley and Bond
(1997)](https://psycnet.apa.org/doiLanding?doi=10.1037%2F1082-989X.2.3.278).
Handling for missing values is provided.

*Please note the newer R package
[srm](https://CRAN.R-project.org/package=srm) which provides SRM
analyses based on structural equation modeling, and supports both
maximum likelihood estimation and least squares estimation.*

The package is described in:

Schönbrodt, F. D., Back, M. D., & Schmukle, S. C. (2012). TripleR: An R
package for social relations analyses based on round-robin designs.
*Behavior Research Methods, 44*, 455–470.
<doi:%5B10.3758/s13428-011-0150-4>\](<https://link.springer.com/article/10.3758/s13428-011-0150-4>).
[Open Access version](https://osf.io/xx267/)

Please cite this paper if you use the package in your publication.

## Installation

The stable version can be installed from
[CRAN](http://cran.r-project.org/web/packages/TripleR/index.html):

    install.packages("TripleR")

The current development version can be installed from this repository:

    install.packages(c("devtools", "plyr", "reshape2", "ggplot2"), dependencies=TRUE)
    library(devtools)
    install_github("nicebread/TripleR", subdir="package")   

# Maintainers

[Felix Schönbrodt](https://www.nicebread.de/)

[Mitja
Back](https://www.uni-muenster.de/PsyIFP/AEBack/members/mitja-back.html)

[Stefan
Schmukle](https://www.lw.uni-leipzig.de/wilhelm-wundt-institut-fuer-psychologie/arbeitsgruppen/persoenlichkeitspsychologie-und-psychologische-diagnostik/team/prof-dr-stefan-schmukle)
