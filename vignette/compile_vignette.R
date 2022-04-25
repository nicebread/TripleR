library(knitr)
knit('TripleR-knitr.Rnw')
pdfres <- system("latexmk -pdf -f -interaction=nonstopmode TripleR-knitr.tex", intern=TRUE, wait=TRUE)
