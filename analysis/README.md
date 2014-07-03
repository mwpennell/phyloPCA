# Statistical and conceptual challenges in the comparative analysis of principal components -- analysis

The purpose of this project was to investigate the effects of different ways of reducing the dimensionality of phylogenetic comparative data on evolutionary inferences. There are two components to this project: a simulation study and an analysis of two empirical datasets.

This project relied on a number of R packages for both the analyses and plotting the data:

* ape
* geiger
* phylolm
* phytools
* MASS
* nmle
* reshape2
* plyr
* foreach
* doMC
* ggplot2
* RColorBrewer
* scales
* lattice

To aid reproducibility we have also integrated our project with [knitr](http://yihui.name/knitr/). This requires the following R packages

* knitr
* markdown
* [sowsear](https://github.com/richfitz/sowsear)

All of these dependencies can be installed using [make](http://www.gnu.org/software/make/)
```
make deps
```

## Analysis and figures
Typing
```
make analysis
```
will generate a Rmd file `phylo-pc.Rmd`, a markdown file `phylo-pc.md` and a html document `phylo-pc.html`. This will also generate the figures used ine manuscript which will be dropped into `output/figs`.

Because the simulations may take a while, we have cached the results in the folder `output/sim-res`.

To regenerate all simulations and results from scratch, type
```
make sim-res
```

## Empirical datasets

To examine the effect of data transformations on real comparative data. The first is morphological measurements of Felidae ("cats"). The data was compiled from [Slater and Van Valkenburgh 2009](http://www.psjournals.org/doi/abs/10.1666/07061.1) and [Sakamoto, et al. 2010](http://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2009.01922.x/full). For this example we used the Carnivora phylogeny from  [Nyakatura and Bininda-Emonds 2012](http://www.biomedcentral.com/1741-7007/10/12). The second dataset was from Cyrinodon fishes ("pupfishes") from a study by [Martin and Wainwright 2011](http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2011.01294.x/full).

We did some processing of the data prior to analyses. The resulting datasets are stored in `output/data`. To re-run the processing steps, type
```
make clean-emp-data
make emp-data
```

If you have any questions please feel free to contact us by email or submit an issue.


