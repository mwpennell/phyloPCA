## Makefile for project

## empirical data
emp-data: output/data/felidae.rds output/data/cyprinodon.rds

output/data/felidae.rds: output/data/felidae.rds
	Rscript make/data-felidae.R

output/data/cyprinodon.rds: output/data/cyprinodon.rds
	Rscript make/data-cyprinodon.R


## produce analysis
analysis: phylo-pc.html

phylo-pc.Rmd: phylo-pc.R
	Rscript -e "library(sowsear); sowsear('$<', 'Rmd')"

phylo-pc.md: phylo-pc.Rmd
	Rscript -e "library(knitr); knit('$<')"

phylo-pc.html: phylo-pc.md
	Rscript -e "library(markdown);\
		opts <- setdiff(markdownHTMLOptions(TRUE), 'base_64_images');\
		markdownToHTML('$<', '$@', options=opts)"

## cleaning simulation results
clean-outputs:
	rm -rf output/sim-res
	rm -rf output/figs
	rm -rf output/data

## install dependencies
deps:
	Rscript make/dependencies.R
