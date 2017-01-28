# Vizection
 Shiny application for single cell transcriptomics studies data visualisation / exploration, quality control and analysis

Developped during my 6 month internship at RIKEN Center for the Life Science Technologies (CLST) in the Division of Genomics Technologies (DGT) in the Genomics Miniaturization Technology Unit. I was supervised by [Charles Plessy](https://github.com/charles-plessy).

## How to install and use

Install from GitHub with:

    devtools::install_github("shamansim/Vizection", upgrade_dependencies = FALSE) # Note the capital 'V'.

Load the library:

    library(vizection) # Note the lowercase 'v'
    
If you do not have any dataset you can generate some:

```
data(iris)
genes <- iris[, 1:4]
libs <- as.data.frame(iris[, "Species"])
colnames(libs) <- c('group')
libs$samplename <- as.character(1:nrow(genes))
rownames(libs) <- libs$samplename 
libs$counts <- rnorm(n = nrow(genes), mean = 1000, sd = 200)
genes <- as.data.frame(t(genes))

```

Start vizection:

    vizection(genes = genes, libs = libs, local = FALSE) # local = TRUE implies host="0.0.0.0" 

Documentation is under way, some features are still not working properly, please be patient :)
