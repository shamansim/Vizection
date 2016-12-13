# Vizection
 Shiny application for single cell transcriptomics studies data visualisation / exploration and quality control

Developped during my 6 month internship at RIKEN CLST in the Division of Genomics Technologies (DGT) in the Genomics Miniaturization Technology Unit. I was supervised by [Charles Plessy](https://github.com/charles-plessy).

## How to install and use

Install from GitHub with:

    devtools::install_github("shamansim/Vizection", upgrade_dependencies = FALSE)) # Note the capital 'V'.

Load the library:

    library(vizection) # Note the lowercase 'v'

Start vizection

    vizection(genes = genes, libs = libs, local = TRUE) # local = TRUE implies host="0.0.0.0" #pas de majuscule

Documentation is under way, please be patient :)
