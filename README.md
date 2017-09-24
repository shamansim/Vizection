# Vizection
 Shiny application for single cell transcriptomics studies data visualisation / exploration, quality control and analysis

Developped during my 6 month internship at RIKEN Center for the Life Science Technologies (CLST) in the Division of Genomics Technologies (DGT) in the Genomics Miniaturization Technology Unit. I was supervised by [Charles Plessy](https://github.com/charles-plessy).

## How to install and use

Install from GitHub with:

    devtools::install_github("shamansim/Vizection", upgrade_dependencies = FALSE) # Note the capital 'V'.

Load the library:

    library(vizection) # Note the lowercase 'v'
    
If you do not have any dataset you can try vizection on some example
data from Bioconductor, for instance the SummarizedExperiment object
provided by the [airway](https://bioconductor.org/packages/airway) package.

    data("airway", package = "airway")
    airway$group <- airway$dex

Start vizection:

    vizection(airway, local = FALSE) # local = TRUE implies host="0.0.0.0" 

Documentation is under way, some features are still not working properly, please be patient :)
