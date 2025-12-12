[![DOI](https://zenodo.org/badge/940970214.svg)](https://doi.org/10.5281/zenodo.15129174)

# Modeling conformity of SAR

Species-area relationships are foundational to ecology and biogeography, yet small islands often deviate from predictable scaling. This high stochasticity in small-island species-area relationships is typically disregarded as unexplained noise or attributed to idiosyncratic extinction rates. Here, we introduce a statistical framework that explicitly incorporates the stochasticity in species-area relationships as a model parameter. Using a global island plant dataset for atolls (378 islands across 19 atolls) – prototypical examples for small-island dynamics – we reveal that rainfall, not cyclone disturbances, predictably modulates stochasticity in species-area relationships. Atolls receiving higher rainfall exhibit increasingly predictable species-area relationships, i.e., lower stochasticity around the species-area curve. Our findings challenge the widely held view that stochasticity in small-island species-area relationships is due to disturbance regimes, instead favouring resource limitation as a cause. Explicitly modelling stochasticity in species-area studies offers a promising novel tool for investigating the environmental factors influencing the conformity to the law of species-area scaling. 

## Workflow

![](SAR_workflow.png)

We compiled island-level biodiversity data for 19 atolls from literature reports, along with cumulative annual rainfall measurements and the frequency of tropical cyclone encounters per atoll.
We then computed the species-area relationship (SAR) following Arrhenius' Power Law, which relates island species richness S to island area A through a log-log linear relationship with intercept c and slope z.
Then, we determined the effect of rainfall and cyclone disturbance frequency on the slope parameter and the residual error in the species-area relationship, using Bayesian Distributional Modelling (i.e., Bayesian Heteroscedastic Modelling).
We tested the hypotheses that higher rainfall amounts reduce the residual error (i.e. the ‘noise’) in species-area relationship, and increase the slope of the relationship.

## Cite

```bibtex
@article{steibl2025sar,
  title={Rainfall increases conformity and strength of species area relationships},
  author={Steibl, Sebastian and Steiger, Simon and Valente, Luís and Russell, James C},
  journal={Ecography},
  year={2025}
  doi={[tbd](https://doi.org/10.1002/ecog.08159)}
}
```

## Installation

### R script

To run the core data analysis:

1. Install the R programming language ([Windows](https://cran.r-project.org/bin/windows/), [MacOS](https://cran.r-project.org/bin/macosx/), [Linux](https://cran.r-project.org/bin/linux/)) (Version 4.4.0 or later)
2. Install the `renv` package by running `install.packages("renv")`
3. Open an R console in the project folder `SAR-conformity`
4. Run `renv::activate()` and restart the R session
5. Run `renv::restore()` to download all dependencies from the lock file

### Julia script

To run the simulation:

1. Open the file `simulation.html` 
2. Follow the installation tutorial by clicking "Edit or run this notebook" in the top right of the window
