# Modeling conformity of SAR

Abstract

## Workflow

![](SAR_workflow.png)

## Cite

```bibtex
@article{steibl2025sar,
  title={Rainfall increases conformity and strength of species area relationships in atolls},
  author={Steibl, Sebastian and Steiger, Simon and Valente, Luís and Russell, James C},
  journal={submitted},
  year={2025}
  doi={tbd}
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
