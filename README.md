# MAGENTA

Individual-based simulation model of malaria epidemiology and genomics.

## Version History

22.02.2017  State pointer passing

16.02.2017  No mosquito/strain version checked 

06.12.2016  Package initialised

### What is this?

*MAGENTA* is an individual-based simulation model of malaria epidemiology and genomics.
*MAGENTA* extends the imperial malaria model by tracking the infection history of 
individuals. With this additional genetic characteristics of the parasite can be 
assessed.

***
> To view the tutorial please click [here](https://github.com/OJWatson/MAGENTA/blob/master/tutorials/MAGENTA_tutorial.md)

***

### Installing *MAGENTA*

To install the development version from github the package [*devtools*](https://github.com/hadley/devtools) is required.

In order to install devtools you need to make sure you have a working development environment:

1. **Windows**: Install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)**. For help on how to install **Rtools** please see the following [guide](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows).

2. **Mac**: Install Xcode from the Mac App Store.

3. **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

Once a working development environment is ready, then devtools can be installed from CRAN:

```r
install.packages("devtools")
library(devtools)
```
Once installed, the package can be installed and loaded using:

```r
devtools::install_github("OJWatson/MAGENTA")
library(MAGENTA)
```

***

#### Asking a question

For bug reports, feature requests, contributions, use github's [issue system.](https://github.com/OJWatson/MAGENTA/issues)