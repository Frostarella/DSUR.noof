# DSUR.noof

## Description
DSUR.noof is a package for R, including (almost) all functions of the DSUR package mentioned in Andy Field's 'Discovering Statistics Using R'. This is by no means an official version of the DSUR package -- just my private collection which I wanted to share.

## Development
This package is a work in progress. It is usable, but the functions are not fully documented yet. Please refer to the book 'Discovering Statistics Using R' when in need of further explanation for a function and be sure to check back regularily for updated versions.
Also, if you find any mistakes of any sort please let me know.

## Installation
The package can be easily installed by using the devtools to install it directly from GitHub. 
Just copy the following commands into your R console and execute them.
```R
package_install("devtools")
library(devtools)
install_github("Frostarella/DSUR.noof")
```

## Credits
* Credit goes to Andy Field, Jeremy Miles, and Zoë Field (Authors of 'Discovering Statistics Using R'), as well as G. Jay Kerns and Antonio Trujillo-Ortiz et al., for the functions in this package.
* This package was made by following several online resources on how to write R packages:
    * Hilary Parker's blog post ([read it here](https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/))
    * Karl Broman's tutorial ([read it here](http://kbroman.org/pkg_primer/))
    * Hadley Wickham's book (only parts of it, though) ([read it here](http://r-pkgs.had.co.nz))
