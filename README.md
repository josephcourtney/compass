COMPASS - Comparative, Objective Measurement of Protein Architectures by Scoring Shifts =============

This software package implements the COMPASS algorithm for the identification of the best structures from a set of models by numerical comparison with unassigned NMR data.

The paper describing this method:
http://dx.doi.org/10.1016/j.str.2015.07.019

[Experimental Protein Structure Verification by Scoring with a Single, Unassigned NMR Spectrum](http://dx.doi.org/10.1016/j.str.2015.07.019)

by

Joseph M. Courtney, Qing Ye, Anna E. Nesbitt, Ming Tang, Marcus D. Tuttle, Eric D. Watt, Kristin M. Nuzzio, Lindsay J. Sperling, Gemma Comellas, Joseph R. Peterson, James H. Morrissey, Chad M. Rienstra

## How to install

This package requires:
* numpy version >= 0.15.0
* [Rosetta](https://www.rosettacommons.org/software) version >= 3.4
* [MODELLER](https://salilab.org/modeller/)

>git clone https://github.com/josephcourtney/compass.git
>cd compass
>python setup.py build
>sudo python setup.py install
