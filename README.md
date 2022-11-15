# oceanwaves

Author: Andrew Reed

## Overview
---

```oceanwaves``` is a package to calculate the wave spectra and summary statistics for surface waves based on measurements from the SeaBird SBE26 pressure sensor operating in waveburst mode. The package includes the following functions:
* ```wavenumber```: Calculate the wavenumber from the wave frequency and water depth using a polynomial approximation following Hunt (1979) as implemented by Voulgaris (1992)
* ```detrendHeight```: Remove very-low frequency (such as tides) trend from a time series to derive a deviation from mean sea-surface-elevation estimate
* ```pressureCorrection```: Correct for depth attenuation of a water surface elevation pressure signal.
* ```waveSpectra```: Calculate wave statistics from the SSH data, corrected for pressure attentuation. This function calculate the wae statistics such as significant wave height, peak period, the average period using two different methods, and the estimates of spectral width from the  sea surface height data measured by a subsurface pressure sensor. 
    
    
    
methods for deriving the estimated sea surface height from pressure, detrending time series of SSH, correcting the depth attenuation of the pressure signal, and calculating the wave spectra and key wave statistics. 

This package is based on the R-language package ```oceanwaves```, developed by Luke Miller, to derive ocean wave height statistics from Ocean Wave Height Loggers, as well as the MatLab ```waves``` package developed by Urs Neumeier. Functions are further modified for the SeaBird SBE26 Wave Burst measurements following the Sea-Bird SBE26 User's Manual.

---
## TODO

There are a number of different methods for calculating wave statistics from instruments not explicitly designed for wave measurements. This package will eventually include methods for deriving wave statistics from:
1. Three Axis Motion Sensor (MOPAK)
2. Direct Covariance Flux (FDCHP)

---
#### References
1. Sea-Bird Scientific (2015). _SBE 26plus Seagauge Wave & Tide Record: User Manual_ (Version 019). Retrieved from [https://www.seabird.com/cms-portals/seabird_com/cms/documents/discontinued-products/manual-sbe26plus.pdf](https://www.seabird.com/cms-portals/seabird_com/cms/documents/discontinued-products/manual-sbe26plus.pdf)
2. Miller, L. (2021, June 2). millerlp/oceanwaves: v0.2.0. [Software]. Zenodo. [http://doi.org/10.5281/zenodo.4925489](http://doi.org/10.5281/zenodo.4925489)
3. Neumeier, U. (2006, Jan. 16). Processing of wave data from pressure sensor. [Software]. [http://neumeier.perso.ch/matlab/waves.html](http://neumeier.perso.ch/matlab/waves.html)

```python

```
