# Tutorial 12: Performing FCS and FCCS

Fluorescence correlation spectroscopy (FCS) and fluorescence cross-correlation spectroscopy (FCCS) is only recommended for images generated with photophysical conversions.
Please download `img0-20000_fs530_T1_I_0.1_0.2.tiff` from University of Alberta Dataverse  DOI: https://doi.org/10.7939/DVN/F3JKZH version 3.0

## 1. FCS
```bash
siliscipy prop --method fcs --file img0-20000_fs530_T1_I_0.1_0.2.tiff --paramfile parameters.dat  
```
The parameter `fcs_tmax` is important.


## 2. FCCS
```bash
siliscopy prop --method fccs --file img0-20000_fs530_T1_I_0.1_0.2.tiff --paramfile parameters.dat
```
The parameter `fcs_tmax` and `fccs_cols` is important.

## 3. Sample plot
```bash
python plot.py
```
