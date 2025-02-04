***************************************************************************************** 
Title: README_profit 
Author: Trevor Mendel (trevor.mendel@anu.edu.au) 
Date of creation: 08/11/2023 
Last modified: 08/11/2023 
Purpose: This file describes the format and creation of the MAGPI_ProfitSersicCat_vXX.csv
catalogues.  
Dependencies: This product is based on v2.2.1 of the reduced cubes, the  manually edited
dilated and undilated segmentation maps, as well as the source catalogues for each field.
Other info: 
v0.0 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502,
1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529) 
v0.1 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502,
1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1533, 1534,
2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2310) 
v0.2 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502,
1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1533, 1534,
2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2310, 1509, 1512, 1524, 1526, 1531, 2309,
2311) 
*****************************************************************************************

Relevant concepts:
-----------------
This file uses Profit (Robotham et al 2017) to model the 2D surface brightness distrubution
assuming a single Sersic component. The fitting is done in the MAGPI i-band (i.e. on a
mock i-band image computed from the MAGPI data cubes) using the scipy.minimize function
and the 'Nelder-Mead' method. The galaxy images fit here are identical to those used in the
creation of the MAGPI_GalfitSersicCat_vXX.csv files. In the case of Profit, we fit only
the galaxy of interest using the dilated ProFound segmentation map. As with the Galfit
fits, the background is held fixed during fitting.  

Data product description:
-------------------------
This catalogue includes outputs for all galaxies in the included fields with secure
redshifts (QOP = 3 or 4). The ordering is such that these files match the stellar mass
measurements in the MAGPI_ProSpectCat_vXX.csv files of the same version. 

Usage notes:
-----------
The parameters are based strictly on the scipy.minimize output.  In the case of complex
morphology (spiral arms, multiple photometric components, etc.) the adopted single Sersic
profile may not provide a good description of the data.  It is therefore recommended to
verify the fidelity of the model (through, say, visual inspection) in the case of complex
systems.

Sersic indices are bounded betwee 0.5 and 8, and galaxies with values at these boundaries
should be approached with caution.

Location of files:
------------------
The catalogue is located at MAGPI/valueadded/Profit_v2.2.0/MAGPI_ProfitSersicCat_vXX.csv

Relevant formatting information:
-------------------------------
MAGPIID					ID of galaxy
z						redshift	
re		                semi-major axis half-light size [arcsec]
n                       Sersic index 
q                       axis ratio (b/a)
mag                     i-band magnitude
pa                      position angle [degrees]
chi2                    reduced chi2 of the fit
