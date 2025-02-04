***************************************************************************************** 
Title: README_galfit 
Author: Trevor Mendel (trevor.mendel@anu.edu.au) 
Date of creation: 08/11/2023 
Last modified: 08/11/2023 
Purpose: This file describes the format and creation of the MAGPI_GalfitSersicCat_vXX.csv
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
This file uses galfit (Peng et al. 2002) to model the 2D surface brightness distrubution
assuming a single Sersic component. The fitting is done in the MAGPI i-band (i.e. on a
mock i-band image computed from the MAGPI data cubes). We model each galaxy as well as any
neighbouring galaxies within 5 magnitudes and a projected separation R_p < 5 (R_galaxy +
R_neighbour). In the case that a neighbour satisfying the R_p cut falls outside of the
postage stamp used for fitting, its central position is held fixed, while all other
parameters (size, sersic index, b/a, pa, magnitude) are left free. All other sources in
the image are masked according to the ProFound segmentation maps. Prior to fitting we
compute a local background and remove it from the image; the background is held fixed
during fitting. When performing the fits we use the reconstructed i-band PSF included in
the MAGPI data cubes. 

Data product description:
-------------------------
This catalogue includes outputs for all galaxies in the included fields with secure
redshifts (QOP = 3 or 4). The ordering is such that these files match the stellar mass
measurements in the MAGPI_ProSpectCat_vXX.csv files of the same version. 

Usage notes:
-----------
The parameters and corresponding uncertainties are based strictly on the GalFit output.
Note also that the reported chi2 values are for the entire scene being modelled
(source+neighbours+background) and so arencertainty for each parameter.  Note that in the
case of complex morphology, (spiral arms, multiple photometric components, etc.) the
adopted single Sersic profile may not provide a good description of the data.  It is
therefore recommended to verify the fidelity of the model (through, say, visual
inspection) in the case of complex systems.

Sersic indices are bounded betwee 0.5 and 8, and galaxies with values at these boundaries
should be approached with caution.

Location of files:
------------------
The catalogue is located at MAGPI/valueadded/Galfit_v2.2.0/MAGPI_GalFitSersicCat_vXX.csv

Relevant formatting information:
-------------------------------
MAGPIID					ID of galaxy
z						redshift	
re		                semi-major axis half-light size [arcsec]
re_err		            1-sigma uncertainty on semi-major axis size [arcsec]
n                       Sersic index 
n_err                   1-sigma uncertainty on Sersic index
q                       axis ratio (b/a)
q_err                   1-sigma uncertainty on axis ratio
mag                     i-band magnitude
mag_err                 1-sigma uncertainty on i-band magnitude
pa                      position angle [degrees]
pa_err                  1-sigma uncertainty on position angle [degrees]
chi2                    reduced chi2 of the fit (note that this is computed over the full image)
