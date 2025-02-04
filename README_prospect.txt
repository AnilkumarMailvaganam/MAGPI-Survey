********************************************************************************
Title: README_prospect
Author: Sabine Bellstedt (sabine.bellstedt@uwa.edu.au)
Date of creation: 13/10/2022
Last modified: 16/09/2024
Purpose: This file describes the format and creation of the MAGPI_ProSpectCat_vXX.csv and MAGPI_ProSpectCat_SFH_vXX.csv catalogues
Dependencies: This product is based on the v2.2.1 (...cloud.../MAGPI/reducedcubes/MAGPI_v2.2/MAGPIXXXX.fits) of the reduced cubes and associated manually edited dilated and undilated segmentation maps (...cloud.../MAGPI/valueadded/analysed_v2.2.1/MAGPIXXXX/segmap/MAGPI1201_manual_segmap.fits, used to generate forced photometry over the u/g/r/i/Z/Y/J/H/Ks bands from GAMA imaging. These photometry are used to conduct ProSpect SED fitting.
Other info:
v0.0 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502, 1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529)
v0.1 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502, 1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1533, 1534, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2310)
v0.2 - includes fields 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1501, 1502, 1503, 1505, 1506, 1507, 1508, 1511, 1522, 1523, 1525, 1527, 1528, 1529, 1530, 1533, 1534, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2310, 1509, 1512, 1524, 1526, 1531, 2309, 2311)
v0.3 - expanded the column list to include the age values and formation epochs in the Property cat, and adding the photometry cat with the modelled ProSpect photometry  in u-Ks bands. Uncertainty computation from the chain has been updated, to more accurately reflect the 1-sigma sampling errors. 
v0.4 - now includes fields 1504, 1513, 1514, 1515, 1516, 1517, 1519, 1521, 1532, 1535, 2312
********************************************************************************

Relevant concepts:
-----------------
The stellar population templates by Bruzual & Charlot (2003) are used in this work. In order to attenuate the stellar light by dust, the dust attenuayion Law by Charlot & Fall (2000) is applied. Finally, the dust emission is conducted according to the Dale et al. (2014) prescription. In this analysis, AGN is not accounted for.  
As described in Bellstedt et al. (2020b), the star formation history is parametrised using the massfunc_snorm_trunc parametrisation, which is a skewed Normal distribution, truncated in the early Universe to ensure that the SFR of each galaxy is 0 at the maximum age of the Universe, assumed for this implementation to be 13.38 Gyr. The mathematical form of this SFH is described in equations 1-5 of Bellstedt et al. (2020b). We fit the following free parameters in this implementation:

mSFR - fitted in log space in the range [-3, 4]
mpeak - fitted in linear space in the range [-(2+t_lb), 13.38-t_lb] in Gyr, where t_lb is the lookback time of the object
mperiod - fitted in log space in the range [log10(0.3), 2]
mskew - fitted in linear space in the range [-0.5, 1]

The truncation period is fixed to 2Gyr. 

** NB ** 
When regenerating the SFH based on the best-fitting parameters, note that the history will always be computed relative to the epoch of observation. Therefore, to convert the SFH to lookback time, always add on the lookback time of the galaxy itself (given in the SFH catalogue as the agemax parameter). 
********

We use the Zfunc_massmap_lin parametrisation to describe the gas-phase metallicity history of each galaxy. The shape of the metallicity evolution is determined by the shape of the stellar mass evolution for each object, with the final target metallicity being a free parameter (Zfinal). This parameter is fitted in logarithmic space, in the range [-4, -1.3]. These limits are set by the limits of the Bruzual & Charlot (2003) stellar population templates that are used for the SED fitting.

Data product description:
-------------------------
This catalogue includes outputs for all galaxies in the included fields with secure redshifts (QOP = 3 or 4), based on u-Ks MAGPI-matched photometry from KiDS and VIKING (as used for the GAMA survey).

Usage notes:
-----------
For each of the parameters presented in this catalogue (Stellar mass and SFR), both the best-fitting and the median parameters are provided, based on the MCMC chain for each galaxy. Furthermore, _16 and _84 values are provided from the chains to provide an estimate of the uncertainty for each parameter. 
Note that the usage case may prescribe whether the best-fitting or the median parameter is preferred, however it is expected that the median parameter will be more useful for most purposes. If SFH parameters are to be used, then it is likely that the best-fitting values will be more useful, as they will directly relate to the best-fitting SFH. 

Location of files:
------------------
The property catalogue is located at MAGPI/valueadded/ProSpect_v2.2.1/MAGPI_ProSpectCat_vXX.csv
The SFH catalogue is located at MAGPI/valueadded/ProSpect_v2.2.1/MAGPI_ProSpectCat_SFH_vXX.csv
The modelled photometry catalogue is located at MAGPI/valueadded/ProSpect_v2.2.1/MAGPI_ProSpectCat_ModelledPhotometry_vXX.csv

Relevant formatting information:
-------------------------------
Columns for the property cat:
MAGPIID					ID of galaxy
z						redshift	
StellarMass_bestfit		Stellar Mass from best-fitting ProSpect step
StellarMass_median		Stellar Mass from median of ProSpect chain
StellarMass_16			Stellar mass at 16th percentile from  ProSpect chain
StellarMass_84			Stellar mass at 84th percentile from  ProSpect chain
SFRburst_bestfit		SFR over the last 100Myr from best-fitting ProSpect step		
SFRburst_median			SFR over the last 100Myr from median of ProSpect chain	
SFRburst_16				SFR over the last 100Myr at 16th percentile from  ProSpect chain
SFRburst_84				SFR over the last 100Myr at 84th percentile from  ProSpect chain
Age_10_bestfit			Age (yr) at which 10% of the stars formed from best-fitting ProSpect step		
Age_10_median			Age (yr) at which 10% of the stars formed from median of ProSpect chain	
Age_10_16				Age (yr) at which 10% of the stars formed at 16th percentile from  ProSpect chain
Age_10_84				Age (yr) at which 10% of the stars formed at 84th percentile from  ProSpect chain
Age_20_bestfit			Age (yr) at which 20% of the stars formed from best-fitting ProSpect step		
Age_20_median			Age (yr) at which 20% of the stars formed from median of ProSpect chain	
Age_20_16				Age (yr) at which 20% of the stars formed at 16th percentile from  ProSpect chain
Age_20_84				Age (yr) at which 20% of the stars formed at 84th percentile from  ProSpect chain
Age_50_bestfit			Age (yr) at which 50% of the stars formed (half-mass age) from best-fitting ProSpect step		
Age_50_median			Age (yr) at which 50% of the stars formed (half-mass age) from median of ProSpect chain	
Age_50_16				Age (yr) at which 50% of the stars formed (half-mass age) at 16th percentile from  ProSpect chain
Age_50_84				Age (yr) at which 50% of the stars formed (half-mass age) at 84th percentile from  ProSpect chain
Age_80_bestfit			Age (yr) at which 80% of the stars formed from best-fitting ProSpect step		
Age_80_median			Age (yr) at which 80% of the stars formed from median of ProSpect chain	
Age_80_16				Age (yr) at which 80% of the stars formed at 16th percentile from  ProSpect chain
Age_80_84				Age (yr) at which 80% of the stars formed at 84th percentile from  ProSpect chain
Age_90_bestfit			Age (yr) at which 90% of the stars formed from best-fitting ProSpect step		
Age_90_median			Age (yr) at which 90% of the stars formed from median of ProSpect chain	
Age_90_16				Age (yr) at which 90% of the stars formed at 16th percentile from  ProSpect chain
Age_90_84				Age (yr) at which 90% of the stars formed at 84th percentile from  ProSpect chain
AgeDelta_10_90_bestfit	Epoch (yr) between 10% and 90% of stars forming from best-fitting ProSpect step		
AgeDelta_10_90_median	Epoch (yr) between 10% and 90% of stars forming from median of ProSpect chain	
AgeDelta_10_90_16		Epoch (yr) between 10% and 90% of stars forming at 16th percentile from  ProSpect chain
AgeDelta_10_90_84		Epoch (yr) between 10% and 90% of stars forming at 84th percentile from  ProSpect chain
AgeDelta_20_80_bestfit	Epoch (yr) between 10% and 90% of stars forming from best-fitting ProSpect step		
AgeDelta_20_80_median	Epoch (yr) between 10% and 90% of stars forming from median of ProSpect chain	
AgeDelta_20_80_16		Epoch (yr) between 10% and 90% of stars forming at 16th percentile from  ProSpect chain
AgeDelta_20_80_84		Epoch (yr) between 10% and 90% of stars forming at 84th percentile from  ProSpect chain
MeanAge_bestfit			Mean age of the galaxy from best-fitting ProSpect step		
MeanAge_median			Mean age of the galaxy from median of ProSpect chain	
MeanAge_16				Mean age of the galaxy at 16th percentile from  ProSpect chain
MeanAge_84				Mean age of the galaxy at 84th percentile from  ProSpect chain

Columns for the SFH cat:
MAGPIID					ID of galaxy
agemax 					Lookback time of the galaxy (in yrs)
mSFR_bestfit			peak SFR
mSFR_median				peak SFR (median of ProSpect chain)
mSFR_16					peak SFR (16th percentile from  ProSpect chain)
mSFR_84					peak SFR (84th percentile from  ProSpect chain)
mpeak_bestfit			age of peak SFR
mpeak_median			age of peak SFR (median of ProSpect chain)
mpeak_16				age of peak SFR (16th percentile from  ProSpect chain)
mpeak_84				age of peak SFR (84th percentile from  ProSpect chain)
mperiod_bestfit			period of skewed Normal SFH
mperiod_median			period of skewed Normal SFH (median of ProSpect chain)
mperiod_16				period of skewed Normal SFH (16th percentile from  ProSpect chain)
mperiod_84				period of skewed Normal SFH (84th percentile from  ProSpect chain)
mskew_bestfit			skewness of skewed Normal SFH
mskew_median			skewness of skewed Normal SFH (median of ProSpect chain)
mskew_16				skewness of skewed Normal SFH (16th percentile from  ProSpect chain)
mskew_84				skewness of skewed Normal SFH (84th percentile from  ProSpect chain)

Columns for the modelled photometry cat:
flux.u_VST_bestfit		ProSpect-modelled u-band, in Jy (best fit)
flux.u_VST_median		ProSpect-modelled u-band, in Jy (median of ProSpect chain)
flux.u_VST_16			ProSpect-modelled u-band, in Jy (16th percentile from  ProSpect chain)
flux.u_VST_84			ProSpect-modelled u-band, in Jy (84th percentile from  ProSpect chain)
flux.g_VST_bestfit		ProSpect-modelled g-band, in Jy (best fit)
flux.g_VST_median		ProSpect-modelled g-band, in Jy (median of ProSpect chain)
flux.g_VST_16			ProSpect-modelled g-band, in Jy (16th percentile from  ProSpect chain)
flux.g_VST_84			ProSpect-modelled g-band, in Jy (84th percentile from  ProSpect chain)
flux.r_VST_bestfit		ProSpect-modelled r-band, in Jy (best fit)
flux.r_VST_median		ProSpect-modelled r-band, in Jy (median of ProSpect chain)
flux.r_VST_16			ProSpect-modelled r-band, in Jy (16th percentile from  ProSpect chain)
flux.r_VST_84			ProSpect-modelled r-band, in Jy (84th percentile from  ProSpect chain)
flux.i_VST_bestfit		ProSpect-modelled i-band, in Jy (best fit)
flux.i_VST_median		ProSpect-modelled i-band, in Jy (median of ProSpect chain)
flux.i_VST_16			ProSpect-modelled i-band, in Jy (16th percentile from  ProSpect chain)
flux.i_VST_84			ProSpect-modelled i-band, in Jy (84th percentile from  ProSpect chain)
flux.Z_VISTA_bestfit	ProSpect-modelled Z-band, in Jy (best fit)
flux.Z_VISTA_median		ProSpect-modelled Z-band, in Jy (median of ProSpect chain)
flux.Z_VISTA_16			ProSpect-modelled Z-band, in Jy (16th percentile from  ProSpect chain)
flux.Z_VISTA_84			ProSpect-modelled Z-band, in Jy (84th percentile from  ProSpect chain)
flux.Y_VISTA_bestfit	ProSpect-modelled Y-band, in Jy (best fit)
flux.Y_VISTA_median		ProSpect-modelled Y-band, in Jy (median of ProSpect chain)
flux.Y_VISTA_16			ProSpect-modelled Y-band, in Jy (16th percentile from  ProSpect chain)
flux.Y_VISTA_84			ProSpect-modelled Y-band, in Jy (84th percentile from  ProSpect chain)
flux.J_VISTA_bestfit	ProSpect-modelled J-band, in Jy (best fit)
flux.J_VISTA_median		ProSpect-modelled J-band, in Jy (median of ProSpect chain)
flux.J_VISTA_16			ProSpect-modelled J-band, in Jy (16th percentile from  ProSpect chain)
flux.J_VISTA_84			ProSpect-modelled J-band, in Jy (84th percentile from  ProSpect chain)
flux.H_VISTA_bestfit	ProSpect-modelled H-band, in Jy (best fit)
flux.H_VISTA_median		ProSpect-modelled H-band, in Jy (median of ProSpect chain)
flux.H_VISTA_16			ProSpect-modelled H-band, in Jy (16th percentile from  ProSpect chain)
flux.H_VISTA_84			ProSpect-modelled H-band, in Jy (84th percentile from  ProSpect chain)
flux.K_VISTA_bestfit	ProSpect-modelled Ks-band, in Jy (best fit)
flux.K_VISTA_median		ProSpect-modelled Ks-band, in Jy (median of ProSpect chain)
flux.K_VISTA_16			ProSpect-modelled Ks-band, in Jy (16th percentile from  ProSpect chain)
flux.K_VISTA_84			ProSpect-modelled Ks-band, in Jy (84th percentile from  ProSpect chain)