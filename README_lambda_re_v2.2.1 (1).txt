*********************************************************************************
Title: README_lambda_Re_v2.2.1
Author: [Caro Derkenne](caro.derkenne@hdr.mq.edu.au)
Date of creation: 23/02/2023
Last modified: 11/01/2024
Purpose: This file describes basic details of how the lambda Re data
    catalogue was created. Please contact the author for additional details or
    if you notice anything awry.

Dependencies:
    [minicubes v2.2.1 (requires login)](https://datacentral.org.au/teamdata/MAGPI/valueadded/analysed_v2.2.1/README_minicubes)
    [profound catalogues (requires login)](https://datacentral.org.au/teamdata/MAGPI/valueadded/analysed_v2.2.1/README_profoundsegmap)
    [mockimages v2.2.1 (requires login)](https://datacentral.org.au/teamdata/MAGPI/valueadded/analysed_v2.2.1/README_mockimages)
    [Galfit catalogue v.2.20 (requires login)](https://datacentral.org.au/teamdata/MAGPI/valueadded/Galfit_v2.2.0/MAGPI_GalfitSersicCat_v0.1.csv)
    [Stellar kinematics v.2.2.1 (requires login)] (https://datacentral.org.au/teamdata/MAGPI/valueadded/StellarKinmeatics_v2.2.1/README_stellar_kinematics)

Other info: 
There was a correction to how the elliptical radii were handled in the code, leading to updated values. The most recent catalogue version is "lambda_Re_seeing_corrected_v1.csv"
*********************************************************************************


***Relevant concepts:***

Object naming follows the minicubes, as per the minicubes README:

Unique 10-digit MAGPI IDs are assigned as a concatenation of the 4 digits FieldID and the 3+3 digits (Xmax,Ymax) position of the brightest pixel in the white-light image. First 2 digits of the FieldID encode the GAMA field number (12=G12, 15=G15, 23=G23), and 3rd and 4th digits encode MAGPI field ID. These IDs may change slightly between releases as the brightest pixel may vary depending on the reduction.

***Data product description:*** 

Lambda Re is defined as (Sigma F_i R_i |V_i|)/(Sigma F_i R_i srqt(V_i^2 + sigma_i^2)), where "Sigma" is a sum from spaxel = 1 to spaxel = N within some aperture, R_i is the galactocentric radius, and F_i is the flux in the ith spaxel. V_i is the observed line-of-sight velocity, and sigma_i is the line-of-sight velocity dispersion. 

V/sigma is the ratio of velocity to velocity dispersion, measured within some aperture, defined as sqrt((Sigma F_i V_i^2)/(Sigma F_i sigma_i^2)). 

For this catalogue, both measurements are taken in a 1 Re elliptical aperture, where the radius of the ith pixel is equal to the semi-major axis of the ellipse upon which it lies (following van de Sande 2017).

These data products have been corrected for seeing effects.

The process for each measurement is as follows: 

1) An Multi-Gaussian Expansion (MGE) is made of the PSF, given in the minicubes fits file PSF extension, using the mgefit software (https://pypi.org/project/mgefit/). The axial ratios of the PSF MGE components are forced to be equal (ie., the PSF must be circular).

2) Each MAGPI object with v2.2.1 Stellar Kinematics has an MGE model constructed. A 20 Re thumbnail is cut from the r-band mock image, where an initial estimate of Re is taken from the ProFound catalogue. Segmentation maps are used to mask out other sources. The model is analytically convolved with the PSF and then compared to data for optimisation. The resulting model is therefore PSF deconvolved. The elliptical effective radius (the semi-major axis of the 1 Re ellipse) and the ellipticity of the 1Re isophote are measured from the best-fitting MGE model. 

3) The elliptical effective radius and ellipticity of the 1Re isophote are used to construct an aperture on the stellar kinematic velocity and velocity dispersion fields, and the flux image constructed from the minicube for the relevant galaxy.

4) Once all potential pixels for making the measurement are located, they are filtered according to the stellar kinematic quality cut recommendation: The velocity and velocity dispersion values must be finite and the dispersion error must satisfy: data['ERR_SIGMA'] < 25. + 0.1 * data['SIGMA']. Only these "valid" pixels are used for calculations. 

5) 2000 realisations of the velocity and velocity dispersion fields were created by drawing from a random normal distribution with sigma-width on each pixel equal to the quoted error for that pixel given in the stellar kinematics fits files. The measurements of v/sigma and Lambda Re are repeated, and the standard deviation of the resulting distribution divided by sqrt(2) is given as the Monte-Carlo error. 


5) The resulting v/sigma and lambda Re measurements are corrected for seeing effects using the Galfit estimate of Sersic index, the resolution of the object (sigma_psf/Re), and the ellipticity, using Kate Harborne's kinematic corrections code: https://github.com/kateharborne/kinematic_corrections.

Galaxies that failed the MGE step do not have an associated lambda Re measurement. This can happen if the mgefit software could not converge to a solution, or if the solution contained NaNs, and occurs only for particularly ill-resolved galaxies. Only galaxies that had zero valid pixels failed the Lambda Re and v/sigma measurement step. Galaxies that did not have Sersic indices failed the seeing correction step.

A flux image with the elliptical aperture marked in red, and the velocity and velocity dispersion pixels used for the calculation, is included for every object ("lambda_r_field_MGE_XXXXXXXXXX.pdf").

The data products consists of an array  (13 columns x 357 rows). Below is a column-ordered description of the array: 

1)  ["MAGPIID"]: MAGPIID following the minicube/standard MAGPI convention. The first four digits give the field identifier, the remaining 6 give the x and y position in the field based on the brightest pixel.
2)  ["obs_lambda_re"] : Lambda Re measured within an elliptical 1 Re aperture, not seeing corrected. 
3)  ["corr_lambda_re"] : Seeing corrected Lambda Re measured with an elliptical 1 Re aperture. 
4)  ["err_obs_lambda_re"] : Monte-Carlo estimate of uncertainty on the Lambda Re measurement.
5)  ["obs_v_on_sigma_obs"]: V/sigma measurement within the elliptical 1 Re aperture, not seeing corrected.
6)  ["corr_v_on_sigma_obs"]: Seeing corrected V/sigma measurement within the elliptical 1 Re aperture. 
7) ["err_v_on"sigma"]: Monte-carlo estimate of uncertainty on the v/sigma measurement. 
8)  ["ellip"]: Ellipticity of the 1 Re isophote of the best-fitting MGE model for the object, PSF deconvolved. Ellipticity  = 1 - b/a, where "b" and "a" are the semi-minor and semi-major axes, respectively.
9)  ["sersic_n"]: Copy of the Sersic n from the Galfit catalogue, required for the seeing correction. 
10)  ["Re_arcsec"]: Elliptical semi-major axis of the 1 Re isophote from the best-fitting MGE model, PSF deconvolved, in arcseconds. 
11) ["sigma_psf_over_re"]: The PSF sigma (= FWHM/~2.355) of the MGE PSF model divided by the semi-major axis of the 1Re isophote of the MGE galaxy model.
12) ["fill_fraction"]: The pixel-area of the 1 Re elliptical aperture divided by the number of **valid** pixels used in the calculations. 
13) ["npix"]: The number of **valid** pixels used in the calculations. 


***Usage notes:***

NaN values have been indicated with -99.99 in the data frame. This indicates where no measurement is available. 

Included data values in the catalogue that are useful as quality indicators are: the sigma of the PSF as a fraction of the circularised effective radius ("sigma_psf_over_re"), the number of pixels used to calculate the lambda Re and v/sigma value ("npix"), and the fill fraction ("fill_fraction"). 

The catalogue includes lambda Re values measured from kinematics with very few valid pixels, or low fill fractions, or that are poorly resolved, (or all three), and these should be filtered from any science sample by using cuts on the above properties. Eg, npix > 50, sigma_psf_over_re < 0.4, and fill_fraction > 0.85. 

The csv file can be read in via (e.g.): 

lambda_re_cat = pandas.read_csv("lambda_Re_seeing_corrected_v1.csv")
keys          = lambda_re_cat.keys() # see all included information
lambda_re     = lambda_re_cat["corr_lambda_re"].values # extracting array for a particular key

***Location of files:*** 

All relevant files (including a copy of this README can be found on the cloud at: 
https://cloud.datacentral.org.au/cloud/MAGPI/valueadded/lambda_Re_v2.2.1

***Relevant formatting information:***

(13 columns x 357 rows) .csv file, columns as given under "***Data product description:*** ".
