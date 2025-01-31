#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 16:57:21 2025

@author: 44185316
"""

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
import csv
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats

from matplotlib import cm as CM
from scipy.stats import spearmanr
import statistics

from matplotlib import pyplot as plt

from astropy.io import fits

import glob

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)



#%%
#MAGPHYS

obsfile = np.genfromtxt("/Users/44185316/Downloads/MAGPI_flux_gg_jy.dat", names=True, dtype=int, usecols=0)



magphys_dir = "/Users/44185316/Downloads/MAGPI_result_inichanged/"

Mstar_magpis = []
Mstar_ustd_magpis = []
Mstar_lstd_magpis = []
SFR_magpis = []
SFR_ustd_magpis = []
SFR_lstd_magpis = []

Mdust_magpis = []
Mdust_ustd_magpis=[]
Mdust_lstd_magpis=[]

chi2_magpis = []



name_mag=[]


tauVism_magpis=[]
tauV_magpis=[]
tauV_ustd_magpis=[]
tauV_lstd_magpis=[]

tauVism_ustd_magpis=[]
tauVism_lstd_magpis=[]

sSFR_magpis=[]
sSFR_ustd_magpis=[]
sSFR_lstd_magpis=[]

for ii in range(len(obsfile)):
    coord = glob.glob(magphys_dir+'%s.fit' % (obsfile[ii][0]))
    if coord == []:
        pass
    else:
        #print(coord)
        read_fit = open(coord[0], 'r')
        lines = read_fit.readlines()
        order = []
        for i, item in enumerate(lines):
            if "#" in item:
                order.append(i)
        
        name_mag.append(obsfile[ii][0])
        chi2 = np.loadtxt(coord[0], skiprows=8, max_rows=1)[-2]
        Mdust = 10**(np.loadtxt(coord[0], skiprows=order[39], max_rows=1)[2])
        Mdust_ustd = 10**(np.loadtxt(coord[0], skiprows=order[39], max_rows=1)[3]) - 10**(np.loadtxt(coord[0], skiprows=order[39], max_rows=1)[2])
        Mdust_lstd = 10**(np.loadtxt(coord[0], skiprows=order[39], max_rows=1)[2]) - 10**(np.loadtxt(coord[0], skiprows=order[39], max_rows=1)[1])
        Tw = np.loadtxt(coord[0], skiprows=order[27], max_rows=1)[2]
        Tw_ustd = np.loadtxt(coord[0], skiprows=order[27], max_rows=1)[3] - np.loadtxt(coord[0], skiprows=order[27], max_rows=1)[2]
        Tw_lstd = np.loadtxt(coord[0], skiprows=order[27], max_rows=1)[2] - np.loadtxt(coord[0], skiprows=order[27], max_rows=1)[1]
        Tc = np.loadtxt(coord[0], skiprows=order[25],max_rows=1)[2]
        Tc_ustd = np.loadtxt(coord[0], skiprows=order[25],max_rows=1)[3] - np.loadtxt(coord[0], skiprows=order[25],max_rows=1)[2]
        Tc_lstd = np.loadtxt(coord[0], skiprows=order[25],max_rows=1)[2] - np.loadtxt(coord[0], skiprows=order[25],max_rows=1)[1]
        tauV = np.loadtxt(coord[0], skiprows=order[17],max_rows=1)[2]
        tauV_ustd = np.loadtxt(coord[0], skiprows=order[17],max_rows=1)[3] - np.loadtxt(coord[0], skiprows=order[17],max_rows=1)[2]
        tauV_lstd = np.loadtxt(coord[0], skiprows=order[17],max_rows=1)[2] - np.loadtxt(coord[0], skiprows=order[17],max_rows=1)[1]
        tauVism = np.loadtxt(coord[0], skiprows=order[37],max_rows=1)[2]
        tauVism_ustd = np.loadtxt(coord[0], skiprows=order[37],max_rows=1)[3] - np.loadtxt(coord[0], skiprows=order[37],max_rows=1)[2]
        tauVism_lstd = np.loadtxt(coord[0], skiprows=order[37],max_rows=1)[2] - np.loadtxt(coord[0], skiprows=order[37],max_rows=1)[1]
        Mstar = 10**np.loadtxt(coord[0], skiprows=order[21],max_rows=1)[2]
        Mstar_ustd = 10**(np.loadtxt(coord[0], skiprows=order[21],max_rows=1)[3]) - 10**(np.loadtxt(coord[0], skiprows=order[21],max_rows=1)[2])
        Mstar_lstd = 10**(np.loadtxt(coord[0], skiprows=order[21],max_rows=1)[2]) - 10**(np.loadtxt(coord[0], skiprows=order[21],max_rows=1)[1])
        SFR = 10**np.loadtxt(coord[0], skiprows=order[41],max_rows=1)[2]
        SFR_ustd =10**np.loadtxt(coord[0], skiprows=order[41],max_rows=1)[3] - 10**np.loadtxt(coord[0], skiprows=order[41],max_rows=1)[2]
        SFR_lstd = 10**np.loadtxt(coord[0], skiprows=order[41],max_rows=1)[2] - 10**np.loadtxt(coord[0], skiprows=order[41],max_rows=1)[1]
        sSFR = np.loadtxt(coord[0], skiprows=order[19],max_rows=1)[2]
        sSFR_ustd =np.loadtxt(coord[0], skiprows=order[19],max_rows=1)[3] - np.loadtxt(coord[0], skiprows=order[19],max_rows=1)[2]
        sSFR_lstd = np.loadtxt(coord[0], skiprows=order[19],max_rows=1)[2] - np.loadtxt(coord[0], skiprows=order[19],max_rows=1)[1]

    

    
        Mstar_magpis.append(Mstar)
        Mstar_ustd_magpis.append(Mstar_ustd)
        Mstar_lstd_magpis.append(Mstar_lstd)


        Mdust_magpis.append(Mdust)
        Mdust_ustd_magpis.append(Mdust_ustd)
        Mdust_lstd_magpis.append(Mdust_lstd)    


    
        SFR_magpis.append(SFR)
        SFR_ustd_magpis.append(SFR_ustd)
        SFR_lstd_magpis.append(SFR_lstd)

        


        chi2_magpis.append(chi2)
        
        
        tauVism_magpis.append(tauVism)
        tauV_magpis.append(tauV)
        tauV_ustd_magpis.append(tauV_ustd)
        tauV_lstd_magpis.append(tauV_lstd)
        
        
        tauVism_ustd_magpis.append(tauVism_ustd)
        tauVism_lstd_magpis.append(tauVism_lstd)
        
        
        sSFR_magpis.append(sSFR)
        sSFR_ustd_magpis.append(sSFR_ustd)
        sSFR_lstd_magpis.append(SFR_lstd)

#%%


name_magpis = name_mag
print(name_magpis)


#%%
Mstar_magpis=np.array(Mstar_magpis)
SFR_magpis=np.array(SFR_magpis)
Mdust_magpis=np.array(Mdust_magpis)
sSFR_magpis=np.array(sSFR_magpis)
sSFR_ustd_magpis=np.array(sSFR_ustd_magpis)
sSFR_lstd_magpis=np.array(sSFR_lstd_magpis)






#%%



# Convert lists to numpy arrays
Mstar = np.array(Mstar_magpis)
Mstar_lerr = np.array(Mstar_lstd_magpis)
Mstar_uerr = np.array(Mstar_ustd_magpis)

Mdust = np.array(Mdust_magpis)
Mdust_lerr = np.array(Mdust_lstd_magpis)
Mdust_uerr = np.array(Mdust_ustd_magpis)

SFR = np.array(SFR_magpis)
SFR_lerr = np.array(SFR_lstd_magpis)
SFR_uerr = np.array(SFR_ustd_magpis)

# Transform to log space
log_Mstar = np.log10(Mstar)
log_Mstar_lerr = Mstar_lerr / (Mstar * np.log(10))
log_Mstar_uerr = Mstar_uerr / (Mstar * np.log(10))

log_Mdust = np.log10(Mdust)
log_Mdust_lerr = Mdust_lerr / (Mdust * np.log(10))
log_Mdust_uerr = Mdust_uerr / (Mdust * np.log(10))

log_SFR = np.log10(SFR)
log_SFR_lerr = SFR_lerr / (SFR * np.log(10))
log_SFR_uerr = SFR_uerr / (SFR * np.log(10))

# Calculate log ratios
MAGPI_log_dust_to_stellar_mass = log_Mdust - log_Mstar
MAGPI_dust_to_stellar_mass_err = np.sqrt(log_Mdust_lerr**2 + log_Mstar_lerr**2)

#both upper and lower
MAGPI_dust_to_stellar_mass_lerr = np.sqrt(log_Mdust_lerr**2 + log_Mstar_lerr**2)
MAGPI_dust_to_stellar_mass_uerr = np.sqrt(log_Mdust_uerr**2 + log_Mstar_uerr**2)


MAGPI_log_sSFR = log_SFR - log_Mstar
MAGPI_sSFR_err = np.sqrt(log_SFR_lerr**2 + log_Mstar_lerr**2)

#both upper and lower
MAGPI_sSFR_lerr = np.sqrt(log_SFR_lerr**2 + log_Mstar_lerr**2)
MAGPI_sSFR_uerr = np.sqrt(log_SFR_uerr**2 + log_Mstar_uerr**2)


#take average of the two
MAGPI_dust_to_stellar_mass_symmetric_error = (MAGPI_dust_to_stellar_mass_lerr + MAGPI_dust_to_stellar_mass_uerr) / 2

MAGPI_sSFR_symmetric_error = (MAGPI_sSFR_lerr + MAGPI_sSFR_uerr) / 2


#%%
namez=[]
z=[]
#MAGPIID,segID,xmax,ymax,RAmax,Decmax,z,QOP,mag_it,mag_rt,R50_it,R90_it,axrat_it,ang_it
#with open("E:\MAGPI_Emission_SEG\MAGPI_master_source_catalogue.csv") as csvfile:
#with open("/Users/44185316/Library/CloudStorage/OneDrive-MacquarieUniversity/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
#with open("E:/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
with open("/Users/44185316/Library/CloudStorage/OneDrive-MacquarieUniversity/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
    csvReader = csv.reader(csvfile, delimiter=',')
    for row in csvReader:
        namez.append(int(row[0])) 
        z.append(float(row[6]))
        
        
redshift_MAGPI=[]
for i in range(len(name_mag)):
    for j in range(len(namez)):
        if name_mag[i]==namez[j]:
            redshift_MAGPI.append(float(z[j]))

print(np.min(redshift_MAGPI))
print(np.max(redshift_MAGPI))





###################################################################################
###################################################################################
########      Emission dust-to-stellar mass ratio     ############################
###################################################################################
###################################################################################



# MAGPIID2=[]
# TOTAL_DUST_MASS_ALL=[]
# TOTAL_DUST_MASSerr_ALL=[]
# VALID_PIX_TOTAL_DUST=[]
# TOTAL_SFR_ALL=[]
# TOTAL_SFRerr=[]
# TOTAL_Attenuation=[]

# #   	MAGPIID	VALID_PIX_DUST	AREA_KPC_10_DUST	AREA_KPC_50_DUST	AREA_KPC_90_DUST	RAD_KPC_10_DUST	RAD_KPC_50_DUST	RAD_KPC_90_DUST	VALID_PIX_SFR	AREA_KPC_10_SFR	AREA_KPC_50_SFR	AREA_KPC_90_SFR	RAD_KPC_10_SFR	RAD_KPC_50_SFR	RAD_KPC_90_SFR	TOTAL_SFR	TOTAL_DUST_MASS	HA_PIXEL	HB_PIXEL	TOTAL_DUST_MASSerr ALL	TOTAL_SFRerr	TOTAL_Attenuation	TOTAL_Attenuationerr	PIXEL_AREA	DUST_MASS_DENSITY	DIST_CM2	FLUX_LUM	HA_SUM
# #     1               2                   3                   4                   5              6             7              8                   9            10               11             12            13                14          15            16               17          18          19           20                      21                22
# #		MAGPIID	VALID_PIX_DUST	AREA_KPC_10_DUST	AREA_KPC_50_DUST	AREA_KPC_90_DUST	RAD_KPC_10_DUST	RAD_KPC_50_DUST	RAD_KPC_90_DUST	VALID_PIX_SFR	AREA_KPC_10_SFR	AREA_KPC_50_SFR	AREA_KPC_90_SFR	RAD_KPC_10_SFR	RAD_KPC_50_SFR	RAD_KPC_90_SFR	TOTAL_SFR	TOTAL_DUST_MASS	HA_PIXEL	HB_PIXEL	TOTAL_DUST_MASSerr ALL	TOTAL_SFRerr	TOTAL_Attenuation	TOTAL_Attenuationerr	PIXEL_AREA	DUST_MASS_DENSITY	DIST_CM2	FLUX_LUM	HA_SUM
# #with open("/Users/44185316/Downloads/MAGPI_full_new_att.csv") as csvfile: #cut is BD=2.86
# #with open("/Users/44185316/Downloads/MAGPI_full_new_att_new_SFR (1).csv") as csvfile: #cut is BD=2.86
# #with open("/Users/44185316/Downloads/MAGPI_full_new_att_new_SFR_BD_att_area_added_flat_all (1).csv")as csvfile: #flat lambda cdm
# #with open("/Users/44185316/Downloads/MAGPI_full_new_att_new_SFR_BD_att_area_added_flat_all2_CF00.csv") as csvfile: #using CFOO
# 

#     csvReader = csv.reader(csvfile, delimiter=',')
#     for row in csvReader:
#         MAGPIID2.append(int(row[1]))
#         TOTAL_DUST_MASS_ALL.append(float(row[17]))
#         TOTAL_DUST_MASSerr_ALL.append(float(row[20]))
#         TOTAL_SFRerr.append(float(row[21]))
#         VALID_PIX_TOTAL_DUST.append(int(row[2]))
#         TOTAL_SFR_ALL.append(float(row[16]))
#         TOTAL_Attenuation.append(float(row[22]))

# MAGPIID2=np.array(MAGPIID2)
# TOTAL_DUST_MASS_ALL=np.array(TOTAL_DUST_MASS_ALL)
# TOTAL_DUST_MASSerr_ALL=np.array(TOTAL_DUST_MASSerr_ALL)
# TOTAL_SFR_ALL=np.array(TOTAL_SFR_ALL)
# TOTAL_SFRerr=np.array(TOTAL_SFRerr)
# TOTAL_Attenuation=np.array(TOTAL_Attenuation)
# print(np.shape(MAGPIID2))


#%%

MAGPIID2=[]
TOTAL_DUST_MASS_ALL=[]
TOTAL_DUST_MASSerr_ALL=[]
VALID_PIX_TOTAL_DUST=[]
TOTAL_SFR_ALL=[]

#Limited header
with open("/Users/44185316/Downloads/MAGPI_full_new_att_new_SFR_BD_att_area_added_flat_all2_CF00_summing_2_86_area.csv")as csvfile: #by calculating sum(ha)/sum(hb)

    csvReader = csv.reader(csvfile, delimiter=',')
    for row in csvReader:
        MAGPIID2.append(int(row[1]))
        TOTAL_DUST_MASS_ALL.append(float(row[3]))
        TOTAL_DUST_MASSerr_ALL.append(0)
        VALID_PIX_TOTAL_DUST.append(float(row[2]))
        TOTAL_SFR_ALL.append(float(row[2]))

MAGPIID2=np.array(MAGPIID2)
TOTAL_DUST_MASS_ALL=np.array(TOTAL_DUST_MASS_ALL)
TOTAL_DUST_MASSerr_ALL=np.array(TOTAL_DUST_MASSerr_ALL)
TOTAL_SFR_ALL=np.array(TOTAL_SFR_ALL)
TOTAL_SFRerr=np.array(TOTAL_SFR_ALL)
TOTAL_Attenuation=np.array(TOTAL_DUST_MASS_ALL)
print(np.shape(MAGPIID2))



#%%
print(np.max(VALID_PIX_TOTAL_DUST))
#%%



find=np.where(TOTAL_DUST_MASS_ALL>0)
TOTAL_DUST_MASS_ALL2=TOTAL_DUST_MASS_ALL[find]
TOTAL_DUST_MASSerr_ALL2=TOTAL_DUST_MASSerr_ALL[find]
TOTAL_SFR_ALL2=TOTAL_SFR_ALL[find]
TOTAL_SFRerr2=TOTAL_SFRerr[find]
MAGPIID3=MAGPIID2[find]
TOTAL_Attenuation2=TOTAL_Attenuation[find]


TOTAL_DUST_MASS_ALL=TOTAL_DUST_MASS_ALL2
TOTAL_DUST_MASSerr_ALL=TOTAL_DUST_MASSerr_ALL2
TOTAL_SFR_ALL=TOTAL_SFR_ALL2
TOTAL_SFRerr=TOTAL_SFRerr2
MAGPIID2=MAGPIID3



#%%

TOTAL_DUST_MASS_ALL = np.array(TOTAL_DUST_MASS_ALL)
TOTAL_DUST_MASSerr_ALL = np.array(TOTAL_DUST_MASSerr_ALL)
TOTAL_SFR_ALL = np.array(TOTAL_SFR_ALL)
TOTAL_SFRerr = np.array(TOTAL_SFRerr)




log_TOTAL_DUST_MASS_ALL2 = np.log10(TOTAL_DUST_MASS_ALL)
log_TOTAL_DUST_MASSerr_ALL2 = TOTAL_DUST_MASSerr_ALL / (TOTAL_DUST_MASS_ALL * np.log(10))



log_TOTAL_SFR_ALL2 = np.log10(TOTAL_SFR_ALL)
log_TOTAL_SFRerr2 = TOTAL_SFRerr / (TOTAL_SFR_ALL * np.log(10))






#%%

############################################################
###################COMBINE BOTH SAMPLE#######################
############################################################3


name_magpis2=[]
Mstar_sed2=[]
Mstar_lerr_sed2=[]
Mstar_uerr_sed2=[]
Mdust_sed2=[]
Mdust_lerr_sed2=[]
Mdust_uerr_sed2=[]
SFR_sed2=[]
SFR_ustd_sed2=[]
SFR_lstd_sed2=[]
chi2_magpis2=[]
TOTAL_DUST_MASS_ALL2=[]
TOTAL_DUST_MASSerr_ALL2=[]
TOTAL_SFR_ALL2=[]
TOTAL_SFRerr2=[]
tauV_magpis2=[]
tauVism_magpis2=[]
tauV_ustd_magpis2=[]
tauV_lstd_magpis2=[]
tauVism_ustd_magpis2=[]
tauVism_lstd_magpis2=[]
Nama2=[]
#chi2_magpis2=[]
for i in range(len(MAGPIID2)):
    for j in range(len(name_magpis)):
        if MAGPIID2[i]==name_magpis[j]:
            #chi2_magpis2.append(chi2_magpis[j])
            
            Nama2.append(MAGPIID2[i])
            
            Mdust_sed2.append(log_Mdust[j])

            
            SFR_sed2.append(log_SFR[j])

            
            Mstar_sed2.append(log_Mstar[j])
            
            
            tauV_magpis2.append(float(tauV_magpis[j]))
            tauVism_magpis2.append(float(tauVism_magpis[j]))
            
            tauV_ustd_magpis2.append(tauV_ustd_magpis[j])
            tauV_lstd_magpis2.append(tauV_lstd_magpis[j])
            
            tauVism_ustd_magpis2.append(tauVism_ustd_magpis[j])
            tauVism_lstd_magpis2.append(tauVism_lstd_magpis[j])

            
 
            TOTAL_DUST_MASS_ALL2.append(log_TOTAL_DUST_MASS_ALL2[i])
            TOTAL_DUST_MASSerr_ALL2.append(log_TOTAL_DUST_MASSerr_ALL2[i])

            TOTAL_SFR_ALL2.append(log_TOTAL_SFR_ALL2[i])
            TOTAL_SFRerr2.append(log_TOTAL_SFRerr2[i])
            
#%%



name_g=[]
z=[]
R50_it=[]
mag_it=[]
mag_rt=[]
axrat_it=[]
#MAGPIID,segID,xmax,ymax,RAmax,Decmax,z,QOP,mag_it,mag_rt,R50_it,R90_it,axrat_it,ang_it
#with open("E:\MAGPI_Emission_SEG\MAGPI_master_source_catalogue.csv") as csvfile:
#with open("/Users/44185316/Library/CloudStorage/OneDrive-MacquarieUniversity/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
#with open("E:/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
with open("/Users/44185316/Library/CloudStorage/OneDrive-MacquarieUniversity/MAGPI_Emission_SEG/MAGPI_master_source_catalogue.csv") as csvfile:
    csvReader = csv.reader(csvfile, delimiter=',')
    for row in csvReader:
        name_g.append(int(row[0])) 
        z.append(float(row[6]))
        R50_it.append(float(row[10]))
        mag_it.append(float(row[8]))
        mag_rt.append(float(row[9]))
        axrat_it.append(float(row[12]))
        
        
        
redshift_MAGPI=[]
R50_it_MAGPI=[]
Nama_seg=[]
axrat_it_MAGPI=[]
for i in range(len(Nama2)):
    for j in range(len(name_g)):
        if Nama2[i]==name_g[j]:
            Nama_seg.append(Nama2[i])
            redshift_MAGPI.append(float(z[j]))
            R50_it_MAGPI.append(float(R50_it[j]))
            axrat_it_MAGPI.append(float(axrat_it[j]))




#%%

from astropy.cosmology import FlatLambdaCDM


radius_kpc5=[]
for i in range(len(R50_it_MAGPI)):
    radius=R50_it_MAGPI[i]
    Redshift=redshift_MAGPI[i]
    #cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) #with curvature

    cosmo=FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725,
                  Neff=3.04, Ob0=None)

    # Conversion factors
    MAGPI_pixel_size = 0.2  # arcsec per pixel
    kpc_per_arcmin = cosmo.kpc_comoving_per_arcmin(Redshift).to_value()
    kpc_per_arcsec = kpc_per_arcmin / 60
    kpc_per_pixel = kpc_per_arcsec * MAGPI_pixel_size
    radius_kpc2=radius*kpc_per_pixel
    radius_kpc5.append(radius_kpc2)



log_radius=np.log10(radius_kpc5)


#%%

Mdust_densityVism = []
Mdust_densityVism_upper = []
Mdust_densityVism_lower = []
# Step 1: Convert tauV and its bounds to attenuation, then to dust density
for i in range(len(tauVism_magpis2)):
    # Convert tauV to attenuation
    attenuation = tauVism_magpis2[i] * 1.086
    attenuation_upper = (tauVism_magpis2[i] + tauVism_ustd_magpis2[i]) * 1.086
    attenuation_lower = (tauVism_magpis2[i] - tauVism_lstd_magpis2[i]) * 1.086
    print(i,attenuation_lower,attenuation_upper,attenuation)

    # Convert attenuation to dust density
    density = attenuation * 5.7e5
    density_upper = attenuation_upper * 5.7e5
    density_lower = attenuation_lower * 5.7e5

    Mdust_densityVism.append(density)
    Mdust_densityVism_upper.append(density_upper)
    Mdust_densityVism_lower.append(density_lower)





#%%
plt.scatter(tauV_magpis2,tauVism_magpis2)
plt.xlabel('TauV')
plt.ylabel('TauVism')
plt.show()
plt.close()



#%%



Mdust_densityV = []
Mdust_densityV_upper = []
Mdust_densityV_lower = []
# Step 1: Convert tauV and its bounds to attenuation, then to dust density
for i in range(len(tauV_magpis2)):
    # Convert tauV to attenuation
    attenuation = tauV_magpis2[i] * 1.086
    attenuation_upper = (tauV_magpis2[i] + tauV_ustd_magpis2[i]) * 1.086
    attenuation_lower = (tauV_magpis2[i] - tauV_lstd_magpis2[i]) * 1.086
    print(i,attenuation_lower,attenuation_upper,attenuation)

    # Convert attenuation to dust density
    density = attenuation * 5.7e5
    density_upper = attenuation_upper * 5.7e5
    density_lower = attenuation_lower * 5.7e5

    Mdust_densityV.append(density)
    Mdust_densityV_upper.append(density_upper)
    Mdust_densityV_lower.append(density_lower)






#%%

# Convert tauV to attenuation
attenuation = np.array(tauV_magpis2) * 1.086
attenuation_upper = (np.array(tauV_magpis2) + np.array(tauV_ustd_magpis2)) * 1.086
attenuation_lower = (np.array(tauV_magpis2) - np.array(tauV_lstd_magpis2)) * 1.086


# Compute error bars
y_errors = [attenuation - attenuation_lower, attenuation_upper - attenuation]

# Create the x-axis positions
x_values = np.arange(len(Mdust_densityV))

# Plot the error bar chart
plt.figure(figsize=(8,6))
plt.errorbar(x_values, attenuation, yerr=y_errors, fmt='o', color='blue', 
             ecolor='black', capsize=5, capthick=1.5, elinewidth=1.5, label=r'$\tau_V$ Attenuation')

# Customize labels and title
plt.xlabel("Index", fontsize=14)
plt.ylabel(r'Attenuation ($\tau_V \times 1.086$)', fontsize=14)
plt.title("Attenuation Values with Error Bars", fontsize=16)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
#plt.yscale('log')
#plt.xscale('log')
# Show the plot
plt.show()
plt.close()
    
    


#%%

# Step 2: Calculate dust mass and its bounds
dust_mass_tauV = []
dust_mass_tauV_upper = []
dust_mass_tauV_lower = []

for i in range(len(radius_kpc5)):
    # Radius in kpc
    semi_major = 10**log_radius[i]  # Convert log radius to linear radius (kpc)
    
    ratio=axrat_it_MAGPI[i]
    
    semi_minor=semi_major*ratio
    
    
    # Dust density
    density = Mdust_densityV[i]
    density_upper = Mdust_densityV_upper[i]
    density_lower = Mdust_densityV_lower[i]

    # Area in kpc² (assuming circular regions)
    area = np.pi * semi_major * semi_minor #Elliptical

    # Calculate dust mass (M_sun)
    mass = density * area
    mass_upper = density_upper * area
    mass_lower = density_lower * area

    dust_mass_tauV.append(mass)
    dust_mass_tauV_upper.append(mass_upper)
    dust_mass_tauV_lower.append(mass_lower)

# Step 3: Convert to log scale
log_dust_mass_tauV = np.log10(dust_mass_tauV)
log_dust_mass_tauV_upper = np.log10(dust_mass_tauV_upper)
log_dust_mass_tauV_lower = np.log10(dust_mass_tauV_lower)

#%%

dust_mass_tauVism = []
dust_mass_tauVism_upper = []
dust_mass_tauVism_lower = []

for i in range(len(radius_kpc5)):
    # Radius in kpc
    semi_major = 10**log_radius[i]  # Convert log radius to linear radius (kpc)
    
    ratio=axrat_it_MAGPI[i]
    
    semi_minor=semi_major*ratio
    
    
    # Dust density
    density = Mdust_densityVism[i]
    density_upper = Mdust_densityVism_upper[i]
    density_lower = Mdust_densityVism_lower[i]

    # Area in kpc² (assuming circular regions)
    area = np.pi * semi_major * semi_minor

    # Calculate dust mass (M_sun)
    mass = density * area
    mass_upper = density_upper * area
    mass_lower = density_lower * area

    dust_mass_tauVism.append(mass)
    dust_mass_tauVism_upper.append(mass_upper)
    dust_mass_tauVism_lower.append(mass_lower)

# Step 3: Convert to log scale
log_dust_mass_tauVism = np.log10(dust_mass_tauVism)
log_dust_mass_tauVism_upper = np.log10(dust_mass_tauVism_upper)
log_dust_mass_tauVism_lower = np.log10(dust_mass_tauVism_lower)


#%%

plt.scatter(log_dust_mass_tauV,log_dust_mass_tauVism)
plt.xlabel('Dust mass(TauV)')
plt.ylabel('Dust mass (TauVism)')
plt.show()
plt.close()

#%%

TOTAL_DUST_MASS_ALL2=np.array(TOTAL_DUST_MASS_ALL2)
log_dust_mass_tauV=np.array(log_dust_mass_tauV)
TOTAL_DUST_MASSerr_ALL2=np.array(TOTAL_DUST_MASSerr_ALL2)



#%%

from scipy.stats import norm


# Calculate the offset
dust_offsetsss = TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauV



sigma_d_sed = (log_dust_mass_tauV_upper - log_dust_mass_tauV_lower) / 2


chi2 = np.sum((TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauV)**2 / (TOTAL_DUST_MASSerr_ALL2**2 + sigma_d_sed**2))
print(chi2)
reduced_chi2 = chi2 / len(TOTAL_DUST_MASS_ALL2)  # Reduced chi-square
print(reduced_chi2)



# Calculate the 16th and 84th percentiles
percentile_16 = np.percentile(dust_offsetsss, 16)
percentile_84 = np.percentile(dust_offsetsss, 84)

# Percentile-based median and error
median = np.median(dust_offsetsss)
percentile_error = (percentile_84 - percentile_16) / 2

# Gaussian median and sigma
gauss_median = np.median(dust_offsetsss)
gauss_sigma = np.std(dust_offsetsss)

# Create figure for stacked plots
plt.figure(figsize=(6, 10))

plt.tight_layout()  # Automatically adjust spacing between subplots

# Fine-tune the layout manually if needed
plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.1, right=0.95, top=0.95, bottom=0.1)


# Top scatter plot with custom tick parameters
ax1 = plt.subplot(2, 1, 1)  # Use subplot to stack scatter plot and histogram vertically
ax1.errorbar(TOTAL_DUST_MASS_ALL2, log_dust_mass_tauV, yerr=[log_dust_mass_tauV_lower,log_dust_mass_tauV_upper], xerr=TOTAL_DUST_MASSerr_ALL2,
             c='k', markeredgecolor='black', markersize=5, fmt='ko', ecolor='grey', elinewidth=1, capsize=3)
x1, y1 = [4, 10], [4, 10]
min_dust1 = 5.7558748556724915
ax1.plot(x1, y1, c='k', linestyle='--')
ax1.axvline(x=min_dust1, color='k', linestyle='-.')
ax1.axhline(y=min_dust1, color='k', linestyle='-.')
ax1.set_xlabel(r'log $M_{\mathrm{d,emission}}$', fontsize=14)
ax1.set_ylabel(r'log $M_{d, SED}$', fontsize=14)
ax1.set_xlim(4, 10)
ax1.set_ylim(4, 10)

# Custom tick parameters for ax1
ax1.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=True, right=True)
ax1.tick_params(which='minor', length=5, width=1, direction="in", top=True, right=True)
ax1.minorticks_on()

# Remove grid from ax1
ax1.grid(False)

# Bottom histogram plot
ax2 = plt.subplot(2, 1, 2)
ax2.hist(dust_offsetsss, bins=20, color='white', edgecolor='black', density=True)
ax2.axvline(0, color='k', linestyle='--', linewidth=2)  # Zero offset line
ax2.axvline(percentile_16, color='r', linestyle='--', linewidth=2)
ax2.axvline(percentile_84, color='r', linestyle='--', linewidth=2)

# Overlay a Gaussian curve with a solid line
mean = np.mean(dust_offsetsss)
std_dev = np.std(dust_offsetsss)
x_range = np.linspace(min(dust_offsetsss), max(dust_offsetsss), 100)
gaussian_curve = norm.pdf(x_range, mean, std_dev)
ax2.plot(x_range, gaussian_curve, 'k-')  # Solid line style

# Annotate both percentile-based and Gaussian values
#ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {(percentile_84 - percentile_16) / 2:.2g} dex\n'
#                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex',
#                    f'Reduced $\chi^2$ = {reduced_chi2:.2g}',
#         transform=ax2.transAxes, fontsize=10, color='k', verticalalignment='top')

ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {percentile_error:.2g} dex\n'
                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex\n'
                    f'Reduced $\chi^2$ = {reduced_chi2:.2f}',
         transform=ax2.transAxes, fontsize=11, color='k', verticalalignment='top')


#ax2.set_xlabel(r'$\Delta logM_{\mathrm{d}}$ (emission - SED)', fontsize=13)
ax2.set_xlabel(r'$\Delta$ log $M_{\mathrm{d, emission - SED}}$ ', fontsize=14)
ax2.set_ylabel('Number', fontsize=14)

# Remove x-axis major and minor ticks
ax2.tick_params(which='major', bottom=False)  # Remove major ticks
ax2.tick_params(which='minor', bottom=False)  # Remove minor ticks

# Customize remaining ticks
ax2.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=False, right=True)
ax2.tick_params(which='minor', length=5, width=1, labelsize=12, direction="in", top=False, right=True)
ax2.minorticks_on()

ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Show the combined plot
plt.tight_layout()
#plt.savefig("/Users/44185316/Desktop/MAGPI Plots/nebular_vs_SED_black5", bbox_inches='tight', dpi=1000)
plt.show()
plt.close()



#%%

from scipy.stats import norm


# Calculate the offset
dust_offsetsss = TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauV



sigma_d_sed = (log_dust_mass_tauV_upper - log_dust_mass_tauV_lower) / 2


chi2 = np.sum((TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauV)**2 / (TOTAL_DUST_MASSerr_ALL2**2 + sigma_d_sed**2))
print(chi2)
reduced_chi2 = chi2 / len(TOTAL_DUST_MASS_ALL2)  # Reduced chi-square
print(reduced_chi2)



# Calculate the 16th and 84th percentiles
percentile_16 = np.percentile(dust_offsetsss, 16)
percentile_84 = np.percentile(dust_offsetsss, 84)

# Percentile-based median and error
median = np.median(dust_offsetsss)
percentile_error = (percentile_84 - percentile_16) / 2

# Gaussian median and sigma
gauss_median = np.median(dust_offsetsss)
gauss_sigma = np.std(dust_offsetsss)

# Create figure for stacked plots
plt.figure(figsize=(6, 10))
plt.title('TauV')
plt.tight_layout()  # Automatically adjust spacing between subplots

# Fine-tune the layout manually if needed
plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.1, right=0.95, top=0.95, bottom=0.1)


# Top scatter plot with custom tick parameters
ax1 = plt.subplot(2, 1, 1)  # Use subplot to stack scatter plot and histogram vertically
ax1.errorbar(TOTAL_DUST_MASS_ALL2, log_dust_mass_tauV, xerr=TOTAL_DUST_MASSerr_ALL2,
             c='k', markeredgecolor='black', markersize=5, fmt='ko', ecolor='grey', elinewidth=1, capsize=3)
x1, y1 = [4, 10], [4, 10]
min_dust1 = 5.7558748556724915
ax1.plot(x1, y1, c='k', linestyle='--')
ax1.axvline(x=min_dust1, color='k', linestyle='-.')
ax1.axhline(y=min_dust1, color='k', linestyle='-.')
ax1.set_xlabel(r'log $M_{\mathrm{d,emission}}$', fontsize=14)
ax1.set_ylabel(r'log $M_{d, SED}$', fontsize=14)
ax1.set_xlim(4, 9.5)
ax1.set_ylim(4, 9.5)

# Custom tick parameters for ax1
ax1.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=True, right=True)
ax1.tick_params(which='minor', length=5, width=1, direction="in", top=True, right=True)
ax1.minorticks_on()

# Remove grid from ax1
ax1.grid(False)

# Bottom histogram plot
ax2 = plt.subplot(2, 1, 2)
ax2.hist(dust_offsetsss, bins=20, color='white', edgecolor='black', density=True)
ax2.axvline(0, color='k', linestyle='--', linewidth=2)  # Zero offset line
ax2.axvline(percentile_16, color='r', linestyle='--', linewidth=2)
ax2.axvline(percentile_84, color='r', linestyle='--', linewidth=2)

# Overlay a Gaussian curve with a solid line
mean = np.mean(dust_offsetsss)
std_dev = np.std(dust_offsetsss)
x_range = np.linspace(min(dust_offsetsss), max(dust_offsetsss), 100)
gaussian_curve = norm.pdf(x_range, mean, std_dev)
ax2.plot(x_range, gaussian_curve, 'k-')  # Solid line style

# Annotate both percentile-based and Gaussian values
#ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {(percentile_84 - percentile_16) / 2:.2g} dex\n'
#                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex',
#                    f'Reduced $\chi^2$ = {reduced_chi2:.2g}',
#         transform=ax2.transAxes, fontsize=10, color='k', verticalalignment='top')

ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {percentile_error:.2g} dex\n'
                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex\n'
                    f'Reduced $\chi^2$ = {reduced_chi2:.2f}',
         transform=ax2.transAxes, fontsize=11, color='k', verticalalignment='top')


#ax2.set_xlabel(r'$\Delta logM_{\mathrm{d}}$ (emission - SED)', fontsize=13)
ax2.set_xlabel(r'$\Delta$ log $M_{\mathrm{d, emission - SED}}$ ', fontsize=14)
ax2.set_ylabel('Number', fontsize=14)

# Remove x-axis major and minor ticks
ax2.tick_params(which='major', bottom=False)  # Remove major ticks
ax2.tick_params(which='minor', bottom=False)  # Remove minor ticks

# Customize remaining ticks
ax2.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=False, right=True)
ax2.tick_params(which='minor', length=5, width=1, labelsize=12, direction="in", top=False, right=True)
ax2.minorticks_on()

ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Show the combined plot
plt.tight_layout()
#plt.savefig("/Users/44185316/Desktop/MAGPI Plots/nebular_vs_SED_black5", bbox_inches='tight', dpi=1000)
plt.show()
plt.close()



#%%


from scipy.stats import norm


# Calculate the offset
dust_offsetsss = TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauVism



sigma_d_sed = (log_dust_mass_tauVism_upper - log_dust_mass_tauVism_lower) / 2


chi2 = np.sum((TOTAL_DUST_MASS_ALL2 - log_dust_mass_tauVism)**2 / (TOTAL_DUST_MASSerr_ALL2**2 + sigma_d_sed**2))
print(chi2)
reduced_chi2 = chi2 / len(TOTAL_DUST_MASS_ALL2)  # Reduced chi-square
print(reduced_chi2)



# Calculate the 16th and 84th percentiles
percentile_16 = np.percentile(dust_offsetsss, 16)
percentile_84 = np.percentile(dust_offsetsss, 84)

# Percentile-based median and error
median = np.median(dust_offsetsss)
percentile_error = (percentile_84 - percentile_16) / 2

# Gaussian median and sigma
gauss_median = np.median(dust_offsetsss)
gauss_sigma = np.std(dust_offsetsss)

# Create figure for stacked plots
plt.figure(figsize=(6, 10))


plt.tight_layout()  # Automatically adjust spacing between subplots

# Fine-tune the layout manually if needed
plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.1, right=0.95, top=0.95, bottom=0.1)
plt.title('TauIsm')

# Top scatter plot with custom tick parameters
ax1 = plt.subplot(2, 1, 1)  # Use subplot to stack scatter plot and histogram vertically

ax1.errorbar(TOTAL_DUST_MASS_ALL2, log_dust_mass_tauVism, xerr=TOTAL_DUST_MASSerr_ALL2,
             c='k', markeredgecolor='black', markersize=5, fmt='ko', ecolor='grey', elinewidth=1, capsize=3)
x1, y1 = [4, 10], [4, 10]
min_dust1 = 5.7558748556724915
ax1.plot(x1, y1, c='k', linestyle='--')
ax1.axvline(x=min_dust1, color='k', linestyle='-.')
ax1.axhline(y=min_dust1, color='k', linestyle='-.')
ax1.set_xlabel(r'log $M_{\mathrm{d,emission}}$', fontsize=14)
ax1.set_ylabel(r'log $M_{d, SED}$', fontsize=14)
ax1.set_xlim(4, 9.5)
ax1.set_ylim(4, 9.5)

# Custom tick parameters for ax1
ax1.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=True, right=True)
ax1.tick_params(which='minor', length=5, width=1, direction="in", top=True, right=True)
ax1.minorticks_on()

# Remove grid from ax1
ax1.grid(False)

# Bottom histogram plot
ax2 = plt.subplot(2, 1, 2)
ax2.hist(dust_offsetsss, bins=20, color='white', edgecolor='black', density=True)
ax2.axvline(0, color='k', linestyle='--', linewidth=2)  # Zero offset line
ax2.axvline(percentile_16, color='r', linestyle='--', linewidth=2)
ax2.axvline(percentile_84, color='r', linestyle='--', linewidth=2)

# Overlay a Gaussian curve with a solid line
mean = np.mean(dust_offsetsss)
std_dev = np.std(dust_offsetsss)
x_range = np.linspace(min(dust_offsetsss), max(dust_offsetsss), 100)
gaussian_curve = norm.pdf(x_range, mean, std_dev)
ax2.plot(x_range, gaussian_curve, 'k-')  # Solid line style

# Annotate both percentile-based and Gaussian values
#ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {(percentile_84 - percentile_16) / 2:.2g} dex\n'
#                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex',
#                    f'Reduced $\chi^2$ = {reduced_chi2:.2g}',
#         transform=ax2.transAxes, fontsize=10, color='k', verticalalignment='top')

ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {percentile_error:.2g} dex\n'
                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex\n'
                    f'Reduced $\chi^2$ = {reduced_chi2:.2f}',
         transform=ax2.transAxes, fontsize=11, color='k', verticalalignment='top')


#ax2.set_xlabel(r'$\Delta logM_{\mathrm{d}}$ (emission - SED)', fontsize=13)
ax2.set_xlabel(r'$\Delta$ log $M_{\mathrm{d, emission - SED}}$ ', fontsize=14)
ax2.set_ylabel('Number', fontsize=14)

# Remove x-axis major and minor ticks
ax2.tick_params(which='major', bottom=False)  # Remove major ticks
ax2.tick_params(which='minor', bottom=False)  # Remove minor ticks

# Customize remaining ticks
ax2.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=False, right=True)
ax2.tick_params(which='minor', length=5, width=1, labelsize=12, direction="in", top=False, right=True)
ax2.minorticks_on()

ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Show the combined plot
plt.tight_layout()
#plt.savefig("/Users/44185316/Desktop/MAGPI Plots/nebular_vs_SED_black5", bbox_inches='tight', dpi=1000)
plt.show()
plt.close()


#%%
log_dust_mass_tau_sum=log_dust_mass_tauV+log_dust_mass_tauVism
#%%


from scipy.stats import norm


# Calculate the offset
dust_offsetsss = TOTAL_DUST_MASS_ALL2 - log_dust_mass_tau_sum



sigma_d_sed = (log_dust_mass_tauVism_upper - log_dust_mass_tauVism_lower) / 2


chi2 = np.sum((TOTAL_DUST_MASS_ALL2 - log_dust_mass_tau_sum)**2 / (TOTAL_DUST_MASSerr_ALL2**2 + sigma_d_sed**2))
print(chi2)
reduced_chi2 = chi2 / len(TOTAL_DUST_MASS_ALL2)  # Reduced chi-square
print(reduced_chi2)



# Calculate the 16th and 84th percentiles
percentile_16 = np.percentile(dust_offsetsss, 16)
percentile_84 = np.percentile(dust_offsetsss, 84)

# Percentile-based median and error
median = np.median(dust_offsetsss)
percentile_error = (percentile_84 - percentile_16) / 2

# Gaussian median and sigma
gauss_median = np.median(dust_offsetsss)
gauss_sigma = np.std(dust_offsetsss)

# Create figure for stacked plots
plt.figure(figsize=(6, 10))


plt.tight_layout()  # Automatically adjust spacing between subplots

# Fine-tune the layout manually if needed
plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.1, right=0.95, top=0.95, bottom=0.1)
plt.title('Tau Sum ')

# Top scatter plot with custom tick parameters
ax1 = plt.subplot(2, 1, 1)  # Use subplot to stack scatter plot and histogram vertically

ax1.errorbar(TOTAL_DUST_MASS_ALL2, log_dust_mass_tau_sum, xerr=TOTAL_DUST_MASSerr_ALL2,
             c='k', markeredgecolor='black', markersize=5, fmt='ko', ecolor='grey', elinewidth=1, capsize=3)
x1, y1 = [4, 10], [4, 10]
min_dust1 = 5.7558748556724915
ax1.plot(x1, y1, c='k', linestyle='--')
ax1.axvline(x=min_dust1, color='k', linestyle='-.')
ax1.axhline(y=min_dust1, color='k', linestyle='-.')
ax1.set_xlabel(r'log $M_{\mathrm{d,emission}}$', fontsize=14)
ax1.set_ylabel(r'log $M_{d, SED}$', fontsize=14)
#ax1.set_xlim(4, 9.5)
#ax1.set_ylim(4, 9.5)

# Custom tick parameters for ax1
ax1.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=True, right=True)
ax1.tick_params(which='minor', length=5, width=1, direction="in", top=True, right=True)
ax1.minorticks_on()

# Remove grid from ax1
ax1.grid(False)

# Bottom histogram plot
ax2 = plt.subplot(2, 1, 2)
ax2.hist(dust_offsetsss, bins=20, color='white', edgecolor='black', density=True)
ax2.axvline(0, color='k', linestyle='--', linewidth=2)  # Zero offset line
ax2.axvline(percentile_16, color='r', linestyle='--', linewidth=2)
ax2.axvline(percentile_84, color='r', linestyle='--', linewidth=2)

# Overlay a Gaussian curve with a solid line
mean = np.mean(dust_offsetsss)
std_dev = np.std(dust_offsetsss)
x_range = np.linspace(min(dust_offsetsss), max(dust_offsetsss), 100)
gaussian_curve = norm.pdf(x_range, mean, std_dev)
ax2.plot(x_range, gaussian_curve, 'k-')  # Solid line style

# Annotate both percentile-based and Gaussian values
#ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {(percentile_84 - percentile_16) / 2:.2g} dex\n'
#                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex',
#                    f'Reduced $\chi^2$ = {reduced_chi2:.2g}',
#         transform=ax2.transAxes, fontsize=10, color='k', verticalalignment='top')

ax2.text(0.02, 0.8, f'Median = {median:.2g} ± {percentile_error:.2g} dex\n'
                    f'Gauss mean = {gauss_median:.2g} ± {gauss_sigma:.2g} dex\n'
                    f'Reduced $\chi^2$ = {reduced_chi2:.2f}',
         transform=ax2.transAxes, fontsize=11, color='k', verticalalignment='top')


#ax2.set_xlabel(r'$\Delta logM_{\mathrm{d}}$ (emission - SED)', fontsize=13)
ax2.set_xlabel(r'$\Delta$ log $M_{\mathrm{d, emission - SED}}$ ', fontsize=14)
ax2.set_ylabel('Number', fontsize=14)

# Remove x-axis major and minor ticks
ax2.tick_params(which='major', bottom=False)  # Remove major ticks
ax2.tick_params(which='minor', bottom=False)  # Remove minor ticks

# Customize remaining ticks
ax2.tick_params(which='major', length=10, width=1.5, labelsize=12, direction="in", top=False, right=True)
ax2.tick_params(which='minor', length=5, width=1, labelsize=12, direction="in", top=False, right=True)
ax2.minorticks_on()

ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Show the combined plot
plt.tight_layout()
#plt.savefig("/Users/44185316/Desktop/MAGPI Plots/nebular_vs_SED_black5", bbox_inches='tight', dpi=1000)
plt.show()
plt.close()



