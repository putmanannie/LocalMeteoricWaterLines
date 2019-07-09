# -*- coding: utf-8 -*-
"""
Created on Mon May 07 16:02:40 2018
This script is intended to take water isotope data and calculate
Local Meteoric Water Lines using orthogonal distance regression (ODR).
A site timeseries must have at least 48 months of data, with at least three 
samples in each 3-month non-dry season. This script also identifies the temporal 
variability in LMWL.

@author: Annie Putman
"""
def lmwlf(B, x):
    #local meteoric water line function
    #for use in ODR regression
    return B[0]*x +B[1]

import os    
import copy
import numpy as np
import pandas as pd
import netCDF4 as nc4
import scipy.odr as odr

#my modules
#this part needs to be changed by the user
suppfileloc = 'C:\Users\u0929173\Documents\Python\WaterIsotopeAnalysis_012016\MeteoricWaterLine\ForSupplementalToPublish'
os.chdir(os.path.join(suppfileloc, 'c:Scripts'))
import retrievedatasets
os.chdir(suppfileloc)


dataloc = """c:Data"""
(t, iso, Pre, Tmp, Lats, 
     Lons, Elev, SiteID, ProjectID) = retrievedatasets.getiso(os, np, nc4, copy, dataloc)
koppen = retrievedatasets.getkoppen(os, pd, dataloc)
    
iso = np.where(iso <=-9000, np.nan, iso)
iso = iso[:, :, 1:3]
mons = t.astype('datetime64[M]').astype(int)%12+1
seas = mons/3
seas = np.where(seas == 4, 0, seas)
ann = t.astype('datetime64[Y]').astype(int)-min(t.astype('datetime64[Y]').astype(int))

#get seasonal average precipitation for data quality threshold
Pre2 = np.where(Pre == -1, 0, Pre)
PreAvgs = np.zeros((Pre.shape[0], 4))
for i in range(3):
    inds = np.where((seas) == i)[0]
    PreAvgs[:, i] = np.average(Pre2[:, inds], axis = 1)

          
#%%    ################# CALCULATE LMWL w ODR ##########################
#initialize the odr model
lmwlmod = odr.Model(lmwlf)
#instantiate with GMWL
beta0 = [8.0, 10.0]
regresults = pd.DataFrame(columns = ["site", "m", "mse", "b", "bse", 
                                     "inv_condnum", "res_var", 
                                     "n", "lat", "lon", "elev", "proj", 
                                     "koppen", "clim"])

for i in range(iso.shape[0]):    
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    mask = np.multiply(~np.isnan(d2H), ~np.isnan(d18O))
    #apply the 48 month threshold
    if len(d2H[mask])>= 48:
        vals, counts = np.unique(seas[mask], return_counts = True)
        #these are the months where we expect at least one sample to be represented
        precexpected = np.where(PreAvgs[i, :] > 10)[0] 
        datthresh = np.zeros(4)
        for j in range(4): 
            if (j in precexpected) and (j in vals):
                if counts[j == vals] >= 3:
                    datthresh[j] = 1
            elif(j not in precexpected):
                datthresh[j] = 1
                    
        if all(datthresh == 1):
            #standard deviations of x vals
            sx = np.ones(len(d18O[mask]))*0.1
            #standard deviations of y vals
            sy = np.ones(len(d18O[mask]))*0.5
            realdata = odr.RealData(d18O[mask], d2H[mask], sx=sx, sy=sy)
            model = odr.ODR(realdata, lmwlmod, beta0 = beta0)
            results = model.run()
            #get koppen climate classification

            minind = np.argmin(np.add(abs(koppen['Lat'].values-Lats[i]), 
                                      abs(koppen['Lon'].values-Lons[i])))
            if ((abs(koppen['Lat'].values-Lats[i])[minind] < 0.75) and 
                (abs(koppen['Lon'].values-Lons[i])[minind] < 0.75)):
                kcls = koppen['Cls'].loc[minind]
            else:
                kcls = 'N'
            regresults.loc[i] = np.array([SiteID[i], results.beta[0], results.sd_beta[0], 
                                         results.beta[1], results.sd_beta[1], 
                                         results.inv_condnum, results.res_var, 
                                         len(d18O[mask]), Lats[i], Lons[i], Elev[i], 
                                         np.unique(ProjectID[i][mask]).astype(int)[0], 
                                         kcls, kcls[0]])

regresults[["m", "mse", "b", "bse", "inv_condnum", "res_var", "n", 
            "lat", "lon", "elev", "proj"]] = regresults[["m", "mse", "b", "bse", 
                                        "inv_condnum", "res_var", "n", "lat", "lon", "elev", 
                                        "proj"]].astype(float)
                                        
regresults.to_csv(os.path.join(dataloc, "LMWL_data_ODR.csv"), index = True)

#%% ######### CALCULATE TEMPORAL VARIABILITY IN LMWL #########################


#number of years of data required 
beta0 = [8.0, 10.0]
lmwlmod = odr.Model(lmwlf)
nyrs = 4
nmons = 48
monspan = 120
#minimum total timespan required for analysis to be useful
timecov = 4   
temporalresults = pd.DataFrame(columns = ["site", "start", "end", "m", "mse", 
                                          "b", "bse", "inv_condnum", "res_var", 
                                          "n", "lat", "lon", "elev", "proj", 
                                          "var", "koppen", "clim"])
allinds = np.array(range(iso.shape[1]))
rows = 0
for i in range(iso.shape[0]): 
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    mask = np.multiply(~np.isnan(d2H), ~np.isnan(d18O)) 
    precexpected = np.where(PreAvgs[i, :] > 10)[0]
    if len(allinds[mask])> 3*nyrs*len(precexpected):
        tcoverage = (max(allinds[mask]) - min(allinds[mask]))/12.0
        #check that there are at least 10 years of data (when prec expected)
        #check that the timespan covers at least 25 years
        if tcoverage >= timecov:
            inds = allinds[mask] 
            for k in range(len(inds)-nmons):
                    indssub = inds[k:(k+nmons)]
                    #make sure the lmwl does not span more than 36 months
                    if (max(indssub) - min(indssub)) <= monspan:
                        vals, counts = np.unique(seas[indssub], return_counts = True)
                        #need to work on this stuff...
                        datthresh = np.zeros(4)
                        for j in range(4): 
                            if (j in precexpected) and (j in vals):
                                if counts[j == vals] >= 3:
                                    datthresh[j] = 1
                            elif(j not in precexpected):
                                datthresh[j] = 1

                    if all(datthresh == 1):
                        sx = np.ones(len(d18O[indssub]))*0.1
                        #standard deviations of y vals
                        sy = np.ones(len(d18O[indssub]))*0.5
                        realdata = odr.RealData(d18O[indssub], d2H[indssub], sx=sx, sy=sy)
                        model = odr.ODR(realdata, lmwlmod, beta0 = beta0)
                        results = model.run()
                        #results = sm.OLS(d2H[indssub], sm.add_constant(d18O[indssub]), missing = 'drop').fit()
                        #get koppen climate classification

                        minind = np.argmin(np.add(abs(koppen['Lat'].values-Lats[i]), 
                                                  abs(koppen['Lon'].values-Lons[i])))
                        if ((abs(koppen['Lat'].values-Lats[i])[minind] < 0.75) 
                            and (abs(koppen['Lon'].values-Lons[i])[minind] < 0.75)):
                            kcls = koppen['Cls'].loc[minind]
                        else:
                            kcls = 'N'
                        temporalresults.loc[rows] = np.array([SiteID[i], t[min(indssub)], 
                                        t[max(indssub)], results.beta[0], 
                                        results.sd_beta[0], results.beta[1], 
                                        results.sd_beta[1], results.inv_condnum, 
                                        results.res_var, len(indssub), 
                                        Lats[i], Lons[i], Elev[i], 
                                        np.unique(ProjectID[i, indssub]).astype(int)[0], 
                                        np.var(d2H[indssub]), kcls, 
                                        kcls[0]])
                        rows = rows+1
                    else:
                        #if it fails when the dataset is the largest it can be, then do not try at
                        #smaller dataset sizes
                        print('could not create a suitable subset')


temporalresults[["m", "mse", "b", "bse", "inv_condnum",  
            "res_var", "n", "lat", "lon", 
            "elev", "proj", "var"]] = temporalresults[["m", "mse", "b", "bse",  
                                        "inv_condnum", "res_var", "n", "lat", 
                                        "lon", "elev","proj", "var"]].astype(float)
                                        
temporalresults.to_csv(os.path.join(dataloc, "LMWL_temporaldata_ODR.csv"), index = True)
#%% ######### CALCULATE SEASONAL D-EXCESS VALUES FOR MECHANISM TESTING ################    
dexcessres = pd.DataFrame(columns = ["site", "DJF_dex", "JJA_dex", "DJF_dex_s", "JJA_dex_s", "DJF_d18O", "JJA_d18O", "DJF_d18O_s", "JJA_d18O_s", "koppen", "clim"])

for i in range(iso.shape[0]):    
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    mask = np.multiply(~np.isnan(d2H), ~np.isnan(d18O))
    if (len(d2H[mask])>= 48) & (76 not in np.unique(ProjectID[i][mask]).astype(int)):
        vals, counts = np.unique(seas[mask], return_counts = True)
        precexpected = np.where(PreAvgs[i, :] > 10)[0] 
        datthresh = np.zeros(4)
        for j in range(4): 
            if (j in precexpected) and (j in vals):
                if counts[j == vals] >= 3:
                    datthresh[j] = 1
            elif(j not in precexpected):
                datthresh[j] = 1
                    
        if all(datthresh == 1):
            djfinds = np.where((mons[mask] == 12) | (mons[mask] <= 2))[0]
            jjainds = np.where((mons[mask] >= 6) & (mons[mask] <= 8))[0]
            #get koppen climate classification
            minind = np.argmin(np.add(abs(koppen['Lat'].values-Lats[i]), abs(koppen['Lon'].values-Lons[i])))
            if abs(koppen['Lat'].values-Lats[i])[minind] < 0.75 and abs(koppen['Lon'].values-Lons[i])[minind] < 0.75:
                kcls = koppen['Cls'].loc[minind]
            else:
                kcls = 'N'
            dex = np.subtract(d2H[mask], 8*d18O[mask])
            djf_dex = np.mean(dex[djfinds])
            jja_dex = np.mean(dex[jjainds])
            djf_dex_s = np.std(dex[djfinds])
            jja_dex_s = np.std(dex[jjainds])
            
            djf_d18O = np.mean(d18O[mask][djfinds])
            jja_d18O = np.mean(d18O[mask][jjainds])
            djf_d18O_s = np.std(d18O[mask][djfinds])
            jja_d18O_s = np.std(d18O[mask][jjainds])

            dexcessres.loc[i] = np.array([SiteID[i], djf_dex, jja_dex, djf_dex_s, 
                                         jja_dex_s, djf_d18O, jja_d18O, djf_d18O_s, 
                                         jja_d18O_s, kcls, kcls[0]])



dexcessres[["DJF_dex", "JJA_dex", "DJF_dex_s", "JJA_dex_s", "DJF_d18O", 
            "JJA_d18O", "DJF_d18O_s", 
            "JJA_d18O_s",]] = dexcessres[["DJF_dex", "JJA_dex", "DJF_dex_s", "JJA_dex_s", 
                          "DJF_d18O", "JJA_d18O", "DJF_d18O_s", "JJA_d18O_s",]].astype(float)
                                        
dexcessres.to_csv(os.path.join(dataloc, "LMWL_data_dexcess.csv"), index = True)
#%%    

        
