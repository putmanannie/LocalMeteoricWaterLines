# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 14:11:38 2018

This script is to calculate a Baysian Ellipse to characterize how many months 
data are required to characterize a site.

variables are d2H and d18O, rarifaction is from MeteoricWaterLine_V1.py

The analysis is based on 
'Comparing isotopic niche widths among and within communities: 
SIBER – Stable Isotope Bayesian Ellipses in R'
https://doi.org/10.1111/j.1365-2656.2011.01806.x 

@author: Annie Putman
"""
def create_ellipse(Point, affinity, center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1)
    ell = affinity.scale(circ, round(lengths[0], 5), round(lengths[1], 5))
    ellr = affinity.rotate(ell, angle)
    return ellr

def ellipse_estimation(pymc, np, X_zscore):
    
    mu = pymc.Normal('mu', 0, 1.0/(10**3), size = 2)
    #set prior for covariance matrix using wishart distribution (following Jackson, 2011)
    ivcov = pymc.Wishart('ivcov', 2, 2*np.eye(2))
    # incorporate observations        
    obs = pymc.MvNormal("observed returns", mu, tau = ivcov, observed = True, value = X_zscore)
    #define a model
    model = pymc.Model([obs, mu, ivcov])
    mcmc = pymc.MCMC(model)
    #run two separate chains for comparison
    for t in range(ntrace):
        mcmc.sample(iter = 6000, burn = 1000, thin = 2)
    #compare the chains
    grconv = pymc.gelman_rubin(mcmc, return_var = True)
    #get mu and ivcov estimates
    ivcov_post = mcmc.trace('ivcov')[:]                      
    mu_all = mcmc.trace('mu')[:].mean(axis = 0)
    
    return grconv, ivcov_post, mu_all, mcmc

def extract_params(np, sp, ivcov_post, mcmc):
        #make a zero-filled array to hold the major and minor axis, and the rotation
    avec = np.zeros(ivcov_post.shape[0])
    bvec = copy.copy(avec)
    thetavec = copy.copy(avec)
    vcov_post_all = np.zeros((ivcov_post.shape[0], 2, 2))
    #have to go through each inverse covariance matrix produced (5000)
    for j in range(ivcov_post.shape[0]):  
        #get a temporary covariance matrix by inverting the inverse covariance matrix          
        vcov_tmp = np.linalg.inv(ivcov_post[j]) 
        #convert from z-score to full size
        vcov_post_all[j] = np.multiply(vcov_tmp, X_cov)
        #eigen decomposition to get the ellipse parameters
        lamda, vec = np.linalg.eig(vcov_post_all[j])
        #major axis, minor axis, rotation
        avec[j] = 2*np.power(lamda[0], 0.5)
        bvec[j] = 2*np.power(lamda[1], 0.5)
        thetavec[j] = np.arcsin(vec[0, 1])
    #check if the values have converged...
    percdiff = np.zeros(3)
    ttest_mu = np.zeros(2)
    stored = len(avec)/2
    percdiff[0] = (np.mean(avec[:stored])- np.mean(avec[stored:]))*100/np.mean(avec)   
    percdiff[1] = (np.mean(bvec[:stored])- np.mean(bvec[stored:]))*100/np.mean(bvec)
    percdiff[2] = (np.mean(thetavec[:stored])- np.mean(thetavec[stored:]))*100/np.mean(thetavec)
    
    ttest_mu[0] = sp.stats.ttest_ind(mcmc.trace('mu', chain = 0)[:, 0], mcmc.trace('mu', chain = 1)[:, 0], equal_var=False)[1]
    ttest_mu[1] = sp.stats.ttest_ind(mcmc.trace('mu', chain = 0)[:, 1], mcmc.trace('mu', chain = 1)[:, 1], equal_var=False)[1] 
    
    return percdiff, ttest_mu, avec, bvec, thetavec 

import pandas as pd
import numpy as np
import pymc
import os
import copy
import netCDF4 as nc4
import scipy as sp
from shapely.geometry.point import Point
from shapely import affinity

suppfileloc = 'C:\Users\u0929173\Documents\Python\WaterIsotopeAnalysis_012016\MeteoricWaterLine\ForSupplemental'
os.chdir(os.path.join(suppfileloc, 'c:Scripts'))
import retrievedatasets
os.chdir(suppfileloc)

dataloc = """c:Data"""
#get processed isotope data
t, iso, Pre, Tmp, Lats, Lons, Elev, SiteID, ProjectID = retrievedatasets.getiso(os, np, nc4, copy, dataloc)
#get koppen climate classifications    
koppen = retrievedatasets.getkoppen(os, pd, dataloc)
#%%  
iso = np.where(iso <=-9000, np.nan, iso)
iso = iso[:, :, 1:3]

mons = t.astype('datetime64[M]').astype(int)%12+1
uniquemons = range(12)
seas = mons/3
seas = np.where(seas == 4, 0, seas)
ann = t.astype('datetime64[Y]').astype(int)-min(t.astype('datetime64[Y]').astype(int))

#Get the length of each record to identify those longer than 120 months
#an array of indicies
allinds = np.array(range(iso.shape[1]))
#zero array for holding record lengths
reclen = np.zeros(iso.shape[0])
for i in range(iso.shape[0]): 
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    mask = np.multiply(~np.isnan(d2H), ~np.isnan(d18O)) 
    if any(mask):
        inds = allinds[mask]       
        reclen[i] = len(inds)
#make a list of the longest records 
longrecs = np.where(reclen >= 120)[0]

#%%
#counter
rows = 0
#empty array for later use
Pre2 = np.where(Pre == -1, 0, Pre)
PreAvgs = np.zeros((Pre.shape[0], 4))
for i in range(4):
   PreAvgs[:, i] = np.average(Pre2[:, seas == i], axis = 1)    
#need to calculate the seasonal average values for precipitation

#set thresholds/constants for use
calcs = 10
ctr = 0
ntrace = 2
#scaling factor for 95% confidence ellipse
s = 5.99146

#initialize the dataframe with no rows
bayesell = pd.DataFrame(columns = ["site", "n", "d18O", "d2H", "a", "b", "theta", "area", "intersect_area", "perc overlap", 'koppen', 'clim'])
#%%
# estimate the ellipses
for i in longrecs: 
    #find the index in the koppen dataframe closest to the site of the the long record
    minind = np.argmin(np.add(abs(koppen['Lat'].values-Lats[i]), abs(koppen['Lon'].values-Lons[i])))
    #if the minimum distance is less than 0.75 degrees in both lat and lon, then assign.
    #we need to do this because no koppen classes are indicated over the ocean, and most
    #islands don't have koppen climate classes,
    if abs(koppen['Lat'].values-Lats[i])[minind] < 0.75 and abs(koppen['Lon'].values-Lons[i])[minind] < 0.75:
        kcls = koppen['Cls'].loc[minind]
    else:
        #this is typically assigned to islands.
        kcls = 'N'    
    
    #this is for readability and to minimize indexing
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    mask = np.multiply(~np.isnan(d2H), ~np.isnan(d18O)) 
    inds = allinds[mask]       
    ctr = ctr+1
    #these are the months where we expect at least one sample to be represented
    #this is so as not to unfairly penalize seasonally dry locations
    precexpected = np.where(PreAvgs[i, :] > 10)[0]
    
    #These are the results for an ellipse calculated from covariance of the full dataset
    # Non- Bayesian, calculated for comparison. code (0)
    X = np.concatenate((d18O[inds, None], d2H[inds, None]), axis = 1)
    mu_det = np.mean(X, axis = 0)
    vcov_det = np.cov(X.T)
    lamda, vec = np.linalg.eig(vcov_det)
    #length of major axis
    a_det = np.sqrt(s*lamda[0])
    #length of minor axis
    b_det = np.sqrt(s*lamda[1])
    #angle relative to x axis
    theta_det = np.arcsin(vec[0, 1])
    #write to the dataframe at index 'rows'
    bayesell.loc[rows] = np.array([SiteID[i], 0,  mu_det[0], mu_det[1], a_det, 
                                  b_det, theta_det, np.nan, np.nan, np.nan, 
                                  kcls, kcls[0]])
    #increment indexer
    rows = rows+1    
    #initialize these comparison parameters
    percdiff = np.zeros(3)+10
    ttest_mu = np.zeros(2)+10
    #calculate z scores for X - the mean and cov are to transform back from z-scores to values
    X_mean = np.mean(X, axis = 0)
    X_cov = np.cov(X.T)
    X_zscore = np.divide(np.subtract(X, np.mean(X, axis = 0)), np.std(X, axis = 0))
    #iterate until convergence - all values within 5% of one another, 
    #mean values are same to p < 0.1
    while all(abs(percdiff) > 5) and all(ttest_mu > 0.1):

        grconv, ivcov_post, mu_all, mcmc = ellipse_estimation(pymc, np, X_zscore)
        percdiff, ttest_mu, avec, bvec, thetavec = extract_params(np, sp, ivcov_post, mcmc)

    #convert from z score mean back to correct scale/units
    mu_all = np.add(mu_all, X_mean)
    #take the average of all 5000 instances of the chain
    a_all = np.mean(avec)
    b_all = np.mean(bvec)
    #angle relative to x axis
    theta_all = np.mean(thetavec)    
    #create an ellipse using this information
    ellr_all = create_ellipse(Point, affinity, (mu_all[0], mu_all[1]), 
                              (a_all, b_all), angle = theta_all*180/np.pi)    
    #record ellipse parameters for full timeseries
    bayesell.loc[rows] = np.array([SiteID[i], len(inds),  mu_all[0], mu_all[1], 
                                  a_all, b_all, theta_all, ellr_all.area, np.nan, 
                                  1, kcls, kcls[0]]) 
    #increment indexer
    rows = rows+1
    
    #save in the interim
    bayesell.to_csv(os.path.join(dataloc, "Bayesell_Zscore.csv"), index = True)  
    
    for k in np.array([18, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120]):
        trycounter = 0
        rounds = 0
        datthresh = np.zeros(4)
        #select double the number of starting points needed in 
        #case some don't work
        startinds = np.random.choice(len(inds), calcs*3, replace = False)

        while not all(datthresh == 1) and (trycounter < 10) and (rounds < min([calcs, (len(inds)-k)+1])):
            
            datthresh = np.zeros(4) 
            #this part could fail bc startinds[rounds]+k could be greater than len(inds)
            if len(inds) > (startinds[rounds]+k): 
                indssub = inds[startinds[rounds]:(startinds[rounds]+k)]
            else:                                
                #wrap the dataset
                endind = k-(len(inds)-startinds[rounds])
                indssub = np.concatenate((inds[:endind], inds[startinds[rounds]:]))
                
            vals, counts = np.unique(seas[indssub], return_counts = True)
            #minor requirement for seasonality
            for j in range(4): 
                if (j in precexpected) and (j in vals):
                    if counts[j == vals] >= 3:
                        datthresh[j] = 1
                elif(j not in precexpected):
                    datthresh[j] = 1
 
            if all(datthresh == 1):
                #from this :https://stackoverflow.com/questions/29752042/bayesian-covariance-prediction-with-pymc
                #--------------------------------------------------------
                #X is the set of samples (write this in once the sub-sampling is in place)                    
                X = np.concatenate((d18O[indssub, None], d2H[indssub, None]), axis = 1) 
                percdiff = np.zeros(3)+10
                ttest_mu = np.zeros(2)+10
                X_mean = np.mean(X, axis = 0)
                X_cov = np.cov(X.T)
                X_zscore = np.divide(np.subtract(X, np.mean(X, axis = 0)), np.std(X, axis = 0))
                attempts = 0
                #these are the vague priors as described in Jackson et al (2011)
                while all(abs(percdiff) > 5) and all(ttest_mu > 0.1) and attempts <=5: 
                    grconv, ivcov_post, mu_post, mcmc = ellipse_estimation(pymc, np, X_zscore)
                    percdiff, ttest_mu, avec, bvec, thetavec = extract_params(np, sp, ivcov_post, mcmc)
                #these are all of the components of an ellipse
#                    lamda, vec = np.linalg.eig(vcov_post_med)
                #length of major axis
                if attempts <= 5:
                    mu_post = np.add(mu_post, X_mean)
                    a = np.mean(avec)
                    b = np.mean(bvec)
                    #angle relative to x axis
                    theta = np.mean(thetavec)
                    ellr = create_ellipse(Point, affinity, (mu_post[0], mu_post[1]), (a, b), angle = theta*180/np.pi)
                    intersect = ellr.intersection(ellr_all)
                    bayesell.loc[rows] = np.array([SiteID[i], k,  mu_post[0], 
                                                  mu_post[1], a, b, theta, ellr.area, 
                                                 intersect.area, intersect.area/ellr.area, 
                                                 kcls, kcls[0]])
                    rows = rows+1
                    rounds = rounds+1
                    
                datthresh = np.zeros(4)
                trycounter = 0
                #count the number of times that a subset has been tried
                
            else:
                #if it fails when the dataset is the largest it can be, then do not try at
                #smaller dataset sizes
                print('could not create a suitable subset')
                trycounter = trycounter+1  
        
        # save in the interim!        
        bayesell.to_csv(os.path.join(dataloc, "Bayesell_Zscore.csv"), index = True) 
        
#final save
bayesell.to_csv(os.path.join(dataloc, "Bayesell_Zscore.csv"), index = True)  

