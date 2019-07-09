# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 14:23:41 2018

@author: u0929173
"""
def lmwlf(B, x):
    #local meteoric water line function
    #for use in ODR regression
    return B[0]*x +B[1]

#%%
import numpy as np
import os
import statsmodels.api as sm
import pandas as pd
import netCDF4 as nc4
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
import matplotlib.cm as cm
import matplotlib as mpl
import scipy as sp
import scipy.odr as odr
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import copy

suppfileloc = 'C:\Users\u0929173\Documents\Python\WaterIsotopeAnalysis_012016\MeteoricWaterLine\ForSupplementalToPublish'
os.chdir(os.path.join(suppfileloc, 'c:Scripts'))
import retrievedatasets
import MakeFigures
os.chdir(suppfileloc)

#%%
dataloc = """c:Data"""
regresults = pd.read_csv(os.path.join(dataloc, "LMWL_data_ODR_Published.csv"), index_col = 0)
usnipmask = regresults['proj'].values != 76.0
timemask = regresults['n'].values >= 48
mask = usnipmask*timemask

bayesell = pd.read_csv(os.path.join(dataloc, "Bayesell_Zscore.csv"), index_col = 0) 
dexcessres = pd.read_csv(os.path.join(dataloc, "LMWL_data_dexcess.csv"), index_col = 0)

uniquesites = np.unique(bayesell['site'].values)
nmons = np.array([18, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120]) 
#%%
(t, iso, Pre, Tmp, Lats, 
     Lons, Elev, SiteID, ProjectID) = retrievedatasets.getiso(os, np, nc4, copy, dataloc)

allinds = np.array(range(iso.shape[1]))
reclen = np.zeros(iso.shape[0])
for i in range(iso.shape[0]): 
    d2H = iso[i, :, 0]
    d18O = iso[i, :, 1]
    #find where both d2H and d18O have data
    lenmask = np.multiply(d2H > -9000, d18O > -9000) 
    if any(lenmask):
        inds = allinds[lenmask]       
        reclen[i] = len(inds)
        
koppen = retrievedatasets.getkoppen(os, pd, dataloc)

modslopes = retrievedatasets.getSWING(os, np, pd, odr, lmwlf, os.path.join(dataloc, 'SWING'), regresults, mask)
#%%
#Figure 1 - How different numbers of samples ellipses compare to full dataset
MakeFigures.Figure1(os, np, sp, mpl, bayesell, nmons, reclen, uniquesites)

#%% Map with Koppen background (pastel) and slope/intercept values in colorscale
colordict = {'A: Tropical \n':'#e8ddee', 'B: Arid \n':'#fff0cc', 'D: Continental\n':'#dae1f1', 
             'Cf: Temperate \n humid':'#daf1da', 'Cs: Temperate \n dry summer': '#ebf1da' , 
             'Cw: Temperate \n dry winter': '#daf1f0', 'E: Polar \n and alpine':'#e6e6e6'}  
koppengrid1 = retrievedatasets.getkoppengrid(np, koppen)
MakeFigures.Figure2(os, np, mpl, Basemap, koppen, koppengrid1, regresults, mask, colordict)

#%%   Latitude v LMWL parameters, grouped by Koppen Classes 
MakeFigures.Figure3(os, np, mpl, Basemap, Rectangle, regresults, mask)
#%% d-excess by climate class and season
joined = regresults.join(dexcessres, lsuffix = '_caller', rsuffix = '_other')
MakeFigures.Figure4(os, np, mpl, Basemap, joined)
#%% Compare data and model by Koppen Climate Class
MakeFigures.Figure6(os, np, mpl, modslopes, regresults, mask)

