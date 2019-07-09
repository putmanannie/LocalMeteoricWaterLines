# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:18:07 2019

@author: u0929173
"""

def getiso(os, np, nc4, copy, fileloc):
    print('Reading NetCDF4 file...')
    f = nc4.Dataset(os.path.join(fileloc, 'Precipitation_WI_LMWL_notembargoed.nc'), 'r', format = 'NETCDF4')
    #f = nc4.Dataset(os.path.join(fileloc, 'Precipitation_WI_4LMWL.nc'), 'r', format = 'NETCDF4')
    start = copy.copy(f.variables['start'][:])
    start = np.datetime64(start[0])
    t = copy.copy(f.variables['time'][:])
    t = start + t.astype('timedelta64[D]')
    
    iso = copy.copy(f.variables['isoarray'][:])
    Pre = copy.copy(f.variables['precip'][:])
    Tmp = copy.copy(f.variables['temper'][:])
    Lats = copy.copy(f.variables['latitudes'][:])
    Lons = copy.copy(f.variables['longitudes'][:])
    Elev = copy.copy(f.variables['elevations'][:])
    SiteID = copy.copy(f.variables['siteid'][:])
    SiteID = SiteID.astype('|S50')         
    ProjectID = copy.copy(f.variables['projects'][:])
    ProjectID = ProjectID.astype('|S6')
    
    f.close()
    return t, iso, Pre, Tmp, Lats, Lons, Elev, SiteID, ProjectID

def getkoppen(os, pd, fileloc):
    #koppen classification key:
    #http://koeppen-geiger.vu-wien.ac.at/shifts.htm
    #Main Climates:
    #   A: Equatorial, B: Arid C: warm temperate D: snow E: Polar
    #Precipitation
    # W: desert S: Steppe f:fully humid s:summer dry w:winter dry m:monsoonal
    #Temperature
    #   h:hot arid k:cold arid a:hot summer b:warm summer c:cool summer 
    #   d:extremely continental F: polar frost T: polar tundra
    
    koppen = pd.read_table(os.path.join(fileloc, '1976-2000_ASCII.txt'), header = 0, delim_whitespace = True)
    return koppen

def getSWING(os, np, pd, odr, lmwlf, swingloc, regresults, mask):
    
    filenamepre =  ['cam2_1x1_',
     'cam5_1x1_',
     'echam5_1x1_',
     'gissE_1x1_',
     'HadAM3_1x1_',
     'isoGSM_1x1_',
     'lmdz_1x1_',
     'miroc_1x1_']
     
    filenameiso = ['d18O.csv', 'd2H.csv']
    
    modslopes = pd.DataFrame(columns = ["lat", "lon", "cam2_m", "cam2_b", "cam5_m", "cam5_b", 
                                        "echam5_m", "echam5_b", "gissE_m", "gissE_b", 
                                        "HadAM3_m", "HadAM3_b", "isoGSM_m", "isoGSM_b", 
                                        "lmdz_m", "lmdz_b", "miroc_m", "miroc_b"], index = regresults.index[mask])
    lmwlmod = odr.Model(lmwlf)
    for i in regresults.index[mask]:
        modslopes['lat'].loc[i] = regresults['lat'].loc[i]
        modslopes['lon'].loc[i] = regresults['lon'].loc[i]
        for j in range(len(filenamepre)):
            d18O_model = pd.read_csv(os.path.join(swingloc, filenamepre[j]+filenameiso[0]), index_col = 1)
            d2H_model = pd.read_csv(os.path.join(swingloc, filenamepre[j]+filenameiso[1]), index_col = 1)
            #check to make sure that the the measurement data index is in the model data index
            if (i in d2H_model.index) and (i in d18O_model.index):           
                #check to make sure there are enough values to estimate a lmwl
                submask = ~np.isnan(d2H_model.loc[i].values[3:])*~np.isnan(d18O_model.loc[i].values[3:])
                nsamples = np.sum(submask)
                #assume that nans are randomly distributed by season
                if nsamples >= 60:
                    d18O = d18O_model.loc[i].values[3:][submask]
                    d2H = d2H_model.loc[i].values[3:][submask]
                    sx = np.ones(len(d18O))*0.1
                    #standard deviations of y vals
                    sy = np.ones(len(d18O))*0.5
                    realdata = odr.RealData(d18O, d2H, sx=sx, sy=sy)
                    model = odr.ODR(realdata, lmwlmod, beta0 = [8.0, 10.0])
                    results = model.run()
                    modslopes[filenamepre[j][:-5]+'_m'].loc[i] = results.beta[0]
                    modslopes[filenamepre[j][:-5]+'_b'].loc[i] = results.beta[1]
                else:
                    modslopes[filenamepre[j][:-5]+'_m'].loc[i] = np.nan
                    modslopes[filenamepre[j][:-5]+'_b'].loc[i] = np.nan                

    modslopes.to_csv(os.path.join('c:Data', 'ModelData_LMWL_ODR.csv'), index = True)
    
    return modslopes

def getkoppengrid(np, koppen):
    codedict = {'A': 0, 'B': 1, 'f': 2, 's':3, 'w': 4, 'D': 5, 'E': 6}
    koppen['clim'] = [i[0] for i in koppen['Cls'].values]
    koppengrid = np.zeros((len(np.unique(koppen['Lon'].values)), len(np.unique(koppen['Lat'].values))))
    koppengrid.fill(np.nan)
    
    for i, latval in zip(range(len(np.unique(koppen['Lat'].values))), np.unique(koppen['Lat'].values)):
        latmatchinds = np.where(koppen['Lat'].values == latval)[0] #this is list index
        #this is grid index, along with i, which is also grid index
        lonmatchinds = [np.where(koppen.loc[j, 'Lon'] == np.unique(koppen['Lon'].values))[0] for j in latmatchinds]
        #print(len(lonmatchinds), len(latmatchinds))
        codes = np.array([j for j in koppen.loc[latmatchinds, 'clim'].values])
        cinds = np.where(codes == 'C')[0]
        if cinds.any():
            codes[cinds] = [j[1] for j in koppen.loc[latmatchinds, 'Cls'].values[cinds]]
        colors = np.array([codedict[j] for j in codes])
        koppengrid[lonmatchinds, i] = colors[:, None]
    
    koppengrid1 = np.ma.masked_where(np.isnan(koppengrid), koppengrid)
    return koppengrid1