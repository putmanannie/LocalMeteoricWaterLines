# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:45:59 2019

@author: u0929173
"""

def Figure1(os, np, sp, mpl, bayesell, nmons, reclen, uniquesites):
    bayesell['union'] = 0
    bayesell['intersect_union'] = 0
    for site in uniquesites:
        siteinds = np.where(bayesell['site'].values == site)[0]
        nmax = np.max(bayesell.loc[siteinds, 'n'])
        print(nmax)
        fullareaind = np.where(bayesell.loc[siteinds, 'n']== nmax)[0]
        fullareaind = fullareaind[0]
            
        fullarea = bayesell.loc[siteinds[fullareaind], 'area']
        subareaind = np.where((bayesell.loc[siteinds, 'n'] > 0)&(bayesell.loc[siteinds, 'n'] < nmax))[0]
        #calculate the union
        bayesell.loc[siteinds[subareaind], 'union'] = np.subtract(bayesell.loc[siteinds[subareaind], 'area'].values, 
            bayesell.loc[siteinds[subareaind], 'intersect_area'].values)+fullarea
        bayesell.loc[siteinds[subareaind], 'intersect_union'] = np.divide(bayesell.loc[siteinds[subareaind], 'intersect_area'].values, 
                bayesell.loc[siteinds[subareaind], 'union'].values)

    #this index gets things that are not the full dataset
    subset = np.where((bayesell['n']>0)&(bayesell['n']<=120))[0]
    quant_all = np.zeros(len(nmons))
    gaussiankde = np.zeros(len(nmons))
    ns = np.arange(0, 1, 0.05)
    
    dataset = []
    for i in range(len(nmons)):
        nsub = np.where(bayesell['n'].values[subset]==nmons[i])[0]
        quant_all[i] = np.percentile(bayesell['intersect_union'].values[subset][nsub], 10)
        dataset.append(bayesell['intersect_union'].values[subset][nsub])
        #get the maximum value of the array
        gaussiankde[i] = sp.stats.gaussian_kde(bayesell['intersect_union'].values[subset][nsub]).evaluate(ns).max()
    
    #scalingfactor so all widths are scaled to their absolute maximum relative to the 
    #maximum value present for all month nubmers tested
    sf = 15/max(gaussiankde)
    mpl.pyplot.figure(figsize=(4.5, 3))
    elarea = mpl.pyplot.violinplot(dataset, nmons, widths = gaussiankde*sf, showmedians = True, showextrema= False)
    ax1 = mpl.pyplot.gca()
    ax2 = ax1.twinx()
    nsites = [len(reclen[reclen>= nmon]) for nmon in nmons]
    ax2.plot(nmons, nsites, linewidth = 1, color = 'k', linestyle = '--', marker = '.')
    
    for ea in elarea['bodies']:
        ea.set_facecolor('#da5353')
        ea.set_edgecolor('k')
        ea.set_linewidths(0.5)
        ea.set_alpha(1)
        
        elarea['cmedians'].set_color('k')
        elarea['cmedians'].set_linewidths(0.5)
    
    
    ax2.spines['top'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_xticks([24, 36, 48, 60, 72, 84, 96, 108, 120])
    mpl.pyplot.xlabel('Substimeseries length (months)')
    ax1.set_ylabel('$A_{Eall} \cap A_{Esub}$ /$A_{Eall} \cup A_{Esub}$')
    ax2.set_ylabel('Number of sites')
    mpl.pyplot.tight_layout()

    mpl.pyplot.savefig(os.path.join('c:Figs', 'Figure1.png'), dpi = 1000)
    
def Figure2(os, np, mpl, Basemap, koppen, koppengrid1, regresults, mask, colordict):
    climdict = {'A: Tropical': 'o', 'B: Arid':'o', 'C: Temperate':'o', 'D: Continental':'o', 'E: Polar':'o', 'N: None':'s'}

    listedcmap2 = mpl.colors.ListedColormap(['#a82424',  '#d74242', '#e06c6c', '#edabab', 
                                             '#bdbdbd', '#737373', '#525252', 
                                              '#252525'], 'segementedcmap', 8) 
                                              
    listedcmap = mpl.colors.ListedColormap(['#e8ddee', '#fff0cc', '#daf1da', 
                                           '#ebf1da', '#daf1f0', '#dae1f1', 
                                           '#e6e6e6'], 'koppenmap', 7)  
    
    fig, ax = mpl.pyplot.subplots(nrows = 2, ncols = 1, figsize = (7, 6))
    mpl.pyplot.sca(ax[0])
    #these are just going to draw on the current axis
    m1 = Basemap(projection = 'robin', lon_0 = 0, resolution = 'c')
    #m1.drawcountries(linewidth = 0.35)
    m1.drawcoastlines(linewidth = 0.5, color = '#404040', zorder = 2)
    #m1.drawmapboundary(fill_color = None)
    #m1.fillcontinents(color = None, zorder = 0)
    xgrd, ygrd = np.meshgrid(np.unique(koppen['Lon'].values), np.unique(koppen['Lat'].values))
    
    m1.pcolormesh(xgrd,ygrd, koppengrid1.T, cmap = listedcmap, latlon = True, zorder = 1)
    plist = []
    for key in np.sort(climdict.keys()):
        subinds = np.where(regresults['clim'].values[mask] == key[0])[0]
        x, y = m1(regresults['lon'].values, regresults['lat'].values)
        colorarray = np.round(regresults['m'].values[mask][subinds]*2.0, 0)/2
        colorarray = np.where(colorarray < 6.5, 6, colorarray)
        colorarray = np.where(colorarray > 9.5, 10, colorarray)
        p = m1.scatter(x[mask][subinds], y[mask][subinds], c = colorarray, 
                   s = 10, vmin = 6, vmax = 10, marker = climdict[key], 
                   linewidths = 0.5, edgecolors = 'k', cmap = listedcmap2, zorder = 10)
        plist.append(p)
    
    cbar = mpl.pyplot.colorbar(orientation = 'vertical',shrink = 0.90)
    cbar.set_label('slope')
    ax[0].annotate('(a)', xy = (0, 1), xytext = (0.02, 0.98), xycoords = 'axes fraction')
    #ax_legend = fig.add_axes([0.35, 0.95, 0.3, 0.03], zorder=3)
    
    mpl.pyplot.sca(ax[1])
    #these are just going to draw on the current axis
    m1 = Basemap(projection = 'robin', lon_0 = 0, resolution = 'c')
    m1.drawcoastlines(linewidth = 0.35, color = '#404040', zorder = 2)
    #m1.fillcontinents(color = '#dadada', zorder = 0)
    m1.pcolormesh(xgrd,ygrd, koppengrid1.T, cmap = listedcmap, latlon = True, zorder = 1)
    for key in climdict.keys():
        subinds = np.where(regresults['clim'].values[mask] == key[0])[0]    
        x, y = m1(regresults['lon'].values, regresults['lat'].values)
        colorarray = np.round(regresults['m'].values[mask][subinds]*2.0, 0)/2
        colorarray = np.where(colorarray < -10, -10, colorarray)
        colorarray = np.where(colorarray > 25, 25, colorarray)
        m1.scatter(x[mask][subinds], y[mask][subinds], c = regresults['b'].values[mask][subinds], 
                   s = 10, vmin = -10, vmax = 30, marker = climdict[key], 
                   linewidths = 0.5, edgecolors = 'k', cmap = listedcmap2, zorder = 10)
    #make the outline of significant regressions bold
    #m1.scatter(x[mask], y[mask], c = regresults['b'].values[mask], s = 10, linewidth = 0.5, vmin = -30, vmax = 20, edgecolor = 'k')
    cbar = mpl.pyplot.colorbar(orientation = 'vertical',shrink = 0.90)
    cbar.set_label('intercept')
    ax[1].annotate('(b)', xy = (0, 1), xytext = (0.02, 0.98), xycoords = 'axes fraction')
    #ax[0].legend(plist[-2:], ['Koppen class assigned', 'Not assigned'], loc = (0.05, 1.05), ncol = 2) 
    indorder = np.argsort(colordict.keys())
    markers = [mpl.pyplot.Line2D([0,0],[0,0],color=colordict.values()[indo], marker='s', linestyle='') for indo in indorder] 
    ax[0].legend(np.array(markers), np.array(sorted(colordict.keys())), 
                      numpoints=1, loc = 'lower left', ncol = 4, fontsize = 6, 
                      bbox_to_anchor=(0.05, 1.02), frameon = False)  

    mpl.pyplot.savefig(os.path.join('c:Figs', 'Figure2.png'), dpi = 1000)

def Figure3(os, np, mpl, Basemap, Rectangle, regresults, mask):
    fig, ax = mpl.pyplot.subplots(nrows = 2, ncols = 1, figsize = (6, 6), sharex = True)
    
    colordict = {'A: Tropical \n':'#693f7b', 'B: Arid \n':'#fbae00', 'D: Continental\n':'#5779c1', 
                 'Cf: Temperate \n humid':'#379537', 'Cs: Temperate \n dry summer':'#9dba45', 'Cw: Temperate \n dry winter': '#338984', 
                 'E: Polar \n and alpine':'dimgray', 'N: No \n classification':'k'} 
    colordict2 = {'A':'#693f7b', 'B':'#fbae00', 'D':'#5779c1', 'E':'dimgray', 'N':'k'}
    colordict3 = {'f':'#379537', 's':'#9dba45', 'w': '#338984'} 
    inc = 10
    latbins = np.arange(-90.0, 90.0+inc, inc)
    #this part is just for plotting
    latctrs = np.array([(latbins[i]+latbins[i+1])/2 for i in range(len(latbins)-1)])
    
    climclass = np.unique(regresults['clim'].values[mask])

    indorder = np.argsort(colordict.keys())
    markers = [mpl.pyplot.Line2D([0,0],[0,0],color=colordict.values()[indo], marker='o', linestyle='') for indo in indorder]
    
    #offset = {'A':-4, 'B':-2, 'C':0, 'D':2, 'E': 4, 'N':2}
    offset = np.arange(-5, 5, 1.666666667)
    width = 1.666666667
    #get all of the rows that are in each lat bin in a loop
    
    for l in range(2):
        if l == 0:
            key = 'm'
            default = 0.2
        else:
            key = 'b'
            default = 2

        for i in range(len(latctrs)):
            latbininds = regresults.index[mask][(regresults['lat'].values[mask] >= latbins[i]) & (regresults['lat'].values[mask] <= latbins[i+1])]
            #group them by koppen class for histogram.
            #if less than 3, no histogram, just plot
            #just scatterplot the N ones
            posind = 0
            for j in climclass:            
                subinds = latbininds[regresults['clim'][latbininds] == j]
                if len(subinds) > 0:
                    if j == 'C':
                        subclass = np.array([regresults['koppen'][k][1] for k in subinds])
                        
                        for k in np.unique(subclass):
                            fclass = regresults[key][subinds].values[subclass == k]                        
                            position = latctrs[i]+offset[posind]
                            llc = (position, min(fclass))
                            hgt = max(fclass) - min(fclass)
                            if hgt < default:
                                hgt = default
                            
                            pat = Rectangle(xy=llc, width = width, height = hgt, 
                               fill = True, facecolor = colordict3[k], zorder = 5)
                            ax[l].add_patch(pat)
                            posind = posind+1 
                       
                    else:                
                        position = latctrs[i]+offset[posind]
                        llc = (position, min(regresults[key][subinds].values))
                        hgt = max(regresults[key][subinds].values) - min(regresults[key][subinds].values)
                        if hgt < default:
                                hgt = default
                        pat = Rectangle(xy = llc, width = width, height = hgt, 
                                fill = True, facecolor = colordict2[j], zorder = 5)
                        ax[l].add_patch(pat)        
                        posind = posind+1 
    
        ax[l].set_xlim([-90, 90])       
    #    ax[l].spines['right'].set_visible(False)
        ax[l].spines['top'].set_visible(False)
        
       
    ax[0].set_ylim([4, 10]) 
    ax[1].set_ylim([-20, 30])    
    ax[0].set_xticks([])    
    ax[1].set_xticks(np.arange(-90, 100, 10))
    ax[1].set_xticklabels(np.arange(-90, 100, 10))
    ax[0].xaxis.grid(True)
    ax[1].xaxis.grid(True) 
    ax[1].set_xlabel('Latitude')
    ax[0].set_ylabel('Slope')
    ax[1].set_ylabel('Intercept')
    ax[0].annotate('(a)', xy = (0, 1), xytext = (0.01, 0.95), xycoords = 'axes fraction')
    ax[1].annotate('(b)', xy = (0, 1), xytext = (0.01, 0.95), xycoords = 'axes fraction')
    ax[0].legend(np.array(markers), np.array(sorted(colordict.keys())), numpoints=1, loc = 'upper left', 
        ncol = 4, fontsize = 8, bbox_to_anchor=(-0.05, 1.35), frameon = False)

    mpl.pyplot.savefig(os.path.join('c:Figs', 'Figure3.png'), dpi = 1000)
    
def Figure4(os, np, mpl, Basemap, joined):

    Asubset = joined[(joined['clim_caller'] == 'A')&(joined['proj'] != '00076')].index     
    Bsubset = joined[(joined['clim_caller'] == 'B')&(joined['proj'] != '00076')].index 
    Cssubset = joined[((joined['koppen_caller'] == 'Csa')|(joined['koppen_caller'] == 'Csb'))&(joined['proj'] != '00076')].index
    Dsubset = joined[(joined['clim_caller'] == 'D')&(joined['proj'] != '00076')].index
    
    warmseason_dex = np.where(joined.loc[Dsubset, 'lat'] < 0, joined.loc[Dsubset, 'DJF_dex'], joined.loc[Dsubset, 'JJA_dex'])
    coldseason_dex = np.where(joined.loc[Dsubset, 'lat'] > 0, joined.loc[Dsubset, 'DJF_dex'], joined.loc[Dsubset, 'JJA_dex'])
    slope = joined.loc[Dsubset, 'm']
    warmseason_dex_e = np.where(joined.loc[Dsubset, 'lat'] < 0, joined.loc[Dsubset, 'DJF_dex_s'], joined.loc[Dsubset, 'JJA_dex_s'])
    coldseason_dex_e = np.where(joined.loc[Dsubset, 'lat'] > 0, joined.loc[Dsubset, 'DJF_dex_s'], joined.loc[Dsubset, 'JJA_dex_s'])
      
    xvals1 = np.where(joined.loc[Bsubset, 'lat'] > 0, joined.loc[Bsubset, 'JJA_dex'], joined.loc[Bsubset, 'DJF_dex'])
    xvals1i = np.where(joined.loc[Bsubset, 'lat'] < 0, joined.loc[Bsubset, 'JJA_dex'], joined.loc[Bsubset, 'DJF_dex'])
    yvals1 = joined.loc[Bsubset, 'm']
    evals1 = np.where(joined.loc[Bsubset, 'lat'] > 0, joined.loc[Bsubset, 'JJA_dex_s'], joined.loc[Bsubset, 'DJF_dex_s'])
    evals1i = np.where(joined.loc[Bsubset, 'lat'] < 0, joined.loc[Bsubset, 'JJA_dex_s'], joined.loc[Bsubset, 'DJF_dex_s'])
    
    xvals2 = np.where(joined.loc[Cssubset, 'lat'] > 0, joined.loc[Cssubset, 'JJA_dex'], joined.loc[Cssubset, 'DJF_dex'])
    xvals2i = np.where(joined.loc[Cssubset, 'lat'] < 0, joined.loc[Cssubset, 'JJA_dex'], joined.loc[Cssubset, 'DJF_dex'])
    yvals2 = joined.loc[Cssubset, 'm']
    evals2 = np.where(joined.loc[Cssubset, 'lat'] > 0, joined.loc[Cssubset, 'JJA_dex_s'], joined.loc[Cssubset, 'DJF_dex_s'])
    evals2i = np.where(joined.loc[Cssubset, 'lat'] < 0, joined.loc[Cssubset, 'JJA_dex_s'], joined.loc[Cssubset, 'DJF_dex_s'])
    
    xvals3 = np.where(joined.loc[Asubset, 'lat'] > 0, joined.loc[Asubset, 'JJA_dex'], joined.loc[Asubset, 'DJF_dex'])
    xvals3i = np.where(joined.loc[Asubset, 'lat'] < 0, joined.loc[Asubset, 'JJA_dex'], joined.loc[Asubset, 'DJF_dex'])
    yvals3 = joined.loc[Asubset, 'm']
    evals3 = np.where(joined.loc[Asubset, 'lat'] > 0, joined.loc[Asubset, 'JJA_dex_s'], joined.loc[Asubset, 'DJF_dex_s'])
    evals3i = np.where(joined.loc[Asubset, 'lat'] < 0, joined.loc[Asubset, 'JJA_dex_s'], joined.loc[Asubset, 'DJF_dex_s'])
    
    lims = np.array([[-1, -4, -14], [9, 5, 9]])
    
    symbolsize = 10
    fig, ax = mpl.pyplot.subplots(2, 2, figsize = (6,7), sharey = True)
    
    h1 = ax[0, 0].scatter(xvals2, yvals2, s = symbolsize, c = '#9dba45', zorder = 1, label = 'Temperate, hot-dry summer')
    ax[0, 0].errorbar(np.concatenate((xvals1, xvals2, xvals3)), 
                        np.concatenate((yvals1, yvals2, yvals3)), 
                        xerr = np.concatenate((evals1, evals2, evals3)), zorder = -2, fmt = '.', c = '#e6e6e6', label = 'standard deviation')
    h3 = ax[0, 0].scatter(xvals1, yvals1, s = symbolsize, c = '#fbae00', zorder = 2, label = 'Arid')
    h4 = ax[0, 0].scatter(xvals3, yvals3, s = symbolsize, c = '#693f7b', zorder = 3, label = 'Tropical')
   
    ax[0, 0].fill_between(lims[:, 0], y1= [4.5, 4.5], y2= [9.5, 9.5], color = '#fbae00', alpha = 0.1, zorder = -1)
    ax[0, 0].fill_between(lims[:, 1], y1= [4.5, 4.5], y2= [9.5, 9.5], color = '#fbae00', alpha = 0.1, zorder = -1)
    ax[0, 0].fill_between(lims[:, 2], y1= [4.5, 4.5], y2= [9.5, 9.5], color = '#fbae00', alpha = 0.1, zorder = -1)
    
    ax[0, 1].scatter(xvals2i, yvals2, s = symbolsize, c = '#9dba45', zorder = 1, label = 'Temperate, hot-dry summer')
    ax[0, 1].errorbar(np.concatenate((xvals1i, xvals2i, xvals3i)), 
                        np.concatenate((yvals1, yvals2, yvals3)), 
                        xerr = np.concatenate((evals1i, evals2i, evals3i)), zorder = -2, fmt = '.', c = '#e6e6e6', label = 'standard deviation')
    ax[0, 1].scatter(xvals1i, yvals1, s = symbolsize, c = '#fbae00', zorder = 2, label = 'Arid')
    ax[0, 1].scatter(xvals3i, yvals3, s = symbolsize, c = '#693f7b', zorder = 3, label = 'Tropical')
    
    h5 = ax[1,0].scatter(warmseason_dex, slope, s = symbolsize, color = '#5779c1', label = 'Continental')
    ax[1,0].errorbar(warmseason_dex, slope, xerr = warmseason_dex_e, fmt = '.', color = '#e6e6e6', zorder = -2)

    ax[1,0].fill_between(lims[:, 0], y1= [4.5, 4.5], y2= [10, 10], color = '#fbae00', alpha = 0.1, zorder = -1)
    ax[1,0].fill_between(lims[:, 1], y1= [4.5, 4.5], y2= [10, 10], color = '#fbae00', alpha = 0.1, zorder = -1)
    ax[1,0].fill_between(lims[:, 2], y1= [4.5, 4.5], y2= [10, 10], color = '#fbae00', alpha = 0.1, zorder = -1) 

    ax[1,1].scatter(coldseason_dex, slope, s = symbolsize, color = '#5779c1')
    ax[1,1].errorbar(coldseason_dex, slope, xerr = coldseason_dex_e, fmt = '.', color = '#e6e6e6', zorder = -2)
    ax[1,1].fill_between([15, 30], y1= [4.5, 4.5], y2= [10, 10], color = '#5779c1', alpha = 0.1, zorder = -1)
    ax[1,1].fill_between([0, 10], y1= [4.5, 4.5], y2= [10, 10], color = '#5779c1', alpha = 0.1, zorder = -1)

    ax[1,1].annotate('Slope', xy = (0.1, 0.5), xytext = (0.02, 0.5), xycoords = 'figure fraction', rotation = 'vertical')
    ax[0,0].annotate(r'(a) ', xy = (0.1, 0.9), xytext = (0.05, 0.9), xycoords = 'axes fraction')
    ax[0,1].annotate(r'(b) ', xy = (0.1, 0.9), xytext = (0.05, 0.9), xycoords = 'axes fraction')  
    ax[1,1].annotate(r'(d) ', xy = (0.1, 0.9), xytext = (0.05, 0.9), xycoords = 'axes fraction')
    ax[1,0].annotate(r'(c) ', xy = (0.1, 0.9), xytext = (0.05, 0.9), xycoords = 'axes fraction')  
    ax[1,1].annotate('snow from \nmixed-phase \nclouds', xy = (0, 4.5), xytext = (0, 5), xycoords = 'data')
    ax[1,1].annotate('lake effect \nsnow', xy = (15, 4.5), xytext = (14, 5), xycoords = 'data')
    ax[1, 0].annotate('evaporated rain', xy = (-14, 4.5), xytext = (-13, 5), xycoords = 'data')
    ax[1, 0].set_xlabel('Mean warm season $d$ ('+u'\u2030'+')')
    ax[1, 1].set_xlabel('Mean cool season $d$ ('+u'\u2030'+')')
    ax[1,0].set_ylim([4.5, 9.5])
    ax[1,1].set_ylim([4.5, 9.5])
    ax[0,0].set_xlim([-15, 25])
    ax[1,0].set_xlim([-15, 25])
    ax[0,1].set_xlim([-5, 25])
    ax[1,1].set_xlim([-5, 25])
    mpl.pyplot.figlegend(handles = [h1, h3, h4, h5], 
                         labels = ['Temperate, hot-dry summer', 'Arid', 
                                   'Tropical', 'Continental'], 
                         bbox_to_anchor = [0.8, 0.98], ncol = 2, frameon = False)
    mpl.pyplot.savefig(os.path.join('c:Figs', 'Figure4.png'), dpi = 500)


def Figure6(os, np, mpl, modslopes, regresults, mask):
    filenamepre =  ['cam2_1x1_',
     'cam5_1x1_',
     'echam5_1x1_',
     'gissE_1x1_',
     'HadAM3_1x1_',
     'isoGSM_1x1_',
     'lmdz_1x1_',
     'miroc_1x1_']
    
    colnames = ['dimgray', 'rosybrown', 'maroon', 'peru', 'goldenrod', 'olive', 'darkolivegreen', 'cadetblue', 'midnightblue', ]
    maskdict = {0:('B', 'B: Arid'), 1:('A', 'A: Tropical'), 2:('D', 'D: Continental')}
    keys={0:('m', '_m'), 1:('b', '_b')}
    letters = [['(a)', '(b)'], ['(c)','(d)'], ['(e)', '(f)']]
    
    fig, ax = mpl.pyplot.subplots(nrows = 3, ncols = 2, figsize = (6, 7), sharex = True)
    
    #count through rows
    for i in range(len(maskdict)):
        climmask = regresults['clim'][mask] == maskdict[i][0]
        for j in range(len(keys)):
            data = regresults[keys[j][0]][mask][climmask].values[:, None]
    
            for k in range(len(filenamepre)):
                y = modslopes[filenamepre[k][:-5]+keys[j][1]][climmask].values[:, None]
                data = np.concatenate((data, y), axis = 1)
            data = np.array(list(data), dtype = np.float)
    
            for k in range(9):
                flierprops = dict(marker='.', markerfacecolor=colnames[k], markersize=6,
                              linestyle='none', markeredgecolor = 'none')
                bplot = ax[i, j].boxplot(data[:, k][~np.isnan(data[:, k])], positions = [k+1], 
                                   patch_artist = True, vert = True, widths = 0.5, flierprops = flierprops)
                bplot['boxes'][0].set_color(colnames[k])
                bplot['medians'][0].set_color('k')
                
            ax[i, j].annotate(letters[i][j], xy = (0, 1), xytext = (0.02, 0.92), xycoords = 'axes fraction')
            ax[i, j].set_xlim([0, 10])
            ax[i, j].spines['right'].set_visible(False)
            ax[i, j].spines['top'].set_visible(False)
            if j == 0:
                ax[i, j].set_ylim([5, 10])
                ax[i, j].set_ylabel('{}'.format(maskdict[i][1]))
            else:
                ax[i, j].set_ylim([-20, 30])
                #ax[i, j].set_ylabel('Intercept - {}'.format(maskdict[i][1]))
            
    fig.text(0.25, 0.98, 'Slope', ha='center')  
    fig.text(0.75, 0.98, 'Intercept', ha='center')
                   
    ax[2, 0].set_xticks(range(1,10,1))
    ax[2, 1].set_xticks(range(1,10,1))
    legendnames = [fname[:-5] for fname in filenamepre]
    legendnames.insert(0, 'data')
    
    ax[2, 0].set_xticklabels(legendnames, rotation = 70)
    ax[2, 1].set_xticklabels(legendnames, rotation = 70)
    
    fig.text(0.5, 0.02, 'Data Source', ha='center')
    mpl.pyplot.tight_layout(pad = 1.1, w_pad = 1, h_pad = 2.0)
    
    mpl.pyplot.savefig(os.path.join('c:Figs', 'Figure6.png'), dpi = 1000)