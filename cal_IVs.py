#============================================#
# functions for internal variability
#============================================#
import numpy as np
import scipy.stats
import xarray as xr
from eofs.xarray import Eof
import cartopy.crs as ccrs
import matplotlib.pyplot as plt



#===== func: calculate ENSO index (ONI) =====#
def cal_nino(sst, anom_method, index="ONI", standardize=False):
    # calculate anomalies
    #anom_method: 1 for ens mean; 0 for single run
    
    #anomaly
    if anom_method == 0:
        sst_a = sst.groupby("time.month") - sst.groupby("time.month").mean("time")
    elif anom_method ==1: 
        sst_a = sst - sst.mean("ens")

    #standardize
    if standardize == True:
        sst_a = sst_a.groupby("time.month") / sst_a.groupby("time.month").std("time")    
        
    #pick key region
    if index == "ONI" or index == "nino3.4":
        reg =[-5, 5, 190, 240] # 5N-5S; 170W-120W
        
    sst_a = sst_a.sel(lat=slice(reg[0], reg[1]), lon=slice(reg[2],reg[3]))
    weights = np.cos(np.deg2rad(sst_a.lat))
    weights.name = "weights"
    nino = sst_a.weighted(weights).mean(("lat", "lon"))
    
    # smooth
    if index == "ONI":
        nino_smooth = nino.rolling(time=3, center=True).mean()
    return nino_smooth

#===== func: calculate PDO index (TPI) =====#
def cal_tpi(sst, anom_method, standardize=False):
    # calculate anomalies
    #anom_method: 1 for ens mean; 0 for single run
    
    #anomaly
    if anom_method == 0:
        sst_a = sst.groupby("time.month") - sst.groupby("time.month").mean("time")
    elif anom_method ==1: 
        sst_a = sst - sst.mean("ens")

    # pick key regions
    t1 = sst_a.sel(lat=slice(25, 45), lon=slice(140, 215))
    t2 = sst_a.sel(lat=slice(-10, 10), lon=slice(170, 270))
    t3 = sst_a.sel(lat=slice(-50, 15), lon=slice(150, 200))
    
    def gmean(t):
        weights = np.cos(np.deg2rad(t.lat))
        weights.name = "weights"
        tmean = t.weighted(weights).mean(("lat", "lon"))
        return tmean
    
    sst1 = gmean(t1)
    sst2 = gmean(t2)
    sst3 = gmean(t3)
    
    tpi = sst2 - (sst1 + sst3) /2.0
    
    # smooth
    wgt = np.array([1,6,19,42,71,96,106,96,71,42,19,6,1])
    wgt = wgt / np.sum(wgt)
    
    tpi_smooth = tpi.copy()
    for i in range(10):
        tpi_smooth[i,:] = np.convolve(tpi[i,:], wgt[::-1], 'same')
        
    return tpi_smooth

#===== func: calculate NAO =====#
def cal_nao(psl0, standardize = False):
    
    psl = psl0.copy()

    # convert longitude from [0, 360] to [-180,180]
    if psl.lon.max()> 200:
        print("converting longitude...")
        psl['_lon_adjusted'] = xr.where(psl.lon > 180, psl.lon - 360, psl.lon)
        
        psl = (
            psl
            .swap_dims({'lon': '_lon_adjusted'})
            .sel(**{'_lon_adjusted': sorted(psl._lon_adjusted)})
            .drop('lon'))

        psl = psl.rename({'_lon_adjusted': 'lon'})
    
    # select NAO region + select JFM 
    psl = psl.sel(time = np.in1d( psl['time.month'], [1,2,3]), lat = slice(20, 90), lon = slice(-90, 40)).groupby('time.year').mean('time')
    psl = psl.rename({'year':'time'}).transpose('time', ...)
    
    # anomalies
    psl_a = psl - psl.mean("ens")
    
    # EOF
    coslat = np.cos(np.deg2rad(psl_a.lat)).clip(0., 1.)
    wgts = np.sqrt(coslat)

    nao = xr.DataArray(
        dims=["ens", "time"],
        coords=dict(
            ens = np.arange(10),
            time= psl_a.time,
        ),
        attrs=dict(
            long_name = 'NAO index'
        ),
    )
    
    eof = xr.DataArray(
        dims=["ens", "lat", "lon"],
        coords=dict(
            ens = np.arange(10),
            lat = psl_a.lat,
            lon = psl_a.lon,
        ),
        attrs=dict(
            long_name = 'NAO pattern'
        ),
    )    
    
    
    for i in range(10):
        temp = psl_a[:, i, :, :]
        solver = Eof(temp)
        nao[i,:] = solver.pcs(npcs = 1, pcscaling=1).squeeze('mode')
        eof[i, :, :] = solver.eofsAsCovariance(neofs=1, pcscaling=1).squeeze('mode')
        
        
    # ensure lower lats get negative pattern    
    for i in range(10):
        temp = eof[i,]
        if temp.sel(lat = 65, lon = -22, method = 'nearest') > 0:
            eof[i, :, :] = eof[i, :, :]*-1.0
            nao[i, : ] = nao[i, :]*-1.0 
    
    # Plot the leading EOF expressed as covariance in the European/Atlantic domain.
    fig, axs = plt.subplots(nrows=5,ncols=2,
                            subplot_kw={'projection': ccrs.Orthographic(central_longitude=-20, central_latitude=60)},
                            figsize=(20,30))
    axs = axs.flatten()
    
    for i in range(10):
        cs = axs[i].contourf(eof.lon, eof.lat, eof[i,], 
                            cmap=plt.cm.RdBu_r, transform=ccrs.PlateCarree())
        axs[i].coastlines()
        
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    return nao, eof

#===== func: calculate AMO =====#
def cal_amo(sst, anom_method, standardize=False):
    # calculate anomalies
    #anom_method: 1 for ens mean; 0 for single run
    
    #anomaly
    if anom_method == 0:
        sst_a = sst.groupby("time.month") - sst.groupby("time.month").mean("time")
    elif anom_method ==1: 
        sst_a = sst - sst.mean("ens")

    # pick key regions
    amo = sst_a.sel(lat=slice(0, 60), lon=slice(280, 360)).mean(("lat", "lon"))
    amo_a = amo - amo.mean("ens")
    
    # smooth
    wgt = np.array([1,6,19,42,71,96,106,96,71,42,19,6,1])
    wgt = wgt / np.sum(wgt)
    
    amo_smooth = amo_a.copy()
    for i in range(10):
        amo_smooth[i,:] = np.convolve(amo_a[i,:], wgt[::-1], 'same')
        
    return amo_smooth
