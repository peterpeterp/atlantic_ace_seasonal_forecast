#%% Imports and environment setup
import os,sys,importlib

sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

sys.path.append(base_path+'seasonal_indicators')
import TC_support ;  TC_support = importlib.reload(TC_support)
cat_colors={0:'lightblue',1:'#ffffcc',2:'#ffe775',3:'#ffc148',4:'#ff8f20',5:'#ff6060'}

result_path = ''
cutoff = [np.datetime64('1979-01-01'),np.datetime64('2019-01-01')]


#%% Load data
csu = pd.read_csv('data/csu_forecast-verifications_all.csv', sep=';', header=1, usecols=range(36,41))
csu_ace_Obs = xa.DataArray(np.array(csu['Observed.6']), coords=dict(time=np.array(csu['Year.6'])), dims=['time'])

nc_TCs = xa.open_dataset(data_path+'/IBTrACS.NA.v04r00.nc')


#%% open an example file to get the time axis
nc_example = xa.open_dataset(data_path+'/ERA5/ERA5_MSLP-SST_1979-2018_grid1x1.nc')
time_axis = nc_example['time']
time_axis_readable = np.array([int(str(dd).split('-')[0]) + (int(str(dd).split('-')[1]))/13. for dd in nc_example['time'].values])
index = xa.DataArray(np.zeros(nc_example['time'].shape), coords=dict(time=nc_example['time']), dims=['time'])

#%% define a polygon in the atlantic and find all TC days within this polygon
reg_atl=Polygon([[-98.44,19.48],[-92.29,16.13],[-81.74,8.41],[-74.18,5.62],[-52.21,-0.53],[-9.67,0.35],[-17.4,30.9],[-98.44,36.32],[-98.44,32.25],[-98.44,19.48]])
reg_atl_path = matplotlib.path.Path([[-98.44,19.48],[-92.29,16.13],[-81.74,8.41],[-74.18,5.62],[-52.21,-0.53],[-9.67,0.35],[-17.4,30.9],[-98.44,36.32],[-98.44,32.25],[-98.44,19.48]])
notna = np.isfinite(nc_TCs['wmo_wind']) & (nc_TCs['time'] > index.time[0])
in_reg = reg_atl_path.contains_points(np.hstack((nc_TCs['lon'].values.flatten()[:,np.newaxis],nc_TCs['lat'].values.flatten()[:,np.newaxis]))).reshape(notna.shape)

#%% calculate the accumulated cyclone energy of TCs within the polygon from above
for wind,date in zip(nc_TCs['wmo_wind'].values[(notna) & (in_reg)], nc_TCs['time'].values[(notna) & (in_reg)]):
    if wind>34 and date>index.time[0]:
        id_ = index.time[(index.time < date)][-1]
        index.loc[id_] += wind**2 * 0.0001

#%% plot some things
plt.close('all')
with PdfPages(result_path+index_name+'.pdf') as pdf:
    #%% region
    fig,ax = plt.subplots(nrows=1,figsize=(3,3),subplot_kw={'projection': ccrs.Orthographic(central_longitude=-50, central_latitude=30)})
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, facecolor='gray'); ax.add_feature(cartopy.feature.OCEAN,facecolor='w')
    ax.add_geometries([reg_atl], ccrs.PlateCarree(), edgecolor='lightgreen',alpha=0.5,facecolor='lightgreen', color='none',linewidth=3,zorder=100)
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

    #%% region and storm genesis locations
    fig,ax = plt.subplots(nrows=1,figsize=(5,2.5),subplot_kw={'projection': ccrs.PlateCarree()})
    ax.coastlines(); ax.set_extent([-100,0,0,40],crs=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, facecolor='gray'); ax.add_feature(cartopy.feature.OCEAN,facecolor='black')
    ax.add_geometries([reg_atl], ccrs.PlateCarree(), edgecolor='lightgreen',alpha=1,facecolor='none', color='none',linewidth=3,zorder=100)
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

    #%% time series of the index
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,2.5))
    ax.plot(time_axis_readable,index.values,color='darkcyan',marker='.',linestyle='',alpha=0.2)
    ax.plot(np.nanmean(time_axis_readable.reshape((int(len(index.time)/12),12)),axis=-1),np.nansum(index.values.reshape((int(len(index.time)/12),12)),axis=-1),color='darkcyan')
    ax.plot(csu_ace_Obs.time.values+0.6, csu_ace_Obs.values, 'k:')
    ax.set_ylabel(index_longname,color='k')
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

    #%% time series of the index
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,2.5))
    ax.plot(time_axis_readable,index.values,color='darkcyan',marker='.',linestyle='',alpha=0.2)
    ax.plot(np.nanmean(time_axis_readable.reshape((int(len(index.time)/12),12)),axis=-1),np.nansum(index.values.reshape((int(len(index.time)/12),12)),axis=-1),color='darkcyan')
    ax.set_ylabel(index_longname,color='k')
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

    # select, detrend etc
    index_raw=index[(index.time >= cutoff[0]) & (index.time < cutoff[1])]

    index_anom = index_raw.copy()
    # detrend in each time-step
    for j in range(time_cycle):
            index_anom.values[j::12] = scipy.signal.detrend(index_anom.values[j::12], axis = 0)

    # for i in range(12):
    #       print(i,np.sum(index_raw.values[i::12]))

    index_ASO_raw = index_raw.copy() *0
    for i in range(4):
            index_ASO_raw.values[6::12] += index_raw.values[6+i::12]

    index_ASO_anom = index_ASO_raw.copy()
    # detrend in each time-step
    for j in range(time_cycle):
            index_ASO_anom.values[j::time_cycle] = scipy.signal.detrend(index_ASO_anom.values[j::time_cycle], axis = 0)

    #%% time series of the index
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,2.5))
    ax.plot(time_axis_readable[6::12],index_ASO_raw.values[6::12],color='darkcyan')
    ax.plot(time_axis_readable[6::12],index_ASO_anom.values[6::12],color='darkmagenta')
    ax.set_ylabel(index_longname,color='k')
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

    #%% time series of the index
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,2.5))
    ax.plot(np.nanmean(time_axis_readable.reshape((int(len(time_axis_readable)/12),12)),axis=-1),np.nansum(index_raw.values.reshape((int(len(time_axis_readable)/12),12)),axis=-1),color='darkcyan')
    ax.plot(csu_ace_Obs.time+0.6, csu_ace_Obs, 'k:')
    ax.plot(time_axis_readable[6::12],index_ASO_raw.values[6::12],color='darkcyan')
    ax.plot(time_axis_readable[6::12],index_ASO_anom.values[6::12],color='darkmagenta')
    ax.set_ylabel(index_longname,color='k')
    plt.tight_layout(); pdf.savefig(transparent=True); plt.close()


#%% save results
xa.Dataset({index_name+'_allMon_raw':index_raw,
                        index_name+'_allMon':index_anom,
                        index_name+'_raw':index_ASO_raw,
                        index_name:index_ASO_anom,
                        }).to_netcdf(result_path+index_name+'.nc')
