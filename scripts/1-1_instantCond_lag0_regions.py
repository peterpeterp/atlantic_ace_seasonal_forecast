# -*- coding: utf-8 -*-
import os,sys,importlib, gc

# load all required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# initialize folders
result_path = 'instantCond/'
if os.path.isdir(result_path)==False:
	os.system('mkdir '+result_path)

def init_output():
	'''
	This function is only called once to initialize the xarray's in which results are stored
	'''
	coords = dict(lat = np.asarray(lat_grid), lon = np.asarray(lon_grid), lag=np.arange(lag,lag+1,1, dtype=np.short), field=sorted(precursor_dict.keys()), ID = np.arange(1,n_keep_save+1,dtype=np.short), time= var_anoms.time)
	empty = np.zeros([len(lat_grid),len(lon_grid), 1, len(sorted(precursor_dict.keys())), len(np.arange(1,n_keep_save+1,dtype=np.short)), len(var_anoms.time)]) * np.nan
	out_fld = xa.Dataset({
		'corr': xa.DataArray(data=empty.copy()[:,:,:,:,0,0], coords={k: coords[k] for k in ('lat','lon','lag','field')}, dims=['lat','lon','lag','field']),
		'corr_sig': xa.DataArray(data=empty.copy()[:,:,:,:,0,0],coords={k: coords[k] for k in ('lat','lon','lag','field')}, dims=['lat','lon','lag','field']),
		'region': xa.DataArray(data=empty.copy()[:,:,:,:,:,0],coords={k: coords[k] for k in ('lat','lon','lag','field', 'ID')}, dims=['lat','lon','lag','field','ID']),
		'clusters': xa.DataArray(data=empty.copy()[:,:,:,:,0,0],coords={k: coords[k] for k in ('lat','lon','lag','field')}, dims=['lat','lon','lag','field']),
		})
	out_ts = xa.Dataset({
		'actor': xa.DataArray(data=empty.copy()[0,0,:,:,:,:], coords={k: coords[k] for k in ('lag','field','ID','time')}, dims=['lag','field','ID','time']),
		'actor_raw': xa.DataArray(data=empty.copy()[0,0,:,:,:,:], coords={k: coords[k] for k in ('lag','field','ID','time')}, dims=['lag','field','ID','time']),
		})
	return out_fld,out_ts

for identifier in ["ERA5_setAll_instant_lag0_vLocal"]:	#"param_lag0",
	# load some parameter settings stored in scripts/params/...
	# all varaibles that are required in this script are defined there
	sys.path.append('scripts/params/'+'_'.join(identifier.split('_')[:1]))
	exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

	first = True

	# load the target index (result of 0_index_ACE.py)
	# 1D monthly detrended anomalies
	target =  xa.open_dataset(index_name+'.nc')[index_name+'_allMon'].loc[cutoff[0]:cutoff[1]].squeeze()

	# go through the precursor_dict and identify regions
	for name,precursor in precursor_dict.items():
		# load variable
		var_in = Dataset(precursor['file'], 'r')
		if 'level' in precursor.keys():
			lat_grid, lon_grid, var_vals, var_anoms, var_standard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff,time_cycle,level=precursor['level'])
		else:
			lat_grid, lon_grid, var_vals, var_anoms, var_standard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff,time_cycle)

		# compute pointwise correlation
		corr, corr_sig = get_correlation(target[RV_indices], var_anoms[RV_indices])

		# redefine significant grid-cells
		sig = corr_sig > np.nanpercentile(corr_sig, percent_sig)

		# optional: only consider grid-cells within a region
		if inReg:
			reg_atl=Polygon([[-98.44,19.48],[-92.29,16.13],[-81.74,8.41],[-74.18,5.62],[-52.21,-0.53],[-9.67,0.35],[-17.4,30.9],[-98.44,36.32],[-98.44,32.25],[-98.44,19.48]])
			reg_atl_path = matplotlib.path.Path([[-98.44,19.48],[-92.29,16.13],[-81.74,8.41],[-74.18,5.62],[-52.21,-0.53],[-9.67,0.35],[-17.4,30.9],[-98.44,36.32],[-98.44,32.25],[-98.44,19.48]])
			lon_mesh, lat_mesh = np.meshgrid(lon_grid,lat_grid)
			in_reg = np.expand_dims(reg_atl_path.contains_points(np.vstack((lon_mesh.flatten(),lat_mesh.flatten())).T),1)

			sig.values[in_reg.squeeze().reshape(sig.shape) == 0] = 1

		# optional: only consider ocean grid-cells
		if useMask:
			oceanMask = xa.open_dataset('data/'+identifier.split('_')[0]+'_oceanMask.nc')['oceanMask']
			oceanMask.lon.values[oceanMask.lon.values > 180] -= 360
			if np.any(oceanMask.lat != sig.lat):
				oceanMask.reindex(lat=list(reversed(oceanMask.lat)))
			oceanMask.values[oceanMask!=1] = 0
			sig.values[np.where(oceanMask==False)] = 1

		# if it hasen't been done yet, initialize the output arrays
		if first:
			out_fld,out_ts = init_output()
			first = False

		# store data
		out_fld['corr'].loc[:,:,lag,name] = corr
		out_fld['corr_sig'].loc[:,:,lag,name] = sig.values
		for dataset in ['out_fld','out_ts']:
			globals()[dataset].to_netcdf(result_path+identifier+'_'+dataset.split('_')[1]+'.nc')

		# cluster grid-cells intor regions (of same sign correlation)
		regions, out_fld['clusters'].loc[:,:,lag,name], cluster_dict = cluster_sig_points_into_regions(out_fld['corr'].loc[:,:,lag,name], out_fld['corr_sig'].loc[:,:,lag,name], lat_grid, lon_grid, cluster_params)

		# create a latitude-weighting array
		lat_weight = np.cos(np.deg2rad(lat_grid))
		lat_weight_array = np.repeat(lat_weight[np.newaxis,:], var_vals.shape[2], 0).T

		print(sorted(cluster_dict.keys()))
		# store regions and time series averaged over the regions
		for reg_id in sorted(cluster_dict.keys()):
			if reg_id-1 < n_keep_save:
				out_fld['region'].loc[:,:,lag,name,reg_id] = regions.loc[:,:,reg_id]
				mask = regions.loc[:,:,reg_id].copy()
				mask.values[mask==0] = np.nan
				lons = out_fld.lon
				out_ts['actor'].loc[lag,name,reg_id,:] = np.nanmean(var_anoms.loc[:,out_fld.lat,lons.values] * np.expand_dims(mask,0) * lat_weight_array, axis=(1,2))
				out_ts['actor_raw'].loc[lag,name,reg_id,:] = np.nanmean(var_vals.loc[:,out_fld.lat,lons.values] * np.expand_dims(mask,0) * lat_weight_array, axis=(1,2))

		for dataset in ['out_fld','out_ts']:
			globals()[dataset].to_netcdf(result_path+identifier+'_'+dataset.split('_')[1]+'.nc')

		gc.collect()

	###
