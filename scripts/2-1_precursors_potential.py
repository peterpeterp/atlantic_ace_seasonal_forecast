# -*- coding: utf-8 -*-
import os,sys,importlib,gc

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# check whether the script is called from outside
# the "identifier" is the name of a setting file in scripts/params/'+dataset+'/potential/ containing all relevant parameters
try:
	identifier = sys.argv[1]
except:
	identifier = "JRA55_vws_lag4_X3_MSLP_largeRegs"

# just a convention: the first part of the identifier denotes the dataset
dataset = identifier.split('_')[0]

# create a folder for the results
if os.path.isdir(tmp_path+'precursor/'+dataset+'/'+identifier)==False:
	os.system('mkdir '+tmp_path+'precursor/'+dataset+'/'+identifier)
	os.system('mkdir '+tmp_path+'precursor/'+dataset+'/'+identifier+'/logs')

# save used parameters in this folder
params_in = open('scripts/params/'+dataset+'/potential/'+identifier+'.py', 'r').read()
out = open(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/'+identifier+'.py', 'w')
out.write(params_in)
out.close()

# load parameters
# all varaibles that are required in this script are defined there
sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

# go through different varaibles
for field_name,precursor in precursor_dict.items():
	# load cliamte varaible in which potential precursors are expected
	if 'level' in precursor.keys():
		lat_grid, lon_grid, fieldVals, fieldAnoms, fieldStandard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff,time_cycle,level=precursor['level'], cut_box=cut_box)
	else:
		lat_grid, lon_grid, fieldVals, fieldAnoms, fieldStandard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff,time_cycle, cut_box=cut_box)

	# initialize output arrays
	coords = dict(lat = np.asarray(lat_grid), lon = np.asarray(lon_grid), lag=np.arange(lag,lag+1,1, dtype=np.short), set_id=np.array(split_sets.index.levels[0]), ID = np.arange(1,n_keep_save+1,dtype=np.short), time= fieldAnoms.time)
	empty = np.zeros([len(lat_grid),len(lon_grid),1,len(split_sets.index.levels[0])]) * np.nan
	potential_field = xa.Dataset({
		'corr': xa.DataArray(data=empty.copy()[:,:,:,:], coords={k: coords[k] for k in ('lat','lon','lag','set_id')}, dims=['lat','lon','lag','set_id']),
		'corr_sig': xa.DataArray(data=empty.copy()[:,:,:,:],coords={k: coords[k] for k in ('lat','lon','lag','set_id')}, dims=['lat','lon','lag','set_id']),
		'corr_random': xa.DataArray(data=empty.copy()[:,:,:,:],coords={k: coords[k] for k in ('lat','lon','lag','set_id')}, dims=['lat','lon','lag','set_id']),
		'clusters': xa.DataArray(data=empty.copy()[:,:,:,:],coords={k: coords[k] for k in ('lat','lon','lag','set_id')}, dims=['lat','lon','lag','set_id']),
		})
	empty = np.zeros([1,len(split_sets.index.levels[0]), len(np.arange(1,n_keep_save+1,dtype=np.short)),len(fieldAnoms.time)]) * np.nan
	potential_ts = xa.Dataset({
		'actor': xa.DataArray(data=empty.copy()[:,:,:,:], coords={k: coords[k] for k in ('lag','set_id','ID','time')}, dims=['lag','set_id','ID','time']),
		'actor_raw': xa.DataArray(data=empty.copy()[:,:,:,:], coords={k: coords[k] for k in ('lag','set_id','ID','time')}, dims=['lag','set_id','ID','time']),
		})

	# go through all training sets
	for set_id in all_years:
		print('treating test year ',set_id)

		# select split set and identify training years
		splitSet = split_sets.loc[set_id]
		years_train = all_years[np.where(splitSet['train'])]

		# optional:
		# for robustness:
		#	training set is spilt again to check which grid cells are really robustly correlated
		if leaveOut_potential != 0:
			corr_sig = fieldAnoms[0].copy() * 0
			corr = corr_sig.copy() * 0

			# shift a subset of the training set leaveOut_potential times to get different sub-training-sets
			for shift in range(leaveOut_potential):
				print('shift ',shift)
				# define new sub-training-set
				years_tmp = years_train[shift:-(leaveOut_potential-shift)]
				target_timeSteps_tmp = [tst for tst in target_indices.time.values if int(str(tst).split('-')[0]) in years_tmp]
				prec_timeSteps_tmp = [tst for tst in prec_indices.time.values if int(str(tst).split('-')[0]) in years_tmp]

				# compute lagged pointwise correlation
				corr_tmp, corr_sig_tmp = get_correlation(target.loc[target_timeSteps_tmp], fieldAnoms.loc[prec_timeSteps_tmp,:,:])
				corr += corr_tmp

				# treat the percent_sig most significant grid-cells as significant
				# add a 1 to grid-cells that aren't significant
				corr_sig += corr_sig_tmp > np.nanpercentile(corr_sig_tmp, percent_sig)
				gc.collect()

			# treat grid-cells as significant if they have been identified as such in all sub-training-sets
			potential_field['corr_sig'].loc[:,:,lag,set_id] = corr_sig != 0
			# check which grid-cells have been identified as significant in some sub-training-sets
			potential_field['corr_random'].loc[:,:,lag,set_id] = (corr_sig != leaveOut_potential) & (corr_sig != 0)
			# compute the mean correlation
			potential_field['corr'].loc[:,:,lag,set_id] = corr / float(leaveOut_potential)

		# if leaveOut_potential == 0
		# no additional robustness check
		else:
			# get relevant time steps for the training-set
			target_timeSteps_train = [tst for tst in target_indices.time.values if int(str(tst).split('-')[0]) in years_train]
			prec_timeSteps_train = [tst for tst in prec_indices.time.values if int(str(tst).split('-')[0]) in years_train]

			# compute lagged pointwise correlation
			corr, corr_sig = get_correlation(target.loc[target_timeSteps_train], fieldAnoms.loc[prec_timeSteps_train,:,:])

			# treat the percent_sig most significant grid-cells as significant
			corr_sig = corr_sig > np.nanpercentile(corr_sig, percent_sig)

			# store results
			potential_field['corr_sig'].loc[:,:,lag,set_id] = corr_sig != 0
			potential_field['corr'].loc[:,:,lag,set_id] = corr

		# cluster grid-cells intor regions (of same sign correlation)
		regions, potential_field['clusters'].loc[:,:,lag,set_id], cluster_dict = cluster_sig_points_into_regions(potential_field['corr'].loc[:,:,lag,set_id], potential_field['corr_sig'].loc[:,:,lag,set_id], lat_grid, lon_grid, cluster_params)

		# create a latitude-weighting array
		lat_weight = np.cos(np.deg2rad(lat_grid))
		lat_weight_array = np.repeat(lat_weight[np.newaxis,:], corr.shape[1], 0).T

		print(sorted(cluster_dict.keys()))
		# store time series averaged over the regions as detrended anomalies and as absolute values
		for reg_id in sorted(cluster_dict.keys()):
			if reg_id-1 < n_keep_save:
				mask = regions.loc[:,:,reg_id].copy()
				mask.values[mask==0] = np.nan
				lons = potential_field.lon
				potential_ts['actor'].loc[lag,set_id,reg_id,:] = np.nanmean(fieldAnoms.loc[:,potential_field.lat,lons] * np.expand_dims(mask,0) * lat_weight_array, axis=(1,2))
				potential_ts['actor_raw'].loc[lag,set_id,reg_id,:] = np.nanmean(fieldVals.loc[:,potential_field.lat,lons] * np.expand_dims(mask,0) * lat_weight_array, axis=(1,2))

		# save results
		for dataset__ in ['potential_field','potential_ts']:
			globals()[dataset__].to_netcdf(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_'+dataset__+'_'+field_name+'.nc')

		gc.collect()

# finally plot a summary figure showing the number of times a grid-cell is part of a potential precursor region
from plot_precursor_reoccurring import plot_reoccurring
plot_reoccurring(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_potential_field_'+field_name+'.nc')








'''
for dataset in ERA5 JRA55; do sbatch job.sh scripts/2-1_precursors_potential.py ${dataset}_sst_lag4_X3_SST_v1; done;
for dataset in ERA5 JRA55; do sbatch job.sh scripts/2-1_precursors_potential.py ${dataset}_sst_lag2_X3_SST; done;
for dataset in ERA5 JRA55; do sbatch job.sh scripts/2-1_precursors_potential.py ${dataset}_vws_lag4_X3_MSLP_wideClusters; done;
for dataset in ERA5 JRA55; do sbatch job.sh scripts/2-1_precursors_potential.py ${dataset}_vws_lag2_X3_SST; done;
'''








#
