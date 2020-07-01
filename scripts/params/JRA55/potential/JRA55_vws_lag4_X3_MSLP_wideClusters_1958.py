# -*- coding: utf-8 -*-
import sys, os, inspect, importlib
import dimarray as da
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

#=====================================================================================
# Potentail Precursor Fields
#=====================================================================================
precursor_dict={
	# 'sst':{'file':data_path+'/SST/COBE_sst_mon_145x288.nc','var':'sst'},
	'mslp':{'file':data_path+'/JRA55/mon_JRA55_002_prmsl.nc','var':'var2'},
}

# regions over which the correlation maps should be calculated, resp. domain which should contain the precursos communities:
cut_box = dict(lat1=60, lat2=-60, lon1=-180, lon2=360)

percent_sig = 7.5
cluster_params = {'eps':1500, 'min_samples':20}

# how many regions to keep per field
n_keep_save = 20

#=====================================================================================
# Time Parameters
#=====================================================================================
lag = 4

cutoff_years = [1958,2018]
cutoff = [np.datetime64(str(cutoff_years[0])+'-01-01'),np.datetime64(str(cutoff_years[1])+'-12-31')]
all_years = np.arange(cutoff_years[0],cutoff_years[1]+1)

RV_indices,time_range_all,lag_step = get_RV_indices(lag,lag,time_cycle,1,len(all_years),start_step)

#=====================================================================================
# Robustness and sets
#=====================================================================================
split_sets = pd.DataFrame(index=[np.array([[ii]*len(all_years) for ii in all_years]).flatten(), np.array(list(all_years)*len(all_years)).flatten()], columns=['train','test'])
for iy,testYear in enumerate(all_years):
	# set all train years true except the test year and the nSplit-1 years before
	split_sets.loc[testYear,'train'] = True
	for i in range(3):
		if iy-i>=0:
			split_sets.loc[testYear].loc[testYear-i, 'train'] = False

	# only the test year is set true
	split_sets.loc[testYear,'test'] = False
	split_sets.loc[testYear].loc[testYear, 'test'] = True

leaveOut_potential = 0

#=====================================================================================
# TARGET
#=====================================================================================
instant_version = 'JRA55_setAll_instant_lag0_vLocal'
target_name = 'vws'

index1D_tmp = xa.open_dataset(tmp_path+'instantCond/JRA55_setAll_instant_lag0_vLocal_ts.nc')['actor'].loc[0,target_name,1,:].loc[cutoff[0]:cutoff[1]]
index1D_time = index1D_tmp.time
index1D_tmp = index1D_tmp.squeeze()

index1D = index1D_tmp.copy() * 0

print('Jittering: adding small random numbers to index1D in months outside of the season')
for month in range(12):
	index1D[month::time_cycle] = np.array([np.random.random()*1e-20 for i in range(int(len(index1D_time)/time_cycle))])
for i in range(1,n_steps):
	index1D[start_step::time_cycle].values += index1D_tmp[start_step+i::time_cycle].values
index1D/=float(n_steps)

for j in range(time_cycle):
	index1D.values[j::time_cycle] = scipy.signal.detrend(index1D.values[j::time_cycle], axis = 0)
for i in range(int(time_cycle)):
	index1D.values[i::int(time_cycle)] = (index1D.values[i::int(time_cycle)] - np.mean(index1D.values[i::int(time_cycle)], axis = 0))

target = index1D

target_indices = target[RV_indices].time

prec_indices = target[RV_indices-lag].time


#
