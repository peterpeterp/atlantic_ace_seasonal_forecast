# -*- coding: utf-8 -*-
import sys, os, inspect, importlib
import dimarray as da
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

#=====================================================================================
# Potentail Precursor Fields
#=====================================================================================
precursor_dict={
	'sst': 'ERA5_sst_lag4_X3_SST_dropSmall',
}

#=====================================================================================
# Time Parameters
#=====================================================================================
lag = 4

cutoff_years = [1979,2018]
cutoff = [np.datetime64(str(cutoff_years[0])+'-01-01'),np.datetime64(str(cutoff_years[1])+'-12-31')]
all_years = np.arange(cutoff_years[0],cutoff_years[1]+1)

RV_indices,time_range_all,lag_step = get_RV_indices(lag,lag,time_cycle,1,len(all_years),start_step)

#=====================================================================================
# Robustness and sets
#=====================================================================================
split_sets = pd.read_pickle('split_sets/setX3.pkl')
leaveOut_tigramite = 0


#=====================================================================================
# Tigramite settings
#=====================================================================================
mode = 'fixed not'
pmci_alpha_list = np.arange(0.001,0.3,0.001)
n_keep = np.arange(30,10,-1)
nPrecursors_options = [1]
pq_matrix_options = ['q_matrix','p_matrix']



#=====================================================================================
# TARGET
#=====================================================================================
instant_version = 'ERA5_setAll_instant_lag0_vLocal'
target_name = 'sst'

index1D_tmp = xa.open_dataset(tmp_path+'instantCond/ERA5_setAll_instant_lag0_vLocal_ts.nc')['actor'].loc[0,target_name,1,:].loc[cutoff[0]:cutoff[1]]
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
