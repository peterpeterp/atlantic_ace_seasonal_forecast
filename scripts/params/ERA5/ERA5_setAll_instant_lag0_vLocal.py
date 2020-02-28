# -*- coding: utf-8 -*-
import sys, os, inspect, importlib
import numpy as np

os.chdir('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

lag = 0

precursor_dict={
	'sst':{'file':data_path+'/ERA5/ERA5_MSLP-SST_1979-2018_grid1x1.nc','var':'var34'},
	'vws':{'file':data_path+'/ERA5/ERA5_VWS_1979-2018_grid1x1.nc','var':'vws'},
}

FDR = True

# regions over which the correlation maps should be calculated, resp. domain which should contain the precursos communities:
la_min,la_max,lo_min,lo_max = -89,89,-180,360
box = [la_min, la_max, lo_min, lo_max]

n_keep_save = 3

percent_sig = 3

useMask = True
inReg = True

cluster_params = {'eps':500, 'min_samples':15}

cutoff_years = [1979,2019]
cutoff = [np.datetime64(str(cutoff_years[0])+'-01-01'),np.datetime64(str(cutoff_years[1])+'-01-01')]

all_years = np.arange(cutoff_years[0],cutoff_years[1],1,np.int)

RV_indices,time_range_all,lag_step = get_RV_indices(lag,lag,time_cycle,n_steps,len(all_years),start_step)
