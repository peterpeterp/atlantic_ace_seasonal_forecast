# -*- coding: utf-8 -*-
import os,sys,importlib,gc

sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# define cutoff
cutoff_years = [1979,2018]
cutoff = [np.datetime64(str(cutoff_years[0])+'-01-01'),np.datetime64(str(cutoff_years[1])+'-12-31')]
all_years = np.arange(cutoff_years[0], cutoff_years[1]+1)

# load the target index to get the time format
target =  xa.open_dataset(index_name+'.nc')[index_name+'_allMon'].loc[cutoff[0]:cutoff[1]].squeeze()
time_axis = target.time.values

# create folder
os.system('mkdir split_sets')

# for different numbers of years left out training - testing sets are created
for nSplit in range(1,5):
	# initialize multiindex dataframe
	set_X = pd.DataFrame(index=[np.array([[ii]*len(all_years) for ii in all_years]).flatten(), np.array(list(all_years)*len(all_years)).flatten()], columns=['train','test'])

	for iy,testYear in enumerate(all_years):
		# set all train years true except the test year and the nSplit-1 years before
		set_X.loc[testYear,'train'] = True
		for i in range(nSplit):
			if iy-i>=0:
				set_X.loc[testYear].loc[testYear-i, 'train'] = False

		# only the test year is set true
		set_X.loc[testYear,'test'] = False
		set_X.loc[testYear].loc[testYear, 'test'] = True

	set_X.to_pickle('split_sets/setX'+str(nSplit)+'.pkl')
















#
