# -*- coding: utf-8 -*-
import sys, os,pickle, inspect, textwrap, importlib, glob, itertools
from netCDF4 import Dataset, num2date
from datetime import datetime, date, timedelta
import xarray as xa
import pandas
from pandas import DataFrame
import pandas as pd
import scipy
from scipy import signal
from scipy.spatial import ConvexHull
from statsmodels.sandbox.stats import multicomp
from string import ascii_lowercase

import seaborn as sns
import numpy
import numpy as np
#from matplotlib import rc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib import colors
setColorMap = matplotlib.colors.ListedColormap(['#03a17b','yellow',"r", "b", 'm', 'orange', 'c', 'darkmagenta', 'yellow','k'])
#rc('text', usetex=True)

# from pyPdf import PdfFileWriter, PdfFileReader
# from fpdf import FPDF

from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant

from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from scipy.spatial import ConvexHull

import sklearn
import cartopy
import cartopy.crs as ccrs
from itertools import cycle
from cartopy.util import add_cyclic_point

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from scipy.optimize import curve_fit
from scipy.stats.distributions import  t

from sklearn import metrics
from haversine import haversine
from sklearn.cluster import DBSCAN

from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import brier_score_loss
from sklearn.utils import resample

from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.calibration import calibration_curve

import pandas as pd

def save_pkl(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pkl(name ):
    with open( name, 'rb') as f:
        return pickle.load(f)

import statsmodels.api as sm
from sklearn import preprocessing
import statsmodels.api as sm

#=====================================================================================
# 0) Parameters which must be specified
#=====================================================================================

try:
	os.chdir('/Users/peterpfleiderer/Projects/tropical_cyclones/')
	base_path = '/Users/peterpfleiderer/Projects/tropical_cyclones/'
	data_path = '/Users/peterpfleiderer/Projects/data/'
	tmp_path = '/Users/peterpfleiderer/Projects/tropical_cyclones/'
	sys.path.append(base_path+'git-packages/tigramite')
except:
	os.chdir('/home/pepflei')
	base_path = '/home/pepflei/Projects/'
	data_path = '/p/projects/tumble/carls/shared_folder/data/'
	tmp_path = '/p/tmp/pepflei/'
	sys.path.append('/p/projects/tumble/carls/shared_folder/git-packages/tigramite')

#import tigramite
# from tigramite import data_processing as pp
# from tigramite import plotting as tp
# from tigramite.pcmci import PCMCI
# from tigramite.models import LinearMediation, Prediction
# from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb

index_name = "ACE"
index_longname = "$ACE_{JASO}$"
os.chdir(base_path+'ACE_seasonal_forecast')
tmp_path += 'ACE_seasonal_forecast' + '/'

# time-cycle of data
# 12 for 7days, 365 for daily etc...
time_cycle = 12
# Start month is index of first relevant time-step
# start_day = time-cycle, if alls values of the years should be considered
start_step = 6 # means start with November (= 10ther entry in array)

n_steps = 4

seas_indices = range(n_steps)

def get_RV_indices(tau_min,tau_max,time_cycle,n_steps,n_years,start_day):
	#=====================================================================================
	# Calculates the indices which are taken for response-variable
	#=====================================================================================
	lag_steps = tau_max - tau_min +1
	time_range_all = [0, time_cycle * n_years]
	RV_indices = []

	# case 1) complete time-series
	# cut-off the first year
	if n_steps == time_cycle:
		RV_indices = range(time_cycle*n_years)[time_cycle:]

	# # case 2) starts with Jan or tau_max is in previous year:
	# # cutoff the first year
	# elif start_day - tau_max < 0: #5- 5 = 0
	# 	for i in range(n_steps):
	# 		a = list(range(time_cycle*n_years)[time_cycle:][start_day +i::time_cycle])
	# 		RV_indices = RV_indices + a

	#case 3) winter overlap DJ
	# cut-off last year except winter
	elif start_day + n_steps >=time_cycle: #5+4 =9
		for i in range(n_steps):
			a = list(range(time_cycle*n_years)[:-(start_day + n_steps -time_cycle)][start_day +i::time_cycle])
			RV_indices = RV_indices + a

	#case 4) all good.
	else:
		for i in range(n_steps):
			a = list(range(time_cycle*n_years)[start_day +i::time_cycle])
			RV_indices = RV_indices + a

	RV_indices.sort()
	return np.array(RV_indices),time_range_all,lag_steps


def load_nc(filename,varname,cutoff,time_cycle,level=None, cut_box=None):
	nc = xa.open_dataset(filename)
	if level is not None:
		var_vals = nc[varname].loc[:,level,:,:].squeeze()
	else:
		var_vals = nc[varname].squeeze()

	var_vals = var_vals.loc[cutoff[0]:cutoff[1]]

	if cut_box is not None:
		var_vals = var_vals.loc[:,cut_box['lat1']:cut_box['lat2'],cut_box['lon1']:cut_box['lon2']]

	lat_grid, lon_grid = var_vals.lat.values, var_vals.lon.values
	lon_grid[lon_grid > 180 ] -= 360
	var_vals.coords['lon'] = lon_grid

	# set monthly values to the first of month
	var_vals.coords['time'] = np.array([np.datetime64(str(tt)[:8]+'01'+str(tt)[10:]) for tt in var_vals.time.values])

	# # # detrend in each time-step
	var_anoms = var_vals.copy()
	nan_IDs = np.isnan(var_vals)
	var_anoms.values[nan_IDs] = -9e+27
	for j in range(time_cycle):
		var_anoms.values[j::time_cycle] = scipy.signal.detrend(var_anoms.values[j::time_cycle], axis = 0)
	#calculate anomalies
	for i in range(int(time_cycle)):
		var_anoms.values[i::int(time_cycle)] = (var_anoms.values[i::int(time_cycle)] - np.mean(var_anoms.values[i::int(time_cycle)], axis = 0))

	var_standard = var_anoms / np.std(var_vals.values,axis=0)

	return lat_grid, lon_grid, var_vals, var_anoms, var_standard, nan_IDs


def get_correlation(index, field):
	"""
	This function calculates the correlation coefficent r  and the the pvalue p for each grid-point of field D with the response-variable di
	"""
	field_2D = np.reshape(field.values, (field.shape[0],-1))
	x = numpy.ma.zeros(field_2D.shape[1])
	corr_di_D = numpy.ma.array(data = x, mask =False)
	sig_di_D = numpy.array(x)

	for i in list(range(field_2D.shape[1])):
		r, p = scipy.stats.pearsonr(index.values,field_2D[:,i])
		corr_di_D[i]= r
		sig_di_D[i]= p

	corr = field[0,:,:].squeeze().copy()
	corr.values = corr_di_D.reshape(corr.shape)
	corr_sig = field[0,:,:].squeeze().copy()
	corr_sig.values = sig_di_D.reshape(corr.shape)
	return corr, corr_sig

def cluster_sig_points_into_regions(corr, sig, lats, lons, cluster_params):
	cells = np.where((sig==0) & (np.isfinite(sig)))

	points = np.array([[lats[la],lons[lo]] for la,lo in zip(cells[0],cells[1])])
	vals = np.array([[np.sign(corr[la,lo])] for la,lo in zip(cells[0],cells[1])])

	distance_vals = metrics.pairwise_distances(vals) * 9999
	distance = metrics.pairwise_distances(points, metric=haversine)
	dbresult = DBSCAN(eps=cluster_params['eps'], min_samples=cluster_params['min_samples'], metric='precomputed').fit(distance + distance_vals)
	labels = dbresult.labels_ + 1
	label_list = sorted(set(labels))
	if 0 in label_list:
		label_list.remove(0)

	cluster_dict = {}
	for reg_id,label in enumerate(np.array(label_list)[np.argsort([sum(labels == label) for label in label_list])[::-1]]):
		reg_id += 1
		tmp_points = points[labels == label]
		cluster_dict[reg_id] = {'count':sum(labels == label),
				  'lons':np.array([pp[1] for pp in tmp_points]),
				  'lats':np.array([pp[0] for pp in tmp_points]),
				  'point_id':np.where(labels == label)}
		cluster_dict[reg_id]['c_lat'] = np.mean(cluster_dict[reg_id]['lats'].copy())
		if len(np.where(cluster_dict[reg_id]['lons']>310)[0]) > 0 and len(np.where(cluster_dict[reg_id]['lons']<50)[0]) > 0:
			cluster_dict[reg_id]['c_lon'] = np.mean(cluster_dict[reg_id]['lons'][cluster_dict[reg_id]['lons']>150].copy())
		if len(np.where(cluster_dict[reg_id]['lons']>150)[0]) > 0 and len(np.where(cluster_dict[reg_id]['lons']<-150)[0]) > 0:
			cluster_dict[reg_id]['c_lon'] = np.mean(cluster_dict[reg_id]['lons'][cluster_dict[reg_id]['lons']>150].copy())
		else:
			cluster_dict[reg_id]['c_lon'] = np.mean(cluster_dict[reg_id]['lons'].copy())

	clusters = corr.copy() * np.nan
	for reg_id,details in cluster_dict.items():
		for lo,la in zip(details['lons'],details['lats']):
			clusters.loc[la,lo] = reg_id


	regions = xa.DataArray(data = np.zeros([corr.shape[0],corr.shape[1],len(label_list)])*np.nan, coords=dict(lat = lats, lon = lons, ID = label_list), dims=['lat','lon','ID'])
	for reg_id in regions.ID:
		regions.loc[:,:,reg_id] = clusters == reg_id

	print(str(len(cluster_dict.keys()))+' regions detected')

	return regions, clusters, cluster_dict


def get_scores(events,predict,probs,lvl):
	tp = (predict) & (events)
	tn = (predict==False) & (events==False)
	fp = (predict) & (events==False)
	fn = (predict==False) & (events)
	tpr = tp.sum() / float(events.sum())
	fpr = fp.sum() / float(np.sum(events==False))
	ACC = accuracy_score(events, predict)

	bs = brier_score_loss(np.array(events,np.int), probs)
	bs_c = np.sum((events - (1-lvl))**2) / np.float(len(events))
	bss = (bs_c-bs)/bs_c
	return {'tpr':tpr, 'fpr':fpr, 'ACC':ACC, 'bss':bss}

def get_roc(original,predicted,thresh):
	events = np.array(original > thresh,np.bool)
	scores = predicted + np.nanmean(predicted)

	fpr, tpr, thresholds = metrics.roc_curve(events,scores)

	ROCA = metrics.auc(fpr, tpr)

	return {'tpr':tpr, 'fpr':fpr, 'ROCA':ROCA}

def roc_bootstrap(y,pred,thresh,N):
	ROCA_real = get_roc(y,pred,thresh)['ROCA']

	ROCA_shu = []
	pred_shu = pred.copy()
	for i in range(N):
		np.random.shuffle(pred_shu)
		ROCA_shu.append(get_roc(y,pred_shu,thresh)['ROCA'])

	return np.sum(np.array(ROCA_shu) > ROCA_real) / float(N)

def plot_ROC(obs,pred,ax,title='', n_boot=1000):
	finite = np.isfinite(obs) & np.isfinite(pred)
	obs, pred = obs[finite], pred[finite]

	#ax.set_title(title)
	legend_elements=[]
	roc_summary = {}
	for lvl,thresh,label,color,lsty in zip([33,50,66],list(np.nanpercentile(obs,[33,50,66])),['$>33^{rd} perc.$','$>median$','$>66^{th} perc.$'],['#1f8c3c','blue','#9842f5'],[':','-','--']):

		tmp_roc = get_roc(obs,pred,thresh)
		ax.plot(tmp_roc['fpr'],tmp_roc['tpr'], color=color, linestyle=lsty)
		pval_ROCA = roc_bootstrap(obs,pred,thresh,n_boot)
		label = label+' ROCA='+str(round(tmp_roc['ROCA'],2))
		if pval_ROCA < 0.05:	label += '**'
		elif pval_ROCA < 0.1:	label += '*'
		legend_elements.append(Line2D([0], [0], color=color, linestyle=lsty, label=label))
		roc_summary[lvl] = tmp_roc
		roc_summary[lvl]['ROCA_pval'] = pval_ROCA
	legend_elements.append(Line2D([0], [0], color='gray', linestyle='-', label='guess'))
	ax.plot([0, 1], [0, 1], color='gray')
	ax.legend(handles=legend_elements, loc='lower right', frameon=True, facecolor='w', ncol=1, framealpha=1, edgecolor='w', fontsize = 9).set_zorder(1)
	ax.set_ylabel('true positive rate')
	ax.set_xlabel('false positive rate')

	#
	# ax = axes[1]
	# ax.axis('off')
	# ax.legend(handles=legend_elements, loc='center', frameon=True, facecolor='w', ncol=1, framealpha=1, edgecolor='w').set_zorder(1)
	# plt.tight_layout(); pdf.savefig(transparent = False); plt.close()

	# fig, ax = plt.subplots(nrows=1, figsize=(3,2))
	# ax.axis('off')
	# ax.legend(handles=legend_elements, loc='center', frameon=True, facecolor='w', ncol=1, framealpha=1, edgecolor='w').set_zorder(1)
	# plt.tight_layout(); pdf.savefig(transparent = True); plt.close()

	return roc_summary
