# -*- coding: utf-8 -*-
import os,sys,importlib,time

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# create a folder for the results
if os.path.isdir(tmp_path+'additional')==False:
	os.system('mkdir -p '+tmp_path+'additional')

identifier = 'causal_JRA55_vws_lag4_X3_MSLP_wideClusters_n2'
dataset = identifier.split('_')[1]
sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

target_raw = xa.open_dataset(index_name+'.nc')[index_name+'_raw'].squeeze()

data = pd.DataFrame(index=range(1958,2019), columns=['train','test'])
data['year'] = range(1958,2019)
data['train'].loc[1979:2018] = True
data['test'].loc[1979:2018] = False
data['train'].loc[1958:1978] = False
data['test'].loc[1958:1978] = True
nc_mslp = xa.open_dataset(tmp_path+'additional/causal_JRA55_vws_lag4_X3_MSLP_wideClusters_n2_fingerprint_ts.nc')
data['mslp_neg'] = nc_mslp['neg'].values[2::12]
data['mslp_pos'] = nc_mslp['pos'].values[2::12]
nc_sst = xa.open_dataset(tmp_path+'additional/causal_JRA55_sst_lag4_X3_SST_dropSmall_n1_fingerprint_ts.nc')
data['sst_pos'] = nc_sst['pos'].values[2::12]
data['ACE'] = target_raw.values[6::12]

predictors = ['mslp_pos','mslp_neg','sst_pos']

clf = linear_model.LinearRegression()
clf.fit(data.loc[(data.train),predictors], data.loc[(data.train),'ACE'])
frc = clf.predict(data.loc[:,predictors])
stdev = np.std(frc - data.loc[:,'ACE'])
data['hindcast'] = frc
data['hindcast_upper'] = frc + stdev
data['hindcast_lower'] = frc - stdev

scipy.stats.spearmanr(data.loc[(data.test),'ACE'],data.loc[(data.test),'hindcast'])
scipy.stats.pearsonr(data.loc[(data.test),'ACE'],data.loc[(data.test),'hindcast'])

plt.close()
with PdfPages(tmp_path + 'additional/JRA55_early_skill.pdf') as pdf:
	fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))

	ax.plot(data.year,data.ACE,'k')

	ax.plot(data.loc[(data.test),'year'],data.loc[(data.test),'hindcast'],'m')
	ax.fill_between(data.loc[(data.test),'year'],data.loc[(data.test),'hindcast_upper'],data.loc[(data.test),'hindcast_lower'], color='m', alpha=0.5)

	ax.plot(data.loc[(data.train),'year'],data.loc[(data.train),'hindcast'],'c')
	ax.fill_between(data.loc[(data.train),'year'],data.loc[(data.train),'hindcast_upper'],data.loc[(data.train),'hindcast_lower'], color='c', alpha=0.5)

	legend_elements=[]
	for set__,color in zip(['train','test'],['c','m']):
		spearmanr = scipy.stats.spearmanr(data.loc[(data[set__]),'ACE'],data.loc[(data[set__]),'hindcast'])
		spearman_str = str(round(spearmanr[0],2))
		if spearmanr[1] < 0.1: spearman_str += '*'
		if spearmanr[1] < 0.05: spearman_str += '*'
		pearsonr = scipy.stats.pearsonr(data.loc[(data[set__]),'ACE'],data.loc[(data[set__]),'hindcast'])
		pearson_str = str(round(pearsonr[0],2))
		if pearsonr[1] < 0.1: pearson_str += '*'
		if pearsonr[1] < 0.05: pearson_str += '*'
		legend_elements.append(Line2D([0], [0], color=color, label='pearson r: '+pearson_str+ ' spearman r: '+spearman_str))

	ax.legend(handles=legend_elements, loc='best', frameon=True, facecolor='w', ncol=1, framealpha=0.6, edgecolor='none', fontsize = 9).set_zorder(1)
	plt.tight_layout(); pdf.savefig(); plt.close()

	fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	plot_ROC(data.loc[(data.test),'ACE'].values, data.loc[(data.test),'hindcast'].values, ax, n_boot = 1000)
	plt.tight_layout(); pdf.savefig(); plt.close()

#
