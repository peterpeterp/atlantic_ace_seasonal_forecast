# -*- coding: utf-8 -*-
import os,sys,importlib,time

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# define all forecasts that should be tested in a prediction_dict
prediction_dict = {}
for dataset in ['ERA5','JRA55']:
	prediction_dict[dataset+ '_lag4_X3_wideClusters_2'] = {'lag':4, 'dataset':dataset,
			'potential_prec_dict': {
				'causal_'+dataset+'_vws_lag4_X3_MSLP_wideClusters_n2': dict(limit_regs= dict(mslp=2)),
				'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
			}}

	prediction_dict[dataset+ '_lag4_X3_wideClusters_2_mslpOnly'] = {'lag':4, 'dataset':dataset,
			'potential_prec_dict': {
				'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
			}}

	prediction_dict[dataset+ '_lag4_X3_sstOnly'] = {'lag':4, 'dataset':dataset,
			'potential_prec_dict': {
				'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
			}}

	prediction_dict[dataset+ '_lag2_X3_v1'] = {'lag':2, 'dataset':dataset,
			'potential_prec_dict': {
				'causal_'+dataset+'_vws_lag2_X3_SST_n2': dict(limit_regs= dict(sst=2)),
				'causal_'+dataset+'_sst_lag2_X3_SST_n1': dict(limit_regs= dict(sst=1)),
			}}

# define statistical models to be tested
model_dict = {
	'lr' : dict(clf = linear_model.LinearRegression(), color = 'm', weighting = None),
	# 'lr-weighted' : dict(clf = linear_model.LinearRegression(), color = 'm', weighting = {66:1,33:1}),
	# 'randomForest' : dict(clf = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100), color = 'orange'),
}

# define classifiers to be tested
classi_dict = {}
# classi_dict['KNN'] = {lvl: dict(clf=KNeighborsClassifier(3), color='b') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['randomForest'] = {lvl: dict(clf=RandomForestClassifier(max_depth=10, n_estimators=100, max_features=1, class_weight = 'balanced'), color='g') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', class_weight = {0:1-lvl,1:lvl}), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', class_weight = {0:1-lvl,1:lvl}), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='liblinear', multi_class='auto', C=1, class_weight = 'balanced'), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial'), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='auto'), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
#classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(C=1, penalty='l2',solver='saga',multi_class='ovr',max_iter=10000, class_weight = 'balanced'), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
#classi_dict['logReg_noWeight'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial'), color='c') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# from sklearn.gaussian_process import GaussianProcessClassifier
# from sklearn.gaussian_process.kernels import RBF
# classi_dict['RBF'] = {lvl: dict(clf=GaussianProcessClassifier(1.0 * RBF(1.0)), color='m') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# from sklearn.naive_bayes import GaussianNB
# classi_dict['NaivBayes'] = {lvl: dict(clf=GaussianNB(), color='c') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['QDA'] = {lvl: dict(clf=QuadraticDiscriminantAnalysis(), color='darkorange') for lvl in [0.25,0.33,0.5,0.66,0.75]}
#classi_dict['SVC'] = {lvl: dict(clf=SVC(kernel='linear', C=10, probability=True,random_state=0), color='darkorange') for lvl in [0.33,0.5,0.66]}

# go through all forecasts to be tested
for prediction,params in prediction_dict.items():

	# forecast_summary = load_pkl('prediction/summary/'+prediction+'_skill.pkl')
	# asdas

	forecast_summary = {}

	# set all entries of the prediction settings as a global variable
	for key,val in params.items():
		globals()[key] = val

	for tmp_identifier,details in potential_prec_dict.items():
		# load causal precursor
		details['causal'] = load_pkl(tmp_path+'precursor/'+dataset+'/'+tmp_identifier+'/'+tmp_identifier+'_causal.pkl')
		details['potential_ts'] = {}

		# load all parameter settings of the causal precursor calculation
		settings = "precursor.%s.%s.logs.%s" % (dataset,tmp_identifier,tmp_identifier)
		exec("import %s; importlib.reload(%s); from %s import *" % tuple([settings]*3))
		# actually only the precursor_dict is needed
		for field_name,prec_identifier in precursor_dict.items():
			# store all potential precursors
			details['potential_ts'][field_name] = xa.open_dataset(tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_ts_'+field_name+'.nc')

	# load the target
	target = xa.open_dataset(index_name+'.nc')[index_name].loc[cutoff[0]:cutoff[1]].squeeze()
	target_raw = xa.open_dataset(index_name+'.nc')[index_name+'_raw'].loc[cutoff[0]:cutoff[1]].squeeze()
	forecast = {}

	# initialize result dicts
	result = {'anomaly':{}, 'absolute':{}}
	result_classi = {'anomaly':{}, 'absolute':{}}
	for style in result.keys():
		for model in model_dict.keys():
			result[style][model] = pd.DataFrame(columns=['year','observed','forecast','upper','lower'])
		for classi,details in classi_dict.items():
			result_classi[style][classi] = {}
			for lvl in details.keys():
				result_classi[style][classi][lvl] = pd.DataFrame(columns=['year','observed','forecast','forecast_proba','upper','lower'])


	for set_id in split_sets.index.levels[0]:
		# select split-set and identify training years
		splitSet = split_sets.loc[set_id]
		years_train = all_years[np.where(splitSet['train'])]

		# initialize a data_frame for the predictors and target
		data = {'anomaly':splitSet.copy(), 'absolute':splitSet.copy()}
		data['anomaly']['target'] = list(target.loc[target_indices].values)
		data['absolute']['target'] = list(target_raw.loc[target_indices].values)

		predictors = []
		for tmp_identifier,details in potential_prec_dict.items():
			for field_name,potential_ts in details['potential_ts'].items():
				n = 0
				for id_ in potential_ts.ID.values:
					if field_name+'_'+str(int(id_)) in details['causal'][set_id]:
						# add causal precursors
						data['anomaly']['pred_'+field_name+'_'+str(int(id_))] = potential_ts['actor'].loc[:,set_id,id_,:].squeeze().loc[prec_indices].values
						data['absolute']['pred_'+field_name+'_'+str(int(id_))] = potential_ts['actor_raw'].loc[:,set_id,id_,:].squeeze().loc[prec_indices].values
						predictors.append('pred_'+field_name+'_'+str(int(id_)))
						n += 1
					if n == details['limit_regs'][field_name]:
						# stop if precusor restriction number is reached
						break

		for style in result.keys():
			y = data[style].loc[(data[style].train),'target']

			# optional: weight years
			if model_dict[model]['weighting'] is not None:
				weights = y.copy() * 0 + 1
				for weight_lvl,weight_add in model_dict[model]['weighting'].items():
					weights[y > np.nanpercentile(y,weight_lvl)] += weight_add
			else:
				weights = y.copy() * 0 + 1

			for model in model_dict.keys():
				# fit each model and save forecasted test year
				clf = model_dict[model]['clf']
				clf.fit(data[style].loc[(data[style].train),predictors], y, sample_weight=weights)
				frc = clf.predict(data[style].loc[(data[style].test),predictors])[0]
				stdev = np.std(clf.predict(data[style].loc[(data[style].train),predictors]) - y)
				result[style][model] = result[style][model].append({'year':set_id,
								'observed':data[style].loc[(data[style].test),'target'].values[0],
								'forecast':frc, 'upper':frc + stdev, 'lower':frc - stdev}, ignore_index=True)

			# do the same for all classifiers
			for classi,details in classi_dict.items():
				for lvl in details.keys():
					clf = details[lvl]['clf']
					clf.fit(data[style].loc[(data[style].train),predictors], np.array(y>np.percentile(y,lvl*100),np.int))
					frc = clf.predict(data[style].loc[(data[style].test),predictors])[0]
					frc_prob = clf.predict_proba(data[style].loc[(data[style].test),predictors])[0]
					result_classi[style][classi][lvl] = result_classi[style][classi][lvl].append({'year':set_id,
									'observed':data[style].loc[(data[style].test),'target'].values[0],
									'forecast':frc, 'forecast_proba':frc_prob[1]}, ignore_index=True)

		# store results
		forecast[set_id] = dict(anomaly=data['anomaly'], absolute=data['absolute'])

	plt.close('all')
	with PdfPages('prediction/'+prediction+'_skill.pdf') as pdf:
		for style,result_loc in result.items():
			print('******** '+style+' ********')
			forecast_summary[style] = {}
			legend_elements=[]
			for model,result_tmp in result_loc.items():
				spearmanr = scipy.stats.spearmanr(result_tmp.observed,result_tmp.forecast)
				spearman_str = str(round(spearmanr[0],2))
				if spearmanr[1] < 0.1: spearman_str += '*'
				if spearmanr[1] < 0.05: spearman_str += '*'
				pearsonr = scipy.stats.pearsonr(result_tmp.observed,result_tmp.forecast)
				pearson_str = str(round(pearsonr[0],2))
				if pearsonr[1] < 0.1: pearson_str += '*'
				if pearsonr[1] < 0.05: pearson_str += '*'
				legend_elements.append(Line2D([0], [0], color='w', label='pearson r: '+pearson_str+ ' spearman r: '+spearman_str))

			fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(8,3))
			ax = axes[0]
			ax.annotate('a', xy=(0.05,0.93), xycoords='axes fraction', ha='center', va='center', fontweight='bold', fontsize=13)
			ax.plot(result_loc['lr'].year, result_loc['lr'].observed, color='k')
			for model,result_tmp in result_loc.items():
				ax.plot(result_tmp.year, result_tmp.forecast, label=model, color = model_dict[model]['color'])
				ax.fill_between(result_tmp.year, result_tmp.lower, result_tmp.upper, alpha= 0.4, color = model_dict[model]['color'])

			ax.set_ylabel(index_longname+' ('+style+')')
			ax.set_xlabel('Year')
			ax.legend(handles=legend_elements, loc='best', frameon=True, facecolor='w', ncol=1, framealpha=0.6, edgecolor='none', fontsize = 9).set_zorder(1)

			axes[1].annotate('b', xy=(0.05,0.93), xycoords='axes fraction', ha='center', va='center', fontweight='bold', fontsize=13)
			forecast_summary[style][model] = plot_ROC(result_tmp.observed, result_tmp.forecast, axes[1], n_boot = 1000)
			plt.tight_layout(); plt.savefig('prediction/'+prediction+'_lrTS_'+style+'.png', dpi=300)
			plt.tight_layout(); pdf.savefig(); plt.close()

			fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
			ax.plot(result_loc['lr'].year, (result_loc['lr'].observed - result_loc['lr'].observed.mean()) / result_loc['lr'].observed.std(), color='k')

			for label,dicc in result_classi[style].items():
				tmp = np.zeros([len(result_loc['lr'].year)])
				for lvl,result_tmp in dicc.items():
					if lvl in [0.33,0.5,0.66]:
						tmp = result_tmp.forecast_proba
						ax.plot(result_tmp.year, (tmp - tmp.mean()) /tmp.std(), color = classi_dict[label][lvl]['color'], label= label +' '+ str(round(scipy.stats.spearmanr(result_tmp.observed,tmp)[0],3)))

			ax.set_ylabel(index_longname+' ('+style+')')
			ax.set_xlabel('Year')
			ax.legend(loc='best', frameon=True, facecolor='w', ncol=1, framealpha=1, edgecolor='w', fontsize = 7).set_zorder(1)
			plt.tight_layout(); plt.savefig('prediction/'+prediction+'_classiTS_'+style+'.png', dpi=300)
			plt.tight_layout(); pdf.savefig(); plt.close()

			print('******** '+style+' ********')
			for model,result_tmp in result_loc.items():
				print('******** '+model+' ********')
				spearmanr = scipy.stats.spearmanr(result_tmp.observed,result_tmp.forecast)
				print('spearman r: '+str(round(spearmanr[0],3))+' pval: '+str(round(spearmanr[1],2)))
				pearsonr = scipy.stats.pearsonr(result_tmp.observed,result_tmp.forecast)
				print('pearson r: '+str(round(pearsonr[0],3))+' pval: '+str(round(pearsonr[1],2)))

				mse = np.sum((result_tmp.forecast - result_tmp.observed)**2) / float(len(result_tmp.observed))
				mse_ref = np.sum((result_tmp.observed - np.mean(result_tmp.observed))**2) / float(len(result_tmp.observed))
				msss = 1- mse / mse_ref
				print('msss: '+str(round(msss*100,1)))

				forecast_summary[style][model]['pearsonr'] = pearsonr
				forecast_summary[style][model]['spearmanr'] = spearmanr
				forecast_summary[style][model]['msss'] = msss


			for label,dicc in result_classi[style].items():
				x,y = [],[]
				print(label)
				forecast_summary[style][label] = {}
				fig,ax = plt.subplots(nrows=1, figsize=(4,3))
				for lvl,result_tmp in dicc.items():
					predict = np.array(result_tmp.forecast,np.bool)
					probs = np.array(result_tmp.forecast_proba)
					events = result_tmp.observed.values > np.percentile(result_tmp.observed,lvl*100)
					forecast_summary[style][label][lvl] = get_scores(events,predict,probs,lvl)

				lvl_names = {0.25:'25% least active seasons', 0.33:'below $33^{rd}$ percentile seasons', 0.5:'above median seasons', 0.66:'above $66^{th}$ percentile seasons', 0.75:'25% most active seasons'}

				fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(7,5))
				for i,lvl in enumerate([0.5,0.33,0.66]):
					ax = axes.flatten()[i]
					ax.set_title(lvl_names[lvl])
					events = result_tmp.observed.values > np.percentile(result_tmp.observed,lvl*100)
					probs = result_tmp.forecast_proba
					lvl_ = lvl
					if lvl == 0.33:
						events = 1 - events
						probs = 1- probs
						lvl_ = 1 - lvl_

					ax.plot([0,1],[(xx + 1-lvl_) *0.5 for xx in [0,1]],color='gray',linestyle='-')
					ax.plot([1-lvl_,1-lvl_],[0,1],color='k',linestyle='-')
					ax.plot([0,1],[1-lvl_,1-lvl_],color='k',linestyle='-')
					ax.plot([0,1],[0,1],color='gray',linestyle='--')
					ax.fill_between([1-lvl_,1],[1,1],[1-lvl_, (1+1-lvl_) *0.5], color='green', alpha=0.5)
					ax.fill_between([1-lvl_,1],[1-lvl_, (1+1-lvl_) *0.5],[1-lvl_, 1-lvl_], color='limegreen', alpha=0.5)
					ax.fill_between([0,1-lvl_],[0,0],[(1-lvl_)*0.5,1-lvl_], color='green', alpha=0.5)
					ax.fill_between([0,1-lvl_],[(1-lvl_)*0.5,1-lvl_],[1-lvl_,1-lvl_], color='limegreen', alpha=0.5)

					weight,bins = np.histogram(probs,bins=np.linspace(0,1,6))
					#weight,bins = np.histogram(probs,bins=np.percentile(probs, np.linspace(0,1,8)*100))
					mean_predicted_value = np.array([probs[(probs>=x1) & (probs<x2)].mean() for x1,x2 in zip(bins[:-1],bins[1:])])
					fraction_of_positives = np.array([events[(probs>=x1) & (probs<x2)].mean() for x1,x2 in zip(bins[:-1],bins[1:])])
					fraction_of_positives, mean_predicted_value = calibration_curve(events, probs, n_bins=5)
					ax.plot(mean_predicted_value, fraction_of_positives, "-")
					scatter = ax.scatter(mean_predicted_value, fraction_of_positives, s=weight*10, label='bin weight', zorder=100)
					scatter.set_clip_on(False)
					ax.set_ylabel('observed frequency'); ax.set_xlabel('forecast probability')
					ax.set_ylim(0,1); ax.set_xlim(0,1)
					ax.annotate(ascii_lowercase[i], xy=(0.05,0.93), xycoords='axes fraction', ha='center', va='center', fontweight='bold', fontsize=13)
					ax.annotate('BSS='+str(round(forecast_summary[style][label][lvl]['bss'],2))+' ', xy=(1,0.05), xycoords='axes fraction', ha='right', va='center', fontweight='bold')

				ax = axes.flatten()[3]
				ax.axis('off')
				legend_elements=[]
				legend_elements.append(scatter)
				legend_elements.append(Line2D([0], [0], color='gray', linestyle='--', label='perfect reliability'))
				legend_elements.append(Patch(facecolor='darkgreen', edgecolor='darkgreen', alpha=0.3, label='better skill than climatological forecast'))
				legend_elements.append(Line2D([0], [0], color='k', linestyle='-', label='climatology'))
				legend_elements.append(Patch(facecolor='limegreen', alpha=0.3, label='better skill than random guessing'))
				ax.legend(handles=legend_elements, loc='center',fontsize=9,ncol=1, frameon=True, facecolor='w', framealpha=1, edgecolor='w').set_zorder(1)
				plt.tight_layout(); plt.savefig('prediction/'+prediction+'_reliability_'+style+'.png', dpi=300)
				plt.tight_layout(); pdf.savefig(); plt.close()

		save_pkl(forecast_summary, 'prediction/summary/'+prediction+'_skill.pkl')

		# forecast_summary = load_pkl('prediction/summary/'+prediction+'_skill.pkl')














#
