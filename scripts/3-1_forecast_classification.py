# -*- coding: utf-8 -*-
import os,sys,importlib,time

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# define all forecasts that should be tested in a prediction_dict
prediction_dict = {}
for dataset in ['ERA5', 'JRA55']:
	prediction_dict[dataset+ '_lag4_X3_wideClusters_2'] = {'lag':4, 'dataset':dataset,
			'potential_prec_dict': {
				'causal_'+dataset+'_vws_lag4_X3_MSLP_wideClusters_n2': dict(limit_regs= dict(mslp=2)),
				'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
			}}
	#
	# prediction_dict[dataset+ '_lag4_X3_wideClusters_2_mslpOnly'] = {'lag':4, 'dataset':dataset,
	# 		'potential_prec_dict': {
	# 			'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
	# 		}}
	#
	# prediction_dict[dataset+ '_lag4_X3_sstOnly'] = {'lag':4, 'dataset':dataset,
	# 		'potential_prec_dict': {
	# 			'causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1': dict(limit_regs= dict(sst=1)),
	# 		}}

	# prediction_dict[dataset+ '_lag2_X3_v1'] = {'lag':2, 'dataset':dataset,
	# 		'potential_prec_dict': {
	# 			'causal_'+dataset+'_vws_lag2_X3_SST_n2': dict(limit_regs= dict(sst=2)),
	# 			'causal_'+dataset+'_sst_lag2_X3_SST_n1': dict(limit_regs= dict(sst=1)),
	# 		}}

# define statistical models to be tested
model_dict = {
	'lr' : dict(clf = linear_model.LinearRegression(), color = 'm', weighting = None),
	# 'lr-weighted' : dict(clf = linear_model.LinearRegression(), color = 'm', weighting = {66:1,33:1}),
	# 'randomForest' : dict(clf = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100), color = 'orange'),
}

# define classifiers to be tested
classi_dict = {}
# classi_dict['KNN'] = {lvl: dict(clf=KNeighborsClassifier(3), color='b') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['randomForest'] = dict(clf=RandomForestClassifier(max_depth=10, n_estimators=100, max_features=1, class_weight = 'balanced'), color='g')
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', class_weight = {0:1-lvl,1:lvl}), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', class_weight = {0:1-lvl,1:lvl}), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
classi_dict['logReg'] = dict(clf=LogisticRegression(random_state=0, solver='liblinear', multi_class='auto', C=1, class_weight = {-1:0.33,0:0.33,1:0.33}), color='r')
#classi_dict['logReg'] = {lvl: dict(clf=LogisticRegression(C=1, penalty='l2',solver='saga',multi_class='ovr',max_iter=10000, class_weight = 'balanced'), color='r') for lvl in [0.25,0.33,0.5,0.66,0.75]}
#classi_dict['logReg_noWeight'] = {lvl: dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial'), color='c') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# from sklearn.gaussian_process import GaussianProcessClassifier
# from sklearn.gaussian_process.kernels import RBF
# classi_dict['RBF'] = {lvl: dict(clf=GaussianProcessClassifier(1.0 * RBF(1.0)), color='m') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# from sklearn.naive_bayes import GaussianNB
# classi_dict['NaivBayes'] = {lvl: dict(clf=GaussianNB(), color='c') for lvl in [0.25,0.33,0.5,0.66,0.75]}
# classi_dict['QDA'] = dict(clf=QuadraticDiscriminantAnalysis(), color='darkorange')
classi_dict['LogReg'] = dict(clf=LogisticRegression(random_state=0, solver='lbfgs', multi_class='auto', class_weight='balanced'), color='darkorange')#, class_weight = {0:10,1:30,2:30,3:30}),

#classi_dict['SVC'] = {lvl: dict(clf=SVC(kernel='linear', C=10, probability=True,random_state=0), color='darkorange') for lvl in [0.33,0.5,0.66]}

event_dict={
	-1: {'thresh':None, 'freq':0.33},
	0: {'thresh':None, 'freq':0.33},
	1: {'thresh':None, 'freq':0.33},
}

# go through all forecasts to be tested
for prediction,params in prediction_dict.items():
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
			result_classi[style][classi] = pd.DataFrame(columns=['year','observed',-1,0,1])


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
			y_ = data[style].loc[(data[style].train),'target'].copy()
			y = data[style].loc[(data[style].train),'target'].copy()
			y[y_ < np.percentile(y_,33)] = -1
			y[(y_ > np.percentile(y_,33)) & (y_ < np.percentile(y_,66))] = 0
			y[y_ > np.percentile(y_,66)] = 1


			# do the same for all classifiers
			for classi,details in classi_dict.items():
				clf = details['clf']
				clf.fit(data[style].loc[(data[style].train),predictors], y)
				frc_prob = clf.predict_proba(data[style].loc[(data[style].test),predictors])[0]
				result_classi[style][classi] = result_classi[style][classi].append({'year':set_id,
								'observed':data[style].loc[(data[style].test),'target'].values[0], -1:frc_prob[0], 0:frc_prob[1], 1:frc_prob[2]}, ignore_index=True)

		# store results
		forecast[set_id] = dict(anomaly=data['anomaly'], absolute=data['absolute'])

	plt.close('all')
	with PdfPages('prediction/'+prediction+'_skill_classi.pdf') as pdf:
		for style,result_loc in result_classi.items():

			for label in result_loc.keys():
				for column in [-1,0,1]:
					fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
					tmp = result_loc[label]['observed']
					ax.plot(result_loc[label].year, (tmp - tmp.mean()) / tmp.std())
					tmp = result_loc[label][column]
					ax.plot(result_loc[label].year, (tmp - tmp.mean()) / tmp.std())

					ax.set_ylabel(index_longname+' ('+style+')')
					ax.set_xlabel('Year')
					ax.legend(loc='best', frameon=True, facecolor='w', ncol=1, framealpha=1, edgecolor='w', fontsize = 7).set_zorder(1)
					plt.tight_layout(); pdf.savefig(); plt.close()

			forecast_summary[style] = {}
			for label,result_tmp in result_classi[style].items():
				x,y = [],[]
				print(label)
				forecast_summary[style][label] = {}
				y_ = result_tmp.observed.copy().values
				y = result_tmp.observed.copy().values
				y[y_ < np.percentile(y_,33)] = -1
				y[(y_ > np.percentile(y_,33)) & (y_ < np.percentile(y_,66))] = 0
				y[y_ > np.percentile(y_,66)] = 1
				for lvl in [-1,0,1]:
					events = np.array(y == lvl, np.int)
					probs = result_tmp[lvl]
					bs = brier_score_loss(events, probs)
					bs_c = np.sum((events - 0.33)**2) / np.float(len(events))
					bss = (bs_c-bs)/bs_c
					forecast_summary[style][label][lvl] = {'bss':bss}
	asdas

				# lvl_names = {0.25:'25% least active seasons', 0.33:'above 33 percentile seasons', 0.5:'above median seasons', 0.66:'above 66 percentile seasons', 0.75:'25% most active seasons'}
				#
				# fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(7,5))
				# for i,lvl in enumerate([0.5,0.33,0.66]):
				# 	ax = axes.flatten()[i]
				# 	ax.set_title(lvl_names[lvl])
				# 	events = result_tmp.observed.values > np.percentile(result_tmp.observed,lvl*100)
				# 	probs = result_tmp.forecast_proba
				# 	ax.plot([0,1],[(xx + 1-lvl) *0.5 for xx in [0,1]],color='gray',linestyle='-')
				# 	ax.plot([1-lvl,1-lvl],[0,1],color='k',linestyle='-')
				# 	ax.plot([0,1],[1-lvl,1-lvl],color='k',linestyle='-')
				# 	ax.plot([0,1],[0,1],color='gray',linestyle='--')
				# 	ax.fill_between([1-lvl,1],[1,1],[1-lvl, (1+1-lvl) *0.5], color='green', alpha=0.5)
				# 	ax.fill_between([1-lvl,1],[1-lvl, (1+1-lvl) *0.5],[1-lvl, 1-lvl], color='limegreen', alpha=0.5)
				# 	ax.fill_between([0,1-lvl],[0,0],[(1-lvl)*0.5,1-lvl], color='green', alpha=0.5)
				# 	ax.fill_between([0,1-lvl],[(1-lvl)*0.5,1-lvl],[1-lvl,1-lvl], color='limegreen', alpha=0.5)
				#
				# 	weight,bins = np.histogram(probs,bins=np.linspace(0,1,6))
				# 	#weight,bins = np.histogram(probs,bins=np.percentile(probs, np.linspace(0,1,8)*100))
				# 	mean_predicted_value = np.array([probs[(probs>=x1) & (probs<x2)].mean() for x1,x2 in zip(bins[:-1],bins[1:])])
				# 	fraction_of_positives = np.array([events[(probs>=x1) & (probs<x2)].mean() for x1,x2 in zip(bins[:-1],bins[1:])])
				# 	fraction_of_positives, mean_predicted_value = calibration_curve(events, probs, n_bins=5)
				# 	ax.plot(mean_predicted_value, fraction_of_positives, "-")
				# 	scatter = ax.scatter(mean_predicted_value, fraction_of_positives, s=weight*10, label='bin weight', zorder=100)
				# 	scatter.set_clip_on(False)
				# 	ax.set_ylabel('observed frequency'); ax.set_xlabel('forecast probability')
				# 	ax.set_ylim(0,1); ax.set_xlim(0,1)
				# 	ax.annotate(ascii_lowercase[i], xy=(0.05,0.93), xycoords='axes fraction', ha='center', va='center', fontweight='bold', fontsize=13)
				# 	ax.annotate('BSS='+str(round(forecast_summary[style][label][lvl]['bss'],2))+' ', xy=(1,0.05), xycoords='axes fraction', ha='right', va='center', fontweight='bold')
				#
				# ax = axes.flatten()[3]
				# ax.axis('off')
				# legend_elements=[]
				# legend_elements.append(scatter)
				# legend_elements.append(Line2D([0], [0], color='gray', linestyle='--', label='prefect reliability'))
				# legend_elements.append(Patch(facecolor='darkgreen', edgecolor='darkgreen', alpha=0.3, label='skill compared to climatology'))
				# legend_elements.append(Line2D([0], [0], color='k', linestyle='-', label='climatology'))
				# legend_elements.append(Patch(facecolor='limegreen', alpha=0.3, label='skill compared to random guessing'))
				# ax.legend(handles=legend_elements, loc='center',fontsize=9,ncol=1, frameon=True, facecolor='w', framealpha=1, edgecolor='w').set_zorder(1)
				# plt.tight_layout(); pdf.savefig(); plt.close()

		# save_pkl(forecast_summary, 'prediction/summary/'+prediction+'_skill.pkl')

		# forecast_summary = load_pkl('prediction/summary/'+prediction+'_skill.pkl')















#
