# -*- coding: utf-8 -*-
import os,sys,importlib,gc

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# check whether the script is called from outside
# the "identifier" is the name of a setting file in scripts/params/'+dataset+'/causal/ containing all relevant parameters
try:
	identifier = sys.argv[1]
except:
	identifier = "causal_ERA5_vws_lag2_X3_SST_n2_2018"

# just a convention: the second part of the identifier denotes the dataset
dataset = identifier.split('_')[1]

# create a folder for the results
if os.path.isdir(tmp_path+'precursor/'+dataset+'/'+identifier)==False:
	os.system('mkdir '+tmp_path+'precursor/'+dataset+'/'+identifier)
	os.system('mkdir '+tmp_path+'precursor/'+dataset+'/'+identifier+'/logs')

# save used parameters in this folder
params_in = open('scripts/params/'+dataset+'/causal/'+identifier+'.py', 'r').read()
out = open(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/'+identifier+'.py', 'w')
out.write(params_in)
out.close()

# load parameters
# all varaibles that are required in this script are defined there
sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

# initialize a dictionary in which the names of causal precursors are saved
# the names refer to names of potential precursors from the files specified in precursor_dict
causal_precurors = {}
causal_results = {}

# open a pdf file to which plots are appended for each training set
with PdfPages(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal_overview.pdf') as pdf:
	for set_id in split_sets.index.levels[0]:
		print('treating test year ',set_id)

		# select a training set and identify tarining years and the relevant time-stamps
		splitSet = split_sets.loc[set_id]
		years_train = all_years[np.where(splitSet['train'])]
		target_timeSteps_train = [tst for tst in target.time.values if int(str(tst).split('-')[0]) in years_train]
		target_train = target.loc[target_timeSteps_train]

		# create a list of var_names and a network_data and network_mask array containing the traget
		var_names = [target_name]
		network_data = target_train.values
		network_mask = np.ones(network_data.shape,dtype='bool')
		# mask all irrelevant time steps
		network_mask[start_step:: time_cycle] = False

		for field_name,prec_identifier in precursor_dict.items():
			# load potential precursor time series
			potential_precursors_ts = xa.open_dataset(tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_ts_'+field_name+'.nc')['actor']

			# select relevant time stamps
			prec_timeSteps_train = [tst for tst in potential_precursors_ts.time.values if int(str(tst).split('-')[0]) in years_train]
			potPrec_train = potential_precursors_ts.loc[lag,set_id,:,prec_timeSteps_train]

			# go through actors
			for i in potPrec_train.ID.values:
				tmp = potPrec_train.loc[i,:].values
				use_region = True
				# check if all elements are finite
				if np.isfinite(tmp).all() == False:
					use_region = False
				# (add other conditions?)
				if use_region:
					# add data to the network_data
					network_data = numpy.column_stack((network_data, tmp))
					# create a mask for the relevant months
					tmp_mask = np.zeros(target_train.shape,dtype='bool')
					tmp_mask[start_step-lag:: time_cycle] = False
					network_mask = numpy.column_stack((network_mask, tmp_mask))
					var_names = np.append(var_names, field_name+'_'+str(i))

		def run_tigramite_(n_act,pmci_alpha):
			'''
			little function running tigramite
			I just put it in a function because it is caaled twice... not really elegant
			'''
			# select subset of the network_data
			network_input = pp.DataFrame(data=network_data[:,:n_act+1], mask=network_mask[:,:n_act+1])

			# do pcmci
			parcorr = ParCorr(significance='analytic',use_mask =True,mask_type='y',verbosity=0)
			pcmci = PCMCI(dataframe=network_input,cond_ind_test=parcorr,var_names=var_names,selected_variables=None,verbosity=0)
			results = pcmci.run_pcmci(tau_min=lag, tau_max=lag, pc_alpha = pmci_alpha)

			# get corrected p-values
			results['q_matrix'] = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], fdr_method='fdr_bh')

			return pcmci,results

		# this is a mode where significance thresholds are gredually increased until precursors are found
		# I implemented this to get more robust results for all training sets
		if mode != 'fixed':

			# n_keep is the number of causal precursors I want to find
			# check whether there are enough potnetial precursors
			if max(n_keep) > len(var_names):
				if min(n_keep) >= len(var_names):
					n_keep_tmp = [len(var_names)]
				else:
					n_keep_tmp = np.arange(len(var_names),min(n_keep),-1)
			else:
				n_keep_tmp = n_keep.copy()

			found = False
			# first try to find desired number of precursors nPrecursors_options[0]
			# if that is not succesful, reduce nPrecursors
			for nPrecursors in nPrecursors_options:
				# first try with the corrected p-values (q_matrix)
				# if that is not succesful, try with p_matrix /!\ not recommended to allow for this option
				for pq_matrix in pq_matrix_options:
					# first try with all precursors
					# if that isn't succesful, reduce this number /!\ not recommended to allow small numbers of n_act
					for n_act in n_keep_tmp:
						print('testing '+str(n_act)+' actors '+pq_matrix)
						# check whether a result can be found with the highest possible significance threshold
						pcmci,results = run_tigramite_(n_act,max(pmci_alpha_list))
						sig = pcmci._return_significant_parents(pq_matrix=results[pq_matrix],val_matrix=results['val_matrix'], alpha_level=max(pmci_alpha_list))
						# if with the highest possible significance threshold the desired number of causal precursors is found (or more)
						if len([tt[0] for tt in sig['parents'][0] if tt[0]!=0]) > 0:
							# gradually increase pmci_alpha starting from a low value
							for pmci_alpha in pmci_alpha_list:
								pcmci,results = run_tigramite_(n_act,pmci_alpha)
								sig = pcmci._return_significant_parents(pq_matrix=results[pq_matrix],val_matrix=results['val_matrix'], alpha_level=pmci_alpha)
								if len([tt[0] for tt in sig['parents'][0] if tt[0]!=0]) >= nPrecursors:
									found = True
									# if the desired number of causal precursors is found (or more) stop all iterations
									break
						if found:
							break
					if found:
						break
				if found:
					break

		# print some results
		if pq_matrix == 'p_matrix':
			pcmci._print_significant_links(p_matrix = results['p_matrix'],val_matrix = results['val_matrix'], alpha_level = pmci_alpha)
			pcmci._print_significant_links(p_matrix = results['p_matrix'],q_matrix = results['q_matrix'],val_matrix = results['val_matrix'], alpha_level = 0.3)
		if pq_matrix == 'q_matrix':
			pcmci._print_significant_links(p_matrix = results['p_matrix'],q_matrix = results['q_matrix'],val_matrix = results['val_matrix'], alpha_level = pmci_alpha)

		# store casual precursor names
		causal_precurors[set_id] = {}
		causal_prec_id = [iii[0] for iii in sig['parents'][0]]
		for par in causal_prec_id:
			causal_precurors[set_id][var_names[par]] = dict(pval=results['p_matrix'][par][0][-1], qval=results['q_matrix'][par][0][-1])

		print(causal_precurors[set_id])

		# save into a pickle (in principle this could also be done once in the end, but it doesn't hurt here)
		save_pkl(causal_precurors, tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pkl')

		# plot the network
		nodes_to_plot = [0] + causal_prec_id
		fig,ax = plt.subplots(nrows=1, figsize=(4,4))
		tp.plot_graph(
		  fig_ax=(fig,ax),
		  val_matrix= results['val_matrix'],
		  link_matrix= sig['link_matrix'],
		  var_names=[nn.upper().replace('_',' ') for nn in var_names],
		  link_colorbar_label='cross-MCI',
		  node_colorbar_label='auto-MCI',
		  show_colorbar=True,
		  arrow_linewidth=20,
		  arrowhead_size=20,
		  node_label_size=10,
		  link_label_fontsize=0,
		  node_size = 10,
		  )
		ax.set_title(set_id)
		pdf.savefig(transparent=True); plt.close()

		# plot only the relevant nodes
		fig,ax = plt.subplots(nrows=1, figsize=(4,4))
		tp.plot_graph(
		  fig_ax=(fig,ax),
		  val_matrix= results['val_matrix'][nodes_to_plot][:,nodes_to_plot],
		  link_matrix= sig['link_matrix'][nodes_to_plot][:,nodes_to_plot],
		  var_names=[nn.upper().replace('_',' ') for nn in var_names[nodes_to_plot]],
		  link_colorbar_label='cross-MCI',
		  node_colorbar_label='auto-MCI',
		  show_colorbar=True,
		  arrow_linewidth=20,
		  arrowhead_size=20,
		  node_label_size=10,
		  link_label_fontsize=0,
		  node_size = 10,
		  )
		ax.set_title(set_id)
		pdf.savefig(transparent=True); plt.close()

		# save into a pickle (in principle this could also be done once in the end, but it doesn't hurt here)
		causal_results[set_id] = {'results':results,'var_names':var_names, 'sig':sig, 'nodes_to_plot':nodes_to_plot}
		save_pkl(causal_results, tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal_extended.pkl')

		# plot all potential precursors
		for field_name,precursor in precursor_dict.items():
			fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(6,3), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0)})
			ax.set_title(str(set_id) + ' potential', fontsize=12,fontweight='bold')
			ax.coastlines(color='gray')
			# plot clusters
			clusters = xa.open_dataset(tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_field_'+field_name+'.nc')['clusters']
			clusters = clusters.loc[:,:,lag,set_id]
			lons = clusters.lon.values; lons[lons < 0] += 360; clusters.coords['lon'] = lons
			to_plot, lon_cyc = add_cyclic_point(clusters.copy(), clusters.lon.values)
			im = ax.contourf(lon_cyc, clusters.lat, to_plot, cmap='jet', transform=ccrs.PlateCarree(), vmin=0, vmax=5, alpha=0.5)
			# add labels
			for id_ in np.unique(clusters.values)[np.isfinite(np.unique(clusters.values))]:
				var_name = field_name+'_'+str(int(id_))
				points = np.where(clusters == int(var_name.split('_')[1]))
				c_lat = np.median(clusters.lat[points[0]])
				c_lon = np.median(clusters.lon[points[1]])
				ax.text(c_lon-10,c_lat+10,var_name, fontsize=8 ,color='k',fontweight='bold' ,transform=ccrs.PlateCarree(), va='bottom', ha='center', rotation=0)
				ax.plot([c_lon,c_lon-10],[c_lat,c_lat+10],color='k',transform=ccrs.PlateCarree())

			# add a ghost-layer such that the whole globe is shown
			to_plot[:,:] = 1
			im = ax.contour(lon_cyc, clusters.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
			plt.tight_layout(); pdf.savefig(); plt.close()

		fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(6,3), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0)})
		ax.set_title(str(set_id) + ' causal', fontsize=12,fontweight='bold')
		ax.coastlines(color='gray')
		for field_name,precursor in precursor_dict.items():
			clusters = xa.open_dataset(tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_field_'+field_name+'.nc')['clusters']
			clusters = clusters.loc[:,:,lag,set_id]
			lons = clusters.lon.values; lons[lons < 0] += 360; clusters.coords['lon'] = lons
			for var_name in causal_precurors[set_id]:
				if var_name.split('_')[0] == field_name and '_' in var_name:
					points = np.where(clusters == int(var_name.split('_')[1]))
					if len(points[0]) > 1:
						c_lat = np.median(clusters.lat[points[0]])
						c_lon = np.median(clusters.lon[points[1]])
						ax.text(c_lon-10,c_lat+10, ' '.join([var_name,'\np:',str(round(causal_precurors[set_id][var_name]['pval'],3)),'\nq:',str(round(causal_precurors[set_id][var_name]['qval'],3))]), fontsize=8 ,color='k',fontweight='bold' ,transform=ccrs.PlateCarree(), va='bottom', ha='center', rotation=0)
						ax.plot([c_lon,c_lon-10],[c_lat,c_lat+10],color='k',transform=ccrs.PlateCarree())
						to_plot = clusters.values.copy()
						to_plot[to_plot != int(var_name.split('_')[1])] = np.nan
						to_plot, lon_cyc = add_cyclic_point(to_plot.copy(), clusters.lon.values)
						im = ax.contourf(lon_cyc, clusters.lat, to_plot, cmap='jet', transform=ccrs.PlateCarree(), vmin=0, vmax=5, alpha=0.5)

		# add a ghost-layer such that the whole globe is shown
		to_plot[:,:] = 1
		im = ax.contour(lon_cyc, clusters.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
		plt.tight_layout(); pdf.savefig(); plt.close()

# finally plot a summary figure showing the number of times a grid-cell is part of a potential precursor region and the number of times a grid-cell is part of a causal precursor region
from plot_precursor_reoccurring import plot_reoccurring
for field_name,prec_identifier in precursor_dict.items():
	plot_reoccurring(in_file_name = tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_field_'+field_name+'.nc',
					plot_file_name = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pdf',
					causal_pkl = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pkl')



'''
for dataset in JRA55 ERA5; do sbatch job.sh scripts/2-2_precursors_causal.py causal_${dataset}_sst_lag4_X3_SST_n1; done;
for dataset in JRA55 ERA5; do sbatch job.sh scripts/2-2_precursors_causal.py causal_${dataset}_sst_lag2_X3_SST_n1; done;
for dataset in JRA55 ERA5; do sbatch job.sh scripts/2-2_precursors_causal.py causal_${dataset}_vws_lag4_X3_MSLP_wideClusters_n2; done;
for dataset in JRA55 ERA5; do sbatch job.sh scripts/2-2_precursors_causal.py causal_${dataset}_vws_lag2_X3_SST_n2; done;


sbatch job.sh scripts/2-2_precursors_causal.py causal_JRA55_vws_lag4_X3_MSLP_smallRegs_n2
sbatch job.sh scripts/2-2_precursors_causal.py causal_ERA5_vws_lag4_X3_MSLP_smallRegs_n2
'''







#
