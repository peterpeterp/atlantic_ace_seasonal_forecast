# -*- coding: utf-8 -*-
import os,sys,importlib,time

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# create a folder for the results
if os.path.isdir(tmp_path+'additional')==False:
	os.system('mkdir -p '+tmp_path+'additional')


for dataset,curoff_years in zip(['ERA5','JRA55'],[[1979,2018],[1958,2018]]):

	for col,identifier in zip(range(2),['causal_'+dataset+'_vws_lag4_X3_MSLP_wideClusters_n2','causal_'+dataset+'_sst_lag4_X3_SST_dropSmall_n1']):


		if 'vws' in identifier:
			fav = '$VWS_{MDR}$'
			letter1 = 'a'
			letter2 = 'c'
		else:
			fav = '$SST_{MDR}$'
			letter1 = 'b'
			letter2 = 'd'


		cutoff_here = [np.datetime64(str(cutoff_years[0])+'-01-01'),np.datetime64(str(cutoff_years[1])+'-12-31')]

		dataset = identifier.split('_')[1]
		sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
		exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

		field_name = list(precursor_dict.keys())[0]
		in_file_name = tmp_path+'precursor/'+dataset+'/'+precursor_dict[field_name]+'/'+precursor_dict[field_name]+'_potential_field_'+field_name+'.nc'

		dataset = in_file_name.split('/')[-3]
		field_name = in_file_name.split('_')[-1].split('.nc')[0]
		out_fld = xa.open_dataset(in_file_name)

		clusters = out_fld['clusters']
		corr = out_fld['corr']

		causal = load_pkl(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pkl')

		############################################################
		# count reoccurences
		############################################################

		fields = {}
		for set_id in clusters.set_id.values:
			tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
			if field_name.upper()+' Potential' not in fields.keys():
				fields[field_name.upper()+' Potential'] = tmp.copy()
				fields[field_name.upper()+' Potential'].values[:] = 0.0
				fields[field_name.upper()+' Causal'] = fields[field_name.upper()+' Potential'].copy()

			for id_ in np.unique(tmp.values)[np.isfinite(np.unique(tmp.values))]:
				tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
				tmp.values[tmp.values!=id_] = 0
				tmp.values[tmp.values==id_] = 1
				reg_grids = np.where(tmp.values)
				direction = np.sign(np.mean(corr.loc[:,:,:,set_id].values.squeeze()*tmp.values))
				#print(id_,direction)
				fields[field_name.upper()+' Potential'].values += tmp.copy() * direction
				if field_name+'_'+str(int(id_)) in causal[set_id].keys():
					fields[field_name.upper()+' Causal'].values += tmp.copy() * direction

		############################################################
		# define a fingerprint as
		############################################################

		fingerprint = fields[field_name.upper()+' Causal'].copy()
		fingerprint.values[(np.abs(fingerprint) < 18)] = np.nan

		cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#0F2080','w','#F5793A'])

		fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,2.2), subplot_kw={'projection': ccrs.Robinson(central_longitude=-50)})
		ax.coastlines(color='gray', zorder=3)
		im = ax.pcolormesh(fingerprint.lon, fingerprint.lat, fingerprint, cmap=cmap, vmin=-len(clusters.set_id), vmax=len(clusters.set_id), transform=ccrs.PlateCarree())
		ax_all = fig.add_axes([0,0,1,1], zorder=0); ax_all.axis('off')
		#ax_all.annotate(field_name+' Precursors', xy=(0.5,0.9), xycoords='axes fraction',fontsize=12,fontweight='bold',backgroundcolor='w', ha='center')
		#ax.set_title(field_name+' Precursors')
		ax.annotate('Reoccurring Robust Precursors\nof '+fav+' in '+field_name.upper(), xy=(0.5,1.05), xycoords='axes fraction',fontsize=8, ha='center')
		#ax.annotate(letter1, xy=(-0.1,1.1), xycoords='axes fraction',fontsize=10,fontweight='bold', ha='center')

		# ghost plot works better than set_extent
		ghost = fingerprint.copy()
		ghost[:] = 1
		ax.contourf(ghost.lon, ghost.lat, ghost, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
		plt.savefig(tmp_path+'additional/' + identifier + '_fingerprint.pdf', bboxes='tight'); plt.close()


		############################################################
		# get time series of this fingerprint
		############################################################

		# create a latitude-weighting array
		lat_weight = np.cos(np.deg2rad(fingerprint.lat.values))
		lat_weight_array = np.repeat(lat_weight[np.newaxis,:], fingerprint.shape[1], 0).T

		sys.path.append(tmp_path+'precursor/'+dataset+'/'+precursor_dict[field_name]+'/logs/')
		exec("import %s; importlib.reload(%s); from %s import *" % tuple([precursor_dict[field_name]]*3))

		precursor = precursor_dict[field_name]

		var_in = Dataset(precursor['file'], 'r')
		if 'level' in precursor.keys():
			lat_grid, lon_grid, var_vals, var_anoms, var_standard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff_here,time_cycle,level=precursor['level'])
		else:
			lat_grid, lon_grid, var_vals, var_anoms, var_standard, nan_IDs = load_nc(precursor['file'],precursor['var'],cutoff_here,time_cycle)

		mask_pos = fingerprint.copy()
		mask_pos.values[mask_pos < 0 ] = np.nan
		mask_pos.values[mask_pos > 0 ] = 1
		fingerprint_pos = var_vals[:,0,0].squeeze()
		fingerprint_pos.values = np.nanmean(var_vals.loc[:,fingerprint.lat,fingerprint.lon] * np.expand_dims(mask_pos,0) * lat_weight_array, axis=(1,2))
		fingerprint_pos.values *= np.nansum(mask_pos) / np.nansum(np.expand_dims(mask_pos,0) * lat_weight_array)

		mask_neg = fingerprint.copy()
		mask_neg.values[mask_neg > 0 ] = np.nan
		mask_neg.values[mask_neg < 0 ] = 1
		fingerprint_neg = var_vals[:,0,0].squeeze()
		fingerprint_neg.values = np.nanmean(var_vals.loc[:,fingerprint.lat,fingerprint.lon] * np.expand_dims(mask_neg,0) * lat_weight_array, axis=(1,2))
		fingerprint_neg.values *= np.nansum(mask_neg) / np.nansum(np.expand_dims(mask_neg,0) * lat_weight_array)

		fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3))
		for y,color in zip([fingerprint_neg,fingerprint_pos],['#0F2080','#F5793A']):
			if np.any(np.isfinite(y)):
				xx = var_vals.loc[cutoff[0]:cutoff[1],:,:].time[2::12]
				ax.plot(var_vals.time[2::12].dt.year.values, y[2::12], color)
				timeX = sm.add_constant(var_vals.loc[cutoff[0]:cutoff[1],:,:].time[2::12].dt.year.values)
				est = sm.OLS(y.loc[cutoff[0]:cutoff[1]][2::12].values, timeX)
				est2 = est.fit()
				print(est2.summary())
				ax.plot([1979,2018],[est2.params[0]+1979*est2.params[1],est2.params[0]+2018*est2.params[1]], color, label='p-value=%s'%(round(est2.pvalues[1],3)))
				if cutoff_here != cutoff:
					timeX = sm.add_constant(var_vals.time[2::12].dt.year.values)
					est = sm.OLS(y[2::12].values, timeX)
					est2 = est.fit()
					print(est2.summary())
					ax.plot([1958,2018],[est2.params[0]+1958*est2.params[1],est2.params[0]+2018*est2.params[1]], color, linestyle='--', label='p-value=%s'%(round(est2.pvalues[1],3)))

		ax.legend(fontsize=7)
		ax.set_ylabel(field_name.upper())
		#ax.annotate(letter2, xy=(-0.2,1.1), xycoords='axes fraction',fontsize=14,fontweight='bold', ha='center')
		plt.tight_layout()
		plt.savefig(tmp_path+'additional/' + identifier + '_trend.pdf', bboxes='tight'); plt.close()

		xa.Dataset({'pos':fingerprint_pos, 'neg':fingerprint_neg}).to_netcdf(tmp_path+'additional/' + identifier + '_fingerprint_ts.nc')


#
