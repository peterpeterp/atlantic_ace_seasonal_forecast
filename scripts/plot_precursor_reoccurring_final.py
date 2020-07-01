# -*- coding: utf-8 -*-
import os,sys,importlib,time

sys.path.append('scripts')

import __helper_init; importlib.reload(__helper_init); from __helper_init import *

now = time.time()


central_longitude = -50
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#0F2080','w','#F5793A'])


list_of_things_to_plot = [
	{'$VWS_{MDR}$':'causal_ERA5_vws_lag2_X3_SST_n2','$SST_{MDR}$':'causal_ERA5_sst_lag2_X3_SST_n1'},
	# {'$VWS_{MDR}$':'causal_JRA55_vws_lag2_X3_SST_n2','$SST_{MDR}$':'causal_JRA55_sst_lag2_X3_SST_n1'},
	#{'$VWS_{MDR}$':'causal_ERA5_vws_lag4_X3_MSLP_wideClusters_n2','$SST_{MDR}$':'causal_ERA5_sst_lag4_X3_SST_dropSmall_n1'},
	# {'$VWS_{MDR}$':'causal_JRA55_vws_lag4_X3_MSLP_wideClusters_n2','$SST_{MDR}$':'causal_JRA55_sst_lag4_X3_SST_dropSmall_n1'},
	# {'$VWS_{MDR}$':'causal_ERA5_vws_lag4_X3_SST_wideClusters_n2','$SST_{MDR}$':'causal_ERA5_sst_lag4_X3_SST_dropSmall_n1'},
	]

for identifier_dict in list_of_things_to_plot:

	fields = {}
	for fav,identifier in identifier_dict.items():
		dataset = identifier.split('_')[1]
		sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
		exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

		for field_name,prec_identifier in precursor_dict.items():
			in_file_name = tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_field_'+field_name+'.nc'
			causal_pkl = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pkl'

			dataset = in_file_name.split('/')[-3]
			field_name = in_file_name.split('_')[-1].split('.nc')[0]
			out_fld = xa.open_dataset(in_file_name)

			clusters = out_fld['clusters']
			corr = out_fld['corr']

			do_causal = False
			if causal_pkl is not None:
				if os.path.isfile(causal_pkl):
					causal = load_pkl(causal_pkl)
					do_causal = True

			clusters.lon.values[clusters.lon.values<0] += 360
			for set_id in clusters.set_id.values:
				tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
				if 'Potential Precursors of '+fav+' in '+field_name.upper() not in fields.keys():
					fields['Potential Precursors of '+fav+' in '+field_name.upper()] = tmp.copy()
					fields['Potential Precursors of '+fav+' in '+field_name.upper()].values[:] = 0.0
					if do_causal:
						fields['Robust Precursors of '+fav+' in '+field_name.upper()] = fields['Potential Precursors of '+fav+' in '+field_name.upper()].copy()

				for id_ in np.unique(tmp.values)[np.isfinite(np.unique(tmp.values))]:
					tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
					tmp.values[tmp.values!=id_] = 0
					tmp.values[tmp.values==id_] = 1
					reg_grids = np.where(tmp.values)
					direction = np.sign(np.mean(corr.loc[:,:,:,set_id].values.squeeze()*tmp.values))
					#print(id_,direction)
					fields['Potential Precursors of '+fav+' in '+field_name.upper()].values += tmp.copy() * direction
					if do_causal:
						if field_name+'_'+str(int(id_)) in causal[set_id].keys():
							fields['Robust Precursors of '+fav+' in '+field_name.upper()].values += tmp.copy() * direction

	plt.close('all')
	fig,axes = plt.subplots(nrows=3, ncols=2, figsize=(8,4*120./180.+1), subplot_kw={'projection': ccrs.Robinson(central_longitude=central_longitude)}, gridspec_kw = {'height_ratios':[4,4,1]})

	plot_dict = {
		axes[0,0] : dict(field_name = [fld for fld in fields.keys() if 'Potential' in fld and '$VWS_{MDR}$' in fld][0], letter='a'),
		axes[1,0] : dict(field_name = [fld for fld in fields.keys() if 'Robust' in fld and '$VWS_{MDR}$' in fld][0], letter='c'),
		axes[0,1] : dict(field_name = [fld for fld in fields.keys() if 'Potential' in fld and '$SST_{MDR}$' in fld][0], letter='b'),
		axes[1,1] : dict(field_name = [fld for fld in fields.keys() if 'Robust' in fld and '$SST_{MDR}$' in fld][0], letter='d'),
	}

	for ax,details in plot_dict.items():

		to_plot = fields[details['field_name']].copy()
		to_plot.values[to_plot.values==0] = np.nan

		ax.coastlines(color='gray', zorder=3)
		# im = ax.pcolormesh(to_plot.lon.values, to_plot.lat.values, to_plot.values, cmap=cmap, vmin=-len(clusters.set_id), vmax=len(clusters.set_id), transform=ccrs.PlateCarree())
		lons = to_plot.lon.values.copy()
		lons[lons<0] += 360
		cyclic_data, cyclic_lons = add_cyclic_point(to_plot.values, coord=lons)
		im = ax.contourf(cyclic_lons, to_plot.lat, cyclic_data, cmap=cmap, levels=np.arange(-len(clusters.set_id),len(clusters.set_id)+1,5), transform=ccrs.PlateCarree())
		ax.annotate(details['field_name'], xy=(0.5,1.05), xycoords='axes fraction',fontsize=8, ha='center')
		ax.annotate(details['letter'], xy=(0,1), xycoords='axes fraction',fontsize=8,fontweight='bold', ha='center')
		#ax.set_title(details['field_name']+' Precursors', fontsize=8)

		# ghost plot works better than set_extent
		to_plot[:] = 1
		ax.pcolormesh(to_plot.lon, to_plot.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)

	for ax in axes[-1,:]: ax.outline_patch.set_edgecolor('white')
	ax_cb = fig.add_axes([0.2,0.13,0.6,0.05], zorder=0)
	cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal',label='identification count')
	cb.set_ticks([-len(clusters.set_id),-30,-20,-10,0,10,20,30,len(clusters.set_id)])
	cb.set_ticklabels([len(clusters.set_id),30,20,10,0,10,20,30,len(clusters.set_id)])
	ax_cb.annotate('negative links', xy=(0.05,0.3), xycoords='axes fraction', ha='left', color='w')
	ax_cb.annotate('positive links', xy=(0.95,0.3), xycoords='axes fraction', ha='right')
	# plt.tight_layout()
	plt.savefig(tmp_path+'precursor/'+dataset+'/'+'_'.join([dataset,'lag'+str(lag),'causal'])+'.png', dpi=300)
	plt.savefig(tmp_path+'precursor/'+dataset+'/'+'_'.join([dataset,'lag'+str(lag),'causal'])+'.pdf')











	#
