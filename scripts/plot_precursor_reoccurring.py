# -*- coding: utf-8 -*-
import os,sys,importlib,time

sys.path.append('scripts')

import __helper_init; importlib.reload(__helper_init); from __helper_init import *


def plot_reoccurring(in_file_name, plot_file_name=None, central_longitude=0, causal_pkl=None):

	if plot_file_name is None:
		plot_file_name = in_file_name.replace('.nc','_reoccurring_'+str(central_longitude)+'.pdf')

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

	fields = {}

	clusters.lon.values[clusters.lon.values<0] += 360
	for set_id in clusters.set_id.values:
		tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
		if field_name.upper()+' Potential' not in fields.keys():
			fields[field_name.upper()+' Potential'] = tmp.copy()
			fields[field_name.upper()+' Potential'].values[:] = 0.0
			if do_causal:
				fields[field_name.upper()+' Causal'] = fields[field_name.upper()+' Potential'].copy()

		for id_ in np.unique(tmp.values)[np.isfinite(np.unique(tmp.values))]:
			tmp = clusters.loc[:,:,:,set_id].squeeze().copy()
			tmp.values[tmp.values!=id_] = 0
			tmp.values[tmp.values==id_] = 1
			reg_grids = np.where(tmp.values)
			direction = np.sign(np.mean(corr.loc[:,:,:,set_id].values.squeeze()*tmp.values))
			#print(id_,direction)
			fields[field_name.upper()+' Potential'].values += tmp.copy() * direction
			if do_causal:
				if field_name+'_'+str(int(id_)) in causal[set_id].keys():
					fields[field_name.upper()+' Causal'].values += tmp.copy() * direction

	cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#02a368','w','#a30250'])

	plt.close('all')
	with PdfPages(plot_file_name) as pdf:
		for field_name,field in fields.items():

			if np.sum(field.values) != 0:

				to_plot = field.copy()
				to_plot.values[to_plot.values==0] = np.nan

				fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(4,1.8*(max(to_plot.lat)-min(to_plot.lat))/180.+1), subplot_kw={'projection': ccrs.Robinson(central_longitude=central_longitude)})
				ax.coastlines(color='gray', zorder=3)
				im = ax.pcolormesh(to_plot.lon, to_plot.lat, to_plot, cmap=cmap, vmin=-len(clusters.set_id), vmax=len(clusters.set_id), transform=ccrs.PlateCarree())
				ax_all = fig.add_axes([0,0,1,1], zorder=0); ax_all.axis('off')
				ax_all.annotate(field_name+' Precursors', xy=(0.5,0.9), xycoords='axes fraction',fontsize=12,fontweight='bold',backgroundcolor='w', ha='center')
				#ax.set_title(field_name+' Precursors')

				# ghost plot works better than set_extent
				to_plot[:] = 1
				ax.contourf(to_plot.lon, to_plot.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
				plt.tight_layout(); pdf.savefig(); plt.close()

				fig,cax = plt.subplots(nrows=1, ncols=1, figsize=(4,1)); cax.axis('off')
				ax_cb = fig.add_axes([0.1,0.44,0.8,0.24], zorder=0)
				cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal',label='identification count')
				cb.set_ticks([-len(clusters.set_id),-30,-20,-10,0,10,20,30,len(clusters.set_id)])
				cb.set_ticklabels([len(clusters.set_id),30,20,10,0,10,20,30,len(clusters.set_id)])
				cax.annotate('negative links', xy=(0.1,0.95), xycoords='axes fraction', ha='left')
				cax.annotate('positive links', xy=(0.9,0.95), xycoords='axes fraction', ha='right')
				plt.savefig('test.pdf')
				plt.tight_layout(); pdf.savefig(transparent=True); plt.close()


if __name__ == '__main__':
	now = time.time()

	identifier = 'causal_ERA5_vws_lag4_X3_MSLP_wideClusters_n2'
	#identifier = 'causal_ERA5_vws_lag2_X3_SST_n2'
	central_longitude = -50

	dataset = identifier.split('_')[1]
	sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
	exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

	for field_name,prec_identifier in precursor_dict.items():
		plot_reoccurring(in_file_name = tmp_path+'precursor/'+dataset+'/'+prec_identifier+'/'+prec_identifier+'_potential_field_'+field_name+'.nc',
						plot_file_name = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pdf',
						causal_pkl = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal.pkl',
						central_longitude=central_longitude)

	#
