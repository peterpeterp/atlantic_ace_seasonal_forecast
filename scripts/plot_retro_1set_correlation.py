# -*- coding: utf-8 -*-
import os,sys,importlib,gc

sys.path.append('scripts')


identifier = "ERA5_vws_lag4_X3_MSLP_wideClusters"
identifier = "ERA5_sst_lag4_X3_SST_dropSmall"
identifier = "ERA5_vws_lag2_X3_SST"
identifier = "ERA5_sst_lag2_X3_SST"
set_id = 2018
central_longitude = -50


import __helper_init; importlib.reload(__helper_init); from __helper_init import *

dataset = identifier.split('_')[0]
sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

field_name = list(precursor_dict.keys())[0]
out_fld = xa.open_dataset(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_potential_field_'+field_name+'.nc')
lons = out_fld.lon.values; lons[lons < 0] += 360; out_fld.coords['lon'] = lons

'''
plot correlation
'''
plt.close('all')
for name,precursor in precursor_dict.items():
	with PdfPages(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_potential_field_'+str(set_id)+'_corr.pdf') as pdf:
		fig,axes = plt.subplots(nrows=2, ncols=1, figsize=(6,3.5*(max(out_fld.lat)-min(out_fld.lat))/180.+1), subplot_kw={'projection': ccrs.Robinson(central_longitude=central_longitude)}, gridspec_kw = {'height_ratios':[7,1]})
		ax = axes[0]
		ax.coastlines(color='gray', zorder=3)
		#ax.set_extent([-180,180,min(out_fld.lat), max(out_fld.lat)], crs=ccrs.PlateCarree())

		lat = out_fld.lat
		to_plot, lon = add_cyclic_point(out_fld['corr'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		im = ax.contourf(lon,lat, to_plot, cmap='RdBu_r', vmin=-1, vmax=1, transform=ccrs.PlateCarree(), zorder=1)
		to_plot, lon = add_cyclic_point(out_fld['corr_sig'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		ax.contour(lon, lat, to_plot, levels=[0.5,2],transform=ccrs.PlateCarree(), zorder=10)

		markers = {}
		all_lons = []
		all_lats = []
		clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
		for id_ in np.unique(clusters.values)[np.isfinite(np.unique(clusters.values))]:
			clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
			clusters.values[clusters.values!=id_] = 0
			clusters.values[clusters.values==id_] = 1
			points = np.where(clusters.values)
			if len(points[0]) > 1:
				c_lat = np.median(out_fld.lat.values[points[0]])
				c_lon = np.median(out_fld.lon.values[points[1]])
				if c_lon > 180:
					c_lon -= 360
				markers[int(id_)] = dict(c_lon=c_lon, c_lat=c_lat)
				all_lons.append(c_lon)
				all_lats.append(c_lat)

		clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
		to_plot, lon = add_cyclic_point(clusters.values, coord=clusters.lon)
		ax.contourf(lon, clusters.lat, to_plot, colors='none', levels=np.arange(0.5,len(all_lats)+1.5), hatchcolor='gray', hatches=['///','---','+++','\\\\']*5, transform=ccrs.PlateCarree())

		temp = np.array(all_lons).argsort()
		ranks = numpy.empty_like(temp)
		ranks[temp] = numpy.arange(len(all_lons))
		for reg_id,details in markers.items():
			details['rank_lon'] = np.linspace(-180,180,len(all_lons))[ranks[np.array(all_lons) == details['c_lon']]][0]
			details['y'] = {False:max(lat)+10, True:min(lat)-10}[details['c_lat']<np.mean(all_lats)]

		for rr in np.linspace(-180,180,len(all_lons)):
			for reg_id,details in markers.items():
				if details['rank_lon'] == rr:
					break
			t = ax.annotate(reg_id,   xy=(details['c_lon'], details['c_lat']), fontsize=12, fontweight='bold', ha='center', va='center', color='w', clip_on=False, xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), zorder=30)
			t.set_bbox(dict(facecolor='k', alpha=0.5, edgecolor='w', boxstyle='round'))

		# ghost plot works better than set_extent
		to_plot, lon = add_cyclic_point(out_fld['corr'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		ax.contourf(lon, out_fld.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0.1)
		ax.annotate('a', xy=(0,1), xycoords='axes fraction',fontsize=15,fontweight='bold',backgroundcolor='w')

		ax = axes[1]
		ax.outline_patch.set_edgecolor('white')
		cax = fig.add_axes([0.2,0.15,0.6,0.05])#; cax.axis('off')
		cb = fig.colorbar(im, cax=cax, orientation='horizontal',label='correlation')
		plt.tight_layout(); pdf.savefig(transparent=True); plt.close()


		######################################################################################################################################################

		fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(6,3.5*(max(out_fld.lat)-min(out_fld.lat))/180.+1), subplot_kw={'projection': ccrs.Robinson(central_longitude=central_longitude)})
		ax.coastlines(color='gray', zorder=3)
		#ax.set_extent([-180,180,min(out_fld.lat), max(out_fld.lat)], crs=ccrs.PlateCarree())

		lat = out_fld.lat
		to_plot, lon = add_cyclic_point(out_fld['corr'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		im = ax.contourf(lon,lat, to_plot, cmap='RdBu_r', vmin=-1, vmax=1, transform=ccrs.PlateCarree(), zorder=1)
		to_plot, lon = add_cyclic_point(out_fld['corr_sig'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		ax.contour(lon, lat, to_plot, levels=[0.5,2],transform=ccrs.PlateCarree(), zorder=10)


		markers = {}
		all_lons = []
		all_lats = []
		clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
		for id_ in np.unique(clusters.values)[np.isfinite(np.unique(clusters.values))]:
			clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
			clusters.values[clusters.values!=id_] = 0
			clusters.values[clusters.values==id_] = 1
			points = np.where(clusters.values)
			if len(points[0]) > 1:
				c_lat = np.median(out_fld.lat.values[points[0]])
				c_lon = np.median(out_fld.lon.values[points[1]])
				if c_lon > 180:
					c_lon -= 360
				markers[int(id_)] = dict(c_lon=c_lon, c_lat=c_lat)
				all_lons.append(c_lon)
				all_lats.append(c_lat)

		clusters = out_fld['clusters'].loc[:,:,lag,set_id].squeeze().copy()
		to_plot, lon = add_cyclic_point(clusters.values, coord=clusters.lon)
		ax.contourf(lon, clusters.lat, to_plot, colors='none', levels=np.arange(0.5,len(all_lats)+1.5), hatchcolor='gray', hatches=['///','---','+++','\\\\']*5, transform=ccrs.PlateCarree())

		temp = np.array(all_lons).argsort()
		ranks = numpy.empty_like(temp)
		ranks[temp] = numpy.arange(len(all_lons))
		for reg_id,details in markers.items():
			details['rank_lon'] = np.linspace(-180,180,len(all_lons))[ranks[np.array(all_lons) == details['c_lon']]][0]
			details['y'] = {False:max(lat)+10, True:min(lat)-10}[details['c_lat']<np.mean(all_lats)]

		for rr in np.linspace(-180,180,len(all_lons)):
			for reg_id,details in markers.items():
				if details['rank_lon'] == rr:
					break
			t = ax.annotate(reg_id,   xy=(details['c_lon'], details['c_lat']), fontsize=12, fontweight='bold', ha='center', va='center', color='w', clip_on=False, xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), zorder=30)
			t.set_bbox(dict(facecolor='k', alpha=0.5, edgecolor='w', boxstyle='round'))

		# ghost plot works better than set_extent
		to_plot, lon = add_cyclic_point(out_fld['corr'].loc[:,:,lag,set_id].values, coord=out_fld.lon)
		ax.contourf(lon, out_fld.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
		plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

		######################################################################################################################################################

		fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(6,3), subplot_kw={'projection': ccrs.Robinson(central_longitude=0)})
		ax.coastlines(color='gray'); ax.set_extent([-180,180,min(out_fld.lat), max(out_fld.lat)], crs=ccrs.PlateCarree())
		cmap_repeat = mpl.colors.LinearSegmentedColormap.from_list("", list(["r", "b", "g", 'm', 'orange', 'c','gray', 'darkmagenta', 'yellow','k']*10)[:len(list(markers.keys()))])
		ax.pcolormesh(lon, clusters.lat, clusters.values, cmap=cmap_repeat, transform=ccrs.PlateCarree())
		for rr in np.linspace(-180,180,len(all_lons)):
			for reg_id,details in markers.items():
				if details['rank_lon'] == rr:
					break
			t = ax.annotate(reg_id,  xytext=(details['rank_lon'],details['y']), xy=(details['c_lon'], details['c_lat']), arrowprops=dict(edgecolor='k', arrowstyle= '->'), fontsize=10, fontweight='bold', ha='center', va='center', color='k', clip_on=False, xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), zorder=30)
			t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='k', boxstyle='round'))

		ax_all = fig.add_axes([0,0,1,1], zorder=0); ax_all.axis('off')
		ax_all.annotate(name.upper(), xy=(0.025,0.05), xycoords='axes fraction',fontsize=15,fontweight='bold',backgroundcolor='w')
		plt.tight_layout(); pdf.savefig(transparent=True); plt.close()


		to_plot, lon = add_cyclic_point(clusters.values, coord=clusters.lon)
		for reg_id,details in markers.items():
			fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(6,3), subplot_kw={'projection': ccrs.Robinson(central_longitude=0)})
			ax.coastlines(color='gray'); ax.set_extent([-180,180,min(out_fld.lat), max(out_fld.lat)], crs=ccrs.PlateCarree())
			ax.contourf(lon, clusters.lat, to_plot, colors=['#03a17b'], levels=[reg_id-0.5,reg_id+0.5], transform=ccrs.PlateCarree())
			ax_all.annotate(name.upper(), xy=(0.025,0.05), xycoords='figure fraction',fontsize=15,fontweight='bold',backgroundcolor='w')
			plt.tight_layout(); pdf.savefig(transparent=True); plt.close()

			fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(3,3), subplot_kw={'projection': ccrs.Robinson(central_longitude=details['c_lon'])})
			ax.coastlines(color='gray');
			ax.set_extent([details['c_lon']-50,details['c_lon']+50,details['c_lat']-50,details['c_lat']+50], crs=ccrs.PlateCarree())
			ax.contourf(lon, clusters.lat, to_plot, colors=['#03a17b'], levels=[reg_id-0.5,reg_id+0.5], transform=ccrs.PlateCarree())
			ax_all.annotate(name.upper(), xy=(0.025,0.05), xycoords='figure fraction',fontsize=15,fontweight='bold',backgroundcolor='w')
			pdf.savefig(bbox_inches='tight', transparent=True); plt.close()







#
