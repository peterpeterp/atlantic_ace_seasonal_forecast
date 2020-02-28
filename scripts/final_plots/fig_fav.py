# -*- coding: utf-8 -*-
import os,sys,importlib

sys.path.append('scripts')

import __helper_init; importlib.reload(__helper_init); from __helper_init import *

for file_name in [
	# tmp_path+"precursor/ERAint_setAll_vws_lag1_v1/vws_causal_fld.nc",
	# tmp_path+"precursor/ERAint_setAll_sst_lag1_v1/sst_causal_fld.nc",
	# tmp_path+"instantCond/data/ERAint_setAll_instant_lag0_v1_fld.nc",
	# tmp_path+"instantCond/data/JRA55_setAll_instant_lag0_vLocal_fld.nc",
	tmp_path+"instantCond/ERA5_setAll_instant_lag0_vLocal_fld.nc",
	]:

	nc = xa.open_dataset(file_name)
	lons = nc.lon.values; lons[lons < 0] += 360; nc.coords['lon'] = lons

	lons, lats = np.meshgrid(nc.lon.values,nc.lat.values)
	lons[lons<0] += 360

	plt.close('all')
	fig,axes = plt.subplots(nrows=1, ncols=3, figsize=(6,4), subplot_kw={'projection': ccrs.Robinson(central_longitude=0)}, gridspec_kw = {'width_ratios':[4,4,0.5]})

	plot_dict = {
		'a': {'name':'$SST_{MDR}$','var':'sst','id':1, 'ax':axes[0]},
		'b': {'name':'$VWS_{MDR}$','var':'vws','id':1, 'ax':axes[1]}
	}

	for letter, details in plot_dict.items():
		region = nc['region'].loc[:,:,:,details['var'],details['id']].squeeze()
		corr = nc['corr'].loc[:,:,:,details['var']].squeeze()
		ax = details['ax']

		lon_mesh, lat_mesh = np.meshgrid(region.lon,region.lat)
		c_lon = np.mean(lon_mesh[region==1])
		c_lat = np.mean(lat_mesh[region==1])

		ax.coastlines(color='gray');
		ax.set_extent([c_lon-50,c_lon+50,c_lat-50,c_lat+50], crs=ccrs.PlateCarree())
		im = ax.pcolormesh(region.lon,region.lat,corr, transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-0.8, vmax=0.8)
		ax.contour(region.lon,region.lat,region, levels=[0.5,1.5], transform=ccrs.PlateCarree())
		ax.annotate(details['name'], xy=(0.0, 0), xycoords='axes fraction', fontsize=15, ha='left', va='bottom')
		ax.annotate(letter, xy=(0.0, 1), xycoords='axes fraction', fontsize=15, fontweight='bold', ha='left', va='bottom')

	axes[-1].outline_patch.set_edgecolor('white')
	ax_cb = fig.add_axes([0.9,0.2,0.05,0.6], zorder=0)
	cb = fig.colorbar(im, cax=ax_cb, orientation='vertical',label='correlation')

	plt.savefig(file_name.replace('.nc','_fig.png'), bbox_inches='tight', transparent=False, dpi=300)
	plt.savefig(file_name.replace('.nc','_fig.pdf'), bbox_inches='tight', transparent=False); plt.close()
#
