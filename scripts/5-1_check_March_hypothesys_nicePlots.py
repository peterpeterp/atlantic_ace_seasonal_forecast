# -*- coding: utf-8 -*-
import os,sys,importlib,time

# load required functions
sys.path.append('scripts')
import __helper_init; importlib.reload(__helper_init); from __helper_init import *

# create a folder for the results
if os.path.isdir(tmp_path+'hypothesys')==False:
	os.system('mkdir -p '+tmp_path+'hypothesys')


mslp = xa.open_dataset(data_path+'/ERA5/ERA5_MSLP-SST_1979-2018_grid1x1.nc')['var151']
sst = xa.open_dataset(data_path+'/ERA5/ERA5_MSLP-SST_1979-2018_grid1x1.nc')['var34']
vws = xa.open_dataset(data_path+'/ERA5/ERA5_VWS_1979-2018_grid1x1.nc')['vws']
u = xa.open_dataset(data_path+'/ERA5/ERA5_U_200-850mbar_1979-2018_grid1x1.nc')['var131']


weights = np.cos(np.deg2rad(mslp.lat))
weights.name = "weights"

var_names = ['diff MSLP']
network_data = mslp.loc[:,-35:-15,50:100].weighted(weights).mean(('lat','lon')) - mslp.loc[:,-35:-15,170:210].weighted(weights).mean(('lat','lon'))
network_mask = np.ones(network_data.shape,dtype='bool')
for mon in range(1,6):
	network_mask[mon:: 12] = False

# var_names.append('mslp_Indian')
# network_data = numpy.column_stack((network_data, mslp.loc[:,-25:-5,50:100].mean(('lat','lon'))))
# tmp_mask = np.ones(network_data[:,0].shape,dtype='bool')
# tmp_mask[3:: 12] = False
# network_mask = numpy.column_stack((network_mask, tmp_mask))
#
# var_names.append('mslp_Pacific')
# network_data = numpy.column_stack((network_data, mslp.loc[:,-25:-5,170:210].mean(('lat','lon'))))
# tmp_mask = np.zeros(network_data[:,0].shape,dtype='bool')
# tmp_mask[3:: 12] = False
# network_mask = numpy.column_stack((network_mask, tmp_mask))

var_names.append('Trade winds\n(850mbar)')
network_data = numpy.column_stack((network_data, -u.loc[:,85000,-15:5,140:200].squeeze().weighted(weights).mean(('lat','lon'))))
tmp_mask = np.ones(network_data[:,0].shape,dtype='bool')
for mon in range(1,6):
	tmp_mask[mon:: 12] = False
network_mask = numpy.column_stack((network_mask, tmp_mask))

var_names.append('SST ElNino3.4')
network_data = numpy.column_stack((network_data, sst.loc[:,-5:5,190:240].weighted(weights).mean(('lat','lon'))))
tmp_mask = np.ones(network_data[:,0].shape,dtype='bool')
for mon in range(1,6):
	tmp_mask[mon:: 12] = False
network_mask = numpy.column_stack((network_mask, tmp_mask))
#
# var_names.append('vws_MDR')
# network_data = numpy.column_stack((network_data, mslp.loc[:,10:20,280:350].weighted(weights).mean(('lat','lon'))))
# tmp_mask = np.ones(network_data[:,0].shape,dtype='bool')
# for mon in range(5,7):
# 	tmp_mask[mon:: 12] = False
# network_mask = numpy.column_stack((network_mask, tmp_mask))

for j in range(12):
	network_data[j::12] = scipy.signal.detrend(network_data[j::12], axis = 0)
#calculate anomalies
for i in range(int(12)):
	network_data[i::int(12)] = (network_data[i::int(12)] - np.mean(network_data[i::int(12)], axis = 0))

# select subset of the network_data
network_input = pp.DataFrame(data=network_data, mask=network_mask)

# do pcmci
# parcorr = ParCorr(significance='analytic',verbosity=0)
parcorr = ParCorr(significance='analytic',use_mask =True,mask_type='y',verbosity=0)
pcmci = PCMCI(dataframe=network_input,cond_ind_test=parcorr,var_names=var_names,selected_variables=None,verbosity=0)
results = pcmci.run_pcmci(tau_min=1, tau_max=4, pc_alpha = None)
# get corrected p-values
results['q_matrix'] = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], fdr_method='fdr_bh')
sig = pcmci._return_significant_parents(pq_matrix=results['q_matrix'],val_matrix=results['val_matrix'], alpha_level=0.1)
pcmci._print_significant_links(p_matrix = results['p_matrix'],q_matrix = results['q_matrix'],val_matrix = results['val_matrix'], alpha_level = 0.2)



plt.close('all')
with PdfPages(tmp_path+'hypothesys/regions_and_network.pdf') as pdf:

	fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4), subplot_kw={'projection': ccrs.Robinson(central_longitude=160)})
	#ax.set_extent([-180,180,-50,50])
	transform = ccrs.PlateCarree()._as_mpl_transform(ax)
	ax.coastlines(color='gray')
	ax.add_geometries([Polygon([(-170,-5),(-170,5),(-120,5),(-120,-5)])], facecolor='none', edgecolor='lightgreen', crs=ccrs.PlateCarree())
	ax.text(-140,0,'ElNino3.4', color='lightgreen', transform=ccrs.PlateCarree(), ha='center', va='center')

	ax.add_geometries([Polygon([(50,-35),(50,-15),(100,-15),(100,-35)])], facecolor='none', edgecolor='m', crs=ccrs.PlateCarree())
	ax.text(75,-25,'MSLP$_{Indian Ocean}$', color='m', transform=ccrs.PlateCarree(), ha='center', va='center')

	ax.add_geometries([Polygon([(170,-35),(170,-15),(210,-15),(210,-35)])], facecolor='none', edgecolor='c', crs=ccrs.PlateCarree())
	ax.text(190,-25,'MSLP$_{Pacific}$', color='c', transform=ccrs.PlateCarree(), ha='center', va='center')

	ax.annotate('', xy=(170,-35), xytext=(100,-35), xycoords=transform, size=20, arrowprops=dict(facecolor='k', ec = 'none',arrowstyle="fancy", connectionstyle="arc3,rad=0.1"))
	ax.annotate('$\Delta$MSLP\n(MSLP$_{Indian Ocean}$ - MSLP$_{Pacific}$)', xy=(135,-50), xycoords=transform, ha='center', va='center')



	ax.add_geometries([Polygon([(140,-15),(140,5),(200,5),(200,-15)])], facecolor='none', edgecolor='b', crs=ccrs.PlateCarree())
	ax.annotate('', xy=(145, -13), xytext=(190,-13), xycoords=transform, size=20, arrowprops=dict(facecolor='b', ec = 'none',arrowstyle="fancy"))
	ax.text(170,-3,'Trade winds\n(850mbar)', color='b', transform=ccrs.PlateCarree(), ha='center', va='center')

	# ax.add_geometries([Polygon([(-80,10),(-80,20),(-20,20),(-20,10)])], facecolor='none', edgecolor='r', crs=ccrs.PlateCarree())
	# ax.text(-60,15,'VWS', color='r', transform=ccrs.PlateCarree(), ha='center', va='center')

	to_plot = mslp.loc[:,-60:30,40:280][0,:,:].copy()
	ax.contourf(to_plot.lon, to_plot.lat, to_plot, cmap='plasma', transform=ccrs.PlateCarree(), alpha=0)
	#plt.tight_layout();
	plt.savefig(tmp_path+'hypothesys/regions.png')
	pdf.savefig(bboxes='tight'); plt.close()


	tp.plot_timeseries(data=network_data, mask=network_mask, use_mask=True, datatime=mslp.time.dt.year, var_names=var_names)
	pdf.savefig(transparent=True); plt.close()

	correlations = pcmci.get_lagged_dependencies(tau_max=6)
	lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations, setup_args={'var_names':var_names,
                                    'x_base':5, 'y_base':.5})
	pdf.savefig(transparent=True); plt.close()

	# plot the network
	fig,ax = plt.subplots(nrows=1, figsize=(4,4))
	tp.plot_graph(
	  fig_ax=(fig,ax),
	  val_matrix= results['val_matrix'],
	  link_matrix= sig['link_matrix'],
	  var_names=[nn +'\n\n' for nn in var_names],
	  link_colorbar_label='cross-MCI',
	  node_colorbar_label='auto-MCI',
	  show_colorbar=True,
	  arrow_linewidth=20,
	  arrowhead_size=20,
	  node_label_size=10,
	  link_label_fontsize=10,
	  node_size = 10,
	  )
	plt.savefig(tmp_path+'hypothesys/network_graph.png')
	pdf.savefig(transparent=True); plt.close()

	# Plot time series graph
	tp.plot_time_series_graph(
	    val_matrix=results['val_matrix'],
	    link_matrix=sig['link_matrix'],
	    var_names=[nn +'\n\n' for nn in var_names],
	    link_colorbar_label='MCI',
	    )
	plt.savefig(tmp_path+'hypothesys/ts_graph.png')
	pdf.savefig(); plt.close()











#
