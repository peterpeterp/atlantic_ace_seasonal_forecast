# -*- coding: utf-8 -*-
import os,sys,importlib,gc

sys.path.append('scripts')


identifier = "causal_ERA5_vws_lag4_X3_MSLP_wideClusters"
identifier = "causal_ERA5_vws_lag2_X3_SST_n2_2018"
set_id = 2018
central_longitude = -50

import __helper_init; importlib.reload(__helper_init); from __helper_init import *

dataset = identifier.split('_')[1]
sys.path.append(tmp_path+'precursor/'+dataset+'/'+identifier+'/logs/')
exec("import %s; importlib.reload(%s); from %s import *" % tuple([identifier]*3))

causal_pkl = load_pkl(tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal_extended.pkl')
for key,item in causal_pkl[set_id].items():
	globals()[key] = item


# plot only the relevant nodes
fig,axes = plt.subplots(nrows=2, figsize=(4,4), gridspec_kw={'height_ratios':[1,7]})
axes[0].axis('off')
axes[0].annotate('b', xy=(0,1), xycoords='axes fraction', fontweight='bold', fontsize=20)
tp.plot_graph(
  fig_ax=(fig,axes[1]),
  val_matrix= results['val_matrix'][nodes_to_plot][:,nodes_to_plot],
  link_matrix= sig['link_matrix'][nodes_to_plot][:,nodes_to_plot],
  var_names=[nn.upper().replace('_',' ')+'\n\n' for nn in var_names[nodes_to_plot]],
  link_colorbar_label='link strength',
  node_colorbar_label='autodependence',
  show_colorbar=True,
  arrow_linewidth=20,
  arrowhead_size=20,
  node_label_size=15,
  link_label_fontsize=0,
  node_size = 10,
  save_name = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal_overview_'+str(set_id)+'.png'
  )

# plot only the relevant nodes
fig,axes = plt.subplots(nrows=2, figsize=(4,4), gridspec_kw={'height_ratios':[1,7]})
axes[0].axis('off')
axes[0].annotate('b', xy=(0,1), xycoords='axes fraction', fontweight='bold', fontsize=20)
tp.plot_graph(
  fig_ax=(fig,axes[1]),
  val_matrix= results['val_matrix'],
  link_matrix= sig['link_matrix'],
  var_names=[nn.upper().replace('_',' ')+'\n\n' for nn in var_names],
  link_colorbar_label='link strength',
  node_colorbar_label='autodependence',
  show_colorbar=True,
  arrow_linewidth=20,
  arrowhead_size=20,
  node_label_size=15,
  link_label_fontsize=0,
  node_size = 10,
  save_name = tmp_path+'precursor/'+dataset+'/'+identifier+'/'+identifier+'_causal_overview_'+str(set_id)+'_full.png'
  )




#
