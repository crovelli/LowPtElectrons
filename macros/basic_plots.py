import numpy as np
import matplotlib
matplotlib.use('Agg')
import uproot
import matplotlib.pyplot as plt
import root_numpy
import rootpy
import rootpy.plotting as rplt
import json
import pandas as pd
from matplotlib import rc
from pdb import set_trace
import os
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from baseline import baseline
import cosmetics

debug = False
print 'Getting the data'
from datasets import dataset_names, tag, get_data_sync, kmeans_weighter, training_selection

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.mkdirs(plots)

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

all_data = {}
for dataset in dataset_names:
   print 'loading', dataset
   all_data[dataset] = pd.DataFrame(
      get_data_sync(
         dataset, 
         ['is_e', 'is_e_not_matched', 'is_other',
          'gen_pt', 'gen_eta', 'trk_pt'
          ]
         )
      )

plt.figure(figsize=[8,8])
for to_plot, nbins in [
   ('gen_pt', 30),
   ('gen_eta', 30),
   ('trk_pt', 30),]:
   plt.clf()
   for dataset, sample in all_data.iteritems():
      electrons = sample[sample.is_e]
      plt.hist(
         electrons[to_plot], bins=nbins, 
         range=cosmetics.ranges[to_plot],
         histtype='step', normed=True,
         label = dataset_names[dataset],
         )
   plt.xlabel(cosmetics.beauty[to_plot])
   plt.ylabel('Fraction')
   plt.legend(loc='best')
   plt.plot()
   plt.savefig('%s/electrons_%s.png' % (plots, to_plot))
   plt.savefig('%s/electrons_%s.pdf' % (plots, to_plot))
   plt.clf()

   plt.clf()
   for dataset, sample in all_data.iteritems():
      electrons = sample[(sample.is_e) & (sample.trk_pt > 0)]
      plt.hist(
         electrons[to_plot], bins=nbins, 
         range=cosmetics.ranges[to_plot],
         histtype='step', normed=True,
         label = dataset_names[dataset],
         )
   plt.xlabel(cosmetics.beauty[to_plot])
   plt.ylabel('Fraction')
   plt.legend(loc='best')
   plt.plot()
   plt.savefig('%s/electrons_withTrk_%s.png' % (plots, to_plot))
   plt.savefig('%s/electrons_withTrk_%s.pdf' % (plots, to_plot))
   plt.clf()


from features import *
features = id_features+new_features
additional = id_additional

multi_dim_branches = ['gsf_ecal_cluster_ematrix', 'ktf_ecal_cluster_ematrix']
dict_data = get_data_sync(
   'all', 
   features+labeling+additional+multi_dim_branches
)
data = pd.DataFrame(
   {i : dict_data[i] for i in dict_data if i not in multi_dim_branches}
)
multi_dim = {i : dict_data[i] for i in dict_data if i in multi_dim_branches}
data_mask = np.invert(data.is_e_not_matched) & training_selection(data)
for key in multi_dim:
   multi_dim[key] = multi_dim[key][data_mask]
data = data[data_mask] #remove non-matched electrons
data['training_out'] = -1
data['log_trkpt'] = np.log10(data.trk_pt)
#convert bools to integers
for c in features:
   if data[c].dtype == np.dtype('bool'):
      data[c] = data[c].astype(int)

#apply pt-eta reweighting
weights = kmeans_weighter(
   data[['log_trkpt', 'trk_eta']],
   '%s/kmeans_all_weighter.plk' % mods
   ) 
data['weight'] = weights*np.invert(data.is_e) + data.is_e
print 'making plots'

for feat in multi_dim_branches:
   vals = {}
   for dataset in [
      {'name' : 'electrons',
       'mask' : data.is_e,
       'weight' : data[data.is_e].weight},
      {'name' : 'tracks',
       'mask' : np.invert(data.is_e),
       'weight' : data[np.invert(data.is_e)].weight},
      ]:
      plt.clf()
      plt.title(feat.replace('_', ' '))
      masked = multi_dim[feat][dataset['mask']]
      sum_val = masked.sum(axis=-1).sum(axis=-1)
      mask = np.invert(sum_val == 0)
      masked = masked[mask]
      sum_val = sum_val[mask]
      masked /= sum_val[:,None,None]
      heatmap = np.average(masked, axis=0, weights=dataset['weight'][mask])
      vals[dataset['name']] = heatmap
      plt.imshow(heatmap, cmap='viridis', interpolation='nearest')
      plt.colorbar()
      plt.savefig('%s/%s_%s.png' % (plots, dataset['name'], feat))
      plt.savefig('%s/%s_%s.pdf' % (plots, dataset['name'], feat))      
      plt.clf()
   #make ratios
   ratio = (vals['electrons']/vals['tracks'])-1
   plt.clf()
   plt.title(feat.replace('_', ' '))
   plt.imshow(ratio, cmap='RdBu', interpolation='nearest', vmin=-1, vmax=1)
   plt.colorbar()
   plt.savefig('%s/ratio_%s_%s.png' % (plots, dataset['name'], feat))
   plt.savefig('%s/ratio_%s_%s.pdf' % (plots, dataset['name'], feat))      
   plt.clf()

exit()
      
for to_plot in features:
   plt.clf()
   electrons = data[data.is_e]
   tracks = data[np.invert(data.is_e) & np.invert(data.is_e_not_matched)]
   plt.hist(
      electrons[to_plot], bins=50, 
      weights=electrons.weight,
      range=cosmetics.ranges.get(to_plot, None),
      histtype='step', normed=True,
      label = 'Electrons',
      )
   plt.hist(
      tracks[to_plot], bins=50, 
      weights=tracks.weight,
      range=cosmetics.ranges.get(to_plot, None),
      histtype='step', normed=True,
      label = 'Tracks',
      )
   plt.xlabel(cosmetics.beauty.get(to_plot, to_plot.replace('_', ' ')))
   plt.ylabel('Fraction')
   plt.legend(loc='best')
   plt.plot()
   plt.savefig('%s/electrons_vs_tracks_%s.png' % (plots, to_plot))
   plt.savefig('%s/electrons_vs_tracks_%s.pdf' % (plots, to_plot))
   plt.clf()
