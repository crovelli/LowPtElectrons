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
from datasets import dataset_names, tag, get_data_sync, kmeans_weighter, training_selection, pre_process_data, target_dataset

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.mkdirs(plots)

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

##all_data = {}
##for dataset in dataset_names:
##   print 'loading', dataset
##   all_data[dataset] = pd.DataFrame(
##      get_data_sync(
##         dataset, 
##         ['is_e', 'is_e_not_matched', 'is_other',
##          'gen_pt', 'gen_eta', 'trk_pt'
##          ]
##         )
##      )
##
##plt.figure(figsize=[8,8])
##for to_plot, nbins in [
##   ('gen_pt', 30),
##   ('gen_eta', 30),
##   ('trk_pt', 30),]:
##   plt.clf()
##   for dataset, sample in all_data.iteritems():
##      electrons = sample[sample.is_e]
##      plt.hist(
##         electrons[to_plot], bins=nbins, 
##         range=cosmetics.ranges[to_plot],
##         histtype='step', normed=True,
##         label = dataset_names[dataset],
##         )
##   plt.xlabel(cosmetics.beauty[to_plot])
##   plt.ylabel('Fraction')
##   plt.legend(loc='best')
##   plt.plot()
##   try : plt.savefig('%s/electrons_%s.png' % (plots, to_plot))
##   except : pass
##   try : plt.savefig('%s/electrons_%s.pdf' % (plots, to_plot))
##   except : pass
##   plt.clf()
##
##   plt.clf()
##   for dataset, sample in all_data.iteritems():
##      electrons = sample[(sample.is_e) & (sample.trk_pt > 0)]
##      plt.hist(
##         electrons[to_plot], bins=nbins, 
##         range=cosmetics.ranges[to_plot],
##         histtype='step', normed=True,
##         label = dataset_names[dataset],
##         )
##   plt.xlabel(cosmetics.beauty[to_plot])
##   plt.ylabel('Fraction')
##   plt.legend(loc='best')
##   plt.plot()
##   try : plt.savefig('%s/electrons_withTrk_%s.png' % (plots, to_plot))
##   except : pass
##   try : plt.savefig('%s/electrons_withTrk_%s.pdf' % (plots, to_plot))
##   except : pass
##   plt.clf()


from features import *
features = id_features+new_features+seed_features+improved_seed_features
features = list(set(features))
additional = id_additional

multi_dim_branches = ['gsf_ecal_cluster_ematrix', 'ktf_ecal_cluster_ematrix']
data, multi_dim = pre_process_data(
   target_dataset,
   features+labeling+additional+multi_dim_branches
)
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
      try : plt.savefig('%s/%s_%s.png' % (plots, dataset['name'], feat))
      except : pass
      try : plt.savefig('%s/%s_%s.pdf' % (plots, dataset['name'], feat))
      except : pass
      plt.clf()
   #make ratios
   ratio = (vals['electrons']/vals['tracks'])-1
   plt.clf()
   plt.title(feat.replace('_', ' '))
   plt.imshow(ratio, cmap='RdBu', interpolation='nearest', vmin=-1, vmax=1)
   plt.colorbar()
   try : plt.savefig('%s/ratio_%s_%s.png' % (plots, dataset['name'], feat))
   except : pass
   try : plt.savefig('%s/ratio_%s_%s.pdf' % (plots, dataset['name'], feat))
   except : pass
   plt.clf()


      
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
   try : plt.savefig('%s/electrons_vs_tracks_%s.png' % (plots, to_plot))
   except : pass
   try : plt.savefig('%s/electrons_vs_tracks_%s.pdf' % (plots, to_plot))
   except : pass
   plt.clf()
