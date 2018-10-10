import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
import os
from matplotlib.colors import LogNorm

parser = ArgumentParser()
parser.add_argument(
   '--nbins', default=400, type=int
)
parser.add_argument(
   '--nthreads', default=10, type=int
)
args = parser.parse_args()
dataset = 'all'

import matplotlib.pyplot as plt
#import ROOT
import uproot
import json
#import rootpy
#import json
import pandas as pd
from matplotlib import rc
from pdb import set_trace
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import get_data, tag
import os

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.makedirs(plots)

data = pd.DataFrame(
   get_data(dataset, ['trk_pt', 'trk_eta', 'is_e', 'is_e_not_matched', 'is_other'])
)
data = data[np.invert(data.is_e_not_matched)] #remove non-matched electrons
#remove things that do not yield tracks
data = data[(data.trk_pt > 0) & (np.abs(data.trk_eta) < 2.4) & (data.trk_pt < 15)]
data['log_trkpt'] = np.log10(data.trk_pt)

overall_scale = data.shape[0]/float(data.is_e.sum())
reweight_feats = ['log_trkpt', 'trk_eta']

from sklearn.cluster import KMeans
clusterizer = KMeans(n_clusters=args.nbins, n_jobs=-2)
clusterizer.fit(data[data.is_e][reweight_feats])

data['cluster'] = clusterizer.predict(data[reweight_feats])
weights = {}
for cluster, group in data.groupby('cluster'):
   weight = group.shape[0]/float(group.is_e.sum())
   weights[cluster] = weight

from sklearn.externals import joblib
joblib.dump(
   clusterizer, 
   '%s/kmeans_%s_weighter.plk' % (mods, dataset), 
   compress=True
)
weights['features'] = reweight_feats
with open('%s/kmeans_%s_weighter.json' % (mods, dataset), 'w') as ww:
   json.dump(weights, ww)

#vectorize(excluded={2})
apply_weight = np.vectorize(lambda x, y: y.get(x), excluded={2})
data['weight'] = data.is_e*apply_weight(data.cluster, weights)+np.invert(data.is_e)

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = .02     # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = data.log_trkpt.min() - 1, data.log_trkpt.max() + 1
y_min, y_max = data.trk_eta.min() - 1  , data.trk_eta.max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Zlin = clusterizer.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
import cosmetics
Z = Zlin.reshape(xx.shape)
plt.figure(figsize=[8, 8])
plt.clf()
plt.imshow(
   Z, interpolation='nearest',
   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
   cmap=plt.cm.Paired,
   aspect='auto', origin='lower')
plt.title('weighting by clustering')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xlabel(cosmetics.beauty['log_trkpt'])
plt.ylabel(cosmetics.beauty['trk_eta'])
plt.plot()
try : plt.savefig('%s/%s_clusters.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clusters.pdf' % (plots, dataset))
except : pass
plt.clf()


Z = apply_weight(Zlin, weights).reshape(xx.shape)
plt.figure(figsize=[10, 8])
plt.imshow(
   Z, interpolation='nearest',
   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
   cmap=plt.cm.inferno,
   aspect='auto', origin='lower', norm=LogNorm())
plt.title('weight')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xlabel(cosmetics.beauty['log_trkpt'])
plt.ylabel(cosmetics.beauty['trk_eta'])
plt.colorbar()
plt.plot()
try : plt.savefig('%s/%s_clusters_weights.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clusters_weights.pdf' % (plots, dataset))
except : pass
plt.clf()


# plot weight distribution
entries, _, _ = plt.hist(
   data.weight, 
   bins=np.logspace(
      np.log(data.weight.min()), 
      np.log(data.weight.max()), 
      100
      ),
   histtype='stepfilled'
)
plt.xlabel('Weight')
plt.ylabel('Occurrency')
plt.legend(loc='best')
plt.ylim(0.5, entries.max()*1.2)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.plot()
try : plt.savefig('%s/%s_clustering_weights.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clustering_weights.pdf' % (plots, dataset))
except : pass
plt.clf()

# plots with unweighted events
for plot in reweight_feats:
   x_range = min(data[data.is_e][plot].min(),
                 data[np.invert(data.is_e)][plot].min()), \
                 max(data[data.is_e][plot].max(),
                     data[np.invert(data.is_e)][plot].max())
   #if plot in cosmetics.ranges: x_range = cosmetics.ranges[plot]
   plt.hist(
      data[data.is_e][plot], bins=50, normed=True,
      histtype='step', label='electrons', range=x_range,
      )
   plt.hist(
      data[np.invert(data.is_e)][plot], bins=50, normed=True,
      histtype='step', label='background', range=x_range,
      )
   plt.legend(loc='best')
   plt.xlabel(plot if plot not in cosmetics.beauty else cosmetics.beauty[plot])
   plt.ylabel('A.U.')
   try : plt.savefig('%s/%s_unweighted_%s.png' % (plots, dataset, plot))
   except : pass
   try : plt.savefig('%s/%s_unweighted_%s.pdf' % (plots, dataset, plot))
   except : pass
   plt.clf()

# plots with weighted events
for plot in reweight_feats:
   x_range = min(data[data.is_e][plot].min(), data[np.invert(data.is_e)][plot].min()), \
      max(data[data.is_e][plot].max(), data[np.invert(data.is_e)][plot].max())
   if plot in cosmetics.ranges: x_range = cosmetics.ranges[plot]
   plt.hist(
      data[data.is_e][plot], bins=50, normed=True, 
      histtype='step', label='electrons', range=x_range, weights=data[data.is_e].weight
      )
   plt.hist(
      data[np.invert(data.is_e)][plot], bins=50, normed=True, 
      histtype='step', label='background', range=x_range, weights=data[np.invert(data.is_e)].weight
      )
   plt.legend(loc='best')
   plt.xlabel(plot if plot not in cosmetics.beauty else cosmetics.beauty[plot])
   plt.ylabel('A.U.')   
   try : plt.savefig('%s/%s_reweight_%s.png' % (plots, dataset, plot))
   except : pass
   try : plt.savefig('%s/%s_reweight_%s.pdf' % (plots, dataset, plot))
   except : pass
   plt.clf()
   
