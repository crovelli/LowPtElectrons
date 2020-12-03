import numpy as np
import matplotlib 
import matplotlib.font_manager

matplotlib.use('Agg')
from argparse import ArgumentParser
import os
from matplotlib.colors import LogNorm

parser = ArgumentParser()
parser.add_argument(
   '--nbins', default=600, type=int
)
parser.add_argument(
   '--nthreads', default=10, type=int
)
parser.add_argument(
   '--test', action='store_true',
)
parser.add_argument(
   '--dataset'
)
args = parser.parse_args()

import matplotlib.pyplot as plt
import uproot
import json
import pandas as pd
from matplotlib import rc
from pdb import set_trace
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
##############rc('text', usetex=True)
from datasets import tag, apply_weight, get_data_sync, target_dataset, HistWeighter, training_selection
import os
dataset = 'test' if args.test else target_dataset
if args.dataset:
   dataset = args.dataset

mods = './models/%s/' % (tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

plots = './plots/%s/' % (tag)
if not os.path.isdir(plots):
   os.makedirs(plots)

print 'Getting dataset "{:s}"...'.format(dataset)
data = pd.DataFrame(
   get_data_sync(dataset, ['gen_pt', 'gen_eta', 'gen_dR',
                           'gsf_mode_pt', 'gsf_mode_eta', 
                           'evt', 
                           'is_e', 'is_e_not_matched', 'is_other', 'is_egamma'])
)
print '...Done'
data = data[np.invert(data.is_egamma)]          # remove EGamma electrons
data = data[np.invert(data.is_e_not_matched)]   # remove non-matched electrons
mask = training_selection(data)                 # in dataset 
data = data[mask]
data['log_gsfmodept'] = np.log10(data.gsf_mode_pt)

data['original_weight'] = 1. #np.invert(data.is_e)*original_weight.get_weight(data.log_trkpt, data.trk_eta)+data.is_e

overall_scale = data.shape[0]/float(data.is_e.sum())
reweight_feats = ['log_gsfmodept', 'gsf_mode_eta']

print 'clustering...'
from sklearn.cluster import KMeans, MiniBatchKMeans
clusterizer = MiniBatchKMeans(n_clusters=args.nbins, batch_size=3000) #n_jobs=3)
clusterizer.fit(data[reweight_feats]) #fit(data[data.is_e][reweight_feats])
global_ratio = float(data.is_e.sum())/np.invert(data.is_e).sum()

data['cluster'] = clusterizer.predict(data[reweight_feats])
counts = {}
weights = {}
for cluster, group in data.groupby('cluster'):
   nbkg = np.invert(group.is_e).sum()
   nsig = group.is_e.sum()
   if not nbkg: RuntimeError('cluster %d has no background events, reduce the number of bins!' % nbkg)
   elif not nsig: RuntimeError('cluster %d has no electrons events, reduce the number of bins!' % nsig)
   weight = float(nsig)/nbkg if nbkg > 0 else 1.
   weights[cluster] = weight
   counts[cluster] = min(nsig,nbkg)
print "Number of eles: ",data.is_e.sum()
print "Number of fakes:",np.invert(data.is_e).sum()

from sklearn.externals import joblib
joblib.dump(
   clusterizer, 
   '%s/kmeans_%s_weighter.pkl' % (mods, dataset), 
   compress=True
)

weights['features'] = reweight_feats
with open('%s/kmeans_%s_weighter.json' % (mods, dataset), 'w') as ww:
   json.dump(weights, ww)
print '...done'
del weights['features']

data['weight'] = np.invert(data.is_e)*apply_weight(data.cluster, weights)+data.is_e

print 'time for plots!'
# Step size of the mesh. Decrease to increase the quality of the VQ.
h = .01     # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = data.log_gsfmodept.min() - 0.3, data.log_gsfmodept.max() + 0.3
y_min, y_max = data.gsf_mode_eta.min() - 0.3  , data.gsf_mode_eta.max() + 0.3
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
##########plt.xlabel(cosmetics.beauty['log_gsfmodept'])
##########plt.ylabel(cosmetics.beauty['gsf_mode_eta'])
plt.plot()
try : plt.savefig('%s/%s_clusters.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clusters.pdf' % (plots, dataset))
except : pass
plt.clf()

from matplotlib.colors import LogNorm
Z = apply_weight(Zlin, weights).reshape(xx.shape)
plt.figure(figsize=[10, 8])
plt.imshow(
   Z, interpolation='nearest',
   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
   cmap=plt.cm.seismic,
   norm=LogNorm(vmin=10**-4, vmax=10**4),
   aspect='auto', origin='lower')
plt.title('weight')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
###########plt.xlabel(cosmetics.beauty['log_gsfmodept'])
###########plt.ylabel(cosmetics.beauty['gsf_mode_eta'])
plt.colorbar()
plt.plot()
try : plt.savefig('%s/%s_clusters_weights.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clusters_weights.pdf' % (plots, dataset))
except : pass
plt.clf()

from matplotlib.colors import LogNorm
Z = apply_weight(Zlin, counts).reshape(xx.shape)
plt.figure(figsize=[10, 8])
plt.imshow(
   Z, interpolation='nearest',
   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
   cmap=plt.cm.seismic,
   norm=LogNorm(vmin=0.1, vmax=max(counts.values())*10.),
   aspect='auto', origin='lower')
plt.title('counts')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
########plt.xlabel(cosmetics.beauty['log_gsfmodept'])
########plt.ylabel(cosmetics.beauty['gsf_mode_eta'])
plt.colorbar()
plt.plot()
try : plt.savefig('%s/%s_clusters_counts.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clusters_weights.pdf' % (plots, dataset))
except : pass
plt.clf()

# plot weight distribution
entries, _, _ = plt.hist(
   data.weight, 
   bins=np.logspace(
      np.log(max(data.weight.min(), 10**-5)), 
      np.log(data.weight.max()*2.), 
      100
      ),
   histtype='stepfilled'
)
plt.xlabel('Weight')
plt.ylabel('Occurrency')
plt.legend(loc='best')
plt.ylim(0.5, entries.max()*10.)
plt.xlim(max(data.weight[data.weight>0.].min()*0.5,10**-3), data.weight.max()*2.)
#plt.xlim(max(data.weights.min(), 10**-4), data.weights.max()*1.2)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.plot()
try : plt.savefig('%s/%s_clustering_weights.png' % (plots, dataset))
except : pass
try : plt.savefig('%s/%s_clustering_weights.pdf' % (plots, dataset))
except : pass
plt.clf()

for plot in reweight_feats+['gsf_mode_pt']:
   x_range = min(data[data.is_e][plot].min(), data[np.invert(data.is_e)][plot].min()), \
      max(data[data.is_e][plot].max(), data[np.invert(data.is_e)][plot].max())
   x_range = cosmetics.ranges.get(plot, x_range)
   for name, weight in [
      ('unweighted', np.ones(data.shape[0])),
      ('reweight', data.weight),
      ('original', data.original_weight)]:
      plt.hist(
         data[data.is_e][plot], bins=100, normed=True,
         histtype='step', label='electrons', range=x_range, weights=weight[data.is_e]
         )
      plt.hist(
         data[np.invert(data.is_e)][plot], bins=100, normed=True,
         histtype='step', label='background', range=x_range, weights=weight[np.invert(data.is_e)]
         )
      plt.legend(loc='best')
      plt.xlabel(plot if plot not in cosmetics.beauty else cosmetics.beauty[plot])
      plt.ylabel('A.U.')   
      try : plt.savefig('%s/%s_%s_%s.png' % (plots, dataset, name, plot))
      except : pass
      try : plt.savefig('%s/%s_%s_%s.pdf' % (plots, dataset, name, plot))
      except : pass
      plt.clf()

#compute separation with a BDT   
from sklearn.ensemble import GradientBoostingClassifier
from datasets import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score

train_bdt, test_bdt = train_test_split(data.head(1000000), 10, 5)
pre_separation = GradientBoostingClassifier(
   n_estimators=50, learning_rate=0.1,
   max_depth=4, random_state=42, verbose=1
)
pre_separation.fit(train_bdt[reweight_feats], train_bdt.is_e)
test_proba = pre_separation.predict_proba(test_bdt[reweight_feats])[:, 1]
roc_pre = roc_curve(test_bdt[['is_e']],  test_proba)[:2]
auc_pre = roc_auc_score(test_bdt[['is_e']],  test_proba)


post_separation = GradientBoostingClassifier(
   n_estimators=50, learning_rate=0.1,
   max_depth=4, random_state=42, verbose=1
)
post_separation.fit(train_bdt[reweight_feats], train_bdt.is_e, train_bdt.weight)
test_proba = post_separation.predict_proba(test_bdt[reweight_feats])[:, 1]
roc_post = roc_curve(test_bdt[['is_e']],  test_proba, sample_weight=test_bdt.weight)[:2]
auc_post = roc_auc_score(test_bdt[['is_e']],  test_proba, sample_weight=test_bdt.weight)

# make plots
plt.clf()
plt.figure(figsize=[8, 8])
plt.plot(*roc_pre, label='Initial separation (%.3f)' % auc_pre)
plt.plot(*roc_post, label='Separation after reweighting (%.3f)' % auc_post)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.plot()
plt.savefig('%s/%s_reweighting.png' % (plots, dataset))
plt.savefig('%s/%s_reweighting.pdf' % (plots, dataset))
plt.clf()

