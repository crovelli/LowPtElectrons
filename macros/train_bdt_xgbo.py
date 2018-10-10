import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from cmsjson import CMSJson
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument(
   'what', choices=['seeding', 'fullseeding', 'id'], 
)
parser.add_argument(
   '--jobtag', default='', type=str
)

parser.add_argument(
   '--test', action='store_true'
)

args = parser.parse_args()
dataset = 'test' if args.test else 'all' 

import matplotlib.pyplot as plt
import uproot
import json
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import get_data, tag, kmeans_weighter, training_selection
import os

if 'CMSSW_BASE' not in os.environ:
   cmssw_path = dir_path = os.path.dirname(os.path.realpath(__file__)).split('src/LowPtElectrons')[0]
   os.environ['CMSSW_BASE'] = cmssw_path

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.makedirs(mods)

opti_dir = '%s/bdt_bo_%s' % (mods, args.what)
if not os.path.isdir(opti_dir):
   os.makedirs(opti_dir)

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.makedirs(plots)

from features import *

if args.what == 'seeding':
   features = seed_features
   additional = seed_additional
elif args.what == 'fullseeding':
   features = fullseed_features
   additional = seed_additional
elif args.what == 'id':
   features = id_features
   additional = id_additional
else:
   raise ValueError()

fields = features+labeling+additional
if 'gsf_pt' not in fields : fields += ['gsf_pt']
data = pd.DataFrame(
   get_data(dataset, fields)
)
data = data[np.invert(data.is_e_not_matched)] #remove non-matched electrons
data = data[training_selection(data)]
data['training_out'] = -1
data['log_trkpt'] = np.log10(data.trk_pt)
#convert bools to integers
for c in features:
   if data[c].dtype == np.dtype('bool'):
      data[c] = data[c].astype(int)

#apply pt-eta reweighting
weights = kmeans_weighter(
   data[['log_trkpt', 'trk_eta']],
   '%s/kmeans_%s_weighter.plk' % (mods, dataset)
   ) 
data['weight'] = weights*data.is_e + np.invert(data.is_e)

#add baseline seeding (for seeding only)
if args.what in ['seeding', 'fullseeding']:
   data['baseline'] = (
      data.preid_trk_ecal_match | 
      (np.invert(data.preid_trk_ecal_match) & data.preid_trkfilter_pass & data.preid_mva_pass)
      )

from sklearn.model_selection import train_test_split
train, test = train_test_split(data, test_size=0.2, random_state=42)

import xgboost as xgb

xgtrain = xgb.DMatrix(
   train[features],
   label = train.is_e,
   weight = train.weight,
)
xgtest  = xgb.DMatrix(
   test[features],
   label = test.is_e,
   weight = test.weight,
)

params = {'eval_metric':'auc',
          'objective'  :'binary:logitraw'}
#model_default = xgb.train(params, xgtrain, num_boost_round=1000)

#Next we try out the xgbo package.
from xgbo import XgboClassifier

xgbo_classifier = XgboClassifier(
   out_dir=opti_dir,
)

xgbo_classifier.optimize(xgtrain, init_points=5, n_iter=50, acq='ei')
#xgbo_classifier.optimize(xgtrain, init_points=0, n_iter=1, acq='ei')

xgbo_classifier.fit(xgtrain, model="default")
xgbo_classifier.fit(xgtrain, model="optimized")

xgbo_classifier.save_model(features, model="default")
xgbo_classifier.save_model(features, model="optimized")

# Now, let's make predictions with all models we got on the testing sample.
preds_default    = model_default.predict(xgtest)
preds_early_stop = xgbo_classifier.predict(xgtest, model="default")
preds_optimized  = xgbo_classifier.predict(xgtest, model="optimized")

xgbo_classifier._bo.summary()

"""
Finally, we want to plot some ROC curves.
"""
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

fpr_default   , tpr_default, _    = roc_curve(test.is_e, preds_default, pos_label=1)
fpr_early_stop, tpr_early_stop, _ = roc_curve(test.is_e, preds_early_stop, pos_label=1)
fpr_optimized , tpr_optimized, _  = roc_curve(test.is_e, preds_optimized, pos_label=1)

auc_default    = roc_auc_score(test.is_e, preds_default   )
auc_early_stop = roc_auc_score(test.is_e, preds_early_stop)
auc_optimized  = roc_auc_score(test.is_e, preds_optimized )

# make plots
plt.figure(figsize=[8, 8])
plt.title('%s training' % args.what)
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')
plt.plot(fpr_default, tpr_default, label='default (AUC: %.3f)' % auc_default)
plt.plot(fpr_early_stop, tpr_early_stop, label='early stop (AUC: %.3f)' % auc_early_stop)
plt.plot(fpr_optimized, tpr_optimized, label='optimized (AUC: %.3f)' % auc_default)
if args.what in ['seeding', 'fullseeding', 'check']:
   eff = float((data.baseline & data.is_e).sum())/data.is_e.sum()
   mistag = float((data.baseline & np.invert(data.is_e)).sum())/np.invert(data.is_e).sum()
   plt.plot([mistag], [eff], 'o', label='baseline', markersize=5)
elif args.what == 'id':
   mva_v1 = roc_curve(test.is_e, test.ele_mvaIdV1)[:2]   
   mva_v2 = roc_curve(test.is_e, test.ele_mvaIdV2)[:2]
   mva_v1_auc = roc_auc_score(test.is_e, test.ele_mvaIdV1)
   mva_v2_auc = roc_auc_score(test.is_e, test.ele_mvaIdV2)
   plt.plot(*mva_v1, label='MVA ID V1 (AUC: %.2f)'  % mva_v1_auc)
   plt.plot(*mva_v2, label='MVA ID V2 (AUC: %.2f)'  % mva_v2_auc)
else:
   raise ValueError()

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
try : plt.savefig('%s/%s_%s_%s_BDT_bayes.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_BDT_bayes.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
try : plt.savefig('%s/%s_%s_%s_log_BDT_bayes.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_log_BDT_bayes.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.clf()

