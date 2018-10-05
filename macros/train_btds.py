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
   '--ntrees', default=5000, type=int
)
parser.add_argument(
   '--depth', default=4, type=int
)
parser.add_argument(
   '--lrate', default=0.1, type=float
)
parser.add_argument(
   '--rstate', default=42, type=int
)
parser.add_argument(
   '--gamma', default=0, type=float
)
parser.add_argument(
   '--min_child_weight', default=1, type=int
)
parser.add_argument(
   '--subsample', default=1, type=float
)
parser.add_argument(
   '--colsample_bytree', default=1, type=float
)
parser.add_argument(
   '--reg_alpha', default=0, type=float
)
parser.add_argument(
   '--reg_lambda', default=1, type=float
)
parser.add_argument(
   '--nthreads', default=8, type=int
)
parser.add_argument(
   '--no_early_stop', action='store_true'
)

args = parser.parse_args()
dataset = 'all' 
#dataset = 'test'

import matplotlib.pyplot as plt
import ROOT
import uproot
import rootpy
import json
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import get_data, tag, kmeans_weighter, training_selection
import os

mods = '%s/src/LowPtElectrons/LowPtElectrons/macros/models/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(mods):
   os.mkdirs(mods)

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.mkdirs(plots)

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

data = pd.DataFrame(
   get_data(dataset, features+labeling+additional)
)
data = data[np.invert(data.is_e_not_matched)] #remove non-matched electrons
#ensure that there is at least the GSF and a track within meaningful boundaries
data = data[training_selection(data)]
data['training_out'] = -1
data['log_trkpt'] = np.log10(data.trk_pt)
#convert bools to integers
for c in features:
   if data[c].dtype == np.dtype('bool'):
      data[c] = data[c].astype(int)


#apply pt-eta reweighting
## from hep_ml.reweight import GBReweighter
## from sklearn.externals import joblib
## reweighter = joblib.load('%s/%s_reweighting.pkl' % (mods, dataset))
## weights = reweighter.predict_weights(data[['trk_pt', 'trk_eta']])
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
train_test, validation = train_test_split(data, test_size=0.2, random_state=42)
train, test = train_test_split(train_test, test_size=0.16, random_state=42)

from hep_ml.reweight import GBReweighter
from sklearn.externals import joblib
import xgboost as xgb
#
# Train BDTs
#
print 'training'
clf = xgb.XGBClassifier(
   #basic stuff
   max_depth=args.depth, learning_rate=args.lrate, n_estimators=args.ntrees,
   #many different ways of regularization
   gamma=args.gamma, min_child_weight=args.min_child_weight, max_delta_step=0, 
   subsample=args.colsample_bytree, colsample_bylevel=1, 
   reg_alpha=args.reg_alpha, reg_lambda=args.reg_lambda, 
   #running settings and weight balancing
   silent=False, nthread=args.nthreads, scale_pos_weight=1, 
)

early_stop_kwargs = {
   'eval_set' : [(test[features].as_matrix(), test.is_e.as_matrix().astype(int))],
   #'sample_weight_eval_set' : [test.weight.as_matrix()], #undefined in this version
   'eval_metric' : 'auc',
   'early_stopping_rounds' : 10
} if not args.no_early_stop else {}

clf.fit(
   train[features].as_matrix(), 
   train.is_e.as_matrix().astype(int), 
   sample_weight=train.weight.as_matrix(),
   **early_stop_kwargs
)

full_model = '%s/%s_%s_%s_BDT.pkl' % (mods, dataset, args.jobtag, args.what)
joblib.dump(clf, full_model, compress=True)

print 'Training done!'

#
# plot performance
#
from sklearn.metrics import roc_curve, roc_auc_score
args_dict = args.__dict__

rocs = {}
for df, name in [
   (train, 'train'),
   (test, 'test'),
   (validation, 'validation')
   ]:
   training_out = clf.predict_proba(df[features].as_matrix())[:, 1]
   rocs[name] = roc_curve(
      df.is_e.as_matrix().astype(int), 
      training_out)[:2]
   args_dict['%s_AUC' % name] = roc_auc_score(df.is_e, training_out)

with open('%s/%s_%s_%s_BDT.json' % (mods, dataset, args.jobtag, args.what), 'w') as info:
   json.dump(args_dict, info)

# make plots
plt.figure(figsize=[8, 8])
plt.title('%s training' % args.what)
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')
plt.plot(*rocs['validation'], label='Retraining (AUC: %.2f)'  % args_dict['validation_AUC'])
if args.what in ['seeding', 'fullseeding']:
   eff = float((data.baseline & data.is_e).sum())/data.is_e.sum()
   mistag = float((data.baseline & np.invert(data.is_e)).sum())/np.invert(data.is_e).sum()
   plt.plot([mistag], [eff], 'o', label='baseline', markersize=5)
elif args.what == 'id':
   mva_v1 = roc_curve(validation.is_e, validation.ele_mvaIdV1)[:2]   
   mva_v2 = roc_curve(validation.is_e, validation.ele_mvaIdV2)[:2]
   mva_v1_auc = roc_auc_score(validation.is_e, validation.ele_mvaIdV1)
   mva_v2_auc = roc_auc_score(validation.is_e, validation.ele_mvaIdV2)
   plt.plot(*mva_v1, label='MVA ID V1 (AUC: %.2f)'  % mva_v1_auc)
   plt.plot(*mva_v2, label='MVA ID V2 (AUC: %.2f)'  % mva_v2_auc)
else:
   raise ValueError()

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
plt.savefig('%s/%s_%s_%s_BDT.png' % (plots, dataset, args.jobtag, args.what))
plt.savefig('%s/%s_%s_%s_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
plt.savefig('%s/%s_%s_%s_log_BDT.png' % (plots, dataset, args.jobtag, args.what))
plt.savefig('%s/%s_%s_%s_log_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
plt.clf()
