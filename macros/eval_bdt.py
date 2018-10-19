import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from cmsjson import CMSJson
from pdb import set_trace
import os

parser = ArgumentParser()
parser.add_argument('what')
parser.add_argument('model')
parser.add_argument('--dataset')

args = parser.parse_args()

import pandas as pd
import json
import matplotlib.pyplot as plt
from glob import glob
if args.model.endswith('.csv'):    
    bo = pd.read_csv(args.model)
    best = bo.target.argmax()
    pars = dict(bo.loc[best])
    del pars['target']
    base = os.path.dirname(args.model)
    
    model = '%s/model_%d.pkl' % (base, best)
    dataset = glob('%s/*.hdf' % base)[0]
    plots = base
else:
    model = args.model
    dataset = args.dataset
    plots = os.dirname(model)
    if not dataset:
        raise RuntimeError('You must specify a dataset if you are not running in Bayesian Optimization mode')

#this should be outsorced
from features import *
features, additional = get_features(args.what)

from sklearn.externals import joblib
import xgboost as xgb
model = joblib.load(model)

def _monkey_patch():
    return model._Booster

if isinstance(model.booster, basestring):
    model.booster = _monkey_patch

test = pd.read_hdf(dataset, key='data')

#
# plot performance
#
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.special import expit
training_out = expit(model.predict_proba(test[features].as_matrix())[:,1])
roc = roc_curve(
   test.is_e.as_matrix().astype(int), 
   training_out)[:2]
auc_score = roc_auc_score(test.is_e, training_out)

# make plots
plt.figure(figsize=[8, 8])
plt.title('%s training' % args.what)
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')
plt.plot(*roc, label='Retraining (AUC: %.2f)'  % auc_score)
if args.what in ['seeding', 'fullseeding']:
   eff = float((test.baseline & test.is_e).sum())/test.is_e.sum()
   mistag = float((test.baseline & np.invert(test.is_e)).sum())/np.invert(test.is_e).sum()
   plt.plot([mistag], [eff], 'o', label='baseline', markersize=5)
elif 'id' in args.what:
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
plt.savefig('%s/test_NN.png' % (plots))
plt.savefig('%s/test_NN.pdf' % (plots))
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
plt.savefig('%s/test_log_NN.png' % (plots))
plt.savefig('%s/test_log_NN.pdf' % (plots))
plt.clf()

jmap = {
    'roc' : [
        list(roc[0]),
        list(roc[1])
        ],
    'auc' : auc_score,
}
        
with open('%s/roc.json' % plots, 'w') as rr:
    rr.write(json.dumps(jmap))

importances = zip(features, model.feature_importances_)
importances.sort(key=lambda x: -x[1])
with open('%s/importances.txt' % plots, 'w') as rr:
    for name, val in importances:
        rr.write('%s\t\t%.5f\n' % (name, val))
