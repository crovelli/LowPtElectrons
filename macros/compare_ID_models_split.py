import numpy as np
import matplotlib
matplotlib.use('Agg')
from cmsjson import CMSJson
from pdb import set_trace
import os
from glob import glob
import pandas as pd
import json
from pprint import pprint
import matplotlib.pyplot as plt
from features import *
from scipy.interpolate import InterpolatedUnivariateSpline
from sklearn.externals import joblib
import xgboost as xgb
from datasets import HistWeighter
from xgbo.xgboost2tmva import convert_model
from itertools import cycle
from sklearn.metrics import roc_curve , roc_auc_score

def get_model(pkl):
    model = joblib.load(pkl)

    def _monkey_patch():
        return model._Booster

    if isinstance(model.booster, basestring):
        model.booster = _monkey_patch
    return model

# test dataset - must be on the sub-range only
test = pd.read_hdf(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan20/BuToKJpsiToee/BDT_base_PT2to5/'
    '/bdt_cmssw_mva_id_base_testdata.hdf', key='data')
test = test[np.invert(test.is_egamma)] 
test = test[np.invert(abs(test.trk_eta)>=2.4)] 
test = test[np.invert(test.trk_pt<0.5)] 
test = test[np.invert(test.trk_pt>15.)] 

# dedicated training
base = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan20/BuToKJpsiToee/BDT_base_PT2to5/'
    '/2020Jan20PT2to5__cmssw_mva_id_base_BDT.pkl')
based_features, _ = get_features('cmssw_mva_id_base')
test['base_out'] = base.predict_proba(test[based_features].as_matrix())[:,1]
test['base_out'].loc[np.isnan(test.base_out)] = -999 
base_roc = roc_curve(
    test.is_e, test.base_out
    )
base_auc = roc_auc_score(test.is_e, test.base_out)
print base_auc

# general training
ecal = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan20/BuToKJpsiToee/BDT_base/'
    '/2020Jan20__cmssw_mva_id_base_BDT.pkl')
ecal_features, _ = get_features('cmssw_mva_id_base')
test['ecal_out'] = ecal.predict_proba(test[ecal_features].as_matrix())[:,1]
test['ecal_out'].loc[np.isnan(test.ecal_out)] = -999 
ecal_roc = roc_curve(
    test.is_e, test.ecal_out
    )
ecal_auc = roc_auc_score(test.is_e, test.ecal_out)
print ecal_auc

# plots
print "Making plots ..."

# ROCs
plt.figure(figsize=[8, 12])
ax = plt.subplot(111)  
box = ax.get_position()   
ax.set_position([box.x0, box.y0, box.width, box.height*0.666]) 

plt.title('Trainings comparison')
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')

plt.plot(base_roc[0][:-1], base_roc[1][:-1], 
         linestyle='solid', 
         color='black', 
         label='Dedicated training (AUC: %.3f)' %base_auc)

plt.plot(ecal_roc[0][:-1], ecal_roc[1][:-1], 
         linestyle='dashed', 
         color='green', 
         label='One classifier (AUC: %.3f)' %ecal_auc)

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
plt.savefig('DedicatedTraining.png')
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
plt.savefig('DedicatedTraining_log.png')
plt.clf()

