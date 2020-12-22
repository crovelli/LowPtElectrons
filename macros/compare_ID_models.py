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
from itertools import cycle
from sklearn.metrics import roc_curve , roc_auc_score

def get_model(pkl):
    model = joblib.load(pkl)

    def _monkey_patch():
        return model._Booster

    if isinstance(model.booster, basestring):
        model.booster = _monkey_patch
    return model

# test dataset
test = pd.read_hdf(
    'models_all/2020Nov28ULALL_Depth10/bdt_cmssw_mva_id_nnclean2_forUL/'
    '/bdt_cmssw_mva_id_nnclean2_forUL_testdata.hdf', key='data')
test = test[np.invert(test.is_egamma)] 
test = test[np.invert(abs(test.gsf_mode_eta)>=2.4)] 
test = test[np.invert(test.gsf_mode_pt<0.5)] 

print test.size

# default variables
base = get_model(
    'models_all/2020Nov28ULALL_Depth10/bdt_cmssw_mva_id_nnclean2_forUL/'
    '/2020Nov28ULALL__cmssw_mva_id_nnclean2_forUL_BDT.pkl')
based_features, _ = get_features('cmssw_mva_id_nnclean2_forUL')
test['base_out'] = base.predict_proba(test[based_features].as_matrix())[:,1]
test['base_out'].loc[np.isnan(test.base_out)] = -999 
base_roc = roc_curve(
    test.is_e, test.base_out
    )
base_auc = roc_auc_score(test.is_e, test.base_out)
print "ROC done"
print base_auc

# updated variables
ecalA = get_model(
    'models_all/2020Nov28ULALL_Depth12/bdt_cmssw_mva_id_nnclean2_forUL/'
    '/2020Nov28ULALL__cmssw_mva_id_nnclean2_forUL_BDT.pkl')
ecalA_features, _ = get_features('cmssw_mva_id_nnclean2_forUL')
test['ecalA_out'] = ecalA.predict_proba(test[ecalA_features].as_matrix())[:,1]
test['ecalA_out'].loc[np.isnan(test.ecalA_out)] = -999 
ecalA_roc = roc_curve(
    test.is_e, test.ecalA_out
)
ecalA_auc = roc_auc_score(test.is_e, test.ecalA_out)
print ecalA_auc

ecalB = get_model(
    'models_all/2020Nov28ULALL_Depth13/bdt_cmssw_mva_id_nnclean2_forUL/'
    '/2020Nov28ULALL__cmssw_mva_id_nnclean2_forUL_BDT.pkl')
ecalB_features, _ = get_features('cmssw_mva_id_nnclean2_forUL')
test['ecalB_out'] = ecalB.predict_proba(test[ecalB_features].as_matrix())[:,1]
test['ecalB_out'].loc[np.isnan(test.ecalB_out)] = -999 
ecalB_roc = roc_curve(
    test.is_e, test.ecalB_out
)
ecalB_auc = roc_auc_score(test.is_e, test.ecalB_out)
print ecalB_auc

## training version in cmssw
cmssw_roc = roc_curve(
     test.is_e, test.ele_mva_value
    )
cmssw_auc = roc_auc_score(test.is_e, test.ele_mva_value)
print cmssw_auc

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
         linestyle='dashed', 
         color='red', 
         label='Depth10 (AUC: %.3f)' %base_auc)

plt.plot(ecalA_roc[0][:-1], ecalA_roc[1][:-1], 
         linestyle='dashed', 
         color='green', 
         label='Depth12 (AUC: %.3f)' %ecalA_auc)

plt.plot(ecalB_roc[0][:-1], ecalB_roc[1][:-1], 
         linestyle='dashed', 
         color='blue', 
         label='Depth13 (AUC: %.3f)' %ecalB_auc)

plt.plot(cmssw_roc[0][:-1], cmssw_roc[1][:-1],
         linestyle='solid', 
         color='black',
         label='B-Parking Id (AUC: %.3f)' %cmssw_auc)

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
plt.ylim(0., 1)
plt.savefig('ROC_comparison.png')
plt.gca().set_xscale('log')
plt.ylim(0., 1)
plt.xlim(1e-4, 1)
plt.savefig('ROC_comparison_log.png')
plt.clf()


print "Making plots2 ..."

# 1dim distribution
plt.title('BDT output')
basesignal = test.base_out.as_matrix()
basesignal = basesignal[test.is_e==1]
basebkg = test.base_out.as_matrix()
basebkg = basebkg[test.is_e==0]
plt.hist(basesignal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(basebkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUTBase_comparison.png')
plt.clf()

# some working points
print ''
jmap = {}
for base_thr, ecal_thr, wpname in [
        #    (1.83 , 2.61, 'T'),
#    (0.76 , 1.75, 'M'),
#    (-0.48, 1.03, 'L'),
#    (1.45 , 2.61, 'T+'),
#    (0.33 , 1.75, 'M+'),
#    (-0.97, 1.03, 'L+'),
#    ]:
#    (1.83 , 1.83, 'T'),
#    (0.60 , 0.60, 'M'),
#    (-0.56, -0.56, 'L'),
    (5. , 5., 'T1'),
    (4.8 , 4.8, 'T2'),
    (4.6 , 4.6, 'T3'),
    (4.4 , 4.4, 'T4'),
    (4.2 , 4.2, 'T5'),
    (4.0 , 4.0, 'T6'),
    (3.8 , 3.8, 'T7'),
    (3.6 , 3.6, 'T8'),
    (3.4 , 3.4, 'T9'),
    (3.2 , 3.2, 'T10'),
    (3.0 , 3.0, 'T11'),
    (2.8 , 2.8, 'T12'),
    (2.6 , 2.6, 'T13'),
    (2.4 , 2.4, 'T14'),
    (2.2 , 2.2, 'T15'),
    (2.0 , 2.0, 'T16'),
    (1.9 , 1.9, 'T17'),
    (1.8 , 1.8, 'T18'),
    (1.7 , 1.7, 'T19'),
    (1.6 , 1.6, 'T20'),
    (1.5 , 1.5, 'T21'),
    (1.4 , 1.4, 'T22'),
    (1.3 , 1.3, 'T23'),
    (1.2 , 1.2, 'T24'),
    (1.1 , 1.1, 'T25'),
    (1.0 , 1.0, 'T26'),
    (0.9 , 0.9, 'T27'),
    (0.8 , 0.8, 'T28'),
    (0.7 , 0.7, 'T29'),
    (0.6 , 0.6, 'T30'),
    (0.5 , 0.5, 'T31'),
    (0.4 , 0.4, 'T32'),
    (0.3 , 0.3, 'T33'),
    (0.2 , 0.2, 'T34'),
    (0.1 , 0.1, 'T35'),
    (0. , 0., 'T36'),
    (-0.1 , -0.1, 'T37'),
    (-0.2 , -0.2, 'T38'),
    (-0.3 , -0.3, 'T39'),
    (-0.4 , -0.4, 'T40'),
    (-0.5 , -0.5, 'T41'),
    (-0.6 , -0.6, 'T42'),
    (-0.7 , -0.7, 'T43'),
    (-0.8 , -0.8, 'T44'),
    (-0.9 , -0.9, 'T45'),
    (-1. , -1., 'T46'),
    ]:
   print 'WP', wpname
   print 'base:'
   test['base_pass'] = test.base_out > base_thr
#   print 'ecal:'
#   test['ecal_pass'] = test.ecal_out > ecal_thr
    
   eff_base = ((test.base_pass & test.is_e).sum()/float(test.is_e.sum()))
   mistag_base = ((test.base_pass & np.invert(test.is_e)).sum()/float(np.invert(test.is_e).sum()))
#   eff_ecal = ((test.ecal_pass & test.is_e).sum()/float(test.is_e.sum()))
#   mistag_ecal = ((test.ecal_pass & np.invert(test.is_e)).sum()/float(np.invert(test.is_e).sum()))

   jmap[wpname] = [mistag_base, eff_base]
   print 'eff (base): %.3f' % eff_base
   print 'mistag (base): %.3f' % mistag_base
#   print 'eff (ecal): %.3f' % eff_ecal
#   print 'mistag (ecal): %.3f' % mistag_ecal

