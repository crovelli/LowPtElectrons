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

# test dataset
test = pd.read_hdf(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan30/BuToKJpsiToee_renorm/bdt_cmssw_mva_id_nnclean2___paramMauro/'
#    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan30/BuToKJpsiToee_renorm/BDT_base/'
    '/bdt_cmssw_mva_id_nnclean2_testdata.hdf', key='data')
test = test[np.invert(test.is_egamma)] 
test = test[np.invert(abs(test.trk_eta)>=2.4)] 
test = test[np.invert(test.trk_pt<0.5)] 
test = test[np.invert(test.trk_pt>15.)] 
print test.size

# default variables
base = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan30/BuToKJpsiToee_renorm/bdt_cmssw_mva_id_nnclean2___paramMauro/'
    '/2020Jan30__cmssw_mva_id_nnclean2_BDT.pkl')
    #'/eos/cms/store/user/crovelli/LowPtEle/ResultsJan30/BuToKJpsiToee_renorm/BDT_base/'
    #'/2020Jan30__cmssw_mva_id_base_BDT.pkl')
based_features, _ = get_features('cmssw_mva_id_nnclean2')
test['base_out'] = base.predict_proba(test[based_features].as_matrix())[:,1]
test['base_out'].loc[np.isnan(test.base_out)] = -999 
base_roc = roc_curve(
    test.is_e, test.base_out
    )
base_auc = roc_auc_score(test.is_e, test.base_out)
print "ROC done"
print base_auc

# updated variables
#ecal = get_model(
#    '/eos/cms/store/user/crovelli/LowPtEle/ResultsJan30/BuToKJpsiToee_renorm/bdt_cmssw_mva_id_nnclean2/'
#    '/2020Jan30__cmssw_mva_id_nnclean2_BDT.pkl')
#ecal_features, _ = get_features('cmssw_mva_id_nnclean2')
#test['ecal_out'] = ecal.predict_proba(test[ecal_features].as_matrix())[:,1]
#test['ecal_out'].loc[np.isnan(test.ecal_out)] = -999 
#ecal_roc = roc_curve(
#    test.is_e, test.ecal_out
#)
#ecal_auc = roc_auc_score(test.is_e, test.ecal_out)
#print ecal_auc

# training version in cmssw
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
         linestyle='solid', 
         color='black', 
#         label='Aug22 variables (AUC: %.3f)' %base_auc)
         label='Extended set (AUC: %.3f)' %base_auc)

#plt.plot(ecal_roc[0][:-1], ecal_roc[1][:-1], 
#         linestyle='dashed', 
#         color='red', 
#         label='Reduced set (AUC: %.3f)' %ecal_auc)

plt.plot(cmssw_roc[0][:-1], cmssw_roc[1][:-1],
         linestyle='dashed', 
         color='red',
         label='MVA ID V2 (AUC: %.3f)' %cmssw_auc)

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
plt.savefig('ROC_comparison.png')
plt.gca().set_xscale('log')
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

plt.title('gen dR')
basesignal2 = test.gen_dR.as_matrix()
basesignal2 = basesignal2[test.is_e==1]
basebkg2 = test.gen_dR.as_matrix()
basebkg2 = basebkg2[test.is_e==0]
plt.hist(basebkg2,    bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.hist(basesignal2, bins=70, color="green",   lw=0, label='signal',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.gca().set_yscale('log')
plt.xlim(0, 0.3)
plt.savefig('GenDr_comparison.png')
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
    (1.83 , 1.83, 'T'),
    (0.60 , 0.60, 'M'),
    (-0.56, -0.56, 'L'),
    ]:
   print 'WP', wpname
   test['base_pass'] = test.base_out > base_thr
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

#   idx = np.abs(biased_roc[0] - mistag).argmin()
#   print 'similar mistag: %.2f\t%.4f\t%.2f' % (biased_roc[1][idx], biased_roc[0][idx], biased_roc[2][idx])
    

# Eff vs eta and pT for some reference WP
