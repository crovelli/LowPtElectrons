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
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    'nn_bo_cmssw_mva_id_base_testdata.hdf', key='data')

test = test[np.invert(test.is_egamma)] 
test = test[np.invert(abs(test.trk_eta)>=2.4)] 
test = test[np.invert(test.trk_pt<0.5)] 
test = test[np.invert(test.trk_pt>15.)] 
print test.size

# default variables
mod0 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_0.pkl')
mod0d_features, _ = get_features('cmssw_mva_id_base')
test['mod0_out'] = mod0.predict_proba(test[mod0d_features].as_matrix())[:,1]
test['mod0_out'].loc[np.isnan(test.mod0_out)] = -999 
mod0_roc = roc_curve(
    test.is_e, test.mod0_out
    )
mod0_auc = roc_auc_score(test.is_e, test.mod0_out)
print mod0_auc

mod1 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_1.pkl')
mod1d_features, _ = get_features('cmssw_mva_id_base')
test['mod1_out'] = mod1.predict_proba(test[mod1d_features].as_matrix())[:,1]
test['mod1_out'].loc[np.isnan(test.mod1_out)] = -999 
mod1_roc = roc_curve(
    test.is_e, test.mod1_out
    )
mod1_auc = roc_auc_score(test.is_e, test.mod1_out)
print mod1_auc

mod2 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_2.pkl')
mod2d_features, _ = get_features('cmssw_mva_id_base')
test['mod2_out'] = mod2.predict_proba(test[mod2d_features].as_matrix())[:,1]
test['mod2_out'].loc[np.isnan(test.mod2_out)] = -999 
mod2_roc = roc_curve(
    test.is_e, test.mod2_out
    )
mod2_auc = roc_auc_score(test.is_e, test.mod2_out)
print mod2_auc

mod3 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_3.pkl')
mod3d_features, _ = get_features('cmssw_mva_id_base')
test['mod3_out'] = mod3.predict_proba(test[mod3d_features].as_matrix())[:,1]
test['mod3_out'].loc[np.isnan(test.mod3_out)] = -999 
mod3_roc = roc_curve(
    test.is_e, test.mod3_out
    )
mod3_auc = roc_auc_score(test.is_e, test.mod3_out)
print mod3_auc

mod4 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_4.pkl')
mod4d_features, _ = get_features('cmssw_mva_id_base')
test['mod4_out'] = mod4.predict_proba(test[mod4d_features].as_matrix())[:,1]
test['mod4_out'].loc[np.isnan(test.mod4_out)] = -999 
mod4_roc = roc_curve(
    test.is_e, test.mod4_out
    )
mod4_auc = roc_auc_score(test.is_e, test.mod4_out)
print mod4_auc

mod5 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_5.pkl')
mod5d_features, _ = get_features('cmssw_mva_id_base')
test['mod5_out'] = mod5.predict_proba(test[mod5d_features].as_matrix())[:,1]
test['mod5_out'].loc[np.isnan(test.mod5_out)] = -999 
mod5_roc = roc_curve(
    test.is_e, test.mod5_out
    )
mod5_auc = roc_auc_score(test.is_e, test.mod5_out)
print mod5_auc

mod6 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_6.pkl')
mod6d_features, _ = get_features('cmssw_mva_id_base')
test['mod6_out'] = mod6.predict_proba(test[mod6d_features].as_matrix())[:,1]
test['mod6_out'].loc[np.isnan(test.mod6_out)] = -999 
mod6_roc = roc_curve(
    test.is_e, test.mod6_out
    )
mod6_auc = roc_auc_score(test.is_e, test.mod6_out)
print mod6_auc

mod7 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_7.pkl')
mod7d_features, _ = get_features('cmssw_mva_id_base')
test['mod7_out'] = mod7.predict_proba(test[mod7d_features].as_matrix())[:,1]
test['mod7_out'].loc[np.isnan(test.mod7_out)] = -999 
mod7_roc = roc_curve(
    test.is_e, test.mod7_out
    )
mod7_auc = roc_auc_score(test.is_e, test.mod7_out)
print mod7_auc

mod8 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_8.pkl')
mod8d_features, _ = get_features('cmssw_mva_id_base')
test['mod8_out'] = mod8.predict_proba(test[mod8d_features].as_matrix())[:,1]
test['mod8_out'].loc[np.isnan(test.mod8_out)] = -999 
mod8_roc = roc_curve(
    test.is_e, test.mod8_out
    )
mod8_auc = roc_auc_score(test.is_e, test.mod8_out)
print mod8_auc

mod9 = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/ResultsDec17__allLowPtEle__dR0d03/bdt_bo_cmssw_mva_id_base/'
    '/model_9.pkl')
mod9d_features, _ = get_features('cmssw_mva_id_base')
test['mod9_out'] = mod9.predict_proba(test[mod9d_features].as_matrix())[:,1]
test['mod9_out'].loc[np.isnan(test.mod9_out)] = -999 
mod9_roc = roc_curve(
    test.is_e, test.mod9_out
    )
mod9_auc = roc_auc_score(test.is_e, test.mod9_out)
print mod9_auc



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

plt.plot(mod0_roc[0][:-1], mod0_roc[1][:-1], 
         linestyle='solid', 
         color='black', 
         label='Low pT, model0' %mod0_auc)

plt.plot(mod1_roc[0][:-1], mod1_roc[1][:-1], 
         linestyle='solid', 
         color='red', 
         label='Low pT, model1' %mod1_auc)

plt.plot(mod2_roc[0][:-1], mod2_roc[1][:-1], 
         linestyle='solid', 
         color='blue', 
         label='Low pT, model2' %mod2_auc)

plt.plot(mod3_roc[0][:-1], mod3_roc[1][:-1], 
         linestyle='solid', 
         color='yellow', 
         label='Low pT, model3' %mod3_auc)

plt.plot(mod4_roc[0][:-1], mod4_roc[1][:-1], 
         linestyle='solid', 
         color='green', 
         label='Low pT, model4' %mod4_auc)

plt.plot(mod5_roc[0][:-1], mod5_roc[1][:-1], 
         linestyle='solid', 
         color='cyan', 
         label='Low pT, model5' %mod5_auc)

plt.plot(mod6_roc[0][:-1], mod6_roc[1][:-1], 
         linestyle='solid', 
         color='violet', 
         label='Low pT, model6' %mod6_auc)

plt.plot(mod7_roc[0][:-1], mod7_roc[1][:-1], 
         linestyle='dashed', 
         color='black', 
         label='Low pT, model7' %mod7_auc)

plt.plot(mod8_roc[0][:-1], mod8_roc[1][:-1], 
         linestyle='dashed', 
         color='red', 
         label='Low pT, model8' %mod8_auc)

plt.plot(mod9_roc[0][:-1], mod9_roc[1][:-1], 
         linestyle='dashed', 
         color='blue', 
         label='Low pT, model9' %mod9_auc)


plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
plt.savefig('ROC_comparison.png')
plt.clf()


print "Making plots2 ..."

# 1dim distributions
plt.title('BDT output')
mod0signal = test.mod0_out.as_matrix()
mod0signal = mod0signal[test.is_e==1]
mod0bkg = test.mod0_out.as_matrix()
mod0bkg = mod0bkg[test.is_e==0]
plt.hist(mod0signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod0bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod0_comparison.png')
plt.clf()

plt.title('BDT output')
mod1signal = test.mod1_out.as_matrix()
mod1signal = mod1signal[test.is_e==1]
mod1bkg = test.mod1_out.as_matrix()
mod1bkg = mod1bkg[test.is_e==0]
plt.hist(mod1signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod1bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod1_comparison.png')
plt.clf()

plt.title('BDT output')
mod2signal = test.mod2_out.as_matrix()
mod2signal = mod2signal[test.is_e==1]
mod2bkg = test.mod2_out.as_matrix()
mod2bkg = mod2bkg[test.is_e==0]
plt.hist(mod2signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod2bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod2_comparison.png')
plt.clf()

plt.title('BDT output')
mod3signal = test.mod3_out.as_matrix()
mod3signal = mod3signal[test.is_e==1]
mod3bkg = test.mod3_out.as_matrix()
mod3bkg = mod3bkg[test.is_e==0]
plt.hist(mod3signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod3bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod3_comparison.png')
plt.clf()

plt.title('BDT output')
mod4signal = test.mod4_out.as_matrix()
mod4signal = mod4signal[test.is_e==1]
mod4bkg = test.mod4_out.as_matrix()
mod4bkg = mod4bkg[test.is_e==0]
plt.hist(mod4signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod4bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod4_comparison.png')
plt.clf()

plt.title('BDT output')
mod5signal = test.mod5_out.as_matrix()
mod5signal = mod5signal[test.is_e==1]
mod5bkg = test.mod5_out.as_matrix()
mod5bkg = mod5bkg[test.is_e==0]
plt.hist(mod5signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod5bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod5_comparison.png')
plt.clf()

plt.title('BDT output')
mod6signal = test.mod6_out.as_matrix()
mod6signal = mod6signal[test.is_e==1]
mod6bkg = test.mod6_out.as_matrix()
mod6bkg = mod6bkg[test.is_e==0]
plt.hist(mod6signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod6bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod6_comparison.png')
plt.clf()

plt.title('BDT output')
mod7signal = test.mod7_out.as_matrix()
mod7signal = mod7signal[test.is_e==1]
mod7bkg = test.mod7_out.as_matrix()
mod7bkg = mod7bkg[test.is_e==0]
plt.hist(mod7signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod7bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod7_comparison.png')
plt.clf()
 
plt.title('BDT output')
mod8signal = test.mod8_out.as_matrix()
mod8signal = mod8signal[test.is_e==1]
mod8bkg = test.mod8_out.as_matrix()
mod8bkg = mod8bkg[test.is_e==0]
plt.hist(mod8signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod8bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod8_comparison.png')
plt.clf()

plt.title('BDT output')
mod9signal = test.mod9_out.as_matrix()
mod9signal = mod9signal[test.is_e==1]
mod9bkg = test.mod9_out.as_matrix()
mod9bkg = mod9bkg[test.is_e==0]
plt.hist(mod9signal, bins=70, color="green", lw=0, label='signal',normed=1,alpha=0.5)
plt.hist(mod9bkg, bins=70, color="skyblue", lw=0, label='bkg',normed=1,alpha=0.5)
plt.show()
plt.legend(loc='best')
plt.savefig('OUT_mod9_comparison.png')
plt.clf()

