import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from cmsjson import CMSJson
from pdb import set_trace
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.font_manager import FontProperties

parser = ArgumentParser()
parser.add_argument(
   'what'
)
parser.add_argument(
   '--test', action='store_true'
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
parser.add_argument(
   '--config'
)
parser.add_argument(
   '--dataset'
)
parser.add_argument(
   '--selection'
)
parser.add_argument(
   '--as_weight'
)
parser.add_argument(
   '--noweight', action='store_true'
)
parser.add_argument(
   '--SW94X', action='store_true'
)
parser.add_argument(
   '--usenomatch', action='store_true'
)
parser.add_argument(
   '--load_model', action='store_true'
)
parser.add_argument(
   '--notraining', action='store_true'
)

args = parser.parse_args()

import json
if args.config:
   #config overrides eveything
   cfg = json.load(open(args.config))
   args.reg_alpha = cfg['reg_alpha'] 
   args.colsample_bytree = cfg['colsample_bytree'] 
   args.lrate = cfg['learning_rate'] 
   args.min_child_weight = cfg['min_child_weight'] 
   args.ntrees = cfg['n_estimators'] 
   args.subsample = cfg['subsample'] 
   args.reg_lambda = cfg['reg_lambda'] 
   args.depth = cfg['max_depth'] 
   args.gamma = cfg['gamma']

import matplotlib.pyplot as plt
import ROOT
import uproot
import rootpy
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import tag, pre_process_data, target_dataset, get_models_dir, train_test_split
import os

dataset = 'test' if args.test else target_dataset
if args.dataset:
   dataset = args.dataset

mods = '%s/bdt_%s' % (get_models_dir(), args.what)
if not os.path.isdir(mods):
   os.makedirs(mods)

plots = '%s/src/LowPtElectrons/LowPtElectrons/macros/plots/%s/' % (os.environ['CMSSW_BASE'], tag)
if not os.path.isdir(plots):
   os.makedirs(plots)

from features import *
features, additional = get_features(args.what)

fields = features+labeling
if args.SW94X and 'seeding' in args.what:
   fields += seed_94X_additional
else:
   fields += additional

#if 'gsf_pt' not in fields : fields += ['gsf_pt'] #@@

if not args.load_model :# dataset.endswith('.hdf'):
   data = pre_process_data(
      dataset, fields, 
      #is_egamma=False, #@@ train using low pT electrons only!
      for_seeding=('seeding' in args.what),
      keep_nonmatch=args.usenomatch
      )
   egamma = data[data.is_egamma]          # EGamma electrons
   data = data[np.invert(data.is_egamma)] # low pT electrons
   
   #data = data.head(1000).copy() #@@

#   # replicate 'nan' values in old ntuples
#   data.replace(-10.,-1.,inplace=True)
#   data.replace(-10,-1,inplace=True)
#   vars = [x for x in data.columns if x.startswith('eid_')]
#   data[vars].replace(-10.,-666.,inplace=True)

#   # replace -10. with -1. for 
#   vars = ["trk_p","trk_chi2red","gsf_chi2red","sc_E","sc_eta","sc_etaWidth",
#           "sc_phiWidth","match_seed_dEta","match_eclu_EoverP","match_SC_EoverP",
#           "match_SC_dEta","match_SC_dPhi","shape_full5x5_sigmaIetaIeta",
#           "shape_full5x5_sigmaIphiIphi","shape_full5x5_HoverE","shape_full5x5_r9",
#           "shape_full5x5_circularity","rho","brem_frac","ele_pt",]
#   data[vars].replace(-10.,-1.,inplace=True)
#   data[["trk_nhits","gsf_nhits"]].replace(-10,-1,inplace=True)

   if args.selection:
      data = data.query(args.selection)

   if args.as_weight:
      data['weight'] = data[args.as_weight]

   if args.noweight:
      data['weight'] = 1
   train_test, test = train_test_split(data, 10, 8)
   train, validation = train_test_split(train_test, 10, 6)
   train.to_hdf(
      '%s/bdt_%s_traindata.hdf' % (mods, args.what),
      'data'
      ) 
   validation.to_hdf(
      '%s/bdt_%s_valdata.hdf' % (mods, args.what),
      'data'
      ) 
   test.to_hdf(
      '%s/bdt_%s_testdata.hdf' % (mods, args.what),
      'data'
      ) 
else:   
   train = pd.read_hdf('%s/bdt_%s_traindata.hdf' % (mods, args.what), 'data')
   validation = pd.read_hdf('%s/bdt_%s_valdata.hdf' % (mods, args.what), 'data')
   test = pd.read_hdf('%s/bdt_%s_testdata.hdf' % (mods, args.what), 'data')
   if args.selection:
      train = train.query(args.selection)
      validation  = validation.query(args.selection)
      test = test.query(args.selection)

   if args.as_weight:
      train['weight'] = train[args.as_weight]
      validation['weight'] = validation[args.as_weight]
      test['weight'] = test[args.as_weight]
   if args.noweight:
      train['weight'] = 1
      validation['weight'] = 1
      test['weight'] = 1
   dataset = os.path.basename(dataset).split('.')[0]

from sklearn.externals import joblib
import xgboost as xgb
#
# Train BDTs
#

clf = None
if args.notraining :
   print 'No training done, no pre-existing model loaded!'
elif not args.load_model :

   print 'Training'
   print 'Input features:\n',features

   clf = xgb.XGBClassifier(
      #basic stuff
      max_depth=args.depth, learning_rate=args.lrate, n_estimators=args.ntrees,
      objective='binary:logitraw',
      #many different ways of regularization
      gamma=args.gamma, min_child_weight=args.min_child_weight, max_delta_step=0, 
      colsample_bytree=args.colsample_bytree, colsample_bylevel=1, subsample=args.subsample,
      reg_alpha=args.reg_alpha, reg_lambda=args.reg_lambda, 
      #running settings and weight balancing
      silent=False, nthread=args.nthreads, scale_pos_weight=1, 
   )
   
   early_stop_kwargs = {
      'eval_set' : [(validation[features].as_matrix(), validation.is_e.as_matrix().astype(int))],
      #'sample_weight_eval_set' : [validation.weight.as_matrix()], #undefined in this version
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

else :
   
   full_model = '%s/%s_%s_%s_BDT.pkl' % (mods, dataset, args.jobtag, args.what)
   clf = joblib.load(full_model)
   print 'Loaded pre-existing model!'

#
# plot performance
#
from sklearn.metrics import roc_curve, roc_auc_score
args_dict = args.__dict__

rocs = {}
if not args.notraining :
   for df, name in [
      ##(train, 'train'),
      ##(validation, 'validation'),
      (test, 'test')
      ]:
      training_out = clf.predict_proba(df[features].as_matrix())[:, 1]
      rocs[name] = roc_curve(
         df.is_e.as_matrix().astype(int), 
         training_out)[:2]
      args_dict['%s_AUC' % name] = roc_auc_score(df.is_e, training_out)

   with open('%s/%s_%s_%s_BDT.json' % (mods, dataset, args.jobtag, args.what), 'w') as info:
      json.dump(args_dict, info)

# make plots
print "Making plots ..."
plt.figure(figsize=[8, 12])
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height*0.666])
plt.title('%s training' % args.what.replace("_"," "))
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')
if not args.notraining : 

   plt.plot(rocs['test'][0][:-1], rocs['test'][1][:-1], 
            linestyle='solid', 
            color='black', 
            label='Low pT, retraining, AUC: %.3f'  % args_dict['test_AUC'])

if args.what in ['seeding', 'fullseeding']:

   eff = float((test.baseline & test.is_e).sum())/test.is_e.sum()
   mistag = float((test.baseline & np.invert(test.is_e)).sum())/np.invert(test.is_e).sum()
   rocs['baseline'] = [[mistag], [eff]]
   plt.plot([mistag], [eff], 'o', label='baseline', markersize=5)   

elif 'id' in args.what:
   #mva_v2 = roc_curve(test.is_e, test.ele_mvaIdV2)[:2]
   #mva_v2_auc = roc_auc_score(test.is_e, test.ele_mvaIdV2)
   #rocs['mva_v2'] = mva_v2
   #plt.plot(*mva_v2, label='MVA ID V2 (AUC: %.2f)'  % mva_v2_auc)

   ## ID (Low pT) ROC curve
   #mva_lowpt = roc_curve(test.is_e, test.ele_mva_value)[:2] #@@ 
   #mva_lowpt_auc = roc_auc_score(test.is_e, test.ele_mva_value) #@@
   #rocs['mva_lowpt'] = mva_lowpt
   #plt.plot(*mva_lowpt, label='MVA ID low pT (AUC: %.2f)'  % mva_lowpt_auc)

   ##########
   #test['preid_bdtoutOR'] = test[['preid_bdtout1','preid_bdtout2']].max(axis=1)

   # Low pT selections
   has_gsf = (test.has_gsf) & (test.gsf_pt>0.5) & (np.abs(test.gsf_eta)<2.4)
   has_ele = has_gsf & test.has_ele & (test.gsf_bdtout1>test.gsf_bdtout1.min())
   has_ele2 = has_ele & (test.ele_pt>2.0)

   # Low pT seed and ID ROC, AxE
   seed_lowpt = roc_curve(test.is_e, test.gsf_bdtout1)
   seed_lowpt_auc = roc_auc_score(test.is_e, test.gsf_bdtout1)
   #rocs['seed_lowpt_axe'] = seed_lowpt[:2]
   plt.plot(seed_lowpt[0][:-1], seed_lowpt[1][:-1], 
            linestyle='solid', 
            color='green', 
            label='Seed, $\mathcal{A}\epsilon$, AUC: %.2f'%seed_lowpt_auc)
   denom = test.is_e; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e); numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='green', markersize=5 )

   # Low pT seed and ID ROC, Eff(pT>0.5)
   seed_lowpt = roc_curve(test[has_ele].is_e, test[has_ele].gsf_bdtout1)
   seed_lowpt_auc = roc_auc_score(test[has_ele].is_e, test[has_ele].gsf_bdtout1)
   #rocs['seed_lowpt_pt0p5'] = seed_lowpt[:2]
   plt.plot(seed_lowpt[0][:-1], seed_lowpt[1][:-1], 
            linestyle='dashed', 
            color='green', 
            label='Seed, prior: $p^{trk}_{T} > 0.5\, GeV$, AUC: %.2f'%seed_lowpt_auc)
   denom = test.is_e & has_ele; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e) & has_ele; numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='green', markersize=5)

   # Low pT seed and ID ROC, Eff(pT>2.0)
   seed_lowpt = roc_curve(test[has_ele2].is_e, test[has_ele2].gsf_bdtout1)
   seed_lowpt_auc = roc_auc_score(test[has_ele2].is_e, test[has_ele2].gsf_bdtout1)
   #rocs['seed_lowpt_pt0p5'] = seed_lowpt[:2]
   plt.plot(seed_lowpt[0][:-1], seed_lowpt[1][:-1], 
            linestyle='dashdot', 
            color='green', 
            label='Seed, prior: $p^{trk}_{T} > 2.0\, GeV$, AUC: %.2f'%seed_lowpt_auc)
   denom = test.is_e & has_ele2; numer = denom & has_ele2
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e) & has_ele2; numer = denom & has_ele2; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='green', markersize=5)

   ##########

   # Low pT selections
   has_gsf = (test.has_gsf) & (test.gsf_pt>0.5) & (np.abs(test.gsf_eta)<2.4)
   has_ele = has_gsf & test.has_ele & (test.ele_mva_value>test.ele_mva_value.min())
   has_ele2 = has_ele & (test.ele_pt>2.0)

   # Low pT ele and ID ROC, AxE
   mva_lowpt = roc_curve(test.is_e, test.ele_mva_value)
   mva_lowpt_auc = roc_auc_score(test.is_e, test.ele_mva_value)
   rocs['mva_lowpt_axe'] = mva_lowpt[:2]
   plt.plot(mva_lowpt[0][:-1], mva_lowpt[1][:-1], 
            linestyle='solid', 
            color='blue', 
            label='Low pT, $\mathcal{A}\epsilon$, AUC: %.2f'%mva_lowpt_auc)
   # Low pT electron: AxE
   denom = test.is_e; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e); numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='blue', markersize=5 )

   # Low pT ele and ID ROC, Eff(pT>0.5)
   mva_lowpt = roc_curve(test[has_ele].is_e, test[has_ele].ele_mva_value)
   mva_lowpt_auc = roc_auc_score(test[has_ele].is_e, test[has_ele].ele_mva_value)
   rocs['mva_lowpt_pt0p5'] = mva_lowpt[:2]
   plt.plot(mva_lowpt[0][:-1], mva_lowpt[1][:-1], 
            linestyle='dashed', 
            color='blue', 
            label='Low pT, prior: $p^{trk}_{T} > 0.5\, GeV$, AUC: %.2f'%mva_lowpt_auc)
   denom = test.is_e & has_ele; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e) & has_ele; numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='blue', markersize=5)

   # Low pT ele and ID ROC, Eff(pT>2.0)
   mva_lowpt = roc_curve(test[has_ele2].is_e, test[has_ele2].ele_mva_value)
   mva_lowpt_auc = roc_auc_score(test[has_ele2].is_e, test[has_ele2].ele_mva_value)
   rocs['mva_lowpt_pt0p5'] = mva_lowpt[:2]
   plt.plot(mva_lowpt[0][:-1], mva_lowpt[1][:-1], 
            linestyle='dashdot', 
            color='blue', 
            label='Low pT, prior: $p^{trk}_{T} > 2.0\, GeV$, AUC: %.2f'%mva_lowpt_auc)
   denom = test.is_e & has_ele2; numer = denom & has_ele2
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(test.is_e) & has_ele2; numer = denom & has_ele2; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='blue', markersize=5)

   ##########

   # Egamma selections
   has_gsf = (egamma.has_gsf) & (egamma.gsf_pt>0.5) & (np.abs(egamma.gsf_eta)<2.4)
   has_ele = has_gsf & egamma.has_ele & (egamma.ele_mva_value>egamma.ele_mva_value.min())
   has_ele2 = has_ele & (egamma.ele_pt>2.0)

   # PF ele and ID ROC, AxE
   mva_egamma = roc_curve(egamma.is_e, egamma.ele_mva_value)
   mva_egamma_auc = roc_auc_score(egamma.is_e, egamma.ele_mva_value)
   rocs['mva_egamma_axe'] = mva_egamma[:2]
   plt.plot(mva_egamma[0][:-1], mva_egamma[1][:-1], 
            linestyle='solid', 
            color='red', 
            label='PF ele, $\mathcal{A}\epsilon$, AUC: %.2f'%mva_egamma_auc)
   denom = egamma.is_e; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(egamma.is_e); numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='red', markersize=5 )

   # PF ele and ID ROC, Eff(pT>0.5)
   mva_egamma = roc_curve(egamma[has_ele].is_e, egamma[has_ele].ele_mva_value)
   mva_egamma_auc = roc_auc_score(egamma[has_ele].is_e, egamma[has_ele].ele_mva_value)
   rocs['mva_egamma_pt0p5'] = mva_egamma[:2]
   plt.plot(mva_egamma[0][:-1], mva_egamma[1][:-1], 
            linestyle='dashed', 
            color='red', 
            label='PF ele, prior: $p^{trk}_{T} > 0.5\, GeV$, AUC: %.2f'%mva_egamma_auc)
   denom = egamma.is_e & has_ele; numer = denom & has_ele
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(egamma.is_e) & has_ele; numer = denom & has_ele; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='red', markersize=5)

   # PF ele and ID ROC, Eff(pT>2.0)
   mva_egamma = roc_curve(egamma[has_ele2].is_e, egamma[has_ele2].ele_mva_value)
   mva_egamma_auc = roc_auc_score(egamma[has_ele2].is_e, egamma[has_ele2].ele_mva_value)
   rocs['mva_egamma_pt0p5'] = mva_egamma[:2]
   plt.plot(mva_egamma[0][:-1], mva_egamma[1][:-1], 
            linestyle='dashdot', 
            color='red', 
            label='PF ele, prior: $p^{trk}_{T} > 2.0\, GeV$, AUC: %.2f'%mva_egamma_auc)
   denom = egamma.is_e & has_ele2; numer = denom & has_ele2
   eff = float(numer.sum()) / float(denom.sum())
   denom = np.invert(egamma.is_e) & has_ele2; numer = denom & has_ele2; 
   mistag = float(numer.sum()) / float(denom.sum())
   plt.plot([mistag], [eff], marker='o', color='red', markersize=5)

   ##########

   # Adapt legend
   def update_prop(handle, orig):
      handle.update_from(orig)
      handle.set_marker("o")
   plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})

else:
   pass #raise ValueError()

for key in rocs:
   fpr, tpr = rocs[key]
   rocs[key] = [list(fpr), list(tpr)]

with open('%s/%s_%s_%s_ROCS.json' % (plots, dataset, args.jobtag, args.what), 'w') as rr:
   rr.write(json.dumps(rocs))

plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='lower left', bbox_to_anchor=(0., 1.1))
#plt.legend(loc='best')
plt.xlim(0., 1)
try : plt.savefig('%s/%s_%s_%s_BDT.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
try : plt.savefig('%s/%s_%s_%s_log_BDT.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_log_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.clf()
