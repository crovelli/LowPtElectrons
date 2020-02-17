import numpy as np
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import uproot
import root_numpy
import rootpy
import rootpy.plotting as rplt
import json
import pandas as pd
from matplotlib import rc
from pdb import set_trace
from features import *
from sklearn.externals import joblib

matplotlib.use('Agg')

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument(
   '--fakes', action='store_true', help='use all tracks'
)
args = parser.parse_args()

def get_model(pkl):
    model = joblib.load(pkl)

    def _monkey_patch():
        return model._Booster

    if isinstance(model.booster, basestring):
        model.booster = _monkey_patch
    return model

# plotting
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plot_type = 'efficiency' if not args.fakes else 'fakerate'

# test dataset
test = pd.read_hdf(
    '/eos/cms/store/user/crovelli/LowPtEle/Results/models/2019Dec4/bdt_cmssw_mva_id_ecalvai/'
    '/bdt_cmssw_mva_id_ecalvai_testdata.hdf', key='data')
test = test[np.invert(test.is_egamma)] 
test = test[np.invert(abs(test.trk_eta)>=2.4)] 
test = test[np.invert(test.trk_pt<0.5)] 
test = test[np.invert(test.trk_pt>15.)] 
print test.size

# dedicated - chiara
#test = test[np.invert(abs(test.ele_eta)>1.5)] 

# model
base = get_model(
    '/eos/cms/store/user/crovelli/LowPtEle/Results/models/2019Dec4/bdt_cmssw_mva_id_ecalvai/'
    '/2019Dec4__cmssw_mva_id_ecalvai_BDT.pkl')
based_features, _ = get_features('cmssw_mva_id_ecalvai')
test['base_out'] = base.predict_proba(test[based_features].as_matrix())[:,1]
test['base_out'].loc[np.isnan(test.base_out)] = -999 
print "model taken"

# utilities
def plot_efficiency(eff, **kw):
   graph = eff.graph   
   effs = [i for i in graph.y()]
   errs = np.array([i for i in graph.yerr()]).transpose()
   xs = [i for i in graph.x()]
   xerr = np.array([i for i in graph.xerr()]).transpose()
   plt.errorbar(xs, effs, yerr=errs, xerr=xerr, **kw)


# Efficiencies
electrons = test[test.is_e == 1]
if args.fakes:
   electrons = test[test.is_e == 0]   

ordered_masks = [
    ('all', None), 
    ('WP~75', electrons.base_out.as_matrix()>0.60),   # chiara
    ('WP~90', electrons.base_out.as_matrix()>-0.56),  # chiara
]
masks = dict(ordered_masks)

histosPt = {}
for name, mask in masks.iteritems():
    hist = rplt.Hist([0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,13.,15.])
    print name
    masked = electrons[mask] if mask is not None else electrons
    root_numpy.fill_hist(hist, masked.ele_pt)
    histosPt[name] = hist
    print 'masks pT done'

histosEta = {}
for name, mask in masks.iteritems():
    # chiara
    hist = rplt.Hist([0.,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.6,1.8,2.,2.2,2.4])
    #hist = rplt.Hist([0.,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5])   # EB
    #hist = rplt.Hist([1.5,1.6,1.8,2.,2.2,2.4])               # EE
    print name
    masked = electrons[mask] if mask is not None else electrons
    root_numpy.fill_hist(hist, abs(masked.ele_eta))
    histosEta[name] = hist
    print 'masks eta done'

print ''    
markersize = 6

efficienciesPt = {}
efficienciesEta = {}

plt.clf()
for passing, _ in ordered_masks:
    if passing == 'all': continue
    print passing
    efficienciesPt[passing] = rplt.Efficiency(histosPt[passing], histosPt['all'])
    plot_efficiency(efficienciesPt[passing], 
                    fmt="o", markersize=markersize, 
                    label=passing, markeredgewidth=0.0)
    plt.legend(loc='lower right')
    plt.xlabel('Ele $p_{T}$')
    plt.ylabel('Efficiency')
    plt.ylim(0.5, 1.)
    plt.grid()
    plt.savefig('EfficiencyVsPt_%s.png' % (passing))
    plt.clf()

plt.clf()
for passing, _ in ordered_masks:
    if passing == 'all': continue
    print passing
    efficienciesEta[passing] = rplt.Efficiency(histosEta[passing], histosEta['all'])
    plot_efficiency(efficienciesEta[passing], 
                    fmt="o", markersize=markersize, 
                    label=passing, markeredgewidth=0.0)
    plt.legend(loc='lower right')
    plt.xlabel('Ele eta')
    plt.ylabel('Efficiency')
    plt.ylim(0.5, 1.)
    plt.grid()
    plt.savefig('EfficiencyVsEta_%s.png' % (passing))
    plt.clf()

