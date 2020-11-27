from glob import glob
#A single place where to bookkeep the dataset file locations
#
tag = '2020Nov27UL17'
posix = '2020Nov27UL17'
target_dataset = '2020Nov27UL17'

import socket
path = ""
if "cern.ch" in socket.gethostname() : path = '/eos/cms/store/cmst3/group/bpark/electron_training'
elif 'cmg-gpu1080' in socket.gethostname() : path = '/eos/cms/store/cmst3/group/bpark/electron_training'
elif "hep.ph.ic.ac.uk" in socket.gethostname() : path = '/vols/cms/bainbrid/BParking/electron_training'
print socket.gethostname()

import os
from pdb import set_trace
all_sets = glob(path+'/%s/*/*root' % posix)

# chiara: questa parte non funziona con i nuovi formati di store (non sono divisi per dataset)
sets = set([os.path.basename(i).split('_')[0].split('Assoc')[0] for i in all_sets])
sets = sorted(list(sets), key=lambda x: -len(x))
input_files = {i : [] for i in sets}
input_files['all'] = all_sets
for inf in all_sets:
   for name in sets:
      if os.path.basename(inf).startswith(name):
         input_files[name].append(inf)
         break
input_files['limited'] = [j for i, j in enumerate(input_files['all']) if i % 2]
# chiara: questa parte non funziona con i nuovi formati di store (non sono divisi per dataset)

# datasets
input_files['2019Dec16'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch3/BuToKJpsiToee_all_onlyLowPt__genDR0d03.root']
input_files['2019Dec17'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch3/BuToKJpsiToee_all.root']
input_files['2020Jan20'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BuToKJpsiToee/BuToKJpsiToeeALL.root']
input_files['2020Jan23'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BuToKee/BuToKeeALL.root']
input_files['2020Jan24'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BdToKstaree/BdToKstareeALL.root']
input_files['2020Jan27'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BsToPhiJPsi_small/BsToPhiJPsi_small.root']
input_files['2020Jan28'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BsToPhiee_small/BsToPhiee_small.root']
input_files['2020Jan30'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch4/BuToKJpsiToee/BuToKJpsiToeeALL__normalizedVariables.root']
input_files['2020Feb24'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL__normalized.root']
input_files['2020Feb25'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/aod/BuToKJpsiToeeAOD__normalized.root']
input_files['2020Jun5']  = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/DoubleElectronGun__normalized.root']
input_files['2020Jun25'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL_withRegression__normalized.root']
input_files['2020Jun30'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/DoubleElectronGun_withRegression__normalized.root']
input_files['2020Jul08'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL_withRegression_largeBprescale__normalized.root']
input_files['2020Jul20'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL_withNewRegression__normalized.root']
input_files['2020Jul26'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL_withNewRegression_lastRound__normalized.root']
input_files['2020Nov27UL17'] = ['/eos/cms/store/user/crovelli/LowPtEle/Batch1_UL17/BuToKJpsiToee_UL17_AOD__prescale0d15.root'] 

#dataset_names = {
#   'BToKee' : r'B $\to$ K ee',
#   #'BToKstee' : r'B $\to$ K* ee',
#   'BToJPsieeK' : r'B $\to$ K J/$\Psi$(ee)',
#   'BToJPsieeK_0' : r'B $\to$ K J/$\Psi$(ee)',
#   #'BToKstJPsiee' : r'B $\to$ K* J/$\Psi$(ee)',
#}

import os
if not os.path.isdir('/tmp/crovelli/plots/%s' % tag):
   os.mkdir('/tmp/crovelli/plots/%s' % tag)

#import concurrent.futures
import multiprocessing
import uproot
import numpy as np

def get_models_dir():
   if 'CMSSW_BASE' not in os.environ:
      cmssw_path = dir_path = os.path.dirname(os.path.realpath(__file__)).split('src/LowPtElectrons')[0]
      os.environ['CMSSW_BASE'] = cmssw_path
   
   mods = '/tmp/crovelli/models/%s/' % (tag)
   if not os.path.isdir(mods):
      os.makedirs(mods)
   print mods
   return mods

def get_data_sync(dataset, columns, nthreads=2*multiprocessing.cpu_count(), exclude={}, path='ntuplizer/tree'):
   if dataset not in input_files:
      raise ValueError('The dataset %s does not exist, I have %s' % (dataset, ', '.join(input_files.keys())))
   print 'getting files from "%s": ' % dataset
   print ' \n'.join(input_files[dataset])
   infiles = [uproot.open(i) for i in input_files[dataset]]
   print 'available branches:\n',infiles[0][path].keys()
   if columns == 'all':
      columns = [i for i in infiles[0][path].keys() if i not in exclude]
   try:
      ret = infiles[0][path].arrays(columns)
   except KeyError as ex:
      print 'Exception! ', ex
      set_trace()
      raise RuntimeError('Failed to open %s properly' % infiles[0])
   for infile in infiles[1:]:
      try:
         arrays = infile[path].arrays(columns)
      except:
         raise RuntimeError('Failed to open %s properly' % infile)         
      for column in columns:
         ret[column] = np.concatenate((ret[column],arrays[column]))
   return ret

from sklearn.cluster import KMeans
from sklearn.externals import joblib
import json
from pdb import set_trace

apply_weight = np.vectorize(lambda x, y: y.get(x), excluded={2})

def kmeans_weighter(features, fname):
   kmeans = joblib.load(fname)
   cluster = kmeans.predict(features)
   str_weights = json.load(open(fname.replace('.pkl', '.json')))   
   weights = {}
   for i in str_weights:
      try:
         weights[int(i)] = str_weights[i]
      except:
         pass
   return apply_weight(cluster, weights)

def training_selection(df,low=0.5,high=15.):
   #'ensures there is a GSF Track and a KTF track within eta/pt boundaries'
   return (df.gsf_mode_pt > low) & (np.abs(df.gsf_mode_eta) < 2.4) & ( (df.gen_dR<=0.03) | (df.gen_dR>=0.1) ) 

import rootpy.plotting as rplt
import root_numpy

class HistWeighter(object):
   def __init__(self, fname):
      values = [[float(i) for i in j.split()] for j in open(fname)]
      vals = np.array(values)
      xs = sorted(list(set(vals[:,0])))
      ys = sorted(list(set(vals[:,1])))
      vals[:,0] += 0.0001
      vals[:,1] += 0.0001
      mask = (vals[:,2] == 0)
      vals[:,2][mask] = 1 #turn zeros into ones
      vals[:,2] = 1/vals[:,2]
      self._hist = rplt.Hist2D(xs, ys)
      root_numpy.fill_hist(self._hist, vals[:,[0,1]], vals[:, 2])

   def _get_weight(self, x, y):
      ix = self._hist.xaxis.FindFixBin(x)
      iy = self._hist.yaxis.FindFixBin(y)
      return self._hist.GetBinContent(ix, iy)

   def get_weight(self, x, y):
      cnt = lambda x, y: self._get_weight(x, y)
      cnt = np.vectorize(cnt)
      return cnt(x, y)

import pandas as pd
import numpy as np
def pre_process_data(dataset, features, for_seeding=False, keep_nonmatch=False):  
   mods = get_models_dir()
   features = list(set(features+['gen_pt', 'gen_eta', 'gen_dR',
                                 'gsf_mode_pt', 'gsf_mode_eta',
                                 'gsf_pt', 'gsf_eta', 
                                 'ele_pt', 'ele_eta', 
                                 'sc_raw_energy','sc_energy',
                                 'evt', 'weight']))

   data_dict = get_data_sync(dataset, features) # path='features/tree')
   if 'is_e_not_matched' not in data_dict:
      data_dict['is_e_not_matched'] = np.zeros(data_dict['gsf_mode_pt'].shape, dtype=bool)
   multi_dim = {}
   for feat in ['gsf_ecal_cluster_ematrix', 'ktf_ecal_cluster_ematrix']:
      if feat in features:
         multi_dim[feat] = data_dict.pop(feat, None)
   data = pd.DataFrame(data_dict)

   ##FIXME
   ##if 'gsf_ecal_cluster_ematrix' in features:
   ##   flattened = pd.DataFrame(multi_dim['gsf_ecal_cluster_ematrix'].reshape(multi_dim.shape[0], -1))
   ##   new_features = ['crystal_%d' % i for i in range(len(flattened.columns))]
   ##   flattened.columns = new_features
   ##   features += new_features
   ##   data = pd.concat([data, flattened], axis=1)

   # hack to rename weight column to prescale column
   data['prescale'] = data.weight
   data['weight'] = np.ones(data.weight.shape)

   #remove non-matched electrons
   if not keep_nonmatch:
      multi_dim = {i : j[np.invert(data.is_e_not_matched)] for i, j in multi_dim.iteritems()}
      data = data[np.invert(data.is_e_not_matched)] 
   else:
      #make the right fraction as is_e_not_matched are fully kept and normal tracks have a 160 prescale
      notmatched = data[data.is_e_not_matched]
      data = data[np.invert(data.is_e_not_matched)]
      mask = np.random.uniform(size=notmatched.shape[0]) < 1./160
      notmatched = notmatched[mask]
      data = pd.concat((data, notmatched))
   # training pre-selection
   mask = training_selection(data)
   multi_dim = {i : j[mask] for i, j in multi_dim.iteritems()}   
   data = data[mask]
   if 'trk_dxy' in data_dict and 'trk_dxy_err' in data_dict:
      sip = data.trk_dxy/data.trk_dxy_err
      sip[np.isinf(sip)] = 0
      data['trk_dxy_sig'] = sip
      inv_sip = data.trk_dxy_err/data.trk_dxy
      inv_sip[np.isinf(inv_sip)] = 0
      data['trk_dxy_sig_inverted'] = inv_sip
   data['training_out'] = -1
   log_gsfmodept = np.log10(data.gsf_mode_pt)
   log_gsfmodept[np.isnan(log_gsfmodept)] = -9999
   data['log_gsfmodept'] = log_gsfmodept
   
   #apply pt-eta reweighting
   ## from hep_ml.reweight import GBReweighter
   ## from sklearn.externals import joblib
   ## reweighter = joblib.load('%s/%s_reweighting.pkl' % (mods, dataset))
   ## weights = reweighter.predict_weights(data[['trk_pt', 'trk_eta']])
   kmeans_model = '%s/kmeans_%s_weighter.pkl' % (mods, dataset)
   if not os.path.isfile(kmeans_model):
      print 'I could not find the appropriate model, using the general instead'
      kmeans_model = '%s/kmeans_%s_weighter.pkl' % (mods, tag)
   weights = kmeans_weighter(
      data[['log_gsfmodept', 'gsf_mode_eta']],
      kmeans_model
      )    
   data['weight'] = weights*np.invert(data.is_e) + data.is_e

   ## original_weight = HistWeighter('../data/fakesWeights.txt')
   ## data['original_weight'] = np.invert(data.is_e)*original_weight.get_weight(data.log_trkpt, data.trk_eta)+data.is_e

   #
   # pre-process data
   #   
   if 'trk_charge' in data.columns:
      for feat in ['ktf_ecal_cluster_dphi', 'ktf_hcal_cluster_dphi', 'preid_trk_ecal_Dphi']:
         if feat in data.columns:
            data[feat] = data[feat]*data['trk_charge']

#      charge = data.trk_charge
#      for feat in ['gsf_ecal_cluster_ematrix', 'ktf_ecal_cluster_ematrix']:
#         if feat in multi_dim:
#            multi_dim[feat][charge == -1] = np.flip(multi_dim[feat][charge == -1], axis=2)

   #add baseline seeding (for seeding only)
   if for_seeding:
      if 'preid_trk_ecal_match' in data.columns:
         data['baseline'] = (
            data.preid_trk_ecal_match | 
            (np.invert(data.preid_trk_ecal_match) & data.preid_trkfilter_pass & data.preid_mva_pass)
            )
      elif 'trk_pass_default_preid' in data.columns:
         data['baseline'] = data.trk_pass_default_preid
      else:
         data['baseline'] = False

   from features import labeling
   #convert bools to integers
   for c in data.columns:
      if data[c].dtype == np.dtype('bool') and c not in labeling:
         data[c] = data[c].astype(int)

   if not multi_dim:
      return data
   else:
      return data, multi_dim


def train_test_split(data, div, thr):
   mask = data.evt % div
   mask = mask < thr
   return data[mask], data[np.invert(mask)]

def reduce_mem_usage(df):
    """ 
    iterate through all the columns of a dataframe and 
    modify the data type to reduce memory usage.        
    """
    start_mem = df.memory_usage().sum() / 1024**2
    print(('Memory usage of dataframe is {:.2f}' 
                     'MB').format(start_mem))

    print 'before'
    print df
    
    for col in df.columns:
        col_type = df[col].dtype
        
        print col

        if col_type != object:
            c_min = df[col].min()
            c_max = df[col].max()
            if str(col_type)[:3] == 'int':
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)  
            else:
                if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                    df[col] = df[col].astype(np.float16)
                elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
        else:
            df[col] = df[col].astype('category')    

    end_mem = df.memory_usage().sum() / 1024**2
    print(('Memory usage after optimization is: {:.2f}' 
                              'MB').format(end_mem))
    print('Decreased by {:.1f}%'.format(100 * (start_mem - end_mem) 
                                             / start_mem))
    
    print 'after'
    print df

    return df
