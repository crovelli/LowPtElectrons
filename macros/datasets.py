from glob import glob
#A single place where to bookkeep the dataset file locations
tag = '2018Sep20'
input_files = {
   'test' : ['/eos/cms/store/cmst3/user/mverzett/BToKee_Pythia/crab_2018Sep20_BToKee_v1AssocByDR/180920_162924/0000/BToKee_assocByDR_10.root']
}

all_sets = []
for dataset, name in [
   ('BToKee_Pythia', 'BToKee'),
   ('BToKstee_Pythia', 'BToKstee'),
   ('Bu_KJPsi_ee_Pythia', 'BToKJPsiee'),
   ('Bd_KstJPsi_ee_Pythia_GEN-SIM_18_07_01', 'BToKstJPsiee'),
   ('Bu_KJPsi_ee_Pythia_GEN-SIM_18_06_4', 'BToKJPsiee')]:
   if name not in input_files: input_files[name] = []
   files = glob('/eos/cms/store/cmst3/user/mverzett/%s/crab_%s_*/*/*/*.root' % (dataset, tag))
   input_files[name] += files
   all_sets += files

input_files['all'] = all_sets

dataset_names = {
   'BToKee' : r'B $\to$ K ee',
   'BToKstee' : r'B $\to$ K* ee',
   'BToKJPsiee' : r'B $\to$ K J/$\Psi$(ee)',
   'BToKstJPsiee' : r'B $\to$ K* J/$\Psi$(ee)',
}

import os
if not os.path.isdir('plots/%s' % tag):
   os.mkdir('plots/%s' % tag)

import concurrent.futures
import multiprocessing
import uproot
import numpy as np

def get_data(dataset, columns, nthreads=2*multiprocessing.cpu_count()):
   thread_pool = concurrent.futures.ThreadPoolExecutor(nthreads)
   if dataset not in input_files:
      raise ValueError('The dataset %s does not exist, I have %s' % (dataset, ', '.join(input_files.keys())))
   infiles = [uproot.open(i) for i in input_files[dataset]]
   ret = None
   arrays = [i['features/tree'].arrays(columns, executor=thread_pool, blocking=False) for i in infiles]
   ret = arrays[0]()
   for arr in arrays[1:]:
      tmp = arr()
      for column in columns:
         ret[column] = np.concatenate((ret[column],tmp[column]))
   return ret

from sklearn.cluster import KMeans
from sklearn.externals import joblib
import json
from pdb import set_trace

apply_weight = np.vectorize(lambda x, y: y.get(x), excluded={2})

def kmeans_weighter(features, fname):
   kmeans = joblib.load(fname)
   cluster = kmeans.predict(features)
   str_weights = json.load(open(fname.replace('.plk', '.json')))   
   weights = {}
   for i in str_weights:
      try:
         weights[int(i)] = str_weights[i]
      except:
         pass
   return apply_weight(cluster, weights)
