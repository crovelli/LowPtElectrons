#A single place where to bookkeep the dataset file locations
tag = '2018Aug23'
input_files = {
   'BToKeeByDR' : '/eos/cms/store/cmst3/group/bpark/low_pt_datasets/BToKee_assocByDR_2018Aug23.root',
   'BToKeeByHits' : '/eos/cms/store/cmst3/group/bpark/low_pt_datasets/BToKee_assocByHits_2018Aug23.root',
   'test' : '/afs/cern.ch/work/m/mverzett/RK94New/src/LowPtElectrons/LowPtElectrons/run/track_features.root',
}

import os
if not os.path.isdir('plots/%s' % tag):
   os.mkdir('plots/%s' % tag)
