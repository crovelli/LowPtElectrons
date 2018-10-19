import os
import datetime

if not 'REQUEST' in os.environ:
   raise RuntimeError('REQUEST must be set as an environmental variable')

request = os.environ['REQUEST']
samples = {
   'BToJPsieeK'    : '/Bu_KJPsi_ee_Pythia_GEN-SIM_18_06_4/tstreble-Bu_KJPsi_ee_Pythia_PUMix_18_06_21-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
   'BToJPsieeK_v1' : '/Bu_KJPsi_ee_Pythia/tstreble-Bu_KJPsi_ee_Pythia_PUMix_18_06_03-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
   'BToJPsieeKStar': '/Bd_KstJPsi_ee_Pythia_GEN-SIM_18_07_01/tstreble-Bd_KstJPsi_ee_PUMix_18_07_02-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
   'BToKStaree' : '/BToKstee_Pythia/tstreble-BToKstee_Pythia_PUMix_18_04_17-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
   'BToKee_v1' : '/BToKee_Pythia/tstreble-BToKee_Pythia_PUMix_18_03_18-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
   'BToKee_v2' : '/BToKee_Pythia/tstreble-BToKee_Pythia_PUMix_18_03_17-c9b9e020b5bce5ee6bee9ef5f38c415a/USER',
}

splitting = {
   'BToJPsieeK' : 2, 
   'BToKStaree' : 5,
   'BToJPsieeK_v1' : 7,
   'BToJPsieeKStar' : 7,
}

association = ''
if request.endswith('AssocByDR'):
   association = 'DR'
elif request.endswith('AssocByHits'):
   association = 'Hits'
else:
   raise RuntimeError('The request can only end in AssocByDR or AssocByHits')

sample = request.split('Assoc')[0]
if sample not in samples:
   raise RuntimeError('Sample not found, I have: %s' % ', '.join(samples.keys()))

name = '%sAssocBy%s' % (sample, association)
dataset = samples[sample]

date = datetime.date.today().strftime('%Y%b%d')
request = '%s_%s' % (date, name)
idx = 2
while True:
   if not os.path.isdir(request):
      break
   request = '%sv%d_%s' % (date, idx, name)
   idx += 1

from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = request

config.section_('JobType')
config.JobType.psetName = '%s/src/LowPtElectrons/LowPtElectrons/run/mc_features.py' % os.environ['CMSSW_BASE']
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = [
   'outname=BToKee_assocBy%s.root' % association, 
   ('hitAssociation=False' if association == 'DR' else 'hitAssociation=True')]
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.unitsPerJob = splitting.get(sample, 20)
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/cmst3/user/mverzett/'
config.Data.splitting = 'FileBased'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.Data.inputDataset = dataset
config.Site.whitelist = ['T2_UK_London_IC']
config.Data.inputDBS = 'phys03'
