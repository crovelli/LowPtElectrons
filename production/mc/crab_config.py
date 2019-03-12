import os
import datetime

if not 'REQUEST' in os.environ:
   raise RuntimeError('REQUEST must be set as an environmental variable')

request = os.environ['REQUEST']
samples = {
   'BToJPsieeK'    : '/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/GEN-SIM-RAW',
   'BToKStaree' : '/BdToKstar_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/GEN-SIM-RAW',
   'BsToJPsieePhi' : '/BsToPhiJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18DR-PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/GEN-SIM-RAW',
}

splitting = {
   'BToJPsieeK'  : 1,
   'BToKStaree'  : 1,
   'BsToJPsieePhi' : 1,
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
   if not os.path.isdir('crab_%s' % request):
      break
   request = '%sv%d_%s' % (date, idx, name)
   idx += 1

from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = request

config.section_('JobType')
config.JobType.psetName = os.path.abspath('%s/../../run/mc_features.py' % os.environ['PWD'])
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = [
   'outname=%s_assocBy%s.root' % (sample, association), 
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
#config.Site.whitelist = ['T2_UK_London_IC']
#config.Data.inputDBS = 'phys03'
config.Site.blacklist = ['T1_DE_KIT', 'T1_UK_RAL']
