from WMCore.Configuration import Configuration
import os
import datetime

date = datetime.date.today().strftime('%Y%b%d')
request = '%s_BToKeeAssocByDR' % date
idx = 2
while True:
   if not os.path.isdir(request):
      break
   request = '%sv%d_BToKeeAssocByDR' % (date, idx)
   idx += 1

config = Configuration()
config.section_('General')
config.General.requestName = request

config.section_('JobType')
config.JobType.psetName = '%s/src/LowPtElectrons/LowPtElectrons/run/mc_features.py' % os.environ['CMSSW_BASE']
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['outname=BToKee_assocByDR.root', 'hitAssociation=False']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/cmst3/user/mverzett/'
config.Data.splitting = 'FileBased'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.Data.inputDataset = '/BToKee_Pythia/tstreble-BToKee_Pythia_PUMix_18_03_17-c9b9e020b5bce5ee6bee9ef5f38c415a/USER'
config.Site.whitelist = ['T2_UK_London_IC']
config.Data.inputDBS = 'phys03'
