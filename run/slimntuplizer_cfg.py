import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST')

default_file = '/store/mc/RunIIAutumn18RECOBParking/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/FA5679D9-1BB8-304C-AA78-14398A3346F4.root' if options.useAOD else 'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/854B1DC0-2F71-694D-A3F5-8DC1CDE1EF18.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/6056E15A-6AFA-6440-8493-4870F53B71EC.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/491322F0-D192-434E-BC85-7D83DF9482EB.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/4134382B-616B-1847-B991-1B1C0BD409BD.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/30000/8493D2EB-E9BC-FB42-82D4-8E61841ECC93.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/30000/E4E2CB52-92AF-7F4C-8C0D-5932410C8D2C.root','root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/30000/81B1F3E9-C4D2-3D49-8D8F-132F748C893E.root'

output_file = 'output_fromAOD_slim.root' if options.useAOD else 'output_fromMINIAOD_slim.root'

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(default_file),
    secondaryFileNames = cms.untracked.vstring()
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(output_file)
    )

process.ntuplizer_seq = cms.Sequence()

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
if options.useAOD : switchOnVIDElectronIdProducer(process,DataFormat.AOD)
else :           switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] :
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.ntuplizer_seq = cms.Sequence()

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if not options.useAOD : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.ntuplizer_seq *= process.lowPtGsfElectronID

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('LowPtElectrons.LowPtElectrons.IDSlimNtuplizer_cfi')

process.ntuplizer_seq *= process.ntuplizer

process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.schedule = cms.Schedule(process.ntuplizer_path)
