import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.setDefault('maxEvents',1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST')

default_file = 'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18RECOBParking/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/FA5679D9-1BB8-304C-AA78-14398A3346F4.root' if options.useAOD else 'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/AC0EFF4C-C75F-274E-A561-A19D20DB2668.root'
#default_file = '/store/mc/RunIIAutumn18RECOBParking/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/FA5679D9-1BB8-304C-AA78-14398A3346F4.root' if options.useAOD else 'root://xrootd-cms.infn.it//store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/70000/AC0EFF4C-C75F-274E-A561-A19D20DB2668.root'

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

process.slimntuplizer_seq = cms.Sequence()

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
if options.useAOD : switchOnVIDElectronIdProducer(process,DataFormat.AOD)
else :           switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] :
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.slimntuplizer_seq = cms.Sequence()

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if not options.useAOD : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.slimntuplizer_seq *= process.lowPtGsfElectronID

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('LowPtElectrons.LowPtElectrons.IDSlimNtuplizer_cfi')
#process.load("Configuration.Geometry.GeometryReco_cff");
#process.load("Configuration.Geometry.GeometryDB_cff");
#process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load("Configuration.StandardSequences.GeometryDB_cff");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
#process.load("Configuration.StandardSequences.Reconstruction_cff");



process.slimntuplizer_seq *= process.slimntuplizer

process.slimntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.slimntuplizer_seq)
process.schedule = cms.Schedule(process.slimntuplizer_path)
