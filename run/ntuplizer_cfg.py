import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output.root')
options.setDefault('maxEvents',-1)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.register('addSkim',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST') 

# AOD: 9558 + 6230 + 1661 = 17449 events
# MINOAOD: 17449 events
default_file = 'root://xrootd-cms.infn.it//store/relval/CMSSW_11_0_0_pre7/RelValProdZEE_13_pmx25ns/AODSIM/PUpmx25ns_110X_mc2017_realistic_v1-v1/20000/97D8818E-768D-324E-95E2-E14B4370E920.root' if options.useAOD else 'root://xrootd-cms.infn.it//store/relval/CMSSW_11_0_0_pre7/RelValProdZEE_13_pmx25ns/MINIAODSIM/PUpmx25ns_110X_mc2017_realistic_v1-v1/20000/9331CE22-9A6B-C44F-B3BB-FDD77A0413D7.root' 

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(default_file),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mc2017_realistic_v1', '')

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
    )

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('output_filtered.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ntuplizer_path'))
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
if options.useAOD is True : switchOnVIDElectronIdProducer(process,DataFormat.AOD)
else :                      switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in [
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
    ] : 
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplizer_seq = cms.Sequence()

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi')

if options.useAOD is False : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.ntuplizer_seq *= process.lowPtGsfElectronID

process.load('LowPtElectrons.LowPtElectrons.IDNtuplizer_cfi')
process.ntuplizer_seq *= process.ntuplizer

process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.output_path = cms.EndPath(process.output)
if options.addSkim is False : process.schedule = cms.Schedule(process.ntuplizer_path)
else : process.schedule = cms.Schedule(process.ntuplizer_path,process.output_path)
