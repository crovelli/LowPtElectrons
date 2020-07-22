import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.setDefault('maxEvents',100)
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST')

default_file = 'root://xrootd-cms.infn.it//store/relval/CMSSW_11_0_0_pre7/RelValProdZEE_13_pmx25ns/AODSIM/PUpmx25ns_110X_mc2017_realistic_v1-v1/20000/97D8818E-768D-324E-95E2-E14B4370E920.root' if options.useAOD else 'root://xrootd-cms.infn.it//store/relval/CMSSW_11_0_0_pre7/RelValProdZEE_13_pmx25ns/MINIAODSIM/PUpmx25ns_110X_mc2017_realistic_v1-v1/20000/9331CE22-9A6B-C44F-B3BB-FDD77A0413D7.root' 

output_file = 'output_fromAOD_slim11x.root' if options.useAOD else 'output_fromMINIAOD_slim11x.root'

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

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi')
if not options.useAOD : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.ntuplizer_seq *= process.lowPtGsfElectronID

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mc2017_realistic_v1')
# this is for the LowPt energy regression
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")), 
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")), 
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))


process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('LowPtElectrons.LowPtElectrons.IDSlimNtuplizer_cfi')

process.ntuplizer_seq *= process.ntuplizer

process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.schedule = cms.Schedule(process.ntuplizer_path)
