import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('inputFiles','file:input.root')
options.setDefault('outputFile','output.root')
options.setDefault('maxEvents',-1)
options.parseArguments()

process = cms.Process('TEST')

process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring(options.inputFiles),
    #fileNames = cms.untracked.vstring(files[:4]),
    #fileNames = cms.untracked.vstring('/afs/cern.ch/user/b/bainbrid/work/public/5-full-chain/20-retrain-id/CMSSW_10_2_13/src/1-miniaod-from-crab'),
    fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18RECOBParking/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_BParking_Bparking_102X_upgrade2018_realistic_v15-v1/110000/968D9C40-A196-C746-B16C-27E90DFC17DB.root'),
    #fileNames = cms.untracked.vstring('/store/data/Run2018A/ParkingBPH1/AOD/22Mar2019-v1/260005/A032FCE0-D492-D94E-9404-EF96EB3A84BB.root'), # data, AOD
    #fileNames = cms.untracked.vstring('/store/data/Run2018A/ParkingBPH1/MINIAOD/22Mar2019-v1/260002/54516928-947E-6140-A489-4E4099A593CF.root'), # data, MINIAOD
    secondaryFileNames = cms.untracked.vstring()
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
    )

process.ntuplizer_seq = cms.Sequence()

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] : 
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.electronMVAVariableHelper.srcMiniAOD = 'slimmedLowPtElectrons'
#process.electronMVAVariableHelper.vertexCollectionMiniAOD = 'offlineSlimmedPrimaryVertices'
#process.electronMVAVariableHelper.conversionsMiniAOD = 'gsfTracksOpenConversions'
#process.electronMVAVariableHelper.beamSpotMiniAOD = 'offlineBeamSpot'
process.ntuplizer_seq *= process.electronMVAVariableHelper

process.electronMVAValueMapProducer.srcMiniAOD = 'slimmedLowPtElectrons'
process.ntuplizer_seq *= process.electronMVAValueMapProducer

process.load('LowPtElectrons.LowPtElectrons.IDFeatures_cfi')
process.ntuplizer_seq *= process.features

process.ntuplizer_path = cms.Path(process.ntuplizer_seq)
process.schedule = cms.Schedule(process.ntuplizer_path)
