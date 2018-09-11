# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECOWithDQM -s RAW2DIGI,L1Reco,RECO --datatier RECO --runUnscheduled --nThreads 4 --data --era Run2_2017 --scenario pp --conditions 94X_dataRun2_PromptLike_v9 --eventcontent RECO --filein file:/afs/cern.ch/work/m/mverzett/public/RAW-RECO_ZElectron-94X_dataRun2_PromptLike_v5_RelVal_doubEG2017B-v1.root --no_exec
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outname', 'track_features.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('ichunk', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.register('nchunks', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.register('fakePrescale', 0.08,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    ""
)
options.register('data', 'RAWMCTest',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('globalTag', '94X_mc2017_realistic_v12',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('hitAssociation', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register('disableAssociation', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register('checkFromB', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register('edout', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register(
   'pick',
   '',
   VarParsing.multiplicity.list,
   VarParsing.varType.string,
   'Pick single events'
)
options.register('matchAtSeeding', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register('skipEvents', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    ""
)
options.setDefault('maxEvents', -1)
options.parseArguments()

from LowPtElectrons.LowPtElectrons.samples import all_samples
#split into even chunks
input_files = all_samples[options.data] if not options.inputFiles else options.inputFiles
n = len(input_files)/options.nchunks
chunks = [input_files[i:i + n] for i in xrange(0, len(input_files), n)]
leftovers = sum(chunks[options.nchunks:], [])
chunks = chunks[:options.nchunks]
for i in range(len(leftovers)):
   chunks[i].append(leftovers[i])


import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DigiDMPreMix_cff')
process.load('SimGeneral.MixingModule.digi_MixPreMix_cfi')
process.load('Configuration.StandardSequences.DataMixerPreMix_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
process.load('Configuration.StandardSequences.DigiToRawDM_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
   "PoolSource",
   fileNames = cms.untracked.vstring(
      chunks[options.ichunk]
      #'/store/data/Run2017F/SingleElectron/RAW/v1/000/306/459/00000/B6410DA6-C8C5-E711-947F-02163E01A45E.root'
      ),
   secondaryFileNames = cms.untracked.vstring(),
   skipEvents=cms.untracked.uint32(options.skipEvents)
   )

if options.pick:
   process.source.eventsToProcess = cms.untracked.VEventRange(options.pick)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECOWithDQM nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#
# Low pT Electron customization
#
process.electronFeatures = cms.Sequence()
process.load('SimTracker/TrackAssociation/trackingParticleRecoTrackAsssociation_cfi')
#process.reconstruction *= process.simSiPixelDigis
process.electronFeatures *= process.tpClusterProducer
process.electronFeatures *= process.quickTrackAssociatorByHits
process.quickTrackAssociatorByHits.useClusterTPAssociation = False
#process.quickTrackAssociatorByHits.associateStrip = False
process.electronFeatures *= process.trackingParticleRecoTrackAsssociation
from SimGeneral.DataMixingModule.customiseForPremixingInput import customiseForPreMixingInput
customiseForPreMixingInput(process)

#hack through the electron code
#make the tracker driven pass-though for pt > 0.5, produce preID
process.trackerDrivenElectronSeeds._TypedParameterizable__type = 'PassThroughTrackSeeds'
process.trackerDrivenElectronSeeds.MinPt = 1. ##@@
#process.trackerDrivenElectronSeeds.PtThresholdSavePreId = cms.untracked.double(0.) ##@@
process.trackerDrivenElectronSeeds.ProducePreId = True
process.trackerDrivenElectronSeeds.matchToGen = cms.bool(options.matchAtSeeding)
process.trackerDrivenElectronSeeds.genParticles = cms.InputTag("genParticles")

#remove ECAL-driven seeds
process.ecalDrivenElectronSeeds._TypedParameterizable__type = 'EmptySeedProducer'
process.ecalDrivenElectronSeedsFromMultiCl._TypedParameterizable__type = 'EmptySeedProducer'
#set association to match against GSF, which is now run on every single track
process.trackingParticleRecoTrackAsssociation.label_tr = 'electronGsfTracks'

# https://github.com/ICBPHCMS/cmssw/blob/CMSSW_9_4_X/TrackingTools/GsfTracking/python/CkfElectronCandidateMaker_cff.py
# total hack - not checked carefully (ie are max/mins set correct/adequate?)
# combined, they do have an effect
#process.TrajectoryFilterForElectrons.chargeSignificance = 0.
#process.TrajectoryFilterForElectrons.minPt = 0.
#process.TrajectoryFilterForElectrons.minHitsMinPt = -999
#process.TrajectoryFilterForElectrons.maxLostHits = 999
#process.TrajectoryFilterForElectrons.maxNumberOfHits = 999
#process.TrajectoryFilterForElectrons.maxConsecLostHits = 999
#process.TrajectoryFilterForElectrons.nSigmaMinPt = 0.
#process.TrajectoryFilterForElectrons.minimumNumberOfHits = -999
#process.TrajectoryFilterForElectrons.maxCCCLostHits = 999

process.GsfElectronFittingSmoother.MinNumberOfHits = 2 #does not change anything
#process.electronTrajectoryCleanerBySharedHits.fractionShared = 0.9
#
# PUT THE NTUPLIZER HERE!
#
process.load('LowPtElectrons.LowPtElectrons.TrackerElectronsFeatures_cfi')
process.features.hitAssociation = options.hitAssociation
process.features.prescaleFakes = options.fakePrescale
process.features.disableAssociation = options.disableAssociation
process.features.checkFromB = options.checkFromB
# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#'94X_dataRun2_PromptLike_v9', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.electronFeatures *= process.features
process.reconstruction_step *= process.electronFeatures
#process.reconstruction_step *= process.features
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.schedule = cms.Schedule(
   process.raw2digi_step, process.reconstruction_step, 
   process.recosim_step, process.eventinterpretaion_step
   )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)

if options.edout:
   process.AODSIMoutput = cms.OutputModule(
      "PoolOutputModule",
      compressionAlgorithm = cms.untracked.string('LZMA'),
      compressionLevel = cms.untracked.int32(4),
      eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
      fileName = cms.untracked.string('file:EDOutput.root'),
      outputCommands = cms.untracked.vstring( 
         'keep *',
         # 'drop *', 
         # "keep *_offlineBeamSpot_*_*",
         # "keep *_genParticles_*_*",
         # "keep *_electronGsfTracks_*_*",
         # "keep *_particleFlowEGamma_*_*",
         # "keep *_gedGsfElectronCores_*_*",
         # "keep *_gedGsfElectrons_*_*",
         # "keep *_trackingParticleRecoTrackAsssociation_*_*",
         # "keep *_electronCkfTrackCandidates_*_*",
         # "keep *_trackerDrivenElectronSeeds_*_*",
         # 'keep *_generalTracks_*_*',
         # 'keep *_generalTracks_*_*',
         # 'keep *_particleFlowEGamma_*_*',
         # 'keep *_mvaElectrons_*_*',
         # 'keep *_particleFlowBlock_*_*',
         # 'keep *_particleFlowSuperClusterECAL_*_*',
         # 'keep *_pfTrackElec_*_*',
         # 'keep *_pfDisplacedTrackerVertex_*_*',
         # 'keep *_pfTrack_*_*',
         # 'keep *_particleFlowClusterECAL_*_*',
         # 'keep *_particleFlowClusterHCAL_*_*',
         # 'keep *_particleFlowClusterHO_*_*',
         # 'keep *_particleFlowClusterHF_*_*',
         # 'keep *_particleFlowClusterPS_*_*',
         )
      )
   process.end = cms.EndPath(process.AODSIMoutput)
   process.schedule.append(process.end)
   
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(False),
      #SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

# Write ntuple to root file called "options.outname" 
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outname))

#process.pfTrackElec.debugGsfCleaning = True
process.pfTrackElec.applyGsfTrackCleaning = False
