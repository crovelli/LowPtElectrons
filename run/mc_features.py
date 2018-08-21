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
options.setDefault('maxEvents', -1)
options.parseArguments()

from LowPtElectrons.LowPtElectrons.samples import all_samples
#split into even chunks
input_files = all_samples[options.data]
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
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      chunks[options.ichunk]
      #'/store/data/Run2017F/SingleElectron/RAW/v1/000/306/459/00000/B6410DA6-C8C5-E711-947F-02163E01A45E.root'
      ),
    secondaryFileNames = cms.untracked.vstring()
)

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
process.trackerDrivenElectronSeeds.MinPt = 0.
process.trackerDrivenElectronSeeds.PtThresholdSavePreId = cms.untracked.double(0.)
process.trackerDrivenElectronSeeds.ProducePreId = True
#remove ECAL-driven seeds
process.ecalDrivenElectronSeeds._TypedParameterizable__type = 'EmptySeedProducer'
process.ecalDrivenElectronSeedsFromMultiCl._TypedParameterizable__type = 'EmptySeedProducer'
#set association to match against GSF, which is now run on every single track
process.trackingParticleRecoTrackAsssociation.label_tr = 'electronGsfTracks'

# https://github.com/ICBPHCMS/cmssw/blob/CMSSW_9_4_X/TrackingTools/GsfTracking/python/CkfElectronCandidateMaker_cff.py
# total hack - not checked carefully (ie are max/mins set correct/adequate?)
# combined, they do have an effect
process.TrajectoryFilterForElectrons.chargeSignificance = 0.
process.TrajectoryFilterForElectrons.minPt = 0.
process.TrajectoryFilterForElectrons.minHitsMinPt = -999
process.TrajectoryFilterForElectrons.maxLostHits = 999
process.TrajectoryFilterForElectrons.maxNumberOfHits = 999
process.TrajectoryFilterForElectrons.maxConsecLostHits = 999
process.TrajectoryFilterForElectrons.nSigmaMinPt = 0.
process.TrajectoryFilterForElectrons.minimumNumberOfHits = -999
process.TrajectoryFilterForElectrons.maxCCCLostHits = 999

process.GsfElectronFittingSmoother.MinNumberOfHits = 2 #does not change anything
#process.electronTrajectoryCleanerBySharedHits.fractionShared = 0.9
#
# PUT THE NTUPLIZER HERE!
#
process.load('LowPtElectrons.LowPtElectrons.TrackerElectronsFeatures_cfi')

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

open('pydump.py','w').write(process.dumpPython())
