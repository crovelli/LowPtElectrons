# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECOWithDQM -s RAW2DIGI,L1Reco,RECO --datatier RECO --runUnscheduled --nThreads 4 --data --era Run2_2017 --scenario pp --conditions 94X_dataRun2_PromptLike_v9 --eventcontent RECO --filein file:/afs/cern.ch/work/m/mverzett/public/RAW-RECO_ZElectron-94X_dataRun2_PromptLike_v5_RelVal_doubEG2017B-v1.root --no_exec

################################################################################
# Command line args

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
options.register('fakesMultiplier', None,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    ""
)
options.register('data', 'RAWMCTest',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('globalTag', '102X_upgrade2018_realistic_v15',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    ""
)
options.register('hitAssociation', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.register('disableAssociation', True,
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
options.register('MVANtuplizer', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    ""
)
options.setDefault('maxEvents', -1)
options.parseArguments()

################################################################################
# Import samples

from LowPtElectrons.LowPtElectrons.samples import all_samples
#split into even chunks
input_files = all_samples[options.data] if not options.inputFiles else options.inputFiles
n = len(input_files)/options.nchunks
chunks = [input_files[i:i + n] for i in xrange(0, len(input_files), n)]
leftovers = sum(chunks[options.nchunks:], [])
chunks = chunks[:options.nchunks]
for i in range(len(leftovers)):
   chunks[i].append(leftovers[i])

################################################################################
# Standard RECO

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2

process = cms.Process('RECO', eras.Run2_2018, eras.bParking, premix_stage2)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.DigiDMPreMix_cff')
# process.load('SimGeneral.MixingModule.digi_MixPreMix_cfi')
# process.load('Configuration.StandardSequences.DataMixerPreMix_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
process.load('Configuration.StandardSequences.DigiToRawDM_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')

################################################################################
# Input files and events 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
   "PoolSource",
   fileNames = cms.untracked.vstring(
      chunks[options.ichunk]
      ),
   secondaryFileNames = cms.untracked.vstring(),
   skipEvents=cms.untracked.uint32(options.skipEvents)
   )

if options.pick:
   process.source.eventsToProcess = cms.untracked.VEventRange(options.pick)

################################################################################
# Options

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(False),
   #SkipEvent = cms.untracked.vstring('ProductNotFound'),
   #numberOfThreads=cms.untracked.uint32(4)
   #options.numberOfStreams=cms.untracked.uint32(0)
)

################################################################################
# Production Info

process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECOWithDQM nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

################################################################################
# Other statements

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
#'94X_dataRun2_PromptLike_v9', '')

################################################################################
# Ntuplizer steps

if options.MVANtuplizer == True : # Use the Egamma-default MVANtuplizer code

   process.load('LowPtElectrons/LowPtElectrons/ElectronMVANtuplizer_cfi')

else : # Use custom Ntuplizer code

   # Track association by hits
   process.load('SimTracker/TrackAssociation/trackingParticleRecoTrackAsssociation_cfi')
   process.quickTrackAssociatorByHits = cms.EDProducer(
      "QuickTrackAssociatorByHitsProducer",
      AbsoluteNumberOfHits = cms.bool(False),
      Cut_RecoToSim = cms.double(0.75),
      PixelHitWeight = cms.double(1.0),
      Purity_SimToReco = cms.double(0.75),
      Quality_SimToReco = cms.double(0.5),
      SimToRecoDenominator = cms.string('reco'),
      ThreeHitTracksAreSpecial = cms.bool(True),
      associatePixel = cms.bool(True),
      associateStrip = cms.bool(True),
      cluster2TPSrc = cms.InputTag("tpClusterProducer"),
      pixelSimLinkSrc = cms.InputTag("mixData","PixelDigiSimLink"),
      stripSimLinkSrc = cms.InputTag("mixData","StripDigiSimLink"),
      useClusterTPAssociation = cms.bool(False)
      )
   ## .useClusterTPAssociation = False
   ## process.quickTrackAssociatorByHits.associatePixel = cms.bool(True)
   ## process.quickTrackAssociatorByHits.associateStrip = cms.bool(True)
   ## process.quickTrackAssociatorByHits.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis")
   ## process.quickTrackAssociatorByHits.stripSimLinkSrc = cms.InputTag("simSiStripDigis")
   process.ntuplizer_seq *= process.tpClusterProducer
   process.ntuplizer_seq *= process.quickTrackAssociatorByHits
   process.ntuplizer_seq *= process.trackingParticleRecoTrackAsssociation
   #from SimGeneral.DataMixingModule.customiseForPremixingInput import customiseForPreMixingInput
   #customiseForPreMixingInput(process)

   # Electron ID ...
   from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
   dataFormat = DataFormat.AOD
   switchOnVIDElectronIdProducer(process, dataFormat)
   my_id_modules = [
      'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
   ]

   # Add them to the VID producer
   for idmod in my_id_modules :
      setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

   # Change inputs
   process.electronMVAVariableHelper.src = 'lowPtGsfElectrons'
   process.electronMVAValueMapProducer.src = 'lowPtGsfElectrons'
   process.ntuplizer_seq *= process.electronMVAVariableHelper
   process.ntuplizer_seq *= process.electronMVAValueMapProducer

   # Make seeding pass through 
   process.lowPtGsfElectronSeeds.PassThrough = True
   # process.lowPtGsfElectronSuperClusters.MaxDeltaR2 = 9999.

   # Custom Ntuplizer code 
   process.load('LowPtElectrons.LowPtElectrons.LowPtGsfElectronsAnalyzer_cfi')
   process.load('LowPtElectrons.LowPtElectrons.TrackerElectronsFeatures_cfi')
   process.features.hitAssociation = options.hitAssociation
   if options.fakesMultiplier : process.features.fakesMultiplier = options.fakesMultiplier
   process.features.disableAssociation = options.disableAssociation
   process.features.checkFromB = options.checkFromB


################################################################################
# Path and EndPath definitions, TFileService, OutputModule

# ReReco and ntuplize
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
if options.MVANtuplizer == True : 
   process.reconstruction_step *= cms.Sequence(process.ntuplizer)
else :
   process.reconstruction_step *= cms.Sequence(process.ntuplizer_seq)
   process.reconstruction_step *= cms.Sequence(process.simple)
   process.reconstruction_step *= cms.Sequence(process.features)
process.eventinterpretaion_step = cms.Path(process.EIsequence)

process.schedule = cms.Schedule(
   process.raw2digi_step,
   process.reconstruction_step,
   process.recosim_step,
   process.eventinterpretaion_step
   )

# Write ntuple to root file called "options.outname"
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outname))

# EDM output
if options.edout:
   process.AODSIMoutput = cms.OutputModule(
      "PoolOutputModule",
      outputCommands = process.AODSIMEventContent.outputCommands,
      #outputCommands = process.AODEventContent.outputCommands,
      #outputCommands = cms.untracked.vstring('keep *',)
      compressionAlgorithm = cms.untracked.string('LZMA'),
      compressionLevel = cms.untracked.int32(4),
      eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
      fileName = cms.untracked.string('file:EDOutput.root'),
      outputCommands = cms.untracked.vstring( 
         'drop *',
         'keep *_lowPt*_*_*',
         'keep *_generalTracks_*_*',
#         "keep *_offlineBeamSpot_*_*",
#         "keep *_genParticles_*_*",
#         'keep *_generalTracks_*_*',
#         "keep *_trackingParticleRecoTrackAsssociation_*_*",
#         "keep *_trackerDrivenElectronSeeds*_*_*",
#         "keep *_electronCkfTrackCandidates*_*_*",
#         "keep *_electronGsfTracks*_*_*",
#         'keep *_pfTrackElec*_*_*',
#         'keep *_pfTrack*_*_*',
#         "keep *_*GsfElectronCores*_*_*",
#         "keep *_*GsfElectrons*_*_*",
#         'keep *_mvaElectrons_*_*',
#         'keep *_reducedEcalRecHitsEB_*_*',
#         # missing futher collections? remove collections below?
#         "keep *_particleFlowEGamma_*_*",
#         'keep *_pfDisplacedTrackerVertex_*_*',
#         'keep *_particleFlowBlock_*_*',
#         'keep *_particleFlowSuperClusterECAL_*_*',
#         'keep *_particleFlowClusterECAL_*_*',
#         'keep *_particleFlowClusterHCAL_*_*',
#         'keep *_particleFlowClusterHO_*_*',
#         'keep *_particleFlowClusterHF_*_*',
#         'keep *_particleFlowClusterPS_*_*',
         )
      )
   process.end = cms.EndPath(process.AODSIMoutput)
   process.schedule.append(process.end)

################################################################################
# Expert settings ...

#process.Timing = cms.Service(
#   "Timing",
#   summaryOnly = cms.untracked.bool(False),
#   useJobReport = cms.untracked.bool(True)
#   )

# PAT stuff
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
   
# Do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)

# Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

open('pydump.py','w').write(process.dumpPython())
