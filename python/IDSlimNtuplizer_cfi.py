import FWCore.ParameterSet.Config as cms

ntuplizer = cms.EDAnalyzer(
    "IDSlimNtuplizer",
    verbose = cms.int32(0),
    checkFromB = cms.bool(True),
    prescale = cms.double(0.01),
    minTrackPt = cms.double(0.5),  
    maxTrackPt = cms.double(15.),  
    maxTrackEta = cms.double(2.4),  
    # Generic collections
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    genParticles = cms.InputTag("genParticles"), # AOD
    prunedGenParticles = cms.InputTag("prunedGenParticles"), # MINIAOD
    ctfTracks = cms.InputTag("generalTracks"),
    packedCands = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    ebRecHits = cms.InputTag('reducedEcalRecHitsEB'),
    eeRecHits = cms.InputTag('reducedEcalRecHitsEE'),
    ebRecHitsEGM = cms.InputTag('reducedEgamma','reducedEBRecHits'),
    eeRecHitsEGM = cms.InputTag('reducedEgamma','reducedEERecHits'),
    # Low pT collections
    gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
    gsfElectrons = cms.InputTag("lowPtGsfElectrons"), # AOD 
    patElectrons = cms.InputTag("slimmedLowPtElectrons"), # MINIAOD
    gsfTrackLinks = cms.InputTag("lowPtGsfToTrackLinks"), # AOD
    packedCandLinks = cms.InputTag("lowPtGsfLinks:packedCandidates"), # mAOD
    lostTrackLinks = cms.InputTag("lowPtGsfLinks:lostTracks"), # mAOD
    mvaUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
    mvaPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
    mvaValueLowPt = cms.InputTag('lowPtGsfElectronID'),
    dEdx1Tag = cms.InputTag('dedxHarmonic2'),
    # EGamma collections
    #eleSeeds = cms.InputTag("lowPtGsfElectronSeeds"),
    #preIdsEcal = cms.InputTag("lowPtGsfElectronSeeds"),
    #preIdsHcal = cms.InputTag("lowPtGsfElectronSeeds:HCAL"),
    #preIdRefs = cms.InputTag("lowPtGsfElectronSeeds"),
    #eleSeedsEGamma = cms.InputTag("electronMergedSeeds"), # AOD   # trackerDrivenElectronSeeds:SeedsForGsf
    gsfTracksEGamma = cms.InputTag("electronGsfTracks"),                   # AOD
    gsfTracksEGamma_MAOD = cms.InputTag("reducedEgamma:reducedGsfTracks"), # MINIAOD
    gsfElectronsEGamma = cms.InputTag("gedGsfElectrons"),                  # AOD
    patElectronsEGamma = cms.InputTag("slimmedElectrons"),                 # MINIAOD
    mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
    mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90') # wp80?
)
