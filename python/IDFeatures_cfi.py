import FWCore.ParameterSet.Config as cms

features = cms.EDAnalyzer(
    "IDFeatures",
    verbose = cms.int32(0),
    checkFromB = cms.bool(True),
    drMax = cms.double(0.1),
    drThreshold = cms.double(0.02),
    prescale = cms.double(0.01),
    minTrackPt = cms.double(0.5),
    # Generic collections
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    genParticles = cms.InputTag("genParticles"), # AOD
    prunedGenParticles = cms.InputTag("prunedGenParticles"), # MAOD
    packedGenParticles = cms.InputTag("packedGenParticles"), # MAOD
    ctfTracks = cms.InputTag("generalTracks"),
    ebRecHits = cms.InputTag('reducedEcalRecHitsEB'),
    eeRecHits = cms.InputTag('reducedEcalRecHitsEE'),
    barrelSuperClusters = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    endcapSuperClusters = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    # Low pT collections
    eleSeeds = cms.InputTag("lowPtGsfElectronSeeds"),
    preIdsEcal = cms.InputTag("lowPtGsfElectronSeeds"),
    preIdsHcal = cms.InputTag("lowPtGsfElectronSeeds:HCAL"),
    preIdRefs = cms.InputTag("lowPtGsfElectronSeeds"),
    gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
    gsfElectrons = cms.InputTag("lowPtGsfElectrons"), # AOD 
    gsfElectrons_MAOD = cms.InputTag("slimmedLowPtElectrons"), # MAOD
    gsfTrackLinks = cms.InputTag("lowPtGsfToTrackLinks"), # AOD
    packedCandLinks = cms.InputTag("lowPtGsfToTrackLinks:packedCandidates"), # mAOD
    lostTrackLinks = cms.InputTag("lowPtGsfToTrackLinks:lostTracks"), # mAOD
    mvaUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
    mvaPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
    mvaValueLowPt = cms.InputTag('lowPtGsfElectronID'),
    #mvaValue = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
    #mvaId = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp80'),
    # EGamma collections
    eleSeedsEGamma = cms.InputTag("electronMergedSeeds"), # AOD   # trackerDrivenElectronSeeds:SeedsForGsf
    gsfTracksEGamma = cms.InputTag("electronGsfTracks"), # AOD
    gsfTracksEGamma_MAOD = cms.InputTag("reducedEgamma:reducedGsfTracks"), # MAOD
    gsfElectronsEGamma = cms.InputTag("gedGsfElectrons"), # AOD
    gsfElectronsEGamma_MAOD = cms.InputTag("slimmedElectrons"), # MAOD
    mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
    mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'),
    )

    #convVtxFitProb = cms.InputTag('electronMVAVariableHelper:convVtxFitProb'),
