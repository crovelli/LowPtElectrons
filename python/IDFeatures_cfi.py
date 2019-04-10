import FWCore.ParameterSet.Config as cms

features = cms.EDAnalyzer(
    "IDFeatures",
    checkFromB = cms.bool(True),
    drMax = cms.double(0.02),
    fakesMultiplier = cms.double(6.),
    # AOD and MINIAOD
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
    MVASeedUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
    MVASeedPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
    MVAIDLowPt = cms.InputTag('lowPtGsfElectronID'),
    MVAIDV2 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values'),
    # AOD only
    egammaGsfTracks = cms.InputTag("electronGsfTracks"),
    electrons = cms.InputTag("lowPtGsfElectrons"),
    egammaElectrons = cms.InputTag("gedGsfElectrons"),
    genParticles = cms.InputTag("genParticles"),
    # MINIAOD only
    egammaGsfTracks_mAOD = cms.InputTag("reducedEgamma:reducedGsfTracks"),
    electrons_mAOD = cms.InputTag("slimmedLowPtElectrons"),
    egammaElectrons_mAOD = cms.InputTag("slimmedElectrons"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    )

    #convVtxFitProb = cms.InputTag('electronMVAVariableHelper:convVtxFitProb'),
