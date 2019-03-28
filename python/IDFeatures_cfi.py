import FWCore.ParameterSet.Config as cms

features = cms.EDAnalyzer(
    "IDFeatures",
    checkFromB = cms.bool(True),
    drMax = cms.double(0.02),
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
    electrons = cms.InputTag("slimmedLowPtElectrons"),
    MVASeedUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
    MVASeedPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
    MVAIDLowPt = cms.InputTag('lowPtGsfElectronID'),
    #MVAIDV2 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values'),
    #convVtxFitProb = cms.InputTag('electronMVAVariableHelper:convVtxFitProb'),
    )
