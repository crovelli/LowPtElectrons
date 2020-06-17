import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP

ntuplizer = cms.EDAnalyzer(
    "IDSlimNtuplizer",
    verbose = cms.int32(0),
    checkFromB = cms.bool(True),
    prescale = cms.double(0.15),
#    prescale = cms.double(0.05),
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
    mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'), # wp80?
    lowPtRegressionConfig = cms.PSet(
      modifierName = cms.string('EGRegressionModifierLPV1'),       
      rhoTag = cms.string('fixedGridRhoFastjetAll'),
      useClosestToCentreSeedCrysDef = cms.bool(False),
      maxRawEnergyForLowPtEBSigma = cms.double(-1), 
      maxRawEnergyForLowPtEESigma = cms.double(1200.),
      eleRegs = cms.PSet(
        ecalOnlyMean = cms.PSet(
            rangeMinLowEt = cms.double(0.2),
            rangeMaxLowEt = cms.double(2.0),
            rangeMinHighEt = cms.double(-1.),
            rangeMaxHighEt = cms.double(3.0),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
            ),
        ecalOnlySigma = cms.PSet(
            rangeMinLowEt = cms.double(0.0002),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMaxHighEt = cms.double(0.5),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
            ),
        epComb = cms.PSet(
            ecalTrkRegressionConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(0.2),
                rangeMaxHighEt = cms.double(2.0),
                lowEtHighEtBoundary = cms.double(20.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_mean'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_mean'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_mean'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_mean'),
                ),
            ecalTrkRegressionUncertConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                lowEtHighEtBoundary = cms.double(20.),  
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_sigma'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_sigma'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_sigma'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_sigma'),
                ),
            maxEcalEnergyForComb=cms.double(200.),
            minEOverPForComb=cms.double(0.025),
            maxEPDiffInSigmaForComb=cms.double(15.),
            maxRelTrkMomErrForComb=cms.double(10.),                
            )
        ),
      phoRegs = regressionModifier106XUL.phoRegs.clone()
    ),
    gsfRegressionConfig = cms.PSet(
      modifierName = cms.string('EGRegressionModifierV3'),       
      rhoTag = cms.string('fixedGridRhoFastjetAll'),
      useClosestToCentreSeedCrysDef = cms.bool(False),
      maxRawEnergyForLowPtEBSigma = cms.double(-1), 
      maxRawEnergyForLowPtEESigma = cms.double(1200.),
      eleRegs = cms.PSet(
        ecalOnlyMean = cms.PSet(
            rangeMinLowEt = cms.double(0.2),
            rangeMaxLowEt = cms.double(2.0),
            rangeMinHighEt = cms.double(-1.),
            rangeMaxHighEt = cms.double(3.0),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(999999.),
            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
            ),
        ecalOnlySigma = cms.PSet(
            rangeMinLowEt = cms.double(0.0002),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMaxHighEt = cms.double(0.5),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(999999.),
            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
            ),
        epComb = cms.PSet(
            ecalTrkRegressionConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(0.2),
                rangeMaxHighEt = cms.double(2.0),
                lowEtHighEtBoundary = cms.double(999999.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
                ),
            ecalTrkRegressionUncertConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                lowEtHighEtBoundary = cms.double(999999.),  
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
                ),
            maxEcalEnergyForComb=cms.double(200.),
            minEOverPForComb=cms.double(0.025),
            maxEPDiffInSigmaForComb=cms.double(15.),
            maxRelTrkMomErrForComb=cms.double(10.),                
            )
        ),
      phoRegs = regressionModifier106XUL.phoRegs.clone()
    )
)
