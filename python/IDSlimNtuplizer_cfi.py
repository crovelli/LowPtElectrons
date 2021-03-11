import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL

ntuplizer = cms.EDAnalyzer(
    "IDSlimNtuplizer",
    verbose = cms.int32(0),
    checkFromB = cms.bool(True),
######    prescale = cms.double(0.15),            # JPsi usual
    prescale = cms.double(0.30),      
#    prescale = cms.double(0.05),            # DoubleEle
    minTrackPt = cms.double(0.5),  
    maxTrackPt = cms.double(15.),  
    maxTrackEta = cms.double(2.4),  
    # Generic collections
    rho = cms.InputTag('fixedGridRhoFastjetAll'),
    genParticles = cms.InputTag("genParticles"), # AOD
    prunedGenParticles = cms.InputTag("prunedGenParticles"), # MINIAOD
    patMuonCollection = cms.InputTag("slimmedMuons"), # MINIAOD  
    recoMuonCollection = cms.InputTag("muons"), # AOD
    # Low pT collections
    gsfElectrons = cms.InputTag("lowPtGsfElectrons"), # AOD 
    patElectrons = cms.InputTag("slimmedLowPtElectrons"), # MINIAOD
    mvaUnbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
    mvaPtbiased = cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased"),
    mvaValueLowPt = cms.InputTag('lowPtGsfElectronID'),
    dEdx1Tag = cms.InputTag('dedxHarmonic2'), 

    lowPtRegressionConfig = cms.PSet(
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
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_mean"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_mean"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_mean"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_mean"),
            ),
        ecalOnlySigma = cms.PSet(
            rangeMinLowEt = cms.double(0.0002),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMaxHighEt = cms.double(0.5),
            forceHighEnergyTrainingIfSaturated = cms.bool(True),
            lowEtHighEtBoundary = cms.double(20.),
            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_sigma"),
            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_sigma"),
            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_sigma"),
            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_sigma"),
            ),
        epComb = cms.PSet(
            ecalTrkRegressionConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.2),
                rangeMaxLowEt = cms.double(2.0),
                rangeMinHighEt = cms.double(0.2),
                rangeMaxHighEt = cms.double(2.0),
                lowEtHighEtBoundary = cms.double(20.),
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_mean'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_mean'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_mean'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_mean'),
                ),
            ecalTrkRegressionUncertConfig = cms.PSet(
                rangeMinLowEt = cms.double(0.0002),
                rangeMaxLowEt = cms.double(0.5),
                rangeMinHighEt = cms.double(0.0002),
                rangeMaxHighEt = cms.double(0.5),
                lowEtHighEtBoundary = cms.double(20.),  
                forceHighEnergyTrainingIfSaturated = cms.bool(False),
                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_sigma'),
                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_sigma'),
                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_sigma'),
                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_sigma'),
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
