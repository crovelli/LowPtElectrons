import FWCore.ParameterSet.Config as cms

mvaVariablesFile = "LowPtElectrons/LowPtElectrons/data/ElectronIDVariables.txt"

ntuplizer = cms.EDAnalyzer('ElectronMVANtuplizer',
        # AOD case
        src                  = cms.InputTag('lowPtGsfElectronsOpen'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        genParticles         = cms.InputTag('genParticles'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedElectrons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        genParticlesMiniAOD  = cms.InputTag('prunedGenParticles'),
        #
        eleMVAs             = cms.untracked.vstring(),
        eleMVALabels        = cms.untracked.vstring(),
        eleMVAValMaps        = cms.untracked.vstring(),
        eleMVAValMapLabels   = cms.untracked.vstring(),
        eleMVACats           = cms.untracked.vstring(),
        eleMVACatLabels      = cms.untracked.vstring(),
        #
        variableDefinition   = cms.string(mvaVariablesFile),
        isMC                 = cms.bool(True),
        deltaR               = cms.double(0.1),
        )
