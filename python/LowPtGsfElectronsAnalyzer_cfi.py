import FWCore.ParameterSet.Config as cms

simple = cms.EDAnalyzer("LowPtGsfElectronsAnalyzer",
                        ECALClusters = cms.InputTag('particleFlowClusterECAL'),
                        HCALClusters = cms.InputTag('particleFlowClusterHCAL'),
                        PFTracks = cms.InputTag("lowPtGsfElePfTracks"),
                        preId = cms.InputTag("lowPtGsfElectronSeeds"),
                        eleSeeds = cms.InputTag("lowPtGsfElectronSeeds"),
                        trkCandidates = cms.InputTag("lowPtGsfEleCkfTrackCandidates"),
                        gsfTracks = cms.InputTag("lowPtGsfEleGsfTracks"),
                        EGammaGsfTracks = cms.InputTag("electronGsfTracks"),
                        PFGsfTracks = cms.InputTag("lowPtGsfElePfGsfTracks"),
                        electronCaloClusters = cms.InputTag("lowPtGsfElectronSuperClusters"),
                        electronSCs = cms.InputTag("lowPtGsfElectronSuperClusters"),
                        electronSCRefs = cms.InputTag("lowPtGsfElectronSuperClusters"),
                        electronCores = cms.InputTag("lowPtGsfElectronCores"),
                        electrons = cms.InputTag("lowPtGsfElectrons"),
                        )
