import FWCore.ParameterSet.Config as cms

features = cms.EDAnalyzer("TrackerElectronsFeatures",
                          isMC = cms.bool(True),
                          prescaleFakes = cms.double(100),
                          beamspot = cms.InputTag("offlineBeamSpot"),
                          genParticles = cms.InputTag("genParticles"),
                          gsfTracks = cms.InputTag("electronGsfTracks"),
                          #gsfTracks = cms.InputTag("electronGsfTracks"),
                          gedElectrons = cms.InputTag("gedGsfElectrons"),
                          preId = cms.InputTag("trackerDrivenElectronSeeds","preid"),
                          association = cms.InputTag("trackingParticleRecoTrackAsssociation"),
                          trkCandidates = cms.InputTag("electronCkfTrackCandidates"),
                          #pileup = cms.InputTag("addPileupInfo"),
                          #vertices = cms.InputTag("offlinePrimaryVertices"),
                          #generalTracks = cms.InputTag("generalTracks"),
                          )
