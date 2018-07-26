/*
	Ugly copy-paste of RecoParticleFlow/PFTracking/plugins/GoodSeedProducer.cc
	to dump ALL the values used by the seeding to create a new training
	the code should be refactored to extract all the configurable parts
	(matching, extraction of features, etc.)
 */

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoParticleFlow/PFTracking/interface/PFGeometry.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"  
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include <fstream>
#include <string>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;

class TrackerElectronsFeatures: public edm::EDAnalyzer {
  typedef TrajectoryStateOnSurface TSOS;
public:
  explicit TrackerElectronsFeatures(const edm::ParameterSet&);
	~TrackerElectronsFeatures() {
	}

private:
	virtual void beginRun(const edm::Run & run,const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	// ----------member data ---------------------------

	//Features trk_features_;
	edm::Service<TFileService> fs_;
	TTree *tree_;	
	bool isMC_;
	double fake_prescale_;

	const edm::EDGetTokenT< vector<reco::PreId> > preid_;	
	const edm::EDGetTokenT< vector<reco::GsfTrack> > gsf_tracks_;
	const edm::EDGetTokenT< vector<reco::GsfElectron> > ged_electrons_;
	//MC Only
	const edm::EDGetTokenT< reco::RecoToSimCollection > association_;
};

TrackerElectronsFeatures::TrackerElectronsFeatures(const ParameterSet& cfg):
	isMC_{cfg.getParameter<bool>("isMC")},
	fake_prescale_{cfg.getParameter<double>("prescaleFakes")},
	preid_{consumes< vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
	gsf_tracks_   {consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
	ged_electrons_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gedElectrons"))},
	association_{consumes< reco::RecoToSimCollection >(cfg.getParameter<edm::InputTag>("association"))}
{
	tree_ = fs_->make<TTree>("tree", "test");
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackerElectronsFeatures::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerElectronsFeatures::beginRun(const edm::Run & run,
													 const EventSetup& es)
{
}

DEFINE_FWK_MODULE(TrackerElectronsFeatures);
