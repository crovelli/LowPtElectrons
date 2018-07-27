/*
Ntuplizer for everything you need to know about tracker-driven electrons
 */

// system include files
#include <memory>

// user include files

//TODO: cleanup all the includes not used

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

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include <fstream>
#include <string>
#include <map>
#include <set>
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
	const edm::EDGetTokenT< reco::GenParticleCollection > gen_particles_;
};

TrackerElectronsFeatures::TrackerElectronsFeatures(const ParameterSet& cfg):
	isMC_{cfg.getParameter<bool>("isMC")},
	fake_prescale_{cfg.getParameter<double>("prescaleFakes")},
	preid_{consumes< vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
	gsf_tracks_   {consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
	ged_electrons_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gedElectrons"))},
	association_{consumes< reco::RecoToSimCollection >(cfg.getParameter<edm::InputTag>("association"))},
	gen_particles_{consumes< reco::GenParticleCollection >(cfg.getParameter<edm::InputTag>("genParticles"))}
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
	edm::Handle< vector<reco::PreId> > preids;
	iEvent.getByToken(preid_, preids);

	edm::Handle< vector<reco::GsfTrack> > gsf_tracks;
	iEvent.getByToken(gsf_tracks_, gsf_tracks);

	edm::Handle< vector<reco::GsfElectron> > ged_electrons;
	iEvent.getByToken(ged_electrons_, ged_electrons);

	edm::Handle<reco::RecoToSimCollection> matching;
	iEvent.getByToken(association_, matching);

	edm::Handle<vector<reco::GenParticle> > gen_particles;
	iEvent.getByToken(gen_particles_, gen_particles);

	assert(gsf_tracks->size() == preids->size()); //this is bound to fail, but better check

	//gsf2ged
	atd::map<GsfTrackRef, GsfElectronRef> gsf2ged;
	for(size_t idx=0; idx < ged_electrons->size(); ++idx){
		reco::GsfElectronRef ele(ged_electrons, idx);
		GsfTrackRef trk = ele->gsfTrack();
		if(gsf2ged.find(trk) != gsf2ged.end()) {
			std::cout << "THIS SHOULD NEVER HAPPEN! Multiple GSFElectrons matched to the same GSFTrack?!" << std::endl;
		} else {
			gsf2ged.insert(std::pair<GsfTrackRef, GsfElectronRef>(GsfTrackRef, GsfElectronRef));
		}
	}

	//Match gen to reco
	std::set<reco::GenParticleRef> electrons_from_B;
	for(size_t idx=0; idx < gen_particles->size(); idx++) {
		reco::GenParticleRef genp(gen_particles, idx);
		if(genp->isLastCopy() && std::abs(genp->pdgId()) == 11 && 
			 genp->numberOfMothers() >= 1 && 
			 genp->mother()->pdgId() > 510 && genp->mother()->pdgId() < 546) { //is coming from a B
			electrons_from_B.insert(genp);
		}
	}

	//fill quantities
	
	
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerElectronsFeatures::beginRun(const edm::Run & run,
													 const EventSetup& es)
{
}

DEFINE_FWK_MODULE(TrackerElectronsFeatures);
