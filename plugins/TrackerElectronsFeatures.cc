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
#include "LowPtElectrons/LowPtElectrons/interface/ElectronNtuple.h"

#include <fstream>
#include <string>
#include <map>
#include <set>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <assert.h>

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
	ElectronNtuple ntuple_;
	bool isMC_;
	double fake_prescale_;

	const edm::EDGetTokenT< vector<reco::PreId> > preid_;	
	const edm::EDGetTokenT< vector<reco::GsfTrack> > gsf_tracks_;
	const edm::EDGetTokenT< edm::View<reco::Track> > gsf_tracks_view_;
	const edm::EDGetTokenT< vector<reco::GsfElectron> > ged_electrons_;
	//MC Only
	const edm::EDGetTokenT< reco::RecoToSimCollection > association_;
	const edm::EDGetTokenT< reco::GenParticleCollection > gen_particles_;
	const edm::EDGetTokenT< reco::BeamSpot > beamspot_;
};

TrackerElectronsFeatures::TrackerElectronsFeatures(const ParameterSet& cfg):
	ntuple_{},
	isMC_{cfg.getParameter<bool>("isMC")},
	fake_prescale_{cfg.getParameter<double>("prescaleFakes")},
	preid_{consumes< vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
	gsf_tracks_   {consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
	gsf_tracks_view_{consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("gsfTracks"))},
	ged_electrons_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gedElectrons"))},
	association_{consumes< reco::RecoToSimCollection >(cfg.getParameter<edm::InputTag>("association"))},
	gen_particles_{consumes< reco::GenParticleCollection >(cfg.getParameter<edm::InputTag>("genParticles"))},
	beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))}
{
	tree_ = fs_->make<TTree>("tree", "test");
	ntuple_.link_tree(tree_);
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

	edm::Handle< edm::View<reco::Track> > gsf_tracks_view;
	iEvent.getByToken(gsf_tracks_view_, gsf_tracks_view);

	edm::Handle< vector<reco::GsfElectron> > ged_electrons;
	iEvent.getByToken(ged_electrons_, ged_electrons);

	edm::Handle<reco::RecoToSimCollection> matching;
	iEvent.getByToken(association_, matching);

	edm::Handle<vector<reco::GenParticle> > gen_particles;
	iEvent.getByToken(gen_particles_, gen_particles);

	edm::Handle< reco::BeamSpot > beamspot;
	iEvent.getByToken(beamspot_, beamspot);

	assert(gsf_tracks->size() == preids->size()); //this is bound to fail, but better check

	//gsf2ged
	std::map<GsfTrackRef, GsfElectronRef> gsf2ged;
	for(size_t idx=0; idx < ged_electrons->size(); ++idx){
		reco::GsfElectronRef ele(ged_electrons, idx);
		GsfTrackRef trk = ele->gsfTrack();
		if(gsf2ged.find(trk) != gsf2ged.end()) {
			std::cout << "THIS SHOULD NEVER HAPPEN! Multiple GSFElectrons matched to the same GSFTrack?!" << std::endl;
		} else {
			gsf2ged.insert(std::pair<reco::GsfTrackRef, reco::GsfElectronRef>(trk, ele));
		}
	}

	//Match gen to reco
	//Find gen electrons
	std::set<reco::GenParticleRef> electrons_from_B;
	for(size_t idx=0; idx < gen_particles->size(); idx++) {
		reco::GenParticleRef genp(gen_particles, idx);
		if(genp->isLastCopy() && std::abs(genp->pdgId()) == 11 && 
			 genp->numberOfMothers() >= 1 && 
			 genp->mother()->pdgId() > 510 && genp->mother()->pdgId() < 546) { //is coming from a B
			electrons_from_B.insert(genp);
		}
	}

	//Match GEN to GSF
	std::map<reco::GenParticleRef, reco::GsfTrackRef> gen2gsf;
	std::vector<reco::GsfTrackRef> other_electrons;
	std::vector<reco::GsfTrackRef> other_tracks;
	for(size_t idx=0; idx<gsf_tracks->size(); ++idx) {
		RefToBase<reco::Track> key(gsf_tracks_view, idx);
		GsfTrackRef trk(gsf_tracks, idx);
		auto match = matching->find(key);
		
		//check matching
		if(match != matching->end()) { 
			auto tracking_particle = match->val.front().first;
			if(std::abs(tracking_particle->pdgId()) == 11) { //is an electron				
				if(tracking_particle->genParticles().size()) { //the electron is Gen-level
					auto genp = tracking_particle->genParticles()[0];
					if(electrons_from_B.find(genp) != electrons_from_B.end()) { //is coming from a B
						if(gen2gsf.find(genp) != gen2gsf.end()) { //is not filled yet
							gen2gsf.insert(std::pair<reco::GenParticleRef, reco::GsfTrackRef>(genp, trk));
						} else { //something is wrong...
							cout << "THIS SHOULD NEVER HAPPEN! Multiple GSF tracks associated to the same GEN!" << endl;
						}
					} else { //is gen level, but not from B
						other_electrons.push_back(trk);
					}
				} else { //it is not a gen-level electron
						other_electrons.push_back(trk);
				}
			} //ends if(std::abs(tracking_particle->pdgId()) == 11)
		} else { //is not matched or not an electron
			other_tracks.push_back(trk);
		}
	} //end loop on GSF tracks	

	//fill gen electron quantities
	for(auto gen : electrons_from_B) {
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.fill_gen(gen);
		ntuple_.is_e();
		const auto& match = gen2gsf.find(gen);
		if(match != gen2gsf.end()) { //matched to GSF TRK
			ntuple_.fill_gsf_trk(match->second, *beamspot);
			ntuple_.fill_preid(
				preids->at(match->second.index()),
				*beamspot
				);
			const auto& ele_match = gsf2ged.find(match->second);
			if(ele_match != gsf2ged.end()) { //matched to GED Electron
				ntuple_.fill_ele(ele_match->second);
			}
		}
	}
	
	//fill other electron quantities
	for(auto& ele : other_electrons) {
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_e_not_matched();
		ntuple_.fill_gsf_trk(ele, *beamspot);
		ntuple_.fill_preid(
			preids->at(ele.index()),
			*beamspot
			);
		const auto& ele_match = gsf2ged.find(ele);
		if(ele_match != gsf2ged.end()) { //matched to GED Electron
			ntuple_.fill_ele(ele_match->second);
		}
	}

	//(prescaled) fill the backgrounds
	//fill other electron quantities
	for(auto& ele : other_tracks) {
		if(gRandom->Rndm() >= fake_prescale_) continue; //the smaller the few tracks are kept
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_e_not_matched();
		ntuple_.fill_gsf_trk(ele, *beamspot);
		ntuple_.fill_preid(
			preids->at(ele.index()),
			*beamspot
			);
		const auto& ele_match = gsf2ged.find(ele);
		if(ele_match != gsf2ged.end()) { //matched to GED Electron
			ntuple_.fill_ele(ele_match->second);
		}
	}
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerElectronsFeatures::beginRun(const edm::Run & run,
													 const EventSetup& es)
{
}

DEFINE_FWK_MODULE(TrackerElectronsFeatures);
