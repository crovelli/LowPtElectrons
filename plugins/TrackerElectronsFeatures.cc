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
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementGsfTrack.h"

#include <algorithm>
#include <numeric>
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
	bool hit_association_;
	double dr_max_; //for DR Matching only
	double fake_prescale_;

	const edm::EDGetTokenT< vector<reco::PreId> > preid_;	
	const edm::EDGetTokenT< vector<reco::GsfTrack> > gsf_tracks_;
	const edm::EDGetTokenT< edm::View<reco::Track> > gsf_tracks_view_;
	const edm::EDGetTokenT< vector<reco::GsfPFRecTrack> > pf_gsf_tracks_;
	const edm::EDGetTokenT< vector<reco::PFBlock> > pfblocks_;
	const edm::EDGetTokenT< reco::PFCandidateCollection > pf_electrons_;
	const edm::EDGetTokenT< vector<reco::GsfElectronCore> > ged_electron_cores_;
	const edm::EDGetTokenT< vector<reco::GsfElectron> > ged_electrons_;
	//MC Only
	const edm::EDGetTokenT< reco::RecoToSimCollection > association_;
	const edm::EDGetTokenT< reco::GenParticleCollection > gen_particles_;
	const edm::EDGetTokenT< reco::BeamSpot > beamspot_;
	const edm::EDGetTokenT< vector<TrackCandidate> > trk_candidates_;
	const edm::EDGetTokenT< edm::View<TrajectorySeed> > ele_seeds_;
	const edm::EDGetTokenT< reco::TrackToTrackingParticleAssociator > associator_;
	const edm::EDGetTokenT< TrackingParticleCollection > tracking_particles_;
};

TrackerElectronsFeatures::TrackerElectronsFeatures(const ParameterSet& cfg):
	ntuple_{},
	isMC_{cfg.getParameter<bool>("isMC")},
	hit_association_{cfg.getParameter<bool>("hitAssociation")},
	dr_max_{cfg.getParameter<double>("drMax")},
	fake_prescale_{cfg.getParameter<double>("prescaleFakes")},
	preid_{consumes< vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
	gsf_tracks_   {consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
	gsf_tracks_view_{consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("gsfTracks"))},
	pf_gsf_tracks_{consumes< vector<reco::GsfPFRecTrack> >(cfg.getParameter<edm::InputTag>("PFGsfTracks"))},
	pfblocks_{consumes< vector<reco::PFBlock> >(cfg.getParameter<edm::InputTag>("PFBlocks"))},
	pf_electrons_{consumes<reco::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("PFElectrons"))},
	ged_electron_cores_{consumes< vector<reco::GsfElectronCore> >(cfg.getParameter<edm::InputTag>("gedElectronCores"))},
	ged_electrons_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gedElectrons"))},
	association_{consumes< reco::RecoToSimCollection >(cfg.getParameter<edm::InputTag>("association"))},
	gen_particles_{consumes< reco::GenParticleCollection >(cfg.getParameter<edm::InputTag>("genParticles"))},
	beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))},
	trk_candidates_{consumes< vector<TrackCandidate> >(cfg.getParameter<edm::InputTag>("trkCandidates"))},
	ele_seeds_{consumes< edm::View<TrajectorySeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))},
	associator_{consumes< reco::TrackToTrackingParticleAssociator >(cfg.getParameter<edm::InputTag>("associator"))},
	tracking_particles_{consumes< TrackingParticleCollection >(cfg.getParameter<edm::InputTag>("trackingParticles"))}
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

	edm::Handle< vector<reco::GsfPFRecTrack> > pf_gsf_tracks;
	iEvent.getByToken(pf_gsf_tracks_, pf_gsf_tracks);

	edm::Handle< vector<reco::PFBlock> > pfblocks;
	iEvent.getByToken(pfblocks_, pfblocks);

	edm::Handle< reco::PFCandidateCollection > pf_electrons;
	iEvent.getByToken(pf_electrons_, pf_electrons);

	edm::Handle< vector<reco::GsfElectronCore> > ged_electron_cores;
	iEvent.getByToken(ged_electron_cores_, ged_electron_cores);

	edm::Handle< vector<reco::GsfElectron> > ged_electrons;
	iEvent.getByToken(ged_electrons_, ged_electrons);

	edm::Handle<reco::RecoToSimCollection> matching;
	iEvent.getByToken(association_, matching);

	edm::Handle<vector<reco::GenParticle> > gen_particles;
	iEvent.getByToken(gen_particles_, gen_particles);

	edm::Handle< reco::BeamSpot > beamspot;
	iEvent.getByToken(beamspot_, beamspot);

	edm::Handle< vector<TrackCandidate> > trk_candidates;
	iEvent.getByToken(trk_candidates_, trk_candidates);

	edm::Handle< edm::View<TrajectorySeed> > ele_seeds;
	iEvent.getByToken(ele_seeds_, ele_seeds);
	
	edm::Handle< reco::TrackToTrackingParticleAssociator > associator;
	iEvent.getByToken(associator_, associator);
	
	edm::Handle< TrackingParticleCollection > tracking_particles;
	iEvent.getByToken(tracking_particles_, tracking_particles);

	int ntrks = 0;
	for ( unsigned int ii = 0; ii < preids->size(); ++ii ) {
	  if ( preids.product()->at(ii).trackRef().isNonnull() ) { ++ntrks; }
	}

	/*std::cout << "DEBUG"
		  << " preids: " << preids->size()
		  << " trk: " << ntrks
			<< " seeds: " << ele_seeds->size()
		  << " cands: " << trk_candidates->size()
		  << " gsf: " <<  gsf_tracks->size()
		  << std::endl;*/

	//assert(gsf_tracks->size() == preids->size()); //this is bound to fail, but better check

	//
	std::set<GsfTrackRef> pfGSFTrks_sources;
	for(const auto& pfgsf : *pf_gsf_tracks) {
		if(!pfgsf.gsfTrackRef().isNull()) pfGSFTrks_sources.insert(pfgsf.gsfTrackRef());
	}

	std::set<GsfTrackRef> pfBlocks_sources;
	std::set<GsfTrackRef> pfBlocksWSC_sources;
	std::set<GsfTrackRef> pfBlocksWECAL_sources;
	for(const auto& block : *pfblocks) {
		bool has_SC = false;
		bool has_ECAL = false;
		for(const auto& element : block.elements()) {
			has_SC |= (element.type() == PFBlockElement::Type::SC);
			has_ECAL |= (element.type() == PFBlockElement::Type::SC) ||
				(element.type() == PFBlockElement::Type::ECAL);
		}

		for(const auto& element : block.elements()) {
			if(element.type() != PFBlockElement::Type::GSF) continue;			
			const PFBlockElementGsfTrack& casted = dynamic_cast<const PFBlockElementGsfTrack&>(element);
			if(!casted.GsftrackRef().isNull()) {
				pfBlocks_sources.insert(casted.GsftrackRef());
				if(has_SC) pfBlocksWSC_sources.insert(casted.GsftrackRef());
				if(has_ECAL) pfBlocksWECAL_sources.insert(casted.GsftrackRef());
			}
		}
	}


	//pfgsfs
	std::set<GsfTrackRef> pfElectrons_sources;
	for(const auto& pf : *pf_electrons) {
		if(!pf.gsfTrackRef().isNull()) pfElectrons_sources.insert(pf.gsfTrackRef());
	}

	//gsf2core
	std::map<GsfTrackRef, GsfElectronCoreRef> gsf2core;
	for(size_t idx=0; idx < ged_electron_cores->size(); ++idx){
		reco::GsfElectronCoreRef core(ged_electron_cores, idx);
		GsfTrackRef trk = core->gsfTrack();
		if(gsf2core.find(trk) != gsf2core.end()) {
			std::cout << "THIS SHOULD NEVER HAPPEN! Multiple GSFElectronCores matched to the same GSFTrack?!" << std::endl;
		} else {
			gsf2core.insert(std::pair<reco::GsfTrackRef, reco::GsfElectronCoreRef>(trk, core));
		}
	}

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

	//match seed to GSF
	std::map<size_t, std::vector<GsfTrackRef> > seed2gsf;
	for(size_t idx=0; idx<gsf_tracks->size(); ++idx) {
		GsfTrackRef trk(gsf_tracks, idx);
		size_t seed_id = trk->seedRef().castTo<ElectronSeedRef>().index();
		if(seed2gsf.find(seed_id) == seed2gsf.end()) {			
			seed2gsf.insert(std::pair<size_t, std::vector<GsfTrackRef> >(seed_id, {trk}));
		} else {
			std::cout << "THIS SHOULD NEVER HAPPEN! Multiple tracks matched to the same seed." << std::endl;
			seed2gsf[seed_id].push_back(trk);
		}
	}

	//Match gen to seed
	//match seeds to tracking particles
	auto reco2sim = associator->associateRecoToSim(ele_seeds, tracking_particles);
	//auto sim2reco = associator->associateSimToReco();

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
	std::map<reco::GenParticleRef, size_t> gen2seed;
	std::vector<size_t> other_electrons;
	std::vector<size_t> other_tracks;
  if(hit_association_) { //maybe refactor into functions
		for(size_t idx=0; idx<ele_seeds->size(); ++idx) {
			RefToBase<TrajectorySeed> key(ele_seeds, idx);
			auto match = reco2sim.find(key);
		
			//check matching
			if(match != reco2sim.end()) { 
				auto tracking_particle = match->val.front().first;
				if(std::abs(tracking_particle->pdgId()) == 11) { //is an electron				
					if(tracking_particle->genParticles().size()) { //the electron is Gen-level
						auto genp = tracking_particle->genParticles()[0];
						if(electrons_from_B.find(genp) != electrons_from_B.end()) { //is coming from a B
							if(gen2seed.find(genp) == gen2seed.end()) { //is not filled yet
								//cout << "match found" << endl;
								gen2seed.insert(std::pair<reco::GenParticleRef, size_t>(genp, idx));
							} else { //something is wrong...
								cout << "THIS SHOULD NEVER HAPPEN! Multiple GSF tracks associated to the same GEN!" << endl;
							}
						} else { //is gen level, but not from B
							other_electrons.push_back(idx);
						}
					} else { //it is not a gen-level electron
						other_electrons.push_back(idx);
					}
				} //ends if(std::abs(tracking_particle->pdgId()) == 11)
			} else { //is not matched or not an electron
				other_tracks.push_back(idx);
			}
		} //end loop on GSF tracks	

	} else { //DR association
		set<size_t> matched;
		for(const auto& gen : electrons_from_B) {
			size_t best_idx = 666;
			double min_dr = dr_max_;
			for(size_t idx=0; idx<preids->size(); ++idx) {
				double dr = reco::deltaR(*(preids->at(idx).trackRef()), *gen);
				if(dr < min_dr) {
					min_dr = dr;
					best_idx = idx;
				}
			}
			
			if(min_dr < dr_max_) {
				gen2seed.insert(std::pair<reco::GenParticleRef, size_t>(gen, best_idx));
				matched.insert(best_idx);
			}
		}

		other_tracks.reserve(preids->size() - matched.size());
		for(size_t i=0; i<preids->size(); i++) {
			if(matched.find(i) == matched.end()) {
				//be extra safe and check hit matching, anyhow we have a LOT of tracks
				RefToBase<TrajectorySeed> key(ele_seeds, i);
				auto match = reco2sim.find(key);
		
				//check matching
				if(match == reco2sim.end() || std::abs(match->val.front().first->pdgId()) != 11) { 
					other_tracks.push_back(i);
				}
			}
		}
	}
	/*cout << "Found " << electrons_from_B.size() << " gen electrons, " 
			 << gen2seed.size() << " matched electron seeds, "
			 << other_electrons.size() << " non-gen electrons, " 
			 << other_tracks.size() << " other tracks" << endl; */

	//fill gen electron quantities
	for(auto gen : electrons_from_B) {
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_e();
		ntuple_.fill_gen(gen);
		const auto& match = gen2seed.find(gen);
		if(match != gen2seed.end()) { //matched to SEED
			const auto& gsf_match = seed2gsf.find(match->second);
			ntuple_.fill_preid(
				preids->at(match->second),
				*beamspot, (gsf_match != seed2gsf.end())
				);
			
			if(gsf_match != seed2gsf.end()) { //matched to GSF Track
				const double gpt = gen->pt();
				auto best = std::min_element(
					gsf_match->second.begin(),
					gsf_match->second.end(),
					[gpt] (const GsfTrackRef& a, const GsfTrackRef& b) {
						return std::abs(a->pt() - gpt)/gpt < std::abs(b->pt() - gpt)/gpt;
					});
				ntuple_.fill_gsf_trk(*best, *beamspot);
				ntuple_.has_ele_core(gsf2core.find(*best) != gsf2core.end());
				ntuple_.has_pfEgamma(pfElectrons_sources.find(*best) != pfElectrons_sources.end());
				ntuple_.has_pfGSFTrk(
					pfGSFTrks_sources.find(*best) != pfGSFTrks_sources.end()
					);
				ntuple_.has_pfBlock( 
					pfBlocks_sources.find(*best) != pfBlocks_sources.end()
					);
				ntuple_.has_pfBlock_with_SC(
					pfBlocksWSC_sources.find(*best) != pfBlocksWSC_sources.end()
          );
				ntuple_.has_pfBlock_with_ECAL(
					pfBlocksWECAL_sources.find(*best) !=pfBlocksWECAL_sources.end()
					);

				const auto& ele_match = gsf2ged.find(*best);
				if(ele_match != gsf2ged.end()) { //matched to GED Electron
					ntuple_.fill_ele(ele_match->second);
				} //matched to GED Electron
			}//matched to GSF Track
		}//matched to SEED
		tree_->Fill();
	}

	
	//fill other electron quantities
	for(size_t& seed_idx : other_electrons) {
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_e_not_matched();
		const auto& gsf_match = seed2gsf.find(seed_idx);
		ntuple_.fill_preid(
			preids->at(seed_idx),
			*beamspot,
			(gsf_match != seed2gsf.end()) ? gsf_match->second.size() : 0
			);

		if(gsf_match != seed2gsf.end()) { //matched to GSF Track
			ntuple_.fill_gsf_trk(gsf_match->second.at(0), *beamspot);
			ntuple_.has_ele_core(gsf2core.find(gsf_match->second.at(0)) != gsf2core.end());
			ntuple_.has_pfEgamma(pfElectrons_sources.find(gsf_match->second.at(0)) != pfElectrons_sources.end());
			ntuple_.has_pfGSFTrk(
				pfGSFTrks_sources.find(gsf_match->second.at(0)) != pfGSFTrks_sources.end()
				);
			ntuple_.has_pfBlock( 
				pfBlocks_sources.find(gsf_match->second.at(0)) != pfBlocks_sources.end()
				);
			ntuple_.has_pfBlock_with_SC(
				pfBlocksWSC_sources.find(gsf_match->second.at(0)) != pfBlocksWSC_sources.end()
				);
			ntuple_.has_pfBlock_with_ECAL(
				pfBlocksWECAL_sources.find(gsf_match->second.at(0)) !=pfBlocksWECAL_sources.end()
				);

			const auto& ele_match = gsf2ged.find(gsf_match->second.at(0));
			if(ele_match != gsf2ged.end()) { //matched to GED Electron
				ntuple_.fill_ele(ele_match->second);
			} //matched to GED Electron
		}//matched to GSF Track

		tree_->Fill();
	}


	//(prescaled) fill the backgrounds
	//fill other electron quantities
	for(size_t& seed_idx : other_tracks) {
		if(gRandom->Rndm() >= fake_prescale_) continue; //the smaller the few tracks are kept
		ntuple_.reset();
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_other();
		const auto& gsf_match = seed2gsf.find(seed_idx);
		ntuple_.fill_preid(
			preids->at(seed_idx),
			*beamspot,
			(gsf_match != seed2gsf.end()) ? gsf_match->second.size() : 0
			);

		if(gsf_match != seed2gsf.end()) { //matched to GSF Track
			ntuple_.fill_gsf_trk(gsf_match->second.at(0), *beamspot);
			ntuple_.has_ele_core(gsf2core.find(gsf_match->second.at(0)) != gsf2core.end());
			ntuple_.has_pfEgamma(pfElectrons_sources.find(gsf_match->second.at(0)) != pfElectrons_sources.end());
			ntuple_.has_pfGSFTrk(
				pfGSFTrks_sources.find(gsf_match->second.at(0)) != pfGSFTrks_sources.end()
				);
			ntuple_.has_pfBlock( 
				pfBlocks_sources.find(gsf_match->second.at(0)) != pfBlocks_sources.end()
				);
			ntuple_.has_pfBlock_with_SC(
				pfBlocksWSC_sources.find(gsf_match->second.at(0)) != pfBlocksWSC_sources.end()
				);
			ntuple_.has_pfBlock_with_ECAL(
				pfBlocksWECAL_sources.find(gsf_match->second.at(0)) !=pfBlocksWECAL_sources.end()
				);
			
			const auto& ele_match = gsf2ged.find(gsf_match->second.at(0));
			if(ele_match != gsf2ged.end()) { //matched to GED Electron
				ntuple_.fill_ele(ele_match->second);
			} //matched to GED Electron
		}//matched to GSF Track

		tree_->Fill();
	}//*/

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerElectronsFeatures::beginRun(const edm::Run & run,
													 const EventSetup& es)
{
}

DEFINE_FWK_MODULE(TrackerElectronsFeatures);
