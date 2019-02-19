/*
Ntuplizer for everything you need to know about tracker-driven electrons
 */

// system include files
#include <memory>
#include <utility>

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
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

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
	PFClusterRef closest_cluster(const reco::PFTrajectoryPoint& point, const edm::Handle<reco::PFClusterCollection>& clusters);

	virtual void beginRun(const edm::Run & run,const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() override {
		std::cout << "Number of electrons stored: " << npassed_[0] << std::endl
							<< "Number of other electrons: " << npassed_[1] << std::endl
							<< "Number of other tracks: " << npassed_[2] << std::endl;
	}

  class Bin {
  public:
    Bin( float logpt, float eta, float weight ) {
      logpt_ = logpt; eta_ = eta; weight_ = weight;
    }
    float logpt_;
    float eta_;
    float weight_;
  };

  typedef unsigned int uint;
  float weight( float logpt, float eta, float weight = -1. );
  std::pair<uint,uint> indices( float logpt, float eta );
  bool empty_weights();
  void print_weights();
	vector<float> get_iso_rings(const GsfElectronRef& ele, const reco::TrackRef& ele_trk, const vector<reco::TrackRef> &tracks);

  std::pair<float,float> printPfBlock( const reco::GenParticleRef gen,
				       const reco::PreIdRef preid,
				       const reco::PFBlockRef block,
				       const reco::GsfTrackRef gsf,
				       const GsfElectronRef* ele );

	// ----------member data ---------------------------

	//Features trk_features_;
	edm::Service<TFileService> fs_;
	TTree *tree_;	
	ElectronNtuple ntuple_;
	bool isMC_;
	bool printPfBlock_;
	bool disable_association_;
	bool check_from_B_;
	bool hit_association_;
	double dr_max_; //for DR Matching only
	double fake_multiplier_;
        std::string fakes_weights_;
        std::vector<float> bins_logpt_;
        std::vector<float> bins_eta_;
        std::vector< std::vector<float> > weights_;
	const edm::EDGetTokenT< double > rho_;	
	const edm::EDGetTokenT< vector<reco::PreId> > preid_;	
	const edm::EDGetTokenT< vector<reco::GsfTrack> > gsf_tracks_;
	const edm::EDGetTokenT< edm::View<reco::Track> > gsf_tracks_view_;
	const edm::EDGetTokenT< vector<int> > pf_gsf_flags_;
	const edm::EDGetTokenT< reco::PFRecTrackCollection > pf_ktf_tracks_;
	const edm::EDGetTokenT< vector<reco::GsfPFRecTrack> > pf_gsf_tracks_;
	const edm::EDGetTokenT< vector<reco::PFBlock> > pfblocks_;
        const edm::EDGetTokenT< reco::PFCandidateCollection > pf_electrons_;
        const edm::EDGetTokenT< vector<reco::SuperCluster> > ged_electron_sc_;
        const edm::EDGetTokenT< vector<reco::CaloCluster> > ged_electron_clu_;
        const edm::EDGetTokenT< edm::ValueMap<reco::SuperClusterRef> > ged_electron_scref_;
	const edm::EDGetTokenT< vector<reco::GsfElectronCore> > ged_electron_cores_;
	const edm::EDGetTokenT< vector<reco::GsfElectron> > ged_electrons_;
	const edm::EDGetTokenT< reco::PFClusterCollection > ecal_clusters_;
	const edm::EDGetTokenT< reco::PFClusterCollection > hcal_clusters_;
	const edm::EDGetTokenT< EcalRecHitCollection > eb_rechits_;
	const edm::EDGetTokenT< EcalRecHitCollection > ee_rechits_;
	const edm::EDGetTokenT< edm::ValueMap<float> > mvaid_v1_;
	const edm::EDGetTokenT< edm::ValueMap<float> > mvaid_v2_;
	const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;

	//MC Only
	const edm::EDGetTokenT< reco::GenParticleCollection > gen_particles_;
	const edm::EDGetTokenT< reco::BeamSpot > beamspot_;
	const edm::EDGetTokenT< vector<TrackCandidate> > trk_candidates_;
	const edm::EDGetTokenT< edm::View<TrajectorySeed> > ele_seeds_;
	const edm::EDGetTokenT< reco::TrackToTrackingParticleAssociator > associator_;
	const edm::EDGetTokenT< TrackingParticleCollection > tracking_particles_;
	unsigned long long int npassed_[3] = {0,0,0};
};

TrackerElectronsFeatures::TrackerElectronsFeatures(const ParameterSet& cfg):
  ntuple_{},
  isMC_{cfg.getParameter<bool>("isMC")},
  printPfBlock_{cfg.getParameter<bool>("printPfBlock")},
  disable_association_{cfg.getParameter<bool>("disableAssociation")},
  check_from_B_{cfg.getParameter<bool>("checkFromB")},
  hit_association_{cfg.getParameter<bool>("hitAssociation")},
  dr_max_{cfg.getParameter<double>("drMax")},
  fake_multiplier_{cfg.getParameter<double>("fakesMultiplier")},
  fakes_weights_{cfg.getParameter<std::string>("fakesWeights")},
  rho_{consumes< double >(cfg.getParameter<edm::InputTag>("rho"))},	
  preid_{consumes< vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
  gsf_tracks_   {consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
  gsf_tracks_view_{consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("gsfTracks"))},
  pf_gsf_flags_{consumes< vector<int> >(cfg.getParameter<edm::InputTag>("PFGsfFlags"))},
  pf_ktf_tracks_{consumes< reco::PFRecTrackCollection >(cfg.getParameter<edm::InputTag>("PFTracks"))},
  pf_gsf_tracks_{consumes< vector<reco::GsfPFRecTrack> >(cfg.getParameter<edm::InputTag>("PFGsfTracks"))},
  pfblocks_{consumes< vector<reco::PFBlock> >(cfg.getParameter<edm::InputTag>("PFBlocks"))},
  pf_electrons_{consumes<reco::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("PFElectrons"))},
  ged_electron_sc_{consumes< std::vector<reco::SuperCluster> >(cfg.getParameter<edm::InputTag>("gedElectronSCs"))},
  ged_electron_clu_{consumes< std::vector<reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("gedElectronCaloClusters"))},
  ged_electron_scref_{consumes< edm::ValueMap<reco::SuperClusterRef> >(cfg.getParameter<edm::InputTag>("gedElectronSCRefs"))},
  ged_electron_cores_{consumes< vector<reco::GsfElectronCore> >(cfg.getParameter<edm::InputTag>("gedElectronCores"))},
  ged_electrons_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gedElectrons"))},
  ecal_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("ECALClusters"))},
  hcal_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("HCALClusters"))},
  eb_rechits_{consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("EBRecHits"))},
  ee_rechits_{consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("EERecHits"))},
  //mvaid_v1_{consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDV1"))},
  //mvaid_v2_{consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDV2"))},
  //convVtxFitProb_{consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb"))},
  gen_particles_{consumes< reco::GenParticleCollection >(cfg.getParameter<edm::InputTag>("genParticles"))},
  beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))},
  trk_candidates_{consumes< vector<TrackCandidate> >(cfg.getParameter<edm::InputTag>("trkCandidates"))},
  ele_seeds_{consumes< edm::View<TrajectorySeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))},
  associator_{mayConsume< reco::TrackToTrackingParticleAssociator >(cfg.getParameter<edm::InputTag>("associator"))},
  tracking_particles_{consumes< TrackingParticleCollection >(cfg.getParameter<edm::InputTag>("trackingParticles"))}
{
	tree_ = fs_->make<TTree>("tree", "test");
	ntuple_.link_tree(tree_);
}


//
// member functions
//
vector<float> 
TrackerElectronsFeatures::get_iso_rings(const GsfElectronRef& ele, const reco::TrackRef& ele_trk, const vector<reco::TrackRef> &tracks) {
	vector<float> ret = {0., 0., 0., 0.};
	for(const auto& trk : tracks) {
		if(ele_trk == trk) continue;
	 	double dr = deltaR(*ele, *trk);
	 	size_t idx = std::floor(dr/0.1);
	 	if(idx < ret.size())
			ret[idx] += trk->pt(); //this is the conflicting line, for some reason
	}
	return ret;
}

PFClusterRef 
TrackerElectronsFeatures::closest_cluster(const reco::PFTrajectoryPoint& point, const edm::Handle<reco::PFClusterCollection>& clusters) {
	PFClusterRef best_ref;
	if(point.isValid()) { 
		float min_dr = 9999.f;
		for(size_t i=0; i<clusters->size(); i++) {
			float dr = reco::deltaR(clusters->at(i), point.positionREP());
			if(dr < min_dr) {
				best_ref = PFClusterRef(clusters, i);
				min_dr = dr;
			}
		}
	}
	return best_ref;
}

// ------------ method called to produce the data  ------------
void
TrackerElectronsFeatures::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	edm::Handle< double > rho;
	iEvent.getByToken(rho_, rho);

	edm::Handle< vector<reco::PreId> > preids;
	iEvent.getByToken(preid_, preids);

	edm::Handle< vector<reco::GsfTrack> > gsf_tracks;
	iEvent.getByToken(gsf_tracks_, gsf_tracks);

	edm::Handle< edm::View<reco::Track> > gsf_tracks_view;
	iEvent.getByToken(gsf_tracks_view_, gsf_tracks_view);

	edm::Handle< reco::PFRecTrackCollection > pf_ktf_tracks;
	iEvent.getByToken(pf_ktf_tracks_, pf_ktf_tracks);

	edm::Handle< vector<reco::GsfPFRecTrack> > pf_gsf_tracks;
	iEvent.getByToken(pf_gsf_tracks_, pf_gsf_tracks);

	edm::Handle< vector<int> > pf_gsf_flags;
	iEvent.getByToken(pf_gsf_flags_, pf_gsf_flags);

	edm::Handle< vector<reco::PFBlock> > pfblocks;
	iEvent.getByToken(pfblocks_, pfblocks);

	edm::Handle< reco::PFCandidateCollection > pf_electrons;
	iEvent.getByToken(pf_electrons_, pf_electrons);

	edm::Handle< vector<reco::SuperCluster> > ged_electron_sc;
	iEvent.getByToken(ged_electron_sc_, ged_electron_sc);

	edm::Handle< vector<reco::CaloCluster> > ged_electron_clu;
	iEvent.getByToken(ged_electron_clu_, ged_electron_clu);

	edm::Handle< edm::ValueMap<reco::SuperClusterRef> > ged_electron_scref;
	iEvent.getByToken(ged_electron_scref_, ged_electron_scref);

	edm::Handle< vector<reco::GsfElectronCore> > ged_electron_cores;
	iEvent.getByToken(ged_electron_cores_, ged_electron_cores);

	edm::Handle< vector<reco::GsfElectron> > ged_electrons;
	iEvent.getByToken(ged_electrons_, ged_electrons);

	edm::Handle<reco::PFClusterCollection> ecal_clusters;
	iEvent.getByToken(ecal_clusters_, ecal_clusters);

	edm::Handle<reco::PFClusterCollection> hcal_clusters;
	iEvent.getByToken(hcal_clusters_, hcal_clusters);

	noZS::EcalClusterLazyTools ecal_tools(iEvent, iSetup, eb_rechits_, ee_rechits_);

	edm::Handle<vector<reco::GenParticle> > gen_particles;
	iEvent.getByToken(gen_particles_, gen_particles);

	edm::Handle< reco::BeamSpot > beamspot;
	iEvent.getByToken(beamspot_, beamspot);

	edm::Handle< vector<TrackCandidate> > trk_candidates;
	iEvent.getByToken(trk_candidates_, trk_candidates);

	edm::Handle< edm::View<TrajectorySeed> > ele_seeds;
	iEvent.getByToken(ele_seeds_, ele_seeds);
	
	edm::Handle< TrackingParticleCollection > tracking_particles;
	iEvent.getByToken(tracking_particles_, tracking_particles);

	//edm::Handle< edm::ValueMap<float> > mvaid_v1;
	//iEvent.getByToken(mvaid_v1_, mvaid_v1);

	//edm::Handle< edm::ValueMap<float> > mvaid_v2;
	//iEvent.getByToken(mvaid_v2_, mvaid_v2);

	//edm::Handle< edm::ValueMap<float> > convVtxFitProb;
	//iEvent.getByToken(convVtxFitProb_, convVtxFitProb);

	int ntrks = 0;
	std::vector<reco::TrackRef> tracks;
	for ( unsigned int ii = 0; ii < preids->size(); ++ii ) {
	  if ( preids.product()->at(ii).trackRef().isNonnull() ) { 
			tracks.push_back(preids.product()->at(ii).trackRef());
			++ntrks; 
		}
	}

	//Match gen to seed
	//match seeds to tracking particles
	reco::RecoToSimCollectionSeed reco2sim;
	if(!disable_association_) {
		edm::Handle< reco::TrackToTrackingParticleAssociator > associator;
		iEvent.getByToken(associator_, associator);	
		reco2sim = associator->associateRecoToSim(ele_seeds, tracking_particles);
	}

	if (0) {
	  std::cout << "[TrackerElectronsFeatures::analyze]" << std::endl
		    << "  pf_ktf_tracks->size(): " << pf_ktf_tracks->size() << std::endl
		    << "  preids->size(): " << preids->size() << std::endl
		    << "  ele_seeds->size(): " << ele_seeds->size() << std::endl
		    << "  trk_candidates->size(): " << trk_candidates->size() << std::endl
		    << "  gsf_tracks->size(): " << gsf_tracks->size() << std::endl
		    << "  pf_gsf_tracks->size(): " << pf_gsf_tracks->size() << std::endl
		    << "  ged_electron_clu->size(): " << ged_electron_clu->size() << std::endl
		    << "  ged_electron_scref->size(): " << ged_electron_scref->size() << std::endl
		    << "  ged_electron_sc->size(): " << ged_electron_sc->size() << std::endl
		    << "  ged_electron_cores->size(): " << ged_electron_cores->size() << std::endl
		    << "  ged_electrons->size(): " << ged_electrons->size() << std::endl
		    << std::endl;
	}

	//assert(gsf_tracks->size() == preids->size()); //this is bound to fail, but better check


	std::map<reco::TrackRef, reco::PFRecTrackRef> trk2pftrk;
	for(size_t i=0; i<pf_ktf_tracks->size(); i++) {
		reco::PFRecTrackRef pftrk(pf_ktf_tracks, i);
		//it = find (myvector.begin(), myvector.end(), 30);
		//if (it != myvector.end())
		if(trk2pftrk.find(pftrk->trackRef()) != trk2pftrk.end()) assert(false);

		trk2pftrk.insert(pair<reco::TrackRef, reco::PFRecTrackRef>(pftrk->trackRef(), pftrk));
	}

	// PFTrack and Cluster association
	std::map<reco::TrackRef, PFClusterRef> ecal_ktf_clusters_map;
	std::map<reco::TrackRef, PFClusterRef> hcal_ktf_clusters_map;
	for(const auto& preid : *preids) {
		reco::TrackRef trk = preid.trackRef();
		reco::PFRecTrackRef pftrk;
		auto match = trk2pftrk.find(trk);
		if(match == trk2pftrk.end()) assert(false);
		pftrk = match->second;
		
		//get closest ECAL PF cluster
		auto position = pftrk->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
		PFClusterRef best_ref = closest_cluster(position, ecal_clusters);
		if(!best_ref.isNull()) {
			ecal_ktf_clusters_map.insert(
				std::pair<reco::TrackRef, PFClusterRef>(trk, best_ref)
				);
		}
	
		//get closest HCAL PF cluster
		auto hcal_position = pftrk->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
		PFClusterRef hcal_ref = closest_cluster(hcal_position, hcal_clusters);
		if(!hcal_ref.isNull()) {
			hcal_ktf_clusters_map.insert(
				std::pair<reco::TrackRef, PFClusterRef>(trk, hcal_ref)
				);
		}
	}

	// PFGSF and Cluster association
	std::map<GsfTrackRef, reco::GsfPFRecTrackRef> pfGSFTrks_sources;	
	std::map<GsfTrackRef, PFClusterRef> ecal_clusters_map;
	std::map<GsfTrackRef, PFClusterRef> hcal_clusters_map;
	for(size_t igsf=0; igsf<pf_gsf_tracks->size(); igsf++) {
		const auto& pfgsf = pf_gsf_tracks->at(igsf);
		reco::GsfPFRecTrackRef pfgsfRef(pf_gsf_tracks, igsf);
		GsfTrackRef gsfRef = pfgsf.gsfTrackRef();
		if(!pfgsf.gsfTrackRef().isNull()) pfGSFTrks_sources.insert(
			std::pair<GsfTrackRef, GsfPFRecTrackRef>(gsfRef, pfgsfRef)
			);
		
		//get closest ECAL PF cluster
		auto position = pfgsf.extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
		PFClusterRef best_ref = closest_cluster(position, ecal_clusters);
		if(!best_ref.isNull()) {
			ecal_clusters_map.insert(
				std::pair<GsfTrackRef, PFClusterRef>(gsfRef, best_ref)
				);
		}
		
		//get closest HCAL PF cluster
		auto hcal_position = pfgsf.extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
		PFClusterRef hcal_ref = closest_cluster(hcal_position, hcal_clusters);
		if(!hcal_ref.isNull()) {
			hcal_clusters_map.insert(
				std::pair<GsfTrackRef, PFClusterRef>(gsfRef, hcal_ref)
				);
		}
			
	}

	std::map<const GsfTrackRef, const PFBlockRef> pfBlocks_sources;
	std::set<GsfTrackRef> pfBlocksWSC_sources;
	std::set<GsfTrackRef> pfBlocksWECAL_sources;
	for ( unsigned int iblock = 0; iblock < pfblocks->size(); ++iblock ) {
          const PFBlockRef blockref(pfblocks,iblock);
	  bool has_SC = false;
	  bool has_ECAL = false;
		for(const auto& element : blockref->elements()) {
			has_SC |= (element.type() == PFBlockElement::Type::SC);
			has_ECAL |= (element.type() == PFBlockElement::Type::SC) ||
				(element.type() == PFBlockElement::Type::ECAL);
		}

		for(const auto& element : blockref->elements()) {
			if(element.type() != PFBlockElement::Type::GSF) continue;			
			const PFBlockElementGsfTrack& casted = dynamic_cast<const PFBlockElementGsfTrack&>(element);
			if(!casted.GsftrackRef().isNull()) {
			  if ( pfBlocks_sources.find(casted.GsftrackRef()) == pfBlocks_sources.end() ) {
			    pfBlocks_sources.insert( std::make_pair<const GsfTrackRef&,const PFBlockRef&>(casted.GsftrackRef(),blockref) );
			  } else { std::cout << "GSF track already in map!!! " << std::endl; }
			  if(has_SC) pfBlocksWSC_sources.insert(casted.GsftrackRef());
			  if(has_ECAL) pfBlocksWECAL_sources.insert(casted.GsftrackRef());
			}
		}
	}


	// pf candidate (electrons)
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
	//auto sim2reco = associator->associateSimToReco();

	//Find gen electrons
	std::set<reco::GenParticleRef> electrons_from_B;
	for(size_t idx=0; idx < gen_particles->size(); idx++) {
		reco::GenParticleRef genp(gen_particles, idx);
		bool is_ele = genp->isLastCopy() && std::abs(genp->pdgId()) == 11;
		bool comes_from_B = genp->numberOfMothers() >= 1 &&
		  std::abs(genp->mother()->pdgId()) > 510 &&
		  std::abs(genp->mother()->pdgId()) < 546;		
		if (!comes_from_B && // check resonant (J/psi) production?
		    genp->numberOfMothers() >= 1 && genp->mother() && // has mother
		    std::abs(genp->mother()->pdgId()) == 443 && // mother is J/psi
		    genp->mother()->numberOfMothers() >= 1 && genp->mother()->mother() && // has grandmother
		    std::abs(genp->mother()->mother()->pdgId()) > 510 &&
		    std::abs(genp->mother()->mother()->pdgId()) < 546 ) { // grandmother is B
		  comes_from_B = true;
		}
		if(is_ele && (comes_from_B || !check_from_B_)) { //is coming from a B
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
			        if ( ele_seeds->empty() ) { continue; }
				RefToBase<TrajectorySeed> key(ele_seeds, i);
				auto match = reco2sim.find(key);
		
				if(!disable_association_ && match != reco2sim.end() && std::abs(match->val.front().first->pdgId()) == 11)
					other_electrons.push_back(i);
				
				//check matching
				if(disable_association_ || match == reco2sim.end() || std::abs(match->val.front().first->pdgId()) != 11) { 
					other_tracks.push_back(i);
				}
			}
		}
	}
	npassed_[0] += electrons_from_B.size();
	npassed_[1] += other_electrons.size();
	npassed_[2] += other_tracks.size();
	/*cout << "Found " << electrons_from_B.size() << " gen electrons, " 
			 << gen2seed.size() << " matched electron seeds, "
			 << other_electrons.size() << " non-gen electrons, " 
			 << other_tracks.size() << " other tracks" << endl; */

	//fill gen electron quantities
  for(auto gen : electrons_from_B) {
    ntuple_.reset();
    ntuple_.set_rho(*rho);
    ntuple_.fill_evt(iEvent.id());
    ntuple_.is_e();
    ntuple_.fill_gen(gen);
    const auto& match = gen2seed.find(gen);
    if(match != gen2seed.end()) { //matched to SEED
      const auto& gsf_match = seed2gsf.find(match->second);
			const reco::PreIdRef preid(preids,match->second);
      const reco::TrackRef& ktf = preid->trackRef();
      ntuple_.fill_preid(
        *preid,
        *beamspot, (gsf_match != seed2gsf.end())
        );

      auto ktf_ecal = ecal_ktf_clusters_map.find(ktf);        
      reco::PFClusterRef ktf_ecal_cluster;
      if(ktf_ecal != ecal_ktf_clusters_map.end()) {
        ktf_ecal_cluster = ktf_ecal->second;
        ntuple_.fill_KTF_ECAL_cluster_info(
          ktf_ecal->second,
          trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax),
          ecal_tools
          );
      }

      auto ktf_hcal = hcal_ktf_clusters_map.find(ktf);
      reco::PFClusterRef ktf_hcal_cluster;
      if(ktf_hcal != hcal_ktf_clusters_map.end()) {
        ktf_hcal_cluster = ktf_hcal->second;
        ntuple_.fill_KTF_HCAL_cluster_info(
          ktf_hcal->second,
          trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
          );
      }

      
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
	//ntuple_.unpack_pfgsf_flags(pf_gsf_flags->at(best->index()));
        ntuple_.has_pfGSFTrk(
          pfGSFTrks_sources.find(*best) != pfGSFTrks_sources.end()
          );


        // cout << "DEBUG " << endl 
        //      << "#GSF " << gsf_tracks->size() << endl
        //      << "#PFGSF " << pf_gsf_tracks->size() << endl
        //      << "MATCHING " << (pfGSFTrks_sources.find(*best) != pfGSFTrks_sources.end()) << endl;
        auto cluster_matching = ecal_clusters_map.find(*best);        
        if(cluster_matching != ecal_clusters_map.end()) {
          auto point = pfGSFTrks_sources[*best]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
          ntuple_.fill_GSF_ECAL_cluster_info(
            cluster_matching->second,  point, ecal_tools
            );
          ntuple_.is_ECAL_cluster_same((ktf_ecal_cluster == cluster_matching->second));
        }

        auto hcal_matching = hcal_clusters_map.find(*best);
        if(hcal_matching != hcal_clusters_map.end()) {
          ntuple_.fill_GSF_HCAL_cluster_info(
            hcal_matching->second,
            pfGSFTrks_sources[*best]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
            );          
          ntuple_.is_HCAL_cluster_same((ktf_hcal_cluster == hcal_matching->second));
        }

				// PFBlock
				std::map<const GsfTrackRef,const PFBlockRef>::const_iterator iter = pfBlocks_sources.find(*best);
				ntuple_.has_pfBlock(iter != pfBlocks_sources.end());
				// PFBlock+SC
        ntuple_.has_pfBlock_with_SC(
          pfBlocksWSC_sources.find(*best) != pfBlocksWSC_sources.end()
          );
				// PFBlock+ECAL
        ntuple_.has_pfBlock_with_ECAL(
          pfBlocksWECAL_sources.find(*best) !=pfBlocksWECAL_sources.end()
          );

				const GsfElectronRef* ele_ref = 0;
        const auto& ele_match = gsf2ged.find(*best);
        if(ele_match != gsf2ged.end()) { //matched to GED Electron
          float id1 = 0.; // (*mvaid_v1)[ele_match->second];
          float id2 = 0.; // (*mvaid_v2)[ele_match->second];
          float conv_vtx_fit_prob = 0.; // (*convVtxFitProb)[ele_match->second];
					vector<float> iso_rings = get_iso_rings(ele_match->second, ktf, tracks);
          ntuple_.fill_ele(ele_match->second, id1, id2, conv_vtx_fit_prob, iso_rings);
					ele_ref = &(ele_match->second);
        } //matched to GED Electron

				std::pair<float,float> block(-1.,-1.);
				if ( iter != pfBlocks_sources.end() ) {
					block = printPfBlock(gen,
															 preid,
															 iter->second,
															 *best,
															 ele_ref);
				}
				ntuple_.has_pfBlock_size( block.first  );
				ntuple_.has_pfBlock_dr( block.second );

      }//matched to GSF Track
    }//matched to SEED
    tree_->Fill();
  }



	//fill other electron quantities
	for(size_t& seed_idx : other_electrons) {
		ntuple_.reset();
		ntuple_.set_rho(*rho);
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_e_not_matched();
		const auto& gsf_match = seed2gsf.find(seed_idx);
		const auto& preid = preids->at(seed_idx);
		const reco::TrackRef& ktf = preid.trackRef();
		ntuple_.fill_preid(
			preid,
			*beamspot,
			(gsf_match != seed2gsf.end()) ? gsf_match->second.size() : 0
			);
		
		auto ktf_ecal = ecal_ktf_clusters_map.find(ktf);				
		reco::PFClusterRef ktf_ecal_cluster;
		if(ktf_ecal != ecal_ktf_clusters_map.end()) {
			ntuple_.fill_KTF_ECAL_cluster_info(
				ktf_ecal->second,
				trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax),
				ecal_tools
				);
			ktf_ecal_cluster = ktf_ecal->second;
		}
		
		auto ktf_hcal = hcal_ktf_clusters_map.find(ktf);
		reco::PFClusterRef ktf_hcal_cluster;
		if(ktf_hcal != hcal_ktf_clusters_map.end()) {
			ntuple_.fill_KTF_HCAL_cluster_info(
				ktf_hcal->second,
				trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
				);
			ktf_hcal_cluster = ktf_hcal->second;
		}

		if(gsf_match != seed2gsf.end()) { //matched to GSF Track
			ntuple_.fill_gsf_trk(gsf_match->second.at(0), *beamspot);
			ntuple_.has_ele_core(gsf2core.find(gsf_match->second.at(0)) != gsf2core.end());
			ntuple_.has_pfEgamma(pfElectrons_sources.find(gsf_match->second.at(0)) != pfElectrons_sources.end());
			//ntuple_.unpack_pfgsf_flags(gsf_match->second.at(0).index());
			ntuple_.has_pfGSFTrk(
				pfGSFTrks_sources.find(gsf_match->second.at(0)) != pfGSFTrks_sources.end()
				);

			auto cluster_matching = ecal_clusters_map.find(gsf_match->second.at(0));				
			if(cluster_matching != ecal_clusters_map.end()) {
				ntuple_.fill_GSF_ECAL_cluster_info(
					cluster_matching->second,
					pfGSFTrks_sources[gsf_match->second.at(0)]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax),
					ecal_tools
					);
				ntuple_.is_ECAL_cluster_same((ktf_ecal_cluster == cluster_matching->second));
			}

			auto hcal_matching = hcal_clusters_map.find(gsf_match->second.at(0));
			if(hcal_matching != hcal_clusters_map.end()) {
				ntuple_.fill_GSF_HCAL_cluster_info(
					hcal_matching->second,
					pfGSFTrks_sources[gsf_match->second.at(0)]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
					);
				ntuple_.is_HCAL_cluster_same((ktf_hcal_cluster == hcal_matching->second));
			}

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
				float id1 = 0.; // (*mvaid_v1)[ele_match->second];
				float id2 = 0.; // (*mvaid_v2)[ele_match->second];
				float conv_vtx_fit_prob = 0.; // (*convVtxFitProb)[ele_match->second];
				vector<float> iso_rings = get_iso_rings(ele_match->second, ktf, tracks);
				ntuple_.fill_ele(ele_match->second, id1, id2, conv_vtx_fit_prob, iso_rings);
			} //matched to GED Electron
		}//matched to GSF Track

		tree_->Fill();
	}


	//@@
	// fill other electron quantities, prescaled according to (pT,eta) PDF and 'fakesMultiplier' configurable
	std::vector<int> trk_indices;
	unsigned int nfakes = 0;
	while ( trk_indices.size() < other_tracks.size() && // stop when all tracks considered
		( nfakes < fake_multiplier_ || fake_multiplier_ < 0 ) ) { // stop when fakesMultiplier is satisfied
	        // pick a track at random
                int index = int( gRandom->Rndm() * other_tracks.size() );
		// consider each track only once
		if ( std::find( trk_indices.begin(),
				trk_indices.end(),
				index ) != trk_indices.end() ) { continue; }
		trk_indices.push_back(index);
		// Obtain seed index from track
		size_t& seed_idx = other_tracks[index];
		// If map is populated, select tracks according to (pT,eta) PDF ...
		if ( !empty_weights() ) {
		  const reco::TrackRef& ktf = preids->at(seed_idx).trackRef();
		  if ( gRandom->Rndm() > weight(log10(ktf->pt()),fabs(ktf->eta())) ) { continue; }
		}
		nfakes++;

		ntuple_.reset();
		ntuple_.set_rho(*rho);
		ntuple_.fill_evt(iEvent.id());
		ntuple_.is_other();
		const auto& gsf_match = seed2gsf.find(seed_idx);
		const auto& preid = preids->at(seed_idx);
		const reco::TrackRef& ktf = preid.trackRef();
		ntuple_.fill_preid(
			preid,
			*beamspot,
			(gsf_match != seed2gsf.end()) ? gsf_match->second.size() : 0
			);
		
		auto ktf_ecal = ecal_ktf_clusters_map.find(ktf);				
		reco::PFClusterRef ktf_ecal_cluster;
		if(ktf_ecal != ecal_ktf_clusters_map.end()) {
			ntuple_.fill_KTF_ECAL_cluster_info(
				ktf_ecal->second,
				trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax),
				ecal_tools
				);
			ktf_ecal_cluster = ktf_ecal->second;
		}
		
		auto ktf_hcal = hcal_ktf_clusters_map.find(ktf);
		reco::PFClusterRef ktf_hcal_cluster;
		if(ktf_hcal != hcal_ktf_clusters_map.end()) {
			ntuple_.fill_KTF_HCAL_cluster_info(
				ktf_hcal->second,
				trk2pftrk[ktf]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
				);
			ktf_hcal_cluster = ktf_hcal->second;
		}

		if(gsf_match != seed2gsf.end()) { //matched to GSF Track
			ntuple_.fill_gsf_trk(gsf_match->second.at(0), *beamspot);
			ntuple_.has_ele_core(gsf2core.find(gsf_match->second.at(0)) != gsf2core.end());
			ntuple_.has_pfEgamma(pfElectrons_sources.find(gsf_match->second.at(0)) != pfElectrons_sources.end());
			//ntuple_.unpack_pfgsf_flags(gsf_match->second.at(0).index());
			ntuple_.has_pfGSFTrk(
				pfGSFTrks_sources.find(gsf_match->second.at(0)) != pfGSFTrks_sources.end()
				);

			auto cluster_matching = ecal_clusters_map.find(gsf_match->second.at(0));				
			if(cluster_matching != ecal_clusters_map.end()) {
				ntuple_.fill_GSF_ECAL_cluster_info(
					cluster_matching->second,
					pfGSFTrks_sources[gsf_match->second.at(0)]->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax),
					ecal_tools
					);
				ntuple_.is_ECAL_cluster_same((ktf_ecal_cluster == cluster_matching->second));
			}

			auto hcal_matching = hcal_clusters_map.find(gsf_match->second.at(0));
			if(hcal_matching != hcal_clusters_map.end()) {
				ntuple_.fill_GSF_HCAL_cluster_info(
					hcal_matching->second,
					pfGSFTrks_sources[gsf_match->second.at(0)]->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance)
					);
				ntuple_.is_HCAL_cluster_same((ktf_hcal_cluster == hcal_matching->second));
			}

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
				float id1 = 0.; // (*mvaid_v1)[ele_match->second];
				float id2 = 0.; // (*mvaid_v2)[ele_match->second];
				float conv_vtx_fit_prob = 0.; // (*convVtxFitProb)[ele_match->second];
				vector<float> iso_rings = get_iso_rings(ele_match->second, ktf, tracks);
				ntuple_.fill_ele(ele_match->second, id1, id2, conv_vtx_fit_prob, iso_rings);
			} //matched to GED Electron
		}//matched to GSF Track

		tree_->Fill();
	}//*/

}

/*******************************************************************************
 *
 */
std::pair<float,float> TrackerElectronsFeatures::printPfBlock( const reco::GenParticleRef gen,
							       const reco::PreIdRef preid,
							       const reco::PFBlockRef block,
							       const reco::GsfTrackRef gsf,
							       const reco::GsfElectronRef* ele )
{

  const edm::OwnVector<reco::PFBlockElement> elements = block->elements();
  float elements_size = elements.size();
  float elements_max_dr = -1.;
  for ( const auto& ielement : elements ) {
    if ( ielement.type() == PFBlockElement::TRACK ) {
      const PFBlockElementTrack& icast = static_cast<const PFBlockElementTrack&>(ielement);
      if ( icast.trackRef().isNonnull() ) {
	for ( const auto& jelement : elements ) {
	  if ( jelement.type() == PFBlockElement::TRACK ) {
	    const PFBlockElementTrack& jcast = static_cast<const PFBlockElementTrack&>(jelement);
	    if ( jcast.trackRef().isNonnull() ) {
	      float dr = deltaR( *(icast.trackRef()), *(jcast.trackRef()) );
	      if ( dr > elements_max_dr ) { elements_max_dr = dr; }
	    }
	  }
	}
      }
    }
  }
  if ( !printPfBlock_ ) { return std::pair<float,float>(elements_size,elements_max_dr); }

  static const std::string types[] = {"NONE",
				      "TRACK",
				      "PS1",
				      "PS2",
				      "ECAL",
				      "HCAL",
				      "GSF",
				      "BREM",
				      "HFEM",
				      "HFHAD",
				      "SC",
				      "HO",
				      "HGCAL",
				      "kNBETypes"};

  std::cout << "[LowPtEleNtuplizer::printPfBlock]" << std::endl
	    << std::setiosflags(std::ios::right)
	    << std::setiosflags(std::ios::fixed)
	    << std::setprecision(0)
	    << "    Index " << block.index()
	    << " #elements " << elements_size
	    << std::setprecision(2)
	    << " max_dr " << elements_max_dr
	    << std::endl
	    << "    ";
  std::vector<int> histo; histo.resize(14,0);
  for ( unsigned int iele = 0; iele < elements_size; ++iele ) { 
    int type = elements[iele].type() > 13 ? 13 : elements[iele].type();
    histo[type]++;
  }
  for ( unsigned int ii = 0 ; ii < histo.size(); ++ii ) {
    if ( histo[ii] > 0 ) { std::cout << types[ii] << "*" << histo[ii] << ", "; }
  }
  std::cout << std::endl;

  std::cout << "    GEN:"
	    << std::setprecision(2)
	    << std::setiosflags(std::ios::right)
	    << std::setiosflags(std::ios::fixed)
	    << " pt "  << std::setw(5) << gen->pt()
	    << " eta " << std::setw(5) << gen->eta()
	    << " phi " << std::setw(5) << gen->phi()
	    << std::endl;
  int best = -1;
  float min_dr = 1.e6;
  for ( unsigned int iele = 0; iele < elements_size; ++iele ) { 
    const reco::PFBlockElement& ele = elements[iele];
    if ( ele.type() == PFBlockElement::Type::GSF ) {
      const reco::PFBlockElementGsfTrack* trk = 
	static_cast<const reco::PFBlockElementGsfTrack*>(&ele);
      if ( trk->GsftrackRef().isNonnull() ) {
	const reco::GsfTrackRef& ref = trk->GsftrackRef();
	float dr = deltaR(*ref,*gen);
	if (dr < min_dr) {
	  best = iele;
	  min_dr = dr;
	}
      }
    }
  }

  if ( preid.isNonnull() ) {
    std::cout << "    TRK:";
    if ( preid->trackRef().isNonnull() ) {
      float dr  = deltaR(*(preid->trackRef()),*gen);
      std::cout << std::setprecision(2)
		<< std::setiosflags(std::ios::right)
		<< std::setiosflags(std::ios::fixed)
		<< " pt "  << std::setw(5) << preid->trackRef()->pt()
		<< " eta " << std::setw(5) << preid->trackRef()->eta()
		<< " phi " << std::setw(5) << preid->trackRef()->phi()
		<< ", dr " << std::setprecision(3) << dr
		<< std::endl;
    } else { std::cout << " (No track match)" << std::endl; }
  }

  if ( gsf.isNonnull() ) {
    std::cout << "    GSF:"
	      << std::setprecision(2)
	      << std::setiosflags(std::ios::right)
	      << std::setiosflags(std::ios::fixed)
	      << " pt "  << std::setw(5) << gsf->pt()
	      << " eta " << std::setw(5) << gsf->eta()
	      << " phi " << std::setw(5) << gsf->phi()
	      << ", dr " << std::setprecision(3) << min_dr
	      << " (LINKED," << std::setprecision(0) << std::setw(0) << gsf.index() << ")" 
	      << std::endl;
  }

  math::XYZVector outer(1.e6,1.e6,1.e6);
  if ( best >= 0 ) {
    const reco::PFBlockElementGsfTrack* element = 
      static_cast<const reco::PFBlockElementGsfTrack*>(&(elements[best]));
    const reco::GsfTrackRef& ref = element->GsftrackRef();
    std::cout << "    GSF:"
	      << std::setprecision(2)
	      << std::setiosflags(std::ios::right)
	      << std::setiosflags(std::ios::fixed)
	      << " pt "  << std::setw(5) << ref->pt()
	      << " eta " << std::setw(5) << ref->eta()
	      << " phi " << std::setw(5) << ref->phi()
	      << ", dr " << std::setprecision(3) << min_dr
	      << " (MIN DR," << std::setprecision(0) << std::setw(0) << ref.index() << ")" 
	      << std::endl;
    const math::XYZVector& inner = ref->innerMomentum();
    float dr_in  = deltaR(inner,*gen);
    std::cout << "    IN: "
	      << std::setprecision(2)
	      << std::setiosflags(std::ios::right)
	      << std::setiosflags(std::ios::fixed)
	      << " pt "  << std::setw(5) << sqrt(inner.Mag2()) * sin(inner.Theta())
	      << " eta " << std::setw(5) << inner.Eta()
	      << " phi " << std::setw(5) << inner.Phi()
	      << ", dr " << std::setprecision(3) << dr_in
	      << std::endl;
    outer = math::XYZVector(ref->outerMomentum());
    float dr_out = deltaR(outer,*gen);
    std::cout << "    OUT:"
	      << std::setprecision(2)
	      << std::setiosflags(std::ios::right)
	      << std::setiosflags(std::ios::fixed)
	      << " pt "  << std::setw(5) << ref->outerPt()
	      << " eta " << std::setw(5) << ref->outerEta()
	      << " phi " << std::setw(5) << ref->outerPhi()
	      << ", dr " << std::setprecision(3) << dr_out
	      << std::endl;
  }

  for ( unsigned int iele = 0; iele < elements_size; ++iele ) { 
    const reco::PFBlockElement& ele = elements[iele];
    if ( ele.type() == PFBlockElement::Type::ECAL ) {
      const reco::PFBlockElementCluster* ecal = 
	static_cast<const reco::PFBlockElementCluster*>(&ele);
      if ( ecal->clusterRef().isNonnull() ) {
	const reco::PFClusterRef& ref = ecal->clusterRef();
	float dr = 1.e6;
	if ( outer.X() < 1.e5 ) { dr = deltaR(ref->positionREP(),outer); }
	std::cout << "    clu:"
		  << std::setprecision(2)
		  << std::setiosflags(std::ios::right)
		  << std::setiosflags(std::ios::fixed)
		  << " et "  << std::setw(5) << ref->energy()
		  << " eta " << std::setw(5) << ref->positionREP().Eta()
		  << " phi " << std::setw(5) << ref->positionREP().Phi()
		  << ", dr " << std::setprecision(3) << ( dr < 1.e5 ? dr : -999. )
		  << " (w.r.t. OUT)"
		  << std::endl;
      }
    }
  }

  if ( ele != 0 ) { 
    const reco::GsfElectronRef ref = *ele;
    float dr = deltaR(*ref,*gen);
    std::cout << "    ELE:"
	      << std::setprecision(2)
	      << std::setiosflags(std::ios::right)
	      << std::setiosflags(std::ios::fixed)
	      << " pt "  << std::setw(5) << ref->p4().pt()
	      << " eta " << std::setw(5) << ref->p4().eta()
	      << " phi " << std::setw(5) << ref->p4().phi()
	      << ", dr " << std::setprecision(3) << dr
	      << " ELECTRON!"
	      << std::endl;
  } else { std::cout << "    ELE: (No match)" << std::endl; }
  std::cout << std::endl;

  return std::pair<float,float>(elements_size,elements_max_dr);

}

// Initialise the weights LUT to filter fake tracks
void TrackerElectronsFeatures::beginRun( const edm::Run & run,
					 const EventSetup& es ) {

  // Determine path to file containing weights used to filter fakes
  std::string path = "";
  try {
    path = edm::FileInPath(fakes_weights_).fullPath();
  }
  catch (...) {
    std::cout << "[TrackerElectronsFeatures::beginRun] "
	      << "Cannot find file: \"" << fakes_weights_ << "\" in path!"
	      << std::endl;
    return;
  }
  ifstream file(path);
  if ( !file.is_open() ) {
    std::cout << "[TrackerElectronsFeatures::beginRun] "
	      << "Cannot open file: \"" << path << "\"!"
	      << std::endl;
    return;
  } else {
    // Parse file and populate Bin container
    std::vector<Bin> bins;
    while( true ) {
      float pt, eta, weight;
      file >> pt >> eta >> weight;
      if (file.eof()) { break; }
      bins.push_back( Bin(pt,eta,weight) );
    }
    // Determine binning in logpt and eta
    std::set<float> bins_logpt;
    std::set<float> bins_eta;
    std::vector<Bin>::const_iterator iter;
    for ( iter = bins.begin(); iter != bins.end(); ++iter ) {
      bins_logpt.insert(iter->logpt_);
      bins_eta.insert(iter->eta_);
    }
    bins_logpt_.clear();
    bins_eta_.clear();
    std::copy( bins_logpt.begin(), bins_logpt.end(), std::back_inserter(bins_logpt_) );
    std::copy( bins_eta.begin(), bins_eta.end(), std::back_inserter(bins_eta_) );
    std::sort( bins_logpt_.begin(), bins_logpt_.begin() );
    std::sort( bins_eta_.begin(), bins_eta_.begin() );
    // Populate LUT containing weights used to filter fakes
    weights_.clear();
    for ( iter = bins.begin(); iter != bins.end(); ++iter ) {
      weight(iter->logpt_,iter->eta_,iter->weight_);
    }
  }
  // Print the LUT
  //print_weights();
}

// Set or get weight for given (logpt,eta) bin
float TrackerElectronsFeatures::weight( float logpt, float eta, float weight ) {
  // Get (ilogpt,ieta) indices
  std::pair<uint,uint> index_pair = indices(logpt,eta);
  // Check dimensions and their lengths, resize if needed
  if ( weights_.size() <= index_pair.first ) {
    weights_.resize(index_pair.first+1,std::vector<float>());
  }
  if ( weights_[index_pair.first].size() <= index_pair.second ) {
    weights_[index_pair.first].resize(index_pair.second+1,0.);
  }
  // Set weight if nonzero
  if ( weight > 0. ) {
    weights_[index_pair.first][index_pair.second] = weight;
//    std::cout << "Fill LUT ... " << weight << " "
//	      << weights_[index_pair.first][index_pair.second]
//	      << std::endl;
  }
  // Return weight
  return weights_[index_pair.first][index_pair.second];
}

// Convert (logpt,eta) floats into (ilogpt,ieta) indices
std::pair<uint,uint> TrackerElectronsFeatures::indices( float logpt, float eta ) {
  // Bisect the binning vectors to determine indices
  float offset = 1.e-6;
  int ilogpt = int( std::upper_bound(bins_logpt_.cbegin(),
				     bins_logpt_.cend(),
				     logpt+offset) - bins_logpt_.cbegin() );
  int ieta = int( std::upper_bound(bins_eta_.cbegin(),
				   bins_eta_.cend(),
				   eta+offset) - bins_eta_.cbegin() );
  // Checks that do nothing ...
  if ( ilogpt < 0 || ilogpt > (int)bins_logpt_.size() ) {;}
  if ( ieta < 0 || ieta > (int)bins_eta_.size() ) {;}
  // Limit the index ranges ...
  ilogpt = ilogpt > 0 ? ilogpt : 0;
  ieta = ieta > 0 ? ieta : 0;
  ilogpt = ilogpt < (int)bins_logpt_.size() ? ilogpt : (int)bins_logpt_.size()-1;
  ieta = ieta < (int)bins_eta_.size() ? ieta : (int)bins_eta_.size()-1;
  return std::pair<uint,uint>(ilogpt,ieta);
}

// Check if weights LUT is empty
bool TrackerElectronsFeatures::empty_weights() {
  return ( bins_logpt_.empty() ||
	   bins_eta_.empty() ||
	   weights_.empty() );
}

// Print (logpt,eta) binning schemes and weights LUT
void TrackerElectronsFeatures::print_weights() {
  if ( empty_weights() ) {
    std::cout << "[TrackerElectronsFeatures::print_weights] Empty vectors: "
	      << bins_logpt_.size() << " "
	      << bins_eta_.empty()  << " "
	      << weights_.size() << " "
	      << std::endl;
    return;
  }
  std::cout << std::setiosflags(std::ios::right)
	    << std::setiosflags(std::ios::fixed)
	    << std::setw(5)
	    << std::setprecision(2);
  std::cout << "log[pT], " << bins_logpt_.size() << " bins: ";
  for ( auto iter : bins_logpt_ ) { std::cout << iter << " "; }
  std::cout << std::endl;
  std::cout << "|eta|, " << bins_eta_.size() << " bins: ";
  for ( auto iter : bins_eta_ ) { std::cout << iter << " "; }
  std::cout << std::endl;
  std::cout << "weights, outer length: " << weights_.size() << " bins, inner lengths: ";
  for ( auto iter : weights_ ) { std::cout << iter.size() << " "; }
  std::cout << std::endl;
}

DEFINE_FWK_MODULE(TrackerElectronsFeatures);
