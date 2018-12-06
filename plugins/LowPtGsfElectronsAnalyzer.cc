#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>

using namespace std;
using namespace reco;

class LowPtGsfElectronsAnalyzer: public edm::EDAnalyzer {

public:
  
  explicit LowPtGsfElectronsAnalyzer( const edm::ParameterSet& );
  
  ~LowPtGsfElectronsAnalyzer() {
  
  }

private:

  virtual void analyze( const edm::Event&, const edm::EventSetup& );

  const edm::EDGetTokenT< reco::PFClusterCollection > ecal_clusters_;
  const edm::EDGetTokenT< reco::PFClusterCollection > hcal_clusters_;
  const edm::EDGetTokenT< reco::PFRecTrackCollection > pf_ktf_tracks_;
  const edm::EDGetTokenT< std::vector<reco::PreId> > preid_;	
  const edm::EDGetTokenT< edm::View<TrajectorySeed> > ele_seeds_;
  const edm::EDGetTokenT< std::vector<TrackCandidate> > trk_candidates_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsf_tracks_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > egamma_gsf_tracks_;
  const edm::EDGetTokenT< std::vector<reco::GsfPFRecTrack> > pf_gsf_tracks_;
  const edm::EDGetTokenT< std::vector<reco::CaloCluster> > electron_clu_;
  const edm::EDGetTokenT< edm::ValueMap<reco::SuperClusterRef> > electron_scref_;
  const edm::EDGetTokenT< std::vector<reco::SuperCluster> > electron_sc_;
  const edm::EDGetTokenT< std::vector<reco::GsfElectronCore> > electron_cores_;
  const edm::EDGetTokenT< std::vector<reco::GsfElectron> > electrons_;

};

LowPtGsfElectronsAnalyzer::LowPtGsfElectronsAnalyzer( const edm::ParameterSet& cfg ) :
  ecal_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("ECALClusters"))},
  hcal_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("HCALClusters"))},
  pf_ktf_tracks_{consumes< reco::PFRecTrackCollection >(cfg.getParameter<edm::InputTag>("PFTracks"))},
  preid_{consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preId"))},	
  ele_seeds_{consumes< edm::View<TrajectorySeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))},
  trk_candidates_{consumes< std::vector<TrackCandidate> >(cfg.getParameter<edm::InputTag>("trkCandidates"))},
  gsf_tracks_{consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
  egamma_gsf_tracks_{consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("EGammaGsfTracks"))}, 
  pf_gsf_tracks_{consumes< std::vector<reco::GsfPFRecTrack> >(cfg.getParameter<edm::InputTag>("PFGsfTracks"))},
  electron_clu_{consumes< std::vector<reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("electronCaloClusters"))},
  electron_scref_{consumes< edm::ValueMap<reco::SuperClusterRef> >(cfg.getParameter<edm::InputTag>("electronSCRefs"))},
  electron_sc_{consumes< std::vector<reco::SuperCluster> >(cfg.getParameter<edm::InputTag>("electronSCs"))},
  electron_cores_{consumes< std::vector<reco::GsfElectronCore> >(cfg.getParameter<edm::InputTag>("electronCores"))},
  electrons_{consumes< std::vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons"))}
{;}

void LowPtGsfElectronsAnalyzer::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFClusterCollection> ecal_clusters;
  try { iEvent.getByToken(ecal_clusters_, ecal_clusters); }
  catch (...) {;}
  
  edm::Handle<reco::PFClusterCollection> hcal_clusters;
  try { iEvent.getByToken(hcal_clusters_, hcal_clusters); }
  catch (...) {;}
  
  edm::Handle< reco::PFRecTrackCollection > pf_ktf_tracks;
  try { iEvent.getByToken(pf_ktf_tracks_, pf_ktf_tracks); }
  catch (...) {;}
  
  edm::Handle< std::vector<reco::PreId> > preids;
  try { iEvent.getByToken(preid_, preids); }
  catch (...) {;}

  edm::Handle< edm::View<TrajectorySeed> > ele_seeds;
  try { iEvent.getByToken(ele_seeds_, ele_seeds); }
  catch (...) {;}
  
  edm::Handle< std::vector<TrackCandidate> > trk_candidates;
  try { iEvent.getByToken(trk_candidates_, trk_candidates); }
  catch (...) {;}
    
  edm::Handle< std::vector<reco::GsfTrack> > gsf_tracks;
  try { iEvent.getByToken(gsf_tracks_, gsf_tracks); }
  catch (...) {}
  
  edm::Handle< std::vector<reco::GsfTrack> > egamma_gsf_tracks;
  try { iEvent.getByToken(egamma_gsf_tracks_, egamma_gsf_tracks); }
  catch (...) {;}
  
  edm::Handle< std::vector<reco::GsfPFRecTrack> > pf_gsf_tracks;
  try { iEvent.getByToken(pf_gsf_tracks_, pf_gsf_tracks); }
  catch (...) {;}
  
  //	int jtrk = 0;
  //	for ( auto trk : *(pf_gsf_tracks.product()) ) {
  //	  std::cout << "pf gsf pt " << itrk
  //		    << " " << trk.trackRef()->pt()
  //		    << std::endl;
  //	  ++jtrk;
  //	}
  
  edm::Handle< std::vector<reco::CaloCluster> > electron_clu;
  try { iEvent.getByToken(electron_clu_, electron_clu); }
  catch (...) {;}
  
  edm::Handle< edm::ValueMap<reco::SuperClusterRef> > electron_scref;
  try { iEvent.getByToken(electron_scref_, electron_scref); }
  catch (...) {;}
  
  edm::Handle< vector<reco::SuperCluster> > electron_sc;
  try { iEvent.getByToken(electron_sc_, electron_sc); }
  catch (...) {;}
  
  edm::Handle< std::vector<reco::GsfElectronCore> > electron_cores;
  try { iEvent.getByToken(electron_cores_, electron_cores); }
  catch (...) {;}

  edm::Handle< std::vector<reco::GsfElectron> > electrons;
  try { iEvent.getByToken(electrons_, electrons); }
  catch (...) {;}
    
  std::cout << "[LowPtGsfElectronsAnalyzer::analyze]" << std::endl
	    << "  ecal_clusters:     " << int( ecal_clusters.isValid() ? ecal_clusters->size() : -1 ) << std::endl
	    << "  hcal_clusters:     " << int( hcal_clusters.isValid() ? hcal_clusters->size() : -1 ) << std::endl
	    << "  pf_ktf_tracks:     " << int( pf_ktf_tracks.isValid() ? pf_ktf_tracks->size() : -1 ) << std::endl
	    << "  preids:            " << int( preids.isValid() ? preids->size() : -1 ) << std::endl
	    << "  ele_seeds:         " << int( ele_seeds.isValid() ? ele_seeds->size() : -1 ) << std::endl
	    << "  trk_candidates:    " << int( trk_candidates.isValid() ? trk_candidates->size() : -1 ) << std::endl
	    << "  gsf_tracks:        " << int( gsf_tracks.isValid() ? gsf_tracks->size() : -1 ) << std::endl
	    << "  egamma_gsf_tracks: " << int( egamma_gsf_tracks.isValid() ? egamma_gsf_tracks->size() : -1 ) << std::endl
	    << "  pf_gsf_tracks:     " << int( pf_gsf_tracks.isValid() ? pf_gsf_tracks->size() : -1 ) << std::endl
	    << "  electron_clu:      " << int( electron_clu.isValid() ? electron_clu->size() : -1 ) << std::endl
	    << "  electron_scref:    " << int( electron_scref.isValid() ? electron_scref->size() : -1 ) << std::endl
	    << "  electron_sc:       " << int( electron_sc.isValid() ? electron_sc->size() : -1 ) << std::endl
	    << "  electron_cores:    " << int( electron_cores.isValid() ? electron_cores->size() : -1 ) << std::endl
	    << "  electrons:         " << int( electrons.isValid() ? electrons->size() : -1 ) << std::endl
	    << std::endl;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LowPtGsfElectronsAnalyzer);
