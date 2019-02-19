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
#include <string>
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

  const edm::EDGetTokenT< reco::PFClusterCollection > ecal_pf_clusters_;
  const edm::EDGetTokenT< reco::PFClusterCollection > hcal_pf_clusters_;
  const edm::EDGetTokenT<EcalRecHitCollection> eb_rechits_;
  const edm::EDGetTokenT<EcalRecHitCollection> ee_rechits_;
  const edm::EDGetTokenT< reco::PFRecTrackCollection > pf_ktf_tracks_;
  const edm::EDGetTokenT< std::vector<reco::PreId> > preids_;
  const edm::EDGetTokenT< edm::ValueMap<reco::PreIdRef> > preids_valuemap_;
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
  //const edm::EDGetTokenT< edm::ValueMap<double> > mvaId_;
  const std::vector<edm::InputTag> mvaIdTags_;
  std::vector< edm::EDGetTokenT< edm::ValueMap<float> > >  mvaIds_;
  const std::vector<edm::InputTag> mvaSeedTags_;
  std::vector< edm::EDGetTokenT< edm::ValueMap<float> > >  mvaSeeds_;
};

LowPtGsfElectronsAnalyzer::LowPtGsfElectronsAnalyzer( const edm::ParameterSet& cfg ) :
  ecal_pf_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("ECALPFClusters"))},
  hcal_pf_clusters_{consumes<reco::PFClusterCollection>(cfg.getParameter<edm::InputTag>("HCALPFClusters"))},
  eb_rechits_{consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("EBRecHits"))},
  ee_rechits_{consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("EERecHits"))},
  pf_ktf_tracks_{consumes< reco::PFRecTrackCollection >(cfg.getParameter<edm::InputTag>("PFTracks"))},
  preids_{consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIds"))},
  preids_valuemap_{consumes< edm::ValueMap<reco::PreIdRef> >(cfg.getParameter<edm::InputTag>("preIdsValueMap"))},
  ele_seeds_{consumes< edm::View<TrajectorySeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))},
  trk_candidates_{consumes< std::vector<TrackCandidate> >(cfg.getParameter<edm::InputTag>("trkCandidates"))},
  gsf_tracks_{consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))}, 
  egamma_gsf_tracks_{consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("EGammaGsfTracks"))}, 
  pf_gsf_tracks_{consumes< std::vector<reco::GsfPFRecTrack> >(cfg.getParameter<edm::InputTag>("PFGsfTracks"))},
  electron_clu_{consumes< std::vector<reco::CaloCluster> >(cfg.getParameter<edm::InputTag>("electronCaloClusters"))},
  electron_scref_{consumes< edm::ValueMap<reco::SuperClusterRef> >(cfg.getParameter<edm::InputTag>("electronSCRefs"))},
  electron_sc_{consumes< std::vector<reco::SuperCluster> >(cfg.getParameter<edm::InputTag>("electronSCs"))},
  electron_cores_{consumes< std::vector<reco::GsfElectronCore> >(cfg.getParameter<edm::InputTag>("electronCores"))},
  electrons_{consumes< std::vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons"))},
  //mvaId_{consumes< edm::ValueMap<double> >(cfg.getParameter<edm::InputTag>("mvaId"))},
  mvaIdTags_(cfg.getParameter< std::vector<edm::InputTag> >("mvaIds")),
  mvaIds_(),
  mvaSeedTags_(cfg.getParameter< std::vector<edm::InputTag> >("mvaSeeds")),
  mvaSeeds_()
{
  for ( const auto& tag : mvaIdTags_ ) { 
    mvaIds_.push_back( consumes< edm::ValueMap<float> >(tag) ); 
  }
  for ( const auto& tag : mvaSeedTags_ ) { 
    mvaSeeds_.push_back( consumes< edm::ValueMap<float> >(tag) ); 
  }
}

void LowPtGsfElectronsAnalyzer::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFClusterCollection> ecal_pf_clusters;
  try { iEvent.getByToken(ecal_pf_clusters_, ecal_pf_clusters); }
  catch (...) {;}
  
  edm::Handle<reco::PFClusterCollection> hcal_pf_clusters;
  try { iEvent.getByToken(hcal_pf_clusters_, hcal_pf_clusters); }
  catch (...) {;}
  
  edm::Handle<EcalRecHitCollection> eb_rechits;
  try { iEvent.getByToken(eb_rechits_, eb_rechits); }
  catch (...) {;}
  
  edm::Handle<EcalRecHitCollection> ee_rechits;
  try { iEvent.getByToken(ee_rechits_, ee_rechits); }
  catch (...) {;}
  
  edm::Handle< reco::PFRecTrackCollection > pf_ktf_tracks;
  try { iEvent.getByToken(pf_ktf_tracks_, pf_ktf_tracks); }
  catch (...) {;}
  
  edm::Handle< std::vector<reco::PreId> > preids;
  try { iEvent.getByToken(preids_, preids); }
  catch (...) {;}

  edm::Handle< edm::ValueMap<reco::PreIdRef> > preids_valuemap;
  try { iEvent.getByToken(preids_valuemap_, preids_valuemap); }
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

//  edm::Handle< edm::ValueMap<float> > mvaId;
//  try { iEvent.getByToken(mvaId_, mvaId); }
//  catch (...) {;}
    
  std::vector< edm::Handle< edm::ValueMap<float> > > mvaIds;
  for ( const auto& token : mvaIds_ ) { 
    edm::Handle< edm::ValueMap<float> > h;
    try { iEvent.getByToken(token, h); }
    catch (...) {;}
    mvaIds.push_back(h);
  }

  std::vector< edm::Handle< edm::ValueMap<float> > > mvaSeeds;
  for ( const auto& token : mvaSeeds_ ) { 
    edm::Handle< edm::ValueMap<float> > h;
    try { iEvent.getByToken(token, h); }
    catch (...) {;}
    mvaSeeds.push_back(h);
  }
    
  std::cout << "[LowPtGsfElectronsAnalyzer::analyze]" << std::endl
	    << "  ecal_pf_clusters:  " << int( ecal_pf_clusters.isValid() ? ecal_pf_clusters->size() : -1 ) << std::endl
	    << "  hcal_pf_clusters:  " << int( hcal_pf_clusters.isValid() ? hcal_pf_clusters->size() : -1 ) << std::endl
	    << "  eb_rechits:        " << int( eb_rechits.isValid() ? eb_rechits->size() : -1 ) << std::endl
	    << "  ee_rechits:        " << int( ee_rechits.isValid() ? ee_rechits->size() : -1 ) << std::endl
	    << "  pf_ktf_tracks:     " << int( pf_ktf_tracks.isValid() ? pf_ktf_tracks->size() : -1 ) << std::endl
	    << "  preids:            " << int( preids.isValid() ? preids->size() : -1 ) << std::endl
	    << "  preids_valuemap:   " << int( preids_valuemap.isValid() ? preids_valuemap->size() : -1 ) << std::endl
	    << "  ele_seeds:         " << int( ele_seeds.isValid() ? ele_seeds->size() : -1 ) << std::endl
	    << "  trk_candidates:    " << int( trk_candidates.isValid() ? trk_candidates->size() : -1 ) << std::endl
	    << "  gsf_tracks:        " << int( gsf_tracks.isValid() ? gsf_tracks->size() : -1 ) << std::endl
	    << "  egamma_gsf_tracks: " << int( egamma_gsf_tracks.isValid() ? egamma_gsf_tracks->size() : -1 ) << std::endl
	    << "  pf_gsf_tracks:     " << int( pf_gsf_tracks.isValid() ? pf_gsf_tracks->size() : -1 ) << std::endl
	    << "  electron_clu:      " << int( electron_clu.isValid() ? electron_clu->size() : -1 ) << std::endl
	    << "  electron_scref:    " << int( electron_scref.isValid() ? electron_scref->size() : -1 ) << std::endl
	    << "  electron_sc:       " << int( electron_sc.isValid() ? electron_sc->size() : -1 ) << std::endl
	    << "  electron_cores:    " << int( electron_cores.isValid() ? electron_cores->size() : -1 ) << std::endl
	    << "  electrons:         " << int( electrons.isValid() ? electrons->size() : -1 ) << std::endl;
  //<< "  mvaId:             " << int( mvaId.isValid() ? mvaId->size() : -1 ) << std::endl;
  for ( unsigned int iter = 0; iter < mvaIds.size(); ++iter ) {
    std::cout << "  mvaId:             " 
	      << int( mvaIds[iter].isValid() ? mvaIds[iter]->size() : -1 ) 
	      << ", ";
    if ( mvaIds[iter].isValid() &&
	 !mvaIds[iter]->empty() &&
	 electrons.isValid() ) {
      reco::GsfElectronRef ele(electrons,0);
      std::cout << "\"" << mvaIdTags_[iter].instance() << "\"";
      if ( ele.isNonnull() ) { std::cout << "(example value: " << float( (*mvaIds[iter])[ele] ) << ")"; }
    }
    std::cout << std::endl;
  }
  for ( unsigned int iter = 0; iter < mvaSeeds.size(); ++iter ) {
    std::cout << "  mvaSeed:           " 
	      << int( mvaSeeds[iter].isValid() ? mvaSeeds[iter]->size() : -1 ) 
	      << ", ";
    if ( mvaSeeds[iter].isValid() && 
	 !mvaSeeds[iter]->empty() &&
	 gsf_tracks.isValid() ) {
      reco::GsfTrackRef gsf(gsf_tracks,0);
      std::cout << "\"" << mvaSeedTags_[iter].instance() << "\"";
      if ( gsf.isNonnull() ) { std::cout << "(example value: " << float( (*mvaSeeds[iter])[gsf] ) << ")"; }
    }
    std::cout << std::endl;
  }
  //
  std::cout << "  all mvaSeeds:       " << std::endl;
  if ( !mvaSeeds.empty() &&
       mvaSeeds[0].isValid() && 
       !mvaSeeds[0]->empty() &&
       gsf_tracks.isValid() ) {
    for ( unsigned int iter = 0; iter < gsf_tracks->size(); ++iter ) {
      reco::GsfTrackRef gsf(gsf_tracks,iter);
      if ( gsf.isNonnull() ) { 
	std::cout << std::setprecision(2) 
		  << std::setw(8)
		  << " pt: " << gsf->pt()
		  << std::setw(8)
		  << " eta: " << gsf->eta()
		  << std::setw(8)
		  << " phi: " << gsf->phi()
		  << std::setw(8)
	  	  << " mva: " << float( (*mvaSeeds[0])[gsf] ) 
		  << std::endl; 
      }
    }
    std::cout << std::endl;
  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LowPtGsfElectronsAnalyzer);
