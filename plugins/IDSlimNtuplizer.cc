#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <iostream> 
#include <boost/core/demangle.hpp>

namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
typedef std::map<unsigned long,int> PdgIds;

class IDSlimNtuplizer : public edm::EDAnalyzer {
  
public:
  
  explicit IDSlimNtuplizer( const edm::ParameterSet& );
  ~IDSlimNtuplizer() {}


  // ---------------------------------------
  // Main methods
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );      

  // Delete collections at the end                                        
  void deleteCollections(); 

  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ); 

  // GEN-based method to provide a sample of "signal" electrons            
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,  
			  float muon_pt = 7., float muon_eta = 1.5 );

  // ---------------------------------------
  // Utility methods
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );

  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  
private:
  
  // Misc  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDSlimNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double prescale_; 
  int isAOD_;
  bool isMC_;
  double minTrackPt_;   
  bool tag_side_muon;
  
  // Generic collections
  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_;       // AOD
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_; // AOD
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > packedCands_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > packedCandsH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > lostTracks_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > lostTracksH_;

  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;
  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;
  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsEGM_;
  edm::Handle<EcalRecHitCollection> ebRecHitsEGMH_;
  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsEGM_;
  edm::Handle<EcalRecHitCollection> eeRecHitsEGMH_;

  noZS::EcalClusterLazyTools *ecalTools_;       

  const edm::EDGetTokenT<edm::ValueMap< reco::DeDxData > > dEdx1Tag_;  // AOD only
  edm::Handle< edm::ValueMap< reco::DeDxData > > dEdx1H_;
  std::vector<const edm::ValueMap<reco::DeDxData>*> v_dEdx_;

  // Low pT collections
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; // AOD and MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;
  const edm::EDGetTokenT< edm::Association<reco::TrackCollection> > gsfTrackLinks_; // AOD
  edm::Handle<edm::Association<reco::TrackCollection> > gsfTrackLinksH_;
  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;
  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_;
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;

  PdgIds pdgids_;  
};

////////////////////////////////////////////////////////////////////////////////
IDSlimNtuplizer::IDSlimNtuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    prescale_(cfg.getParameter<double>("prescale")),  
    isAOD_(-1),
    isMC_(true),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
    beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
    beamspotH_(),
    genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_(),
    packedCands_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCands"))),
    packedCandsH_(),
    lostTracks_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("lostTracks"))),
    lostTracksH_(),
    ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
    ebRecHitsH_(),
    eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
    eeRecHitsH_(),
    ebRecHitsEGM_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHitsEGM"))),
    ebRecHitsEGMH_(),
    eeRecHitsEGM_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHitsEGM"))),
    eeRecHitsEGMH_(),
    dEdx1Tag_(consumes< edm::ValueMap<reco::DeDxData> >(cfg.getParameter<edm::InputTag>("dEdx1Tag"))),
    // Low pT collections
    gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
    gsfTracksH_(),
    gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
    patElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
    gsfElectronsH_(),
    gsfTrackLinks_(consumes< edm::Association<reco::TrackCollection> >(cfg.getParameter<edm::InputTag>("gsfTrackLinks"))),
    gsfTrackLinksH_(),
    packedCandLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("packedCandLinks"))),
    packedCandLinksH_(),
    lostTrackLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("lostTrackLinks"))),
    lostTrackLinksH_(),
    mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
    mvaUnbiasedH_(),
    mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
    mvaPtbiasedH_(),
    mvaValueLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPt"))),
    mvaValueLowPtH_(),
    pdgids_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "Verbosity level: "<< verbose_ << std::endl;
  }

////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDSlimNtuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) { }

////////////////////////////////////////////////////////////////////////////////
//
void IDSlimNtuplizer::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // Reset ntuple
  ntuple_.reset();

  // Update all handles - MUST be called every event! 
  readCollections(event,setup);
  
  // Gen level electrons from B                       
  std::set<reco::GenParticlePtr> signal_electrons;
  if (isMC_) signalElectrons(signal_electrons);          
  if (!tag_side_muon) return;

  // Loop over low-pT electrons                  
  for( size_t electronlooper = 0; electronlooper < gsfElectronsH_->size(); electronlooper++ ) {

    // ---------------------------------
    // General event info - we save 1 entry per electron
    ntuple_.is_mc_  = isMC_;
    ntuple_.is_aod_ = isAOD_;
    ntuple_.fill_evt( event.id() );
    ntuple_.set_rho( *rhoH_ );
    ntuple_.is_egamma_ = false;
    ntuple_.weight_ = 1.;
  
    // Low pT electrons
    const reco::GsfElectronPtr ele(gsfElectronsH_, electronlooper);
    
    // filter candidates: there must be a trk, a gsf track and a SC linked to the electron (should be always the case). 
    // Ele, trk and gsf_trk must have pT>0.5
    // trk must be high purity
    reco::GsfTrackPtr gsf = edm::refToPtr(ele->gsfTrack());    
    reco::TrackPtr trk;
    reco::SuperClusterRef sc = ele->superCluster();
    if ( !validPtr(gsf) )     continue;     
    if ( !gsfToTrk(gsf,trk) ) continue;
    if ( sc.isNull() )        continue; 
    if ( ele->pt() < minTrackPt_ ) continue;
    if ( gsf->pt() < minTrackPt_ ) continue;
    if ( trk->pt() < minTrackPt_ ) continue;
    if (!trk->quality( reco::TrackBase::qualityByName("highPurity") ) ) continue; 

    // Work on the electron candidate
    TVector3 eleTV3(0,0,0);
    eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());

    
    // ---------------------------------
    // Signal or fake electron, using gen-level info (-999 means nothing found with dR<0.05 )
    float dRGenMin=999.;
    reco::GenParticlePtr theGenParticle;
    for ( auto sig : signal_electrons ) {      
      TVector3 genTV3(0,0,0);
      genTV3.SetPtEtaPhi(sig->pt(), sig->eta(), sig->phi());
      float dR = eleTV3.DeltaR(genTV3);
      if (dR<dRGenMin) { 
	theGenParticle = sig;
	dRGenMin=dR;
      }
    }
    if (dRGenMin<0.05) {
      ntuple_.fill_gen( theGenParticle ); 
      ntuple_.gen_dR_ = dRGenMin;
      ntuple_.is_e_ = true;
      ntuple_.is_other_ = false;
      ntuple_.gen_tag_side_=tag_side_muon;
    } else { 
      ntuple_.fill_gen_default(); 
      ntuple_.gen_dR_ = dRGenMin;
      ntuple_.is_e_ = false;
      ntuple_.is_other_ = true;
      ntuple_.gen_tag_side_=tag_side_muon;
    }

    // prescale fake electrons 
    if (dRGenMin>=0.05) {
      if ( gRandom->Rndm() < prescale_  ) continue;
      ntuple_.weight_ = prescale_;
    }  

    // ---------------------------------
    // Electron ID: dirty hack as ID is not in Event nor embedded in pat::Electron
    float mva_value = -999.;
    int mva_id = -999;
    if ( mvaValueLowPtH_.isValid() && 
	 mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
      mva_value = mvaValueLowPtH_->get( ele.key() );
    } else {
      std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
    }

    // ---------------------------------
    // Fill ele info 
    ntuple_.fill_ele( ele, mva_value, mva_id, -999, *rhoH_ );

    // ---------------------------------
    // Supercluster linked to electron
    if(isAOD_) ntuple_.fill_supercluster(ele, ecalTools_);
    ntuple_.fill_supercluster_miniAOD(ele);  

    // ---------------------------------
    // GSF track linked to electron
    ntuple_.fill_gsf( gsf, *beamspotH_ );
    TVector3 gsfTV3(0,0,0);
    gsfTV3.SetPtEtaPhi(gsf->ptMode(), gsf->etaMode(), gsf->phiMode()); 
    ntuple_.gsf_dr_ = eleTV3.DeltaR(gsfTV3);  
    float unbiasedSeedBdt_ = (*mvaUnbiasedH_)[gsf];
    float ptbiasedSeedBdt_ = (*mvaPtbiasedH_)[gsf];
    ntuple_.fill_bdt( unbiasedSeedBdt_, ptbiasedSeedBdt_ );
    // check track extra x tangenti


    // ---------------------------------
    // KTF track linked to electron
    ntuple_.fill_trk (trk, *beamspotH_ );   
    TVector3 trkTV3(0,0,0);
    trkTV3.SetPtEtaPhi(trk->pt(), trk->eta(), trk->phi());  
    ntuple_.trk_dr_ = eleTV3.DeltaR(trkTV3);  
    PdgIds::const_iterator pos = pdgids_.find(trk.key());
    if ( pos != pdgids_.end() ) { ntuple_.pdg_id_ = pos->second; }
    if ( isAOD_ == 1 ) { 
      v_dEdx_.clear();
      v_dEdx_.push_back(dEdx1H_.product());
      ntuple_.fill_trk_dEdx( trk, v_dEdx_ );
    } else 
      ntuple_.fill_trk_dEdx_default();
    

    tree_->Fill();

  } // electron looper

  // Delete
  deleteCollections();
}

////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(patElectrons_,gsfElectronsH_);
      if ( gsfElectronsH_.isValid() ) { 
	isAOD_ = 0;
	std::cout << "File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << " failed to find a standard AOD or miniAOD particle collection " 
	  << std::endl;
      }
    }
  } else if ( isAOD_ == 1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
  } else if ( isAOD_ == 0 ) {
    event.getByToken(patElectrons_,gsfElectronsH_);
  } else {
    throw cms::Exception(" Invalid value for isAOD: ") 
      << isAOD_ 
      << std::endl;
  }

  // Generic collections 
  event.getByToken(rho_, rhoH_);
  event.getByToken(beamspot_, beamspotH_);
  
  // GEN particles
  if ( isMC_ ) {
    if ( isAOD_ == 1 ) { 
      event.getByToken(genParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
      event.getByToken(prunedGenParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
      }
    }
  }

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
  }

  // RecHits 
  if ( isAOD_ == 1 ) {   
    event.getByToken(ebRecHits_, ebRecHitsH_);
    event.getByToken(eeRecHits_, eeRecHitsH_);
    if (!ebRecHitsH_.isValid()) std::cout << "rechits EB not valid" << std::endl;
    if (!eeRecHitsH_.isValid()) std::cout << "rechits EE not valid" << std::endl;
    if (ebRecHitsH_.isValid() && eeRecHitsH_.isValid()) ecalTools_ = new noZS::EcalClusterLazyTools(event, setup, ebRecHits_, eeRecHits_);
  }
  if ( isAOD_ == 0 ) {   
    event.getByToken(ebRecHitsEGM_, ebRecHitsEGMH_);
    event.getByToken(eeRecHitsEGM_, eeRecHitsEGMH_);
    if (!ebRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EB not valid" << std::endl;
    if (!eeRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EE not valid" << std::endl;
    if (ebRecHitsEGMH_.isValid() && eeRecHitsEGMH_.isValid()) ecalTools_ = new noZS::EcalClusterLazyTools(event, setup, ebRecHitsEGM_, eeRecHitsEGM_);
  }


  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  if ( isAOD_ == 1 ) { 
    event.getByToken(gsfTrackLinks_, gsfTrackLinksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCandLinks_, packedCandLinksH_); 
    event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  }
  
  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);

  // dEdx
  if ( isAOD_ == 1 ) event.getByToken(dEdx1Tag_, dEdx1H_);
}

void IDSlimNtuplizer::deleteCollections( ) {

  delete ecalTools_;
}

// Gen-level electons from B
void IDSlimNtuplizer::signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ) {

  signal_electrons.clear();
  std::set<reco::GenParticlePtr> electrons_from_B;
  genElectronsFromB(electrons_from_B);
  for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
}

// Gen-level electons from B 
void IDSlimNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
				     float muon_pt, float muon_eta ) {   
  
  electrons_from_B.clear();
  tag_side_muon = false;
  
  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
    
    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "ERROR! GenParticlePtr:"
		<< " gen.isNull(): " << gen.isNull()
		<< " gen.isAvailable(): " << gen.isAvailable()
		<< std::endl;
      continue;
    }
    
    // Last copy of GEN electron 
    bool is_ele = std::abs(gen->pdgId()) == 11 && gen->isLastCopy(); //@@ not a method of Candidate
    
    // Does GEN ele comes from B decay?
    bool non_resonant = gen->numberOfMothers() >= 1 && gen->mother() &&   // has mother
      std::abs(gen->mother()->pdgId()) > 510 &&                           // mother is B
      std::abs(gen->mother()->pdgId()) < 546;                             // mother is B
    bool resonant = gen->numberOfMothers() >= 1 && gen->mother() &&       // has mother
      std::abs(gen->mother()->pdgId()) == 443 &&                          // mother is J/psi
      gen->mother()->numberOfMothers() >= 1 && gen->mother()->mother() && // has grandmother
      std::abs(gen->mother()->mother()->pdgId()) > 510 &&                 // grandmother is B
      std::abs(gen->mother()->mother()->pdgId()) < 546;                   // grandmother is B
    
    //  Check for tag side muon
    tag_side_muon |= ( std::abs(gen->pdgId()) == 13 && gen->isLastCopy() 
		       && 
		       ( ( gen->numberOfMothers() >= 1 &&
			   gen->mother() &&
			   std::abs(gen->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->pdgId()) < 546 && 
			   gen->mother()->pt() > muon_pt && 
			   std::abs(gen->mother()->eta()) < muon_eta ) 
			 ||
			 ( gen->numberOfMothers() >= 1 && 
			   gen->mother() &&
			   gen->mother()->numberOfMothers() >= 1 && 
			   gen->mother()->mother() &&
			   std::abs(gen->mother()->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->mother()->pdgId()) < 546 && 
			   gen->mother()->mother()->pt() > muon_pt &&
			   std::abs(gen->mother()->mother()->eta()) < muon_eta ) ) );
    
    // is coming from a B
    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
      electrons_from_B.insert(gen);
      if ( verbose_ > 0 ) {
	std::cout << "electronsFromB: "
		  << " #signal_electrons: " << electrons_from_B.size()
		  << " resonant? " << resonant
		  << " non resonant? " << non_resonant
		  << " tag-side muon? " << tag_side_muon
		  << std::endl;
      }
    }
    
  } // genParticles loop

  // We don't want this for now
  //// if ( !tag_side_muon ) { electrons_from_B.clear(); }
}

////////////////////////////////////////////////////////////////////////////////
template <typename T> 
bool IDSlimNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

////////////////////////////////////////////////////////////////////////////////
bool IDSlimNtuplizer::gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  edm::RefToBase<TrajectorySeed> traj;
  if ( gsf->extra().isNonnull() && gsf->extra().isAvailable() ) { 
    traj = gsf->seedRef(); 
  } else {
    if ( verbose_ > 2 ) { // TrackExtra are not stored by default in MINIAOD
      std::cout << "ERROR: TrackExtra:" 
		<< " gsf->extra().isNull(): " << gsf->extra().isNull()
		<< " gsf->extra().isAvailable(): " << gsf->extra().isAvailable()
		<< std::endl; 
    }
    return false;
  }
  if ( traj.isNull() || !traj.isAvailable() ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR: TrajectorySeedRef:" 
		<< " traj.isNull(): " << traj.isNull()
		<< " traj.isAvailable(): " << traj.isAvailable()
		<< std::endl; 
    }
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool IDSlimNtuplizer::seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  if ( !seed->isTrackerDriven() ) {
    if ( verbose_ > 3 ) {
      std::cout << "INFO! ElectronSeedPtr:"
		<< " seed->isTrackerDriven(): " << seed->isTrackerDriven()
		<< std::endl;
    }
  }
  trk = edm::refToPtr(seed->ctfTrack());
  if ( !validPtr(trk) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "INFO! TrackPtr:"
		<< " trk.isNull(): " << trk.isNull()
		<< " trk.isAvailable(): " << trk.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool IDSlimNtuplizer::gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk ) {   

  // Attempt to navigate via Seed (and TrackExtra) to Track
  reco::ElectronSeedPtr seed;
  if ( gsfToSeed(gsf,seed) && seedToTrk(seed,trk) ) { return true; }

  // ... if above fails (e.g. TrackExtra missing), attempt to use ...
  if ( isAOD_ == 1 ) {
    // ... track Association in AOD
    reco::GsfTrackRef gsf_ref(gsfTracksH_,(unsigned long)gsf.key());
    reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
    trk = edm::refToPtr(trk_ref);
    if ( validPtr(trk) ) { return true; }
  } else if ( isAOD_ == 0 ) {
    // ... "packedCand" Associations in MINIAOD
    reco::GsfTrackRef gsf_ref(gsfTracksH_,(unsigned long)gsf.key());
    pat::PackedCandidateRef packed_ref = (*packedCandLinksH_)[gsf_ref];
    if ( packed_ref.isAvailable() && 
	 packed_ref.isNonnull() && 
	 packed_ref->bestTrack() != nullptr ) { 
      trk = reco::TrackPtr(packed_ref->bestTrack(),packed_ref.key());
      if ( validPtr(trk) ) { return true; }
    }
    // ... "lostTrack" Associations in MINIAOD
    pat::PackedCandidateRef lost_ref = (*lostTrackLinksH_)[gsf_ref];
    if ( lost_ref.isAvailable() && 
	 lost_ref.isNonnull() && 
	 lost_ref->bestTrack() != nullptr ) { 
      trk = reco::TrackPtr(lost_ref->bestTrack(),lost_ref.key());
      if ( validPtr(trk) ) { return true; }
    }
  }

  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDSlimNtuplizer);
