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
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <boost/core/demangle.hpp>

namespace reco { typedef edm::Ref<CaloClusterCollection> CaloClusterRef; }
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<PreId> PreIdPtr; }

//constexpr int NEG_INT = -10;
//constexpr float NEG_FLOAT = -10.;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <class T1, class T2> 
class DeltaR2 {
public:
  DeltaR2( edm::Ptr<T1> ptr1, edm::Ptr<T2> ptr2, double dr2 ) {
    ptr1_ = ptr1; ptr2_ = ptr2; dr2_ = dr2;
  };
  edm::Ptr<T1> ptr1_;
  edm::Ptr<T2> ptr2_;
  double dr2_ = IDNtuple::NEG_FLOAT;
  static bool compare_by_dr2( const DeltaR2<T1,T2>& a, const DeltaR2<T1,T2>& b ) {
    return a.dr2_ < b.dr2_;
  };
};

template <class T1, class T2> 
std::ostream& operator<< ( std::ostream& out, const DeltaR2<T1,T2>& obj ) {
  out << "Class type:    " << boost::core::demangle( typeid(obj).name() ) << "\n"
      << "  Ptr1 type:   " << boost::core::demangle( typeid(obj.ptr1_).name() ) << "\n"
      << "  Ptr2 type:   " << boost::core::demangle( typeid(obj.ptr2_).name() ) << "\n"
      << "  Ptr1 id/key: " << obj.ptr1_.id() << "/" << obj.ptr1_.key() << "\n"
      << "  Ptr2 id/key: " << obj.ptr2_.id() << "/" << obj.ptr2_.key() << "\n"
      << "  dR(1,2):     " << ( obj.dr2_ >= 0. ? sqrt(obj.dr2_) : IDNtuple::NEG_FLOAT );
  return out;
};

typedef DeltaR2<reco::Candidate,reco::Track> SigToTrkR2;
typedef DeltaR2<reco::Candidate,reco::GsfTrack> SigToGsfR2;
typedef DeltaR2<reco::Candidate,reco::GsfElectron> SigToEleR2;
typedef DeltaR2<reco::Track,reco::GsfTrack> TrkToGsfR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfElectron> GsfToEleR2;
typedef DeltaR2<reco::Track,reco::GsfElectron> TrkToEleR2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
class ElectronChain {

public :

  explicit ElectronChain() {;}
  ~ElectronChain() {;}

public:

  bool is_e_ = false;
  bool is_mc_ = false;
  bool is_egamma_ = false;

  // "Signal electron" info 
  reco::CandidatePtr sig_;

  // Track info
  reco::TrackPtr trk_;
  float trk_dr_ = IDNtuple::NEG_FLOAT;
  bool trk_match_ = false;

  // CaloCluster info
  reco::CaloClusterPtr calo_;

  // ElectronSeed info
  reco::ElectronSeedPtr seed_;
  bool seed_tracker_driven_ = false;
  bool seed_ecal_driven_ = false;
  
  // PreId info
  reco::PreIdPtr preid_ecal_;
  reco::PreIdPtr preid_hcal_;

  // Seed BDTs
  float unbiased_ = IDNtuple::NEG_FLOAT;
  float ptbiased_ = IDNtuple::NEG_FLOAT;

  // GsfTrack info
  reco::GsfTrackPtr gsf_;
  float gsf_dr_ = IDNtuple::NEG_FLOAT;
  bool gsf_match_ = false;

  // GsfElectron info
  reco::GsfElectronPtr ele_;
  float ele_dr_ = IDNtuple::NEG_FLOAT;
  bool ele_match_ = false;

};

template <typename T> 
void key_id( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  ss << "id/key: " 
     << ptr.id() << ", " 
     << std::setw(4) << int( ptr.key() );// << "\n";
}

template <typename T> 
void pt_eta_phi( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << " pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->pt() << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();// << "\n";
}

void pt_eta_phi( const reco::CaloClusterPtr& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << " Et/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->energy() / std::cosh(ptr->eta()) << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();// << "\n";
}

std::ostream& operator<< ( std::ostream& out, const ElectronChain& obj ) {
  std::stringstream ss;
  ss << "ElectronChain:"
     << " is_egamma: " << obj.is_egamma_
     << " is_e: " << obj.is_e_
     << " is_mc: " << obj.is_mc_;
  ss << "\n SIG:   "; key_id(obj.sig_,ss); pt_eta_phi(obj.sig_,ss);
  ss << "\n TRK:   "; key_id(obj.trk_,ss); pt_eta_phi(obj.trk_,ss);
  ss << "\n CALO:  "; key_id(obj.calo_,ss); pt_eta_phi(obj.calo_,ss);
  ss << "\n GSF:   "; key_id(obj.gsf_,ss); pt_eta_phi(obj.gsf_,ss);
  ss << "\n ELE:   "; key_id(obj.ele_,ss); pt_eta_phi(obj.ele_,ss);
  ss << "\n SEED:  "; key_id(obj.seed_,ss);
  ss << ", TRK driven: " << obj.seed_tracker_driven_
     << ", ECAL driven: " << obj.seed_ecal_driven_;
  //ss << "\n PREID: "; key_id(obj.preid_ecal_,ss);
  ss << "\n BDTs:  " << "unbiased: " << std::setprecision(4) << obj.unbiased_ 
     << ", ptbiased: " << std::setprecision(4) << obj.ptbiased_;
  ss << "\n MATCH: "
     << "trk: " << obj.trk_match_ << "/" << std::setprecision(4) << obj.trk_dr_ 
     << ", gsf: " << obj.gsf_match_ << "/" << std::setprecision(4) << obj.gsf_dr_ 
     << ", ele: " << obj.ele_match_ << "/" << std::setprecision(4) << obj.ele_dr_;
  out << ss.str();
  return out;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
class IDNtuplizer : public edm::EDAnalyzer {
  
public:
  
  explicit IDNtuplizer( const edm::ParameterSet& );
  ~IDNtuplizer() {}
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
  ///////////////
  // Main methods
  ///////////////

  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );

  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::CandidatePtr>& signal_electrons );

  // GEN-based method to provide a sample of "signal" electrons
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B, 
			  float muon_pt = 7., float muon_eta = 1.5 );

  // Links "signal" electrons to reconstructed objects
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       edm::Handle< std::vector<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand );
  
  // Top-level method that creates ElectronChain objects for EGamma electrons
  void pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		    std::vector<SigToTrkR2>& sig2trk,
		    std::vector<SigToTrkR2>& other_trk );
  
  // Top-level method that creates ElectronChain objects for low pT electrons
  void lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector<SigToTrkR2>& sig2trk,
		       std::vector<SigToTrkR2>& other_trk );

  // Lower-level method that creates ElectronChain objects for low pT electrons
  void signal( bool is_egamma,
	       std::set<reco::CandidatePtr>& signal_electrons,
	       std::vector<SigToTrkR2>& sig2trk,
	       std::vector<SigToGsfR2>& sig2gsf,
	       std::vector<SigToEleR2>& sig2ele,
	       std::vector<TrkToGsfR2>& trk2gsf,
	       std::vector<TrkToEleR2>& trk2ele,
	       std::vector<GsfToEleR2>& gsf2ele);
  
  // Lower-level method that creates ElectronChain objects for low pT electrons
  void fakes( bool is_egamma,
	      std::vector<SigToTrkR2>& other_trk,
	      std::vector<SigToGsfR2>& other_gsf,
	      std::vector<TrkToGsfR2>& trk2gsf,
	      std::vector<TrkToEleR2>& trk2ele,
	      std::vector<GsfToEleR2>& gsf2ele );
  
  // Fills tree per ElectronChain object
  void fill( const edm::Event& event, const edm::EventSetup& setup );

  //////////////////
  // Utility methods
  //////////////////

  template <typename T> bool validPtr( edm::Ptr<T>& ptr );
  
  template <typename T> 
  bool filterCand( edm::Ptr<T>& cand );
  bool filterCand( edm::Ptr<reco::Track>& trk );
  
  template <typename T1, typename T2> 
  float deltaR2( T1& cand1, T2& cand2 ); // templated, used by sigToCandLinks
  float deltaR2( reco::TrackPtr& trk, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::CandidatePtr& sig, reco::GsfTrackPtr& gsf ); // overload
  
  template <typename T>
  void matchElectron( reco::CandidatePtr& sig,
		      std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		      edm::Ptr<T>& cand, float& dr2, bool& match ); // pass by ref
  
  void trkToGsfLinks( edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::vector<TrkToGsfR2>& trk2gsf );

  void trkToEleLinks( edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      edm::Handle< std::vector<reco::GsfElectron> >& gsfElectron,
		      std::vector<TrkToEleR2>& trk2gsf );
  
  void gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< std::vector<reco::GsfElectron> >& gsfElectrons,
		      std::vector<GsfToEleR2>& gsf2ele );
  
  bool eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  bool seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo );
  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk );
  bool eleToTrk( reco::GsfElectronPtr& ele, reco::TrackPtr& trk );
  
private:
  
  // Misc
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double dr_max_; // Max DeltaR value considered
  double dr_threshold_; // Threshold for DeltaR matching
  double prescale_;
  int isAOD_;
  bool isMC_;
  double minTrackPt_;
  
  // Generic collections

  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticles_; // AOD
  edm::Handle< std::vector<reco::GenParticle> > genParticlesH_;

//  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
//  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;
//
//  const edm::EDGetTokenT< edm::View<reco::Candidate> > packedGenParticles_; // MINIAOD
//  edm::Handle< edm::View<reco::Candidate> > packedGenParticlesH_;

  const edm::EDGetTokenT< std::vector<reco::Track> > ctfTracks_; // AOD
  edm::Handle< std::vector<reco::Track> > ctfTracksH_;

  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;

  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;

  //noZS::EcalClusterLazyTools ecalTools_;
  
  const edm::EDGetTokenT<reco::SuperClusterCollection> barrelSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> barrelSCsH_;

  const edm::EDGetTokenT<reco::SuperClusterCollection> endcapSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> endcapSCsH_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> packedCands_; // MINIAOD
  edm::Handle<pat::PackedCandidateCollection> packedCandsH_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracks_; // MINIAOD
  edm::Handle<pat::PackedCandidateCollection> lostTracksH_;

  // Low pT collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeeds_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsEcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsEcalH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsHcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsHcalH_;

  const edm::EDGetTokenT< edm::ValueMap<reco::PreIdRef> > preIdRefs_; // AOD
  edm::Handle< edm::ValueMap<reco::PreIdRef> > preIdRefsH_;

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; // AOD and MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;

  const edm::EDGetTokenT< std::vector<reco::GsfElectron> > gsfElectrons_; // AOD
  edm::Handle< std::vector<reco::GsfElectron> > gsfElectronsH_;

  const edm::EDGetTokenT< std::vector<pat::Electron> > patElectrons_; // MINIAOD
  edm::Handle< std::vector<pat::Electron> > patElectronsH_;

  const edm::EDGetTokenT< edm::Association<reco::TrackCollection> > gsfTrackLinks_; // AOD
  edm::Handle<edm::Association<reco::TrackCollection> > gsfTrackLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; //  on the fly?
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;

//  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValue_; // on the fly?
//  edm::Handle< edm::ValueMap<float> > mvaValueH_;

//  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaId_; // on the fly?
//  edm::Handle< edm::ValueMap<bool> > mvaIdH_;

  // EGamma collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeedsEGamma_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsEGammaH_; // AOD

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_; // AOD
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_MAOD_; // MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksEGammaH_;

  const edm::EDGetTokenT< std::vector<reco::GsfElectron> > gsfElectronsEGamma_; // AOD
  edm::Handle< std::vector<reco::GsfElectron> > gsfElectronsEGammaH_;

  const edm::EDGetTokenT< std::vector<pat::Electron> > patElectronsEGamma_; // MINIAOD
  edm::Handle< std::vector<pat::Electron> > patElectronsEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaIdEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<bool> > mvaIdEGammaH_;

  // Conversions

  //@@ const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;

  std::vector<ElectronChain> chains_;
  
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizer::IDNtuplizer( const edm::ParameterSet& cfg ) :
  tree_(0),
  ntuple_{},
  verbose_(cfg.getParameter<int>("verbose")),
  check_from_B_(cfg.getParameter<bool>("checkFromB")),
  dr_max_(cfg.getParameter<double>("drMax")),
  dr_threshold_(cfg.getParameter<double>("drThreshold")),
  prescale_(cfg.getParameter<double>("prescale")),
  isAOD_(-1),
  isMC_(true),
  minTrackPt_(cfg.getParameter<double>("minTrackPt")),
  // Generic collections
  rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
  rhoH_(),
  beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
  beamspotH_(),
  genParticles_(consumes< std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
  genParticlesH_(),
  //genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
  //prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
  //genParticlesH_(),
  //packedGenParticles_(consumes< edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("packedGenParticles"))),
  //packedGenParticlesH_(),
  ctfTracks_(consumes< std::vector<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
  ctfTracksH_(),
  ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
  ebRecHitsH_(),
  eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
  eeRecHitsH_(),
  //ecalTools_(),
  barrelSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("barrelSuperClusters"))),
  barrelSCsH_(),
  endcapSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("endcapSuperClusters"))),
  endcapSCsH_(),
  packedCands_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
  packedCandsH_(),
  lostTracks_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
  lostTracksH_(),
  // Low pT collections
  eleSeeds_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))),
  eleSeedsH_(),
  preIdsEcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsEcal"))),
  preIdsEcalH_(),
  preIdsHcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsHcal"))),
  preIdsHcalH_(),
  preIdRefs_(consumes< edm::ValueMap<reco::PreIdRef> >(cfg.getParameter<edm::InputTag>("preIdRefs"))),
  preIdRefsH_(),
  gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  gsfTracksH_(),
  gsfElectrons_(consumes< std::vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
  gsfElectronsH_(),
  patElectrons_(consumes< std::vector<pat::Electron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
  patElectronsH_(),
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
  //mvaValue_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValue"))),
  //mvaValueH_(),
  //mvaId_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("mvaId"))),
  //mvaIdH_(),
  // EGamma collections
  eleSeedsEGamma_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeedsEGamma"))),
  eleSeedsEGammaH_(),
  gsfTracksEGamma_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma"))),
  gsfTracksEGamma_MAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma_MAOD"))),
  gsfTracksEGammaH_(),
  gsfElectronsEGamma_(consumes< std::vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma"))),
  gsfElectronsEGammaH_(),
  patElectronsEGamma_(consumes< std::vector<pat::Electron> >(cfg.getParameter<edm::InputTag>("patElectronsEGamma"))),
  patElectronsEGammaH_(),
  mvaValueEGamma_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGamma"))),
  mvaValueEGammaH_(),
  mvaIdEGamma_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("mvaIdEGamma"))),
  mvaIdEGammaH_(),
  // Conversions
  //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
  chains_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "Verbosity level: "<< verbose_ << std::endl;
  }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDNtuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) {
  //@@ ?
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // Update all handles - MUST be called every event!
  readCollections(event,setup);
  
  // Clear ElectronChain vector and populate
  chains_.clear();

  //HERE HERE HERE;
  //signal_electrons; check all udpated code: gen->sig;
  //wrap getElctronsFromB with method that handles both MC and data;

  // Identify signal electrons, from (GEN) MC or data control regions
  std::set<reco::CandidatePtr> signal_electrons;
  signalElectrons(signal_electrons);

  // Match "signal" electrons to generalTracks
  std::vector<SigToTrkR2> sig2trk;
  std::vector<SigToTrkR2> other_trk;
  sigToCandLinks<reco::Track>( signal_electrons, ctfTracksH_, sig2trk, other_trk );
  if ( verbose_ > 2 ) {
    std::cout << "sig2trk.size(): " << sig2trk.size() << std::endl;
    for ( auto iter : sig2trk ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Populate ElectronChain objects using PF electrons
  pfElectrons( signal_electrons, sig2trk, other_trk );

  // Populate ElectronChain objects using low-pT electrons
  lowPtElectrons( signal_electrons, sig2trk, other_trk );
  
  // Print ElectronChain objects
  if ( verbose_ > 0 ) {
    for ( auto iter : chains_ ) { std::cout << iter << std::endl; }
  }
  
  // Fill ntuple
  fill(event,setup);

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_); // std::vector<reco::GsfElectron>
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(patElectrons_,patElectronsH_); // std::vector<pat::Electron>
      if ( patElectronsH_.isValid() ) { 
	isAOD_ = 0;
	std::cout << "File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << " failed to find a standard AOD or miniAOD particle collection " 
	  << std::endl;
      }
    }
  } else if ( isAOD_ == 1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_); // std::vector<reco::GsfElectron>
  } else if ( isAOD_ == 0 ) {
    event.getByToken(patElectrons_,patElectronsH_); // std::vector<pat::Electron>
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
//      event.getByToken(prunedGenParticles_, genParticlesH_);
//      if ( !(genParticlesH_.isValid()) ) { 
//	isMC_ = false;
//	std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
//      }
    }
  }
//  if ( isMC_ && isAOD_ == 0 ) { event.getByToken(packedGenParticles_, packedGenParticlesH_); }

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
  }

  // RecHits and SuperClusters
  if ( isAOD_ == 1 ) { 
    event.getByToken(ebRecHits_, ebRecHitsH_);
    event.getByToken(eeRecHits_, eeRecHitsH_);
    event.getByToken(barrelSCs_, barrelSCsH_);
    event.getByToken(endcapSCs_, endcapSCsH_);
    //ecalTools_ = noZS::EcalClusterLazyTools(event, setup, ebRecHitsH_, eeRecHitsH_);
  }

  // ElectronSeeds and PreIds
  if ( isAOD_ == 1 ) { 
    event.getByToken(eleSeeds_, eleSeedsH_); 
    event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); 
    event.getByToken(preIdsEcal_, preIdsEcalH_); 
    event.getByToken(preIdsHcal_, preIdsHcalH_); 
    event.getByToken(preIdRefs_, preIdRefsH_); 
  }

  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  if      ( isAOD_ == 1 ) { 
    event.getByToken(gsfTrackLinks_, gsfTrackLinksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCandLinks_, packedCandLinksH_); 
    event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  }
  
  // EGamma collections 
  if      ( isAOD_ == 1 ) { event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfTracksEGamma_, gsfTracksEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(gsfTracksEGamma_MAOD_, gsfTracksEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfElectronsEGamma_, gsfElectronsEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(patElectronsEGamma_, patElectronsEGammaH_); }

  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
  //event.getByToken(mvaValue_, mvaValueH_);
  //event.getByToken(mvaId_, mvaIdH_);
  event.getByToken(mvaValueEGamma_, mvaValueEGammaH_);
  event.getByToken(mvaIdEGamma_, mvaIdEGammaH_);

  // Conversions
  //event.getByToken(convVtxFitProb_, convVtxFitProbH_);

  //////////

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::signalElectrons( std::set<reco::CandidatePtr>& signal_electrons ) {

  signal_electrons.clear();

  if ( isMC_ ) {
    // Identify "signal" (GEN) electrons from B decays
    std::set<reco::GenParticlePtr> electrons_from_B;
    genElectronsFromB(electrons_from_B);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
  } else { // data-driven
    // Identify "signal" electrons from data control regions
    //@@ FOR NOW, A HACK ...
    std::set<reco::GenParticlePtr> electrons_from_B;
    genElectronsFromB(electrons_from_B);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
				     float muon_pt, float muon_eta ) {
  
  electrons_from_B.clear();
  bool tag_side_muon = false;
  
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
    bool non_resonant = 
      gen->numberOfMothers() >= 1 && gen->mother() &&                     // has mother
      std::abs(gen->mother()->pdgId()) > 510 &&                           // mother is B
      std::abs(gen->mother()->pdgId()) < 546;                             // mother is B
    bool resonant = 
      gen->numberOfMothers() >= 1 && gen->mother() &&                     // has mother
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
//      std::cout << "electronsFromB: "
//		<< " #signal_electrons: " << electrons_from_B.size()
//		<< " resonant? " << resonant
//		<< " non resonant? " << non_resonant
//		<< " tag-side muon? " << tag_side_muon
//		<< std::endl;
    }
    
  } // genParticles loop
  
  if ( !tag_side_muon ) { electrons_from_B.clear(); }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  edm::Handle< std::vector<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand ) {
  
  sig2cand.clear();
  other_cand.clear();
  
  // DeltaR2 matches for all combinations of signal electrons and reco::Candidate
  std::vector< DeltaR2<reco::Candidate,T> > sig2cand_all;
  sig2cand_all.reserve( signal_electrons.size()*candidates->size() );
  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
    edm::Ptr<T> cand(candidates,icand);
    if ( !validPtr(cand) ) {
      std::cout << "ERROR! CandidatePtr:"
		<< " cand.isNull(): " << cand.isNull()
		<< " cand.isAvailable(): " << cand.isAvailable()
		<< std::endl;
    }
    if ( !filterCand<T>(cand) ) { continue; }
    for ( auto sig : signal_electrons ) {
      sig2cand_all.emplace_back( sig, 
				 cand, 
				 deltaR2< reco::CandidatePtr, edm::Ptr<T> >(sig,cand) );
    }
    // Note: matched candidates are removed below
    if ( prescale_ < 0. || ( gRandom->Rndm() < prescale_ ) ) { 
      reco::CandidatePtr null;
      other_cand.emplace_back( null, cand, IDNtuple::NEG_FLOAT );
    }
  }
  
  // Sort by DeltaR2!!
  std::sort( sig2cand_all.begin(), 
	     sig2cand_all.end(), 
	     DeltaR2<reco::Candidate,T>::compare_by_dr2 );
  
  // Select best matches according to DeltaR2 metric
  sig2cand.clear();
  for ( auto iter : sig2cand_all ) {
    auto sig = std::find_if( sig2cand.begin(), 
			     sig2cand.end(), 
			     [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
			       return dr2.ptr1_ == iter.ptr1_; 
			     }
			     );
    auto cand = std::find_if( sig2cand.begin(), 
			      sig2cand.end(), 
			      [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				return dr2.ptr2_ == iter.ptr2_; 
			      }
			      );
    if ( sig == sig2cand.end() &&
	 cand == sig2cand.end() && 
	 iter.dr2_ < dr_max_*dr_max_ ) {
      if ( sig2cand.size() < signal_electrons.size() ) { sig2cand.push_back(iter); }
      else { break; } // found unique match for all signal electrons
    }
  }
  
  // Remove matched candidates
  for ( auto iter : sig2cand ) {
    auto match = std::find_if( other_cand.begin(), 
			       other_cand.end(), 
			       [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				 return dr2.ptr2_ == iter.ptr2_; 
			       }
			       );
    if ( match != other_cand.end() ) { other_cand.erase(match); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
			       std::vector<SigToTrkR2>& sig2trk,
			       std::vector<SigToTrkR2>& other_trk ) {
  
  // Match "signal" electrons to GsfTracks
  std::vector<SigToGsfR2> sig2gsf;
  std::vector<SigToGsfR2> other_gsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, gsfTracksEGammaH_, 
				  sig2gsf, other_gsf );
  if ( verbose_ > 2 ) {
    std::cout << "sig2gsf.size(): " << sig2gsf.size() << std::endl;
    for ( auto iter : sig2gsf ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Match "signal" electrons to GsfElectrons
  std::vector<SigToEleR2> sig2ele;
  std::vector<SigToEleR2> other_ele;
  sigToCandLinks<reco::GsfElectron>( signal_electrons, gsfElectronsEGammaH_, 
				     sig2ele, other_ele );
  if ( verbose_ > 2 ) {
    std::cout << "sig2ele.size(): " << sig2ele.size() << std::endl;
    for ( auto iter : sig2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Match Tracks to GsfTracks
  std::vector<TrkToGsfR2> trk2gsf;
  trkToGsfLinks( ctfTracksH_, gsfTracksEGammaH_, trk2gsf );
  if ( verbose_ > 2 ) {
    std::cout << "trk2gsf.size(): " << trk2gsf.size() << std::endl;
    for ( auto iter : trk2gsf ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Match Tracks to GsfElectrons
  std::vector<TrkToEleR2> trk2ele;
  trkToEleLinks( ctfTracksH_, gsfElectronsEGammaH_, trk2ele );
  if ( verbose_ > 2 ) {
    std::cout << "trk2ele.size(): " << trk2ele.size() << std::endl;
    for ( auto iter : trk2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Match GsfTracks to GsfElectrons
  std::vector<GsfToEleR2> gsf2ele;
  gsfToEleLinks( gsfTracksEGammaH_, gsfElectronsEGammaH_, gsf2ele );
  if ( verbose_ > 2 ) {
    std::cout << "gsf2ele.size(): " << gsf2ele.size() << std::endl;
    for ( auto iter : gsf2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }

  // Populate ElectronChain objects starting with signal electrons
  signal( true,
	  signal_electrons, 
	  sig2trk,
	  sig2gsf,
	  sig2ele,
	  trk2gsf, 
	  trk2ele,
	  gsf2ele );
  
  // Populate ElectronChain objects with fake electrons
  fakes( true,
	 other_trk,
	 other_gsf, 
	 trk2gsf, 
	 trk2ele,
	 gsf2ele );
  
  if ( verbose_ > 1 ) {
    std::cout << "pfElectrons:" << " isMC: " << isMC_ << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " sig2trk.size(): " << sig2trk.size() << std::endl
	      << " sig2gsf.size(): " << sig2gsf.size() << std::endl
	      << " sig2ele.size(): " << sig2ele.size() << std::endl
	      << " ctfTracksH_->size(): " << ctfTracksH_->size() << std::endl
	      << " other_trk.size(): " << other_trk.size() << std::endl
	      << " trk2gsf.size(): " << trk2gsf.size() << std::endl
	      << " trk2ele.size(): " << trk2ele.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " other_gsf.size(): " << other_gsf.size() << std::endl
	      << " gsf2ele.size(): " << gsf2ele.size() << std::endl
	      << " gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl;
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
				  std::vector<SigToTrkR2>& sig2trk,
				  std::vector<SigToTrkR2>& other_trk ) {
  
  // Match "signal" electrons to GsfTracks
  std::vector<SigToGsfR2> sig2gsf;
  std::vector<SigToGsfR2> other_gsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, gsfTracksH_, 
				  sig2gsf, other_gsf );
  if ( verbose_ > 2 ) {
    std::cout << "sig2gsf.size(): " << sig2gsf.size() << std::endl;
    for ( auto iter : sig2gsf ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }

  // Match "signal" electrons to GsfElectrons
  std::vector<SigToEleR2> sig2ele;
  std::vector<SigToEleR2> other_ele;
  sigToCandLinks<reco::GsfElectron>( signal_electrons, gsfElectronsH_, 
				     sig2ele, other_ele );
  if ( verbose_ > 2 ) {
    std::cout << "sig2ele.size(): " << sig2ele.size() << std::endl;
    for ( auto iter : sig2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }

  // Match Tracks to GsfTracks
  std::vector<TrkToGsfR2> trk2gsf;
  trkToGsfLinks( ctfTracksH_, gsfTracksH_, trk2gsf );
  if ( verbose_ > 2 ) {
    std::cout << "trk2gsf.size(): " << trk2gsf.size() << std::endl;
    for ( auto iter : trk2gsf ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }
  
  // Match Tracks to GsfElectrons
  std::vector<TrkToEleR2> trk2ele;
  trkToEleLinks( ctfTracksH_, gsfElectronsH_, trk2ele );
  if ( verbose_ > 2 ) {
    std::cout << "trk2ele.size(): " << trk2ele.size() << std::endl;
    for ( auto iter : trk2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }

  // Match GsfTracks to GsfElectrons
  std::vector<GsfToEleR2> gsf2ele;
  gsfToEleLinks( gsfTracksH_, gsfElectronsH_, gsf2ele );
  if ( verbose_ > 2 ) {
    std::cout << "gsf2ele.size(): " << gsf2ele.size() << std::endl;
    for ( auto iter : gsf2ele ) { std::cout << iter << std::endl; }
    std::cout << std::endl;
  }

  if ( verbose_ > 1 ) {
    std::cout << "lowPtElectrons" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " sig2trk.size(): " << sig2trk.size() << std::endl
	      << " sig2gsf.size(): " << sig2gsf.size() << std::endl
	      << " sig2ele.size(): " << sig2ele.size() << std::endl
	      << " ctfTracksH_->size(): " << ctfTracksH_->size() << std::endl
	      << " other_trk.size(): " << other_trk.size() << std::endl
	      << " trk2gsf.size(): " << trk2gsf.size() << std::endl
	      << " trk2ele.size(): " << trk2ele.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " other_gsf.size(): " << other_gsf.size() << std::endl
	      << " gsf2ele.size(): " << gsf2ele.size() << std::endl
	      << " gsfElectronsH_->size(): " << gsfElectronsH_->size() << std::endl;
  }

  // Populate ElectronChain objects starting with signal electrons
  signal( false,
	  signal_electrons, 
	  sig2trk,
	  sig2gsf,
	  sig2ele,
	  trk2gsf, 
	  trk2ele,
	  gsf2ele );
  
  // Populate ElectronChain objects with fake electrons  
  fakes( false,
	 other_trk,
	 other_gsf,
	 trk2gsf,
	 trk2ele,
	 gsf2ele );
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::signal( bool is_egamma,
			  std::set<reco::CandidatePtr>& signal_electrons,
			  std::vector<SigToTrkR2>& sig2trk,
			  std::vector<SigToGsfR2>& sig2gsf,
			  std::vector<SigToEleR2>& sig2ele,
			  std::vector<TrkToGsfR2>& trk2gsf,
			  std::vector<TrkToEleR2>& trk2ele,
			  std::vector<GsfToEleR2>& gsf2ele ) {
  
  // Iterate through "signal" electrons and initialise ElectronChain objects
  for ( auto sig : signal_electrons ) {
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_e_ = true;
    chain.is_egamma_ = is_egamma;
    chain.is_mc_ = isMC_;
    chain.sig_ = sig;
    //matchElectron<reco::GsfElectron>(sig,sig2ele,chain.ele_,chain.ele_dr_,chain.ele_match_);
    matchElectron<reco::GsfTrack>(sig,sig2gsf,chain.gsf_,chain.gsf_dr_,chain.gsf_match_);
    matchElectron<reco::Track>(sig,sig2trk,chain.trk_,chain.trk_dr_,chain.trk_match_);
    
    // Update ElectronChain 
    if ( chain.gsf_match_ ) {

      // Update Seed BDTs
      if ( !chain.is_egamma_ ) { 
	chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
	chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
      }
      
      // Update (override?) GsfElectron info
      auto match_ele = std::find_if( gsf2ele.begin(), 
				     gsf2ele.end(), 
				     [chain](const GsfToEleR2& dr2) { 
				       return dr2.ptr1_ == chain.gsf_; 
				     }
				     );
      if ( match_ele != gsf2ele.end() && validPtr(match_ele->ptr2_) ) { 
	chain.ele_ = match_ele->ptr2_; 
	chain.ele_match_ = true;;
	chain.ele_dr_ = deltaR2(chain.sig_,chain.ele_);
      } else {
	chain.ele_ = reco::GsfElectronPtr();
	chain.ele_match_ = false;
	chain.ele_dr_ = IDNtuple::NEG_FLOAT;
      }
      
      // Update ElectronSeed info
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) {
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
	
//	// Update PreId info
//	if ( !chain.is_egamma_ ) { // && isAOD_ == 1
//	  chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
//	  chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
//	}

	// Update Track info
	reco::TrackPtr trk;
	if ( seedToTrk(seed,trk) ) { 
	  chain.trk_ = trk; 
	  chain.trk_match_ = true;
	  chain.trk_dr_ = deltaR2(chain.sig_,chain.trk_);
	} else {
	  chain.trk_ = reco::TrackPtr(); 
	  chain.trk_match_ = false;
	  chain.trk_dr_ = IDNtuple::NEG_FLOAT;
	}

	// Update CaloCluster info
	reco::CaloClusterPtr calo;
	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
	
	// If ECAL-driven _only_, update "surrogate" Track info
	if ( !seed->isTrackerDriven() && seed->isEcalDriven() ) {
	  auto match_gsf = std::find_if( trk2gsf.begin(), 
					 trk2gsf.end(), 
					 [chain](const TrkToGsfR2& dr2) { 
					   return dr2.ptr2_ == chain.gsf_; 
					 }
					 );
	  if ( match_gsf != trk2gsf.end() && validPtr(match_gsf->ptr1_) ) { 
	    chain.trk_ = match_gsf->ptr1_; 
	    chain.trk_match_ = true;
	    chain.trk_dr_ = deltaR2(chain.sig_,chain.trk_);
	  } else {
	    chain.trk_ = reco::TrackPtr(); 
	    chain.trk_match_ = false;
	    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
	  }
	}
      }

    } // if ( chain.gsf_match_ )

  } // for ( auto sig : signal_electrons )
    
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::fakes( bool is_egamma,
			 std::vector<SigToTrkR2>& other_trk,
			 std::vector<SigToGsfR2>& other_gsf,
			 std::vector<TrkToGsfR2>& trk2gsf,
			 std::vector<TrkToEleR2>& trk2ele,
			 std::vector<GsfToEleR2>& gsf2ele ) {
  
  // Iterate through generalTracks
  for ( auto iter : other_trk ) {

    if ( !validPtr(iter.ptr2_) ) { continue; }

    // Create ElectronChain object and start to populate
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_e_ = false;
    chain.is_egamma_ = is_egamma;
    chain.is_mc_ = isMC_;
    
    // Store Track info
    chain.trk_ = iter.ptr2_;
    chain.trk_match_ = true;
    chain.trk_dr_ = IDNtuple::NEG_FLOAT;
    
//@@ Match to GsfTracks only (to prevent instances of matching to unconnected GsfElectron)
//    // Find matched GsfElectron match (via GsfTrack/ElectronSeed or "surrogate" track)
//    auto match_ele = std::find_if( trk2ele.begin(), 
//				   trk2ele.end(), 
//				   [trk](const TrkToEleR2& dr2) { 
//				     return dr2.ptr1_ == trk;
//				   }
//				   );
//    if ( match_ele != trk2ele.end() ) { 
//      if ( validPtr(match_ele->ptr2_) ) {
//	chain.ele_ = match_ele->ptr2_;
//	chain.ele_dr_ = match_ele->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_ele->dr2_);
//	chain.ele_match_ = ( chain.ele_dr_ >= 0. );// && ( chain.ele_dr_ < dr_threshold_ );
//      }
//    } else {
//      if ( verbose_ > 3 ) {
//	std::cout << "INFO! Cannot match TrackPtr to GsfElectronPtr:"
//		  << " TrackPtr: " << trk.id() << "/" << trk.key()
//		  << std::endl;
//      }
//    }
    
    // Find matched GsfTrack match, either via ElectronSeed or just deltaR
    auto match_gsf = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [chain](const TrkToGsfR2& dr2) { 
				     return dr2.ptr1_ == chain.trk_;
				   }
				   );
    if ( match_gsf != trk2gsf.end() && validPtr(match_gsf->ptr2_) ) {
      chain.gsf_ = match_gsf->ptr2_;
      chain.gsf_match_ = true;
      chain.gsf_dr_ = match_gsf->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match_gsf->dr2_);
    } else {
      chain.gsf_ = reco::GsfTrackPtr();
      chain.gsf_match_ = false;
      chain.gsf_dr_ = IDNtuple::NEG_FLOAT;
      if ( verbose_ > 3 ) {
	std::cout << "INFO! Cannot match TrackPtr to GsfElectronPtr:"
		  << " TrackPtr: " << chain.trk_.id() << "/" << chain.trk_.key()
		  << " GsfTrackPtr: " << match_gsf->ptr2_.id() << "/" << match_gsf->ptr2_.key()
		  << std::endl;
      }
    }

    // Check GsfTrackPtr, then update ElectronChain info
    if ( !validPtr(chain.gsf_) ) { continue; } 

    // Update Seed BDTs
    if ( !chain.is_egamma_ ) { 
      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
    }
    
    // Update (override?) GsfElectron info (should be identical to that above?!)
    auto match = std::find_if( gsf2ele.begin(), 
			       gsf2ele.end(), 
			       [chain](const GsfToEleR2& dr2) { 
				 return dr2.ptr1_ == chain.gsf_; 
			       }
			       );
    if ( match != gsf2ele.end() && validPtr(match->ptr2_) ) {
      if ( validPtr(chain.ele_) && chain.ele_ != match->ptr2_ ) {
	std::cout << "ERROR! GsfElectrons are not consistent!"
		  << " matched ele id/key: " << match->ptr2_.id() << "/" << match->ptr2_.key()
		  << " chain.ele id/key: " << chain.ele_.id() << "/" << chain.ele_.key()
		  << std::endl;
      }
      chain.ele_ = match->ptr2_;
      chain.ele_match_ = true;
      chain.ele_dr_ = deltaR2(chain.trk_,chain.ele_);
    } else {
      chain.ele_ = reco::GsfElectronPtr();
      chain.ele_match_ = false;
      chain.ele_dr_ = IDNtuple::NEG_FLOAT;
    }

    // Check GsfElectron and GsfTrack are consistent
    reco::GsfTrackPtr gsf;
    if ( validPtr(chain.ele_) && 
	 eleToGsf(chain.ele_,gsf) && 
	 gsf != chain.gsf_ ) { 
      std::cout << "ERROR! GsfElectron and GsfTrack are not consistent!"
		<< " GsfElectron id/key: " << chain.ele_.id() << "/" << chain.ele_.key()
		<< " GsfTrack id/key: " << chain.gsf_.id() << "/" << chain.gsf_.key()
		<< " matched gsf id/key: " << gsf.id() << "/" << gsf.key() 
		<< std::endl;
    }
    
    // Add ElectronSeed, CaloCluster, check Track
    if ( validPtr(chain.gsf_) ) { 
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) {
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
      }
      
      reco::CaloClusterPtr calo;
      if ( seedToCalo(seed,calo) ) { 
	chain.calo_ = calo; 
      }
      if ( verbose_ > 3 ) {
	reco::TrackPtr trk;
	if ( seedToTrk(seed,trk) && ( trk != chain.trk_ ) && chain.gsf_dr_ < 0. ) { 
	  std::cout << "WARNING! Track is not consistent!"
		    << " Track id/key: " << chain.trk_.id() << "/" << chain.trk_.key()
		    << " matched trk id/key: " << trk.id() << "/" << trk.key() 
		    << ". Likely due to multiple GsfTracks originating from a"
		    << " common ElectronSeed, coupled with the use of gsfToSeed()."
		    << std::endl;
	}
      }
      
    }
    
  } // for ( auto trk : other_trk )
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::fill( const edm::Event& event,
			const edm::EventSetup& setup ) {
  
  for ( auto chain : chains_ ) {
    
    // Init tree here
    ntuple_.reset();

    // Event stuff
    ntuple_.fill_evt( event.id() );
    ntuple_.set_rho( *rhoH_ );

    // EGamma of low pT
    ntuple_.is_egamma( chain.is_egamma_ );

    // Truth label
    ntuple_.is_e( chain.is_e_ );
    ntuple_.is_other( !chain.is_e_ );
    
    // Set background weight
    if ( !chain.is_e_ ) { ntuple_.set_weight( prescale_ ); }

    // "signal" info
    if ( validPtr(chain.sig_) ) {
      ntuple_.fill_gen( chain.sig_ ); // Add sig info
    }

    // Track info
    if ( validPtr(chain.trk_) ) {
      ntuple_.has_trk( chain.trk_match_ );
      ntuple_.fill_trk( chain.trk_, *beamspotH_ );
    }

    // ElectronSeed
    if ( validPtr(chain.seed_) ) {
      ntuple_.has_seed(true);
      ntuple_.fill_seed( chain.seed_tracker_driven_, 
			 chain.seed_ecal_driven_ );
    }

    // PreId
//    noZS::EcalClusterLazyTools ecal_tools(event, setup, ebRecHits_, eeRecHits_);
//    ntuple_.fill_preid( *chain.preid_ecal_,
//			*chain.preid_hcal_,
//			*beamspotH_,
//			*rhoH_,
//			ecal_tools );
    
    // GsfTrack info
    if ( validPtr(chain.gsf_) ) {
      ntuple_.has_gsf( chain.gsf_match_ );
      ntuple_.fill_gsf( chain.gsf_, *beamspotH_ );
      ntuple_.gsf_dr( chain.gsf_dr_ );
      ntuple_.gsf_dr_mode( chain.gsf_dr_ );
      ntuple_.fill_bdt( chain.unbiased_, chain.ptbiased_ );
    }

    // GsfElectron info
    if ( validPtr(chain.ele_) ) {

      ntuple_.has_ele( chain.ele_match_ );

      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float mva_value = -999.;
      int mva_id = -999;
      if ( !chain.is_egamma_ ) {
	if ( mvaValueLowPtH_.isValid() && 
	     mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
	  mva_value = mvaValueLowPtH_->get( chain.ele_.key() );
	} else {
	  std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
      } else {
	if ( mvaValueEGammaH_.isValid() && 
	     mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_value = mvaValueEGammaH_->get( chain.ele_.key() );
	} else {
	  std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
	if ( mvaIdEGammaH_.isValid() && 
	     mvaIdEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_id = mvaIdEGammaH_->get( chain.ele_.key() ) ? 1 : 0;
	} else {
	  std::cout << "ERROR! Issue matching MVA ID to GsfElectrons!" << std::endl;
	}
      }
      
      //@@ dirty hack as is not in Event nor embedded in pat::Electron
      float conv_vtx_fit_prob = -999.;
      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
      //  conv_vtx_fit_prob = convVtxFitProb->get( chain.ele_.key() );
      //}
      
      ntuple_.fill_ele( chain.ele_, mva_value, mva_id, conv_vtx_fit_prob, *rhoH_ );

    }

    //ntuple_.fill_supercluster(chain.ele_);
    
    tree_->Fill(); 
    
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
bool IDNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
bool IDNtuplizer::filterCand( edm::Ptr<T>& cand ) { 
  if ( cand->pt() < minTrackPt_ ) { return false; }
  return true;
}

bool IDNtuplizer::filterCand( edm::Ptr<reco::Track>& trk ) { 
  if ( !(trk->quality(reco::TrackBase::qualityByName("highPurity"))) ) { return false; }
  if ( trk->pt() < minTrackPt_ ) { return false; }
  return true; 
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <typename T1, typename T2>
float IDNtuplizer::deltaR2( T1& cand1, T2& cand2 ) {
  return reco::deltaR2( cand1->eta(), 
			cand1->phi(),
			cand2->eta(),
			cand2->phi() ); 
}

float IDNtuplizer::deltaR2( reco::CandidatePtr& sig,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( sig->eta(), 
			sig->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

float IDNtuplizer::deltaR2( reco::TrackPtr& trk,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( trk->eta(), 
			trk->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::matchElectron( reco::CandidatePtr& sig,
				 std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				 edm::Ptr<T>& cand_ptr, 
				 float& cand_dr, 
				 bool& cand_match ) { 
  
  // Identify GEN electron matches to candidate
  auto match = std::find_if( sig2cand.begin(), 
			     sig2cand.end(), 
			     [sig](const DeltaR2<reco::Candidate,T>& dr2) { 
			       return dr2.ptr1_ == sig; 
			     }
			     );
  if ( match != sig2cand.end() ) {
    cand_ptr = match->ptr2_; // pass by ref
    if ( validPtr(cand_ptr) ) {
      if ( deltaR2(sig,cand_ptr) >= 0. ) { 
	cand_dr = sqrt(deltaR2(sig,cand_ptr)); // pass by ref
      }
      cand_match = ( cand_dr >= 0. ) && ( cand_dr < dr_threshold_ ); // pass by ref
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::trkToGsfLinks( edm::Handle< std::vector<reco::Track> >& ctfTracks,
				 edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 std::vector<TrkToGsfR2>& trk2gsf ) {
  
  // If not, use deltaR match between GsfTrack and "surrogate" Track
  
  std::vector<TrkToGsfR2> trk2gsf_all;
  trk2gsf_all.reserve( ctfTracks->size()*gsfTracks->size() );
  for ( size_t itrk = 0; itrk < ctfTracks->size(); ++itrk ) {
    reco::TrackPtr trk(ctfTracks, itrk);
    for ( size_t igsf = 0; igsf < gsfTracks->size(); ++igsf ) {
      reco::GsfTrackPtr gsf(gsfTracks, igsf);
      reco::TrackPtr ptr;
      if ( gsfToTrk(gsf,ptr) && ptr == trk ) { 
	trk2gsf_all.emplace_back( trk, gsf, IDNtuple::NEG_FLOAT ); // match via ElectronSeed
      } else { 
	trk2gsf_all.emplace_back( trk, gsf, deltaR2(trk,gsf) ); // match via deltaR
      }
    }
  }
  
  // Sort by DeltaR2!! (Matches via ElectronSeed have deltaR2 of IDNtuple::NEG_FLOAT)
  std::sort( trk2gsf_all.begin(), 
	     trk2gsf_all.end(), 
	     TrkToGsfR2::compare_by_dr2 );
  
  // Select best matches according to DeltaR2 metric (i.e. keeps matches via ElectronSeed)
  trk2gsf.clear();
  for ( auto iter : trk2gsf_all ) {
    auto trk = std::find_if( trk2gsf.begin(), 
			     trk2gsf.end(), 
			     [iter](const TrkToGsfR2& dr2) { 
			       return dr2.ptr1_ == iter.ptr1_; 
			     }
			     );
    auto gsf = std::find_if( trk2gsf.begin(), 
			     trk2gsf.end(), 
			     [iter](const TrkToGsfR2& dr2) { 
			       return dr2.ptr2_ == iter.ptr2_; 
			     }
			     );
    if ( trk == trk2gsf.end() && // Check trk not already used
	 gsf == trk2gsf.end() ) { // Check gsf not already used
      if ( iter.dr2_ < 0. ) {
	// For matches via ElectronSeed, update to correct deltaR2
	trk2gsf.emplace_back(iter.ptr1_,
			     iter.ptr2_,
			     deltaR2(iter.ptr1_,iter.ptr2_)); 
      } else {
	// For matches using "surrogate" Tracks, update to IDNtuple::NEG_FLOAT
	trk2gsf.emplace_back(iter.ptr1_,
			     iter.ptr2_,
			     IDNtuple::NEG_FLOAT); 
      }
    }
    if ( trk2gsf.size() >= ctfTracks->size() ) { break; } // found unique match for all tracks
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::trkToEleLinks( edm::Handle< std::vector<reco::Track> >& ctfTracks,
				 edm::Handle< std::vector<reco::GsfElectron> >& gsfElectrons,
				 std::vector<TrkToEleR2>& trk2ele ) {
  
  // If not, use deltaR match between GsfElectron and "surrogate" Track
  
  std::vector<TrkToEleR2> trk2ele_all;
  trk2ele_all.reserve( ctfTracks->size()*gsfElectrons->size() );
  for ( size_t itrk = 0; itrk < ctfTracks->size(); ++itrk ) {
    reco::TrackPtr trk(ctfTracks, itrk);
    for ( size_t iele = 0; iele < gsfElectrons->size(); ++iele ) {
      reco::GsfElectronPtr ele(gsfElectrons, iele);
      reco::TrackPtr ptr;
      if ( eleToTrk(ele,ptr) && ptr == trk ) { 
	trk2ele_all.emplace_back( trk, ele, IDNtuple::NEG_FLOAT ); // match via ElectronSeed
      } else { 
	trk2ele_all.emplace_back( trk, ele, deltaR2(trk,ele) ); // match via deltaR
      }
    }
  }
  
  // Sort by DeltaR2!! (Matches via ElectronSeed have deltaR2 of IDNtuple::NEG_FLOAT)
  std::sort( trk2ele_all.begin(), 
	     trk2ele_all.end(), 
	     TrkToEleR2::compare_by_dr2 );
  
  // Select best matches according to DeltaR2 metric (i.e. keeps matches via ElectronSeed)
  trk2ele.clear();
  for ( auto iter : trk2ele_all ) {
    auto trk = std::find_if( trk2ele.begin(), 
			     trk2ele.end(), 
			     [iter](const TrkToEleR2& dr2) { 
			       return dr2.ptr1_ == iter.ptr1_; 
			     }
			     );
    auto ele = std::find_if( trk2ele.begin(), 
			     trk2ele.end(), 
			     [iter](const TrkToEleR2& dr2) { 
			       return dr2.ptr2_ == iter.ptr2_; 
			     }
			     );
    if ( trk == trk2ele.end() && // Check trk not already used
	 ele == trk2ele.end() ) { // Check ele not already used
      if ( iter.dr2_ < 0. ) {
	// For matches via ElectronSeed, update to correct deltaR2
	trk2ele.emplace_back(iter.ptr1_,
			     iter.ptr2_,
			     deltaR2(iter.ptr1_,iter.ptr1_)); 
      } else {
	// For matches using "surrogate" Tracks, update to IDNtuple::NEG_FLOAT
	trk2ele.emplace_back(iter.ptr1_,
			     iter.ptr2_,
			     IDNtuple::NEG_FLOAT); 
      }
    }
    if ( trk2ele.size() >= ctfTracks->size() ) { break; } // found unique match for all tracks
  }

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 
void IDNtuplizer::gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 const edm::Handle< std::vector<reco::GsfElectron> >& gsfElectrons,
				 std::vector<GsfToEleR2>& gsf2ele ) {
  
  gsf2ele.clear();
  
  // Iterate through GsfElectrons
  for ( size_t idx = 0; idx < gsfElectrons->size(); ++idx ) {
    reco::GsfElectronPtr ele(gsfElectrons, idx);
    if ( !validPtr(ele) ) { 
      std::cout << "TEST" << std::endl;
      continue; 
    } //@@ shouldn't ever happen?!
    
    // Retrieve GsfTrack
    reco::GsfTrackPtr gsf;
    if ( !eleToGsf(ele,gsf) ) { continue; } // GsfElectron not found
    
    // Check if collections match
    if ( gsf.id() != gsfTracks.id() ) { 
      std::cout << "ERROR! GsfTrack collections do not match!:"
		<< " gsf.id(): " << gsf.id()
		<< " gsfTracks.id(): " << gsfTracks.id()
		<< std::endl;
      continue;
    } 
    
    // Check if already stored in map
    auto match = std::find_if( gsf2ele.begin(), 
			       gsf2ele.end(), 
			       [gsf](const GsfToEleR2& dr2) { 
				 return dr2.ptr1_ == gsf; 
			       }
			       );
    if ( match == gsf2ele.end() ) {
      gsf2ele.emplace_back( gsf, ele, IDNtuple::NEG_FLOAT); // do not set deltaR2
    } else { std::cout << "ERROR! GsfTrackPtr is already mapped a GsfElectronPtr!" << std::endl; }
    
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf ) {
  if ( !validPtr(ele) ) {
    std::cout << "ERROR! GsfElectronPtr:"
	      << " ele.isNull(): " << ele.isNull()
	      << " ele.isAvailable(): " << ele.isAvailable()
	      << std::endl;
    return false;
  }
  gsf = edm::refToPtr(ele->gsfTrack());
  if ( !validPtr(gsf) ) {
    std::cout << "ERROR! GsfTrackPtr:"
	      << " gsf.isNull(): " << gsf.isNull()
	      << " gsf.isAvailable(): " << gsf.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    std::cout << "ERROR! GsfTrackPtr:"
	      << " gsf.isNull(): " << gsf.isNull()
	      << " gsf.isAvailable(): " << gsf.isAvailable()
	      << std::endl;
    return false;
  }
  edm::RefToBase<TrajectorySeed> traj = gsf->seedRef();
  if ( traj.isNull() || !traj.isAvailable() ) { 
    std::cout << "ERROR: TrajectorySeedRef:" 
	      << " traj.isNull(): " << traj.isNull()
	      << " traj.isAvailable(): " << traj.isAvailable()
	      << std::endl; 
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    std::cout << "ERROR! ElectronSeedPtr:"
	      << " seed.isNull(): " << seed.isNull()
	      << " seed.isAvailable(): " << seed.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    std::cout << "ERROR! ElectronSeedPtr:"
	      << " seed.isNull(): " << seed.isNull()
	      << " seed.isAvailable(): " << seed.isAvailable()
	      << std::endl;
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
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo ) {
  if ( !validPtr(seed) ) { 
    std::cout << "ERROR! ElectronSeedPtr:"
	      << " seed.isNull(): " << seed.isNull()
	      << " seed.isAvailable(): " << seed.isAvailable()
	      << std::endl;
    return false;
  }

  edm::RefToBase<reco::CaloCluster> base = seed->caloCluster();
  if ( base.isNull() || !base.isAvailable() ) { 
    if ( verbose_ > 3 ) {
      std::cout << "INFO! edm::RefToBase<reco::CaloCluster>:"
		<< " base.isNull(): " << base.isNull()
		<< " base.isAvailable(): " << base.isAvailable()
		<< std::endl;
    }
    return false;
  }
  calo = edm::refToPtr(seed->caloCluster().castTo<reco::CaloClusterRef>());
  if ( !validPtr(calo) ) { 
    std::cout << "ERROR! CaloClusterPtr:"
	      << " calo.isNull(): " << calo.isNull()
	      << " calo.isAvailable(): " << calo.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk ) {
  reco::ElectronSeedPtr seed;
  if ( gsfToSeed(gsf,seed) && seedToTrk(seed,trk) ) { return true; }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::eleToTrk( reco::GsfElectronPtr& ele, reco::TrackPtr& trk ) {
  reco::GsfTrackPtr gsf;
  reco::ElectronSeedPtr seed;
  if ( eleToGsf(ele,gsf) && 
       gsfToSeed(gsf,seed) && 
       seedToTrk(seed,trk) ) { return true; }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDNtuplizer);




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//

//  void signalElectrons1( bool is_egamma,
//			 std::set<reco::CandidatePtr>& signal_electrons,
//			 std::vector<SigToTrkR2>& sig2trk,
//			 std::vector<SigToGsfR2>& sig2gsf,
//			 std::vector<TrkToGsfR2>& trk2gsf,
//			 std::vector<GsfToEleR2>& gsf2ele );

//  signalElectrons1( false,
//		    signal_electrons, 
//		    sig2trk,
//		    sig2gsf,
//		    trk2gsf, 
//		    gsf2ele );

//void IDNtuplizer::signalElectrons1( bool is_egamma,
//				    std::set<reco::CandidatePtr>& signal_electrons,
//				    std::vector<SigToTrkR2>& sig2trk,
//				    std::vector<SigToGsfR2>& sig2gsf,
//				    std::vector<TrkToGsfR2>& trk2gsf,
//				    std::vector<GsfToEleR2>& gsf2ele ) {
//  
//  // Iterate through "signal" electrons
//  for ( auto sig : signal_electrons ) {
//    
//    // Create ElectronChain object and start to populate
//    chains_.push_back(ElectronChain());
//    ElectronChain& chain = chains_.back();
//    chain.is_egamma_ = true;
//    chain.is_e_ = true;
//    chain.is_mc_ = isMC_;
//    if ( isMC_ ) { chain.sig_ = sig; }
//    
//    // Find first instance of "signal" electron match (with smallest DeltaR2, as sorted) 
//    auto match = std::find_if( sig2gsf.begin(), 
//			       sig2gsf.end(), 
//			       [sig](const SigToGsfR2& dr2) { 
//				 return dr2.ptr1_ == sig; 
//			       }
//			       );
//
//    // Find match: "signal" electron to GsfElectron 
//    if ( match != sig2gsf.end() ) {
//      
//      // Create and check GsfElectronPtr
//      reco::GsfTrackPtr gsf = match->ptr2_;
//      if ( gsf.isNonnull() && gsf.isAvailable() ) {
//	
//	// Store gsfTrack info
//	chain.gsf_ = gsf;
//	double dr2 = reco::deltaR2(gsf->etaMode(),
//				   gsf->phiMode(),
//				   sig->eta(), 
//				   sig->phi());
//	chain.gsf_dr_ = dr2 < 0. ? IDNtuple::NEG_FLOAT : sqrt(dr2);
//	chain.gsf_match_ = ( chain.gsf_dr_ >= 0. ) && ( chain.gsf_dr_ < dr_threshold_ );
//	
//	// Check if GsfTrack matched to GsfElectron
//	auto match = std::find_if( gsf2ele.begin(), 
//				   gsf2ele.end(), 
//				   [gsf](const GsfToEleR2& dr2) { 
//				     return dr2.ptr1_ == gsf; 
//				   }
//				   );
//
//	// Check match and store electron info
//	if ( match != gsf2ele.end() ) {
//	  reco::GsfElectronPtr ele = match->ptr2_;
//	  chain.ele_ = ele;
//	  chain.ele_dr_ = match->dr2_ < 0. ? IDNtuple::NEG_FLOAT : sqrt(match->dr2_);
//	  chain.ele_match_ = ( chain.ele_dr_ >= 0. ) && ( chain.ele_dr_ < dr_threshold_ );
//	}
//
//	// Create and check ElectronSeedPtr
//	reco::ElectronSeedPtr seed;
//	if ( gsfToSeed(gsf,seed) ) {
//	  
//	  // Store ElectronSeed info
//	  chain.seed_ = seed;
//	  chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	  chain.seed_ecal_driven_ = seed->isEcalDriven();
//	  
//	  // Create, check, and store TrackPtr info
//	  reco::TrackPtr trk;
//	  if ( seedToTrk(seed,trk) ) { 
//	    chain.trk_ = trk; 
//	    chain.trk_match_ = ( chain.trk_dr_ >= 0. ) && ( chain.trk_dr_ < dr_threshold_ );
//	  }
//	  
//	  // Create, check, and store CaloClusterPtr info
//	  reco::CaloClusterPtr calo;
//	  if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
//	  
//	  // If ECAL-driven _only_, then match to "surrogate" track
//	  if ( !seed->isTrackerDriven() && seed->isEcalDriven() ) {
//	    //std::cout << "Surrogate tracks!" << std::endl; // !!! SURROGATE TRACKS !!!
//	    auto match = std::find_if( trk2gsf.begin(), 
//				       trk2gsf.end(), 
//				       [gsf](const TrkToGsfR2& dr2) { 
//					 return dr2.ptr2_ == gsf; 
//				       }
//				       );
//	    if ( match != trk2gsf.end() ) { 
//	      chain.trk_ = match->ptr1_; 
//	      chain.trk_match_ = ( chain.trk_dr_ >= 0. ) && ( chain.trk_dr_ < dr_threshold_ );
//	    }
//	  }
//	  
//	} // if ( gsfToSeed(gsf,seed) )
//      } // if ( gsf.isNonnull() && gsf.isAvailable() )
//    } // if ( match != sig2gsf.end() )
//    else { 
//      
////      // Find first instance of "signal" electron match (with smallest DeltaR2, as sorted) 
////      auto match = std::find_if( sig2ele.begin(), 
////      				 sig2ele.end(), 
////      				 [sig](const DeltaR2& dr2) { return dr2.idx1_ == sig.key(); }
////      				 );
////      // Find match: "signal" electron to GsfElectron 
////      if ( match != sig2ele.end() ) {
////      }
//      
//    }
//
//  } // for ( auto iter : signal_electrons )
//  
//}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//

//void matchGsfTrack( TrkToGsfR2& match, ElectronChain& chain );
  
//void IDNtuplizer::matchGsfTrack( TrkToGsfR2& match, 
//				 ElectronChain& chain ) {
//    
//  // Create and check GsfTrackPtr
//  reco::GsfTrackPtr gsf = match.ptr2_;
//  if ( validPtr(gsf) ) {
//    
//    // Store gsfTrack info
//    chain.gsf_ = gsf;
//    chain.gsf_match_ = true;
//    
//    // If dr2 > 0. then match is made via ElectronSeed (see trkToGsf() method)
//    if ( match.dr2_ >= 0. ) {
//      
//      // Create and check ElectronSeedPtr
//      reco::ElectronSeedPtr seed;
//      if ( gsfToSeed(gsf,seed) ) {
//	
//	// Store ElectronSeed info
//	chain.seed_ = seed;
//	chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	chain.seed_ecal_driven_ = seed->isEcalDriven();
//	
//	// Create, check, and store CaloClusterPtr info
//	reco::CaloClusterPtr calo;
//	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
//	
//      }
//      
//    } else { // made is made via deltaR2 (see trkToGsf() method)
//      
//      // Store ElectronSeed info
//      chain.seed_ecal_driven_ = true;
//      
//    }
//    
//  } // if ( validPtr(gsf) )
//  else {
//    std::cout << "ERROR! GsfTrackPtr:"
//	      << " GsfTrackPtr: " << gsf.id() << "/" << gsf.key()
//	      << std::endl;
//  }
//
//}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//void IDNtuplizer::matchGsfElectron( TrkToEleR2& match, 
//				    ElectronChain& chain ) {
//    
//  // Create and check GsfElectronPtr
//  reco::GsfElectronPtr gsf = match.ptr2_;
//  if ( gsf.isNonnull() && gsf.isAvailable() ) {
//    
//    // Store gsfElectron info
//    chain.ele_ = ele;
//    chain.ele_match_ = true;
//    
//    reco::GsfTrackPtr gsf;
//    if ( !eleToGsf(ele,gsf) ) { continue; }
//    
//    // Store gsfElectron info
//    chain.gsf_ = gsf;
//    chain.gsf_match_ = true;
//    
//    // If dr2 > 0. then match is made via ElectronSeed (see trkToGsf() method)
//    if ( match.dr2_ >= 0. ) {
//      
//      // Create and check ElectronSeedPtr
//      reco::ElectronSeedPtr seed;
//      if ( gsfToSeed(gsf,seed) ) {
//	
//	// Store ElectronSeed info
//	chain.seed_ = seed;
//	chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	chain.seed_ecal_driven_ = seed->isEcalDriven();
//	
//	// Create, check, and store CaloClusterPtr info
//	reco::CaloClusterPtr calo;
//	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; }
//	
//      }
//      
//    } else { // made is made via deltaR2 (see trkToGsf() method)
//      
//      // Store ElectronSeed info
//      chain.seed_ecal_driven_ = true;
//      
//    }
//    
//  } // if ( gsf.isNonnull() && gsf.isAvailable() )
//  else {
//    std::cout << "ERROR! GsfElectronPtr:"
//	      << " GsfElectronPtr: " << gsf.id() << "/" << gsf.key()
//	      << std::endl;
//  }
//
//}
