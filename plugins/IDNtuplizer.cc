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
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDFilter.h" // EDAnalyzer.h
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"
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
#include <algorithm>
#include <random>

namespace reco { typedef edm::Ref<CaloClusterCollection> CaloClusterRef; }
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<PreId> PreIdPtr; }
namespace pat { typedef edm::Ptr<PackedCandidate> PackedCandidatePtr; }
typedef std::map<unsigned long,int> PdgIds;
typedef std::pair<float,float> TagMuon;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
template <class T1, class T2> 
class DeltaR2 {
public:
  DeltaR2( const edm::Ptr<T1>& obj1, const edm::Ptr<T2>& obj2, const double dr2 ) {
    obj1_ = obj1; obj2_ = obj2; dr2_ = dr2;
  };
  edm::Ptr<T1> obj1_;
  edm::Ptr<T2> obj2_;
  double dr2_ = IDNtuple::NEG_FLOATSQ; // Use ^2 because this is deltaR^2
  static bool compare_by_dr2( const DeltaR2<T1,T2>& a, const DeltaR2<T1,T2>& b ) {
    return a.dr2_ < b.dr2_;
  };
};

template <class T1, class T2> 
std::ostream& operator<< ( std::ostream& out, const DeltaR2<T1,T2>& obj ) {
  out << "Class type:        " << boost::core::demangle( typeid(obj).name() ) << "\n"
      << "  Obj1 type:       " << boost::core::demangle( typeid(obj.obj1_).name() ) << "\n"
      << "  Obj2 type:       " << boost::core::demangle( typeid(obj.obj2_).name() ) << "\n";
  out << "  Obj1 id/key:     "; 
  if ( obj.obj1_.isNull() || !obj.obj1_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj1_.id() << "/" << obj.obj1_.key(); }
  out << "\n";
  out << "  Obj2 id/key:     ";
  if ( obj.obj2_.isNull() || !obj.obj2_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj2_.id() << "/" << obj.obj2_.key(); }
  out << "\n";
  out << "  Obj1 pt/eta/phi: ";
  if ( obj.obj1_.isNull() || !obj.obj1_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj1_->pt()  << "/" 
	     << obj.obj1_->eta() << "/"
	     << obj.obj1_->phi(); } 
  out  << "\n";
  out << "  Obj2 pt/eta/phi: ";
  if ( obj.obj2_.isNull() || !obj.obj2_.isAvailable() ) { out << "InvalidKey"; }
  else { out << obj.obj2_->pt()  << "/" 
	     << obj.obj2_->eta() << "/"
	     << obj.obj2_->phi(); } 
  out  << "\n";
  out << "  dR2(1,2):        " << obj.dr2_;
  return out;
};

typedef DeltaR2<reco::Candidate,reco::Track> SigToTrkDR2;
typedef DeltaR2<reco::Candidate,reco::GsfTrack> SigToGsfDR2;
typedef DeltaR2<reco::Candidate,reco::GsfElectron> SigToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfTrack> TrkToGsfDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfElectron> GsfToEleDR2;
typedef DeltaR2<reco::Track,reco::GsfElectron> TrkToEleDR2;
typedef DeltaR2<reco::GsfTrack,reco::GsfTrack> GsfToGsfDR2;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
class ElectronChain {

public :

  explicit ElectronChain() {;}
  ~ElectronChain() {;}

public:
  
  int is_mc_ = -1;
  int is_aod_ = -1;

  bool is_e_ = false;
  bool is_egamma_ = false;

  float tag_pt_ = IDNtuple::NEG_FLOAT;
  float tag_eta_ = IDNtuple::NEG_FLOAT;

  // "Signal electron" info 
  reco::CandidatePtr sig_;

  // Track info
  reco::TrackPtr trk_;
  float trk_dr_ = IDNtuple::NEG_FLOAT;
  bool trk_match_ = false;
  int pdg_id_ = 0;

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

  // PF GsfTrack info
  reco::GsfTrackPtr pfgsf_;
  float pfgsf_dr_ = IDNtuple::NEG_FLOAT;
  bool pfgsf_match_ = false;

  // GsfElectron info
  reco::GsfElectronPtr ele_;
  float ele_dr_ = IDNtuple::NEG_FLOAT;
  bool ele_match_ = false;
  float id_ = IDNtuple::NEG_FLOAT;

  // Inner GSF track P4: defines reference coords in eta-phi space 

  float gsf_ref_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_ref_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_ref_R_ = IDNtuple::NEG_FLOAT;
  float gsf_ref_p_ = IDNtuple::NEG_FLOAT;
  float gsf_ref_pt_ = IDNtuple::NEG_FLOAT;

  // GEN channel (2 points: inner P4 + inner P4 extrapolated to ECAL)

  float gen_inner_eta_ = IDNtuple::NEG_FLOAT;
  float gen_inner_phi_ = IDNtuple::NEG_FLOAT;
  float gen_inner_R_ = IDNtuple::NEG_FLOAT;
  float gen_inner_p_ = IDNtuple::NEG_FLOAT;
  float gen_inner_pt_ = IDNtuple::NEG_FLOAT;

  float gen_proj_eta_ = IDNtuple::NEG_FLOAT;
  float gen_proj_phi_ = IDNtuple::NEG_FLOAT;
  float gen_proj_R_ = IDNtuple::NEG_FLOAT;

  // GSF channel (3 points: inner P4 + inner P4 extrapolated to ECAL + outer P4 at ECAL surface)

  float gsf_inner_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_inner_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_inner_R_ = IDNtuple::NEG_FLOAT;
  float gsf_inner_p_ = IDNtuple::NEG_FLOAT;
  float gsf_inner_pt_ = IDNtuple::NEG_FLOAT;

  float gsf_proj_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_proj_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_proj_R_ = IDNtuple::NEG_FLOAT;
  float gsf_proj_p_ = IDNtuple::NEG_FLOAT;

  float gsf_atcalo_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_atcalo_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_atcalo_R_ = IDNtuple::NEG_FLOAT;
  float gsf_atcalo_p_ = IDNtuple::NEG_FLOAT;

  // Cluster constituents of SC
  std::vector<float> clu_eta_;
  std::vector<float> clu_phi_;
  std::vector<float> clu_e_;
  std::vector<int> clu_nhit_;

  // PF candidates
  std::vector<float> pf_eta_;
  std::vector<float> pf_phi_;
  std::vector<float> pf_p_;
  std::vector<int> pf_pdgid_;
  std::vector<int> pf_matched_;
  std::vector<int> pf_lost_;

};

template <typename T> 
void key_id( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  ss << "id/key: ";
  if ( ptr.isNull() || !ptr.isAvailable() ) { ss << "InvalidKey"; }
  else { 
    std::stringstream tmp; tmp << ptr.id();
    ss << std::setw(6) << tmp.str() << "/"
       << std::setw(3) << int( ptr.key() ); }
}

void pt_eta_phi( const edm::Ptr<reco::GsfTrack>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->ptMode() << ", " 
     << std::setw(4) << ptr->etaMode() << ", " 
     << std::setw(4) << ptr->phiMode();
}

template <typename T> 
void pt_eta_phi( const edm::Ptr<T>& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", pt/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->pt() << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();
}

void pt_eta_phi( const reco::CaloClusterPtr& ptr, std::stringstream& ss ) {
  if ( ptr.isNull() || !ptr.isAvailable() ) { return; }
  ss << ", Et/eta/phi: " 
     << std::fixed
     << std::setprecision(2) 
     << std::setw(5) << ptr->energy() / std::cosh(ptr->eta()) << ", " 
     << std::setw(4) << ptr->eta() << ", " 
     << std::setw(4) << ptr->phi();
}

std::ostream& operator<< ( std::ostream& out, const ElectronChain& obj ) {
  std::stringstream ss;
  ss << "ElectronChain:"
     << " is_egamma: " << obj.is_egamma_
     << " is_e: " << obj.is_e_
     << " is_mc: " << obj.is_mc_
     << " is_aod: " << obj.is_aod_
     << " tag pt/eta: " 
     << std::setprecision(2) << obj.tag_pt_ << "/"
     << std::setprecision(2) << obj.tag_eta_;
  ss << "\n SIG:   "; key_id(obj.sig_,ss); pt_eta_phi(obj.sig_,ss);
  ss << "\n TRK:   "; key_id(obj.trk_,ss); pt_eta_phi(obj.trk_,ss); 
  ss << ", PdgId: " << obj.pdg_id_;
  ss << "\n CALO:  "; key_id(obj.calo_,ss); pt_eta_phi(obj.calo_,ss);
  ss << "\n GSF:   "; key_id(obj.gsf_,ss); pt_eta_phi(obj.gsf_,ss);
  ss << "\n PFGSF: "; key_id(obj.pfgsf_,ss); pt_eta_phi(obj.pfgsf_,ss);
  ss << "\n ELE:   "; key_id(obj.ele_,ss); pt_eta_phi(obj.ele_,ss);
  ss << "\n SEED:  "; key_id(obj.seed_,ss);
  ss << ", TRK driven: " << obj.seed_tracker_driven_
     << ", ECAL driven: " << obj.seed_ecal_driven_;
  //ss << "\n PREID: "; key_id(obj.preid_ecal_,ss);
  ss << "\n BDTs:  " << "unbiased: " << std::setprecision(4) << obj.unbiased_ 
     << ", ptbiased: " << std::setprecision(4) << obj.ptbiased_
     << ", ID: " << std::setprecision(4) << obj.id_;
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
class IDNtuplizer : public edm::EDFilter { // edm::EDAnalyzer
  
public:
  
  explicit IDNtuplizer( const edm::ParameterSet& );
  ~IDNtuplizer();
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual bool filter( edm::Event&, const edm::EventSetup& ) override; // analyze(const,const)
  
  ///////////////
  // Main methods
  ///////////////

  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );
  
  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::CandidatePtr>& signal_electrons,
			std::set<reco::CandidatePtr>& tag_side_muons );
  
  // GEN-based method to provide a sample of "signal" electrons
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B, 
			  std::set<reco::GenParticlePtr>& gen_muons, 
			  float tag_muon_pt_threshold = 5.,
			  float tag_muon_eta_threshold = 2.5 );
  
  // 
  void createChains( std::set<reco::CandidatePtr>& signal_electrons,
		     std::set<reco::CandidatePtr>& tag_side_muons );

  // Method that creates ElectronChain objects for both low-pT and PF electron signal candiates
  void signal( std::set<reco::CandidatePtr>& signal_electrons,
	       std::set<reco::CandidatePtr>& tag_side_muons,
	       std::vector<SigToTrkDR2>& sig2trk,
	       std::vector<SigToGsfDR2>& sig2gsf,
	       std::vector<SigToGsfDR2>& sig2pfgsf,
	       std::vector<SigToTrkDR2>& other_trk,
	       std::vector<TrkToGsfDR2>& trk2gsf,
	       std::vector<TrkToGsfDR2>& trk2pfgsf,
	       std::vector<GsfToGsfDR2>& gsf2pfgsf,
	       std::vector<GsfToEleDR2>& gsf2ele,
	       std::vector<GsfToEleDR2>& pfgsf2ele );

  // Method that creates ElectronChain objects for both low-pT and PF electron bkgd candidates
  void bkgd( std::set<reco::CandidatePtr>& signal_electrons,
	     std::set<reco::CandidatePtr>& tag_side_muons,
	     std::vector<SigToTrkDR2>& sig2trk,
	     std::vector<SigToGsfDR2>& sig2gsf,
	     std::vector<SigToGsfDR2>& sig2pfgsf,
	     std::vector<SigToTrkDR2>& other_trk,
	     std::vector<TrkToGsfDR2>& trk2gsf,
	     std::vector<TrkToGsfDR2>& trk2pfgsf,
	     std::vector<GsfToGsfDR2>& gsf2pfgsf,
	     std::vector<GsfToEleDR2>& gsf2ele,
	     std::vector<GsfToEleDR2>& pfgsf2ele );
  
  // Fills tree per ElectronChain object
  void fill( const edm::Event& event, const edm::EventSetup& setup );

  //////////////////
  // Utility methods
  //////////////////

  // Extracts TrackPtrs to common vector, obtained from AOD or mAOD collections via Handles
  void extractTrackPtrs();
  
  // Various navigation methods (high- to low-level objects)
  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk, bool is_egamma = false );
  bool eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  bool seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo );

  // Various navigation maps (low- to high-level objects)
  void trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
		      edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::vector<TrkToGsfDR2>& trk2gsf,
		      bool is_egamma = false );
  void gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
		      std::vector<GsfToEleDR2>& gsf2ele );
  void gsfToPfGsfLinks( edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
			edm::Handle< std::vector<reco::GsfTrack> >& gsfTracksEGamma,
			std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Links "signal" electrons to reconstructed objects
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector< edm::Ptr<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Wraps method above to allow use of Handle<View<T>>
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       edm::Handle< edm::View<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Wraps method above to allow use of Handle<vector<T>>
  template <typename T> 
  void sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
		       edm::Handle< std::vector<T> >& candidates,
		       std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
		       std::vector< DeltaR2<reco::Candidate,T> >& other_cand, 
		       bool append = false );

  // Return by reference the 'cand' that is matched to 'sig' in 'sig2cand' map
  template <typename T>
  void match( reco::CandidatePtr& sig,
	      std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
	      edm::Ptr<T>& cand, float& dr, bool& match ); // pass by ref

  // Wrap reco::deltaR2 method for various types, specifically GsfTrack uses eta and phi 'Mode'
  template <typename T1, typename T2> 
  float deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ); // used by sigToCandLinks
  float deltaR2( reco::TrackPtr& trk, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::GsfTrackPtr& pfgsf, reco::GsfTrackPtr& gsf ); // overload
  float deltaR2( reco::CandidatePtr& sig, reco::GsfTrackPtr& gsf ); // overload

  // Filter track candidates by quality flag, simple pass for all other cands
  template <typename T> 
  bool filterCand( edm::Ptr<T>& cand );
  bool filterCand( edm::Ptr<reco::Track>& trk );
  bool filterCand( edm::Ptr<reco::GsfTrack>& gsf );

  // Check is Ptr is valid and available
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );


  // Electron "images"
  typedef math::XYZVector Vector;
  typedef math::XYZPoint Point;
  BaseParticlePropagator extrapolate_track( const Vector& mom, const Point& pos, int charge,
					    int& reach_ECAL, GlobalPoint& pos_ECAL,
					    int& reach_HCAL, GlobalPoint& pos_HCAL,
					    int& reach_EXIT, GlobalPoint& pos_EXIT );
  //void debug_image( const edm::Event& event, const edm::EventSetup& setup );
  void build_image( const edm::Event& event, const edm::EventSetup& setup );

  inline float adj_eta( float eta, float ref_eta, int charge = 0 ) { 
    return eta-ref_eta; 
    //return eta; //@@ DON'T MODIFY ETA !!!
  }
  inline float adj_phi( float phi, float ref_phi, int charge = 0 ) { 
    return charge == 0 ? reco::deltaPhi(phi,ref_phi) : float(charge)*reco::deltaPhi(phi,ref_phi); 
    //return charge == 0 ? phi : float(charge)*phi; //@@ DON'T SUBTRACT REF_PHI !!!
    //return phi; //@@ DON'T MODIFY PHI !!!
  }
  
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

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_; // AOD
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;
 
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > packedCands_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > packedCandsH_;
 
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > lostTracks_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > lostTracksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > pfToPackedCands_; // MINIAOD
  edm::Handle< edm::Association<pat::PackedCandidateCollection> > pfToPackedCandsH_;
  
  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;

  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;

  //noZS::EcalClusterLazyTools ecalTools_;
  
  const edm::EDGetTokenT<reco::SuperClusterCollection> barrelSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> barrelSCsH_;

  const edm::EDGetTokenT<reco::SuperClusterCollection> endcapSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> endcapSCsH_;

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

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;

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

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth10_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth10H_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPtDepth15_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtDepth15H_;

  // EGamma collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeedsEGamma_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsEGammaH_; // AOD

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_; // AOD
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_MAOD_; // MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksEGammaH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectronsEGamma_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectronsEGamma_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGammaRetrained_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaRetrainedH_;

  // Conversions

  //@@ const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;

  std::vector<ElectronChain> chains_;

  std::vector<reco::TrackPtr> tracks_;
  PdgIds pdgids_;

  //////////////////
  // Obsolete methods
  //////////////////

  // Top-level method that creates ElectronChain objects for low-pT electrons
  void lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		       std::vector<SigToTrkDR2>& sig2trk,
		       std::vector<SigToTrkDR2>& other_trk,
		       std::vector<SigToGsfDR2>& sig2gsf,
		       std::vector<TrkToGsfDR2>& trk2gsf,
		       std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for low-pT electron signal candiates
  void lowPtElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			      std::vector<SigToTrkDR2>& sig2trk,
			      std::vector<SigToTrkDR2>& other_trk,
			      std::vector<SigToGsfDR2>& sig2gsf,
			      std::vector<GsfToGsfDR2>& gsf2pfgsf,
			      std::vector<GsfToEleDR2>& gsf2ele );

  // Method that creates ElectronChain objects for low-pT electron fake candidates
  void lowPtElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			     std::vector<TrkToGsfDR2>& trk2gsf,
			     std::vector<GsfToGsfDR2>& gsf2pfgsf,
			     std::vector<GsfToEleDR2>& gsf2ele );
  
  // Top-level method that creates ElectronChain objects for EGamma electrons
  void pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
		    std::vector<SigToTrkDR2>& sig2trk,
		    std::vector<SigToTrkDR2>& other_trk,
		    std::vector<SigToGsfDR2>& sig2gsf,
		    std::vector<TrkToGsfDR2>& trk2gsf,
		    std::vector<GsfToGsfDR2>& gsf2pfgsf );

  // Method that creates ElectronChain objects for EGamma signal candidates
  void pfElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
			   std::vector<SigToTrkDR2>& sig2trk,
			   std::vector<SigToTrkDR2>& other_trk,
			   std::vector<SigToGsfDR2>& sig2gsf,
			   std::vector<SigToGsfDR2>& sig2pfgsf,
			   std::vector<TrkToGsfDR2>& trk2pfgsf,
			   std::vector<GsfToGsfDR2>& gsf2pfgsf,
			   std::vector<GsfToEleDR2>& pfgsf2ele );
  
  // Method that creates ElectronChain objects for EGamma fake candidates
  void pfElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
			  std::vector<TrkToGsfDR2>& trk2gsf,
			  std::vector<TrkToGsfDR2>& trk2pfgsf,
			  std::vector<GsfToGsfDR2>& gsf2pfgsf,
			  std::vector<GsfToEleDR2>& pfgsf2ele );
  std::vector<int> pf_pdgids_;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizer::~IDNtuplizer() {
  std::cout << "pf_pdgids: ";
  for ( auto id : pf_pdgids_ ) { std::cout << id << " "; }
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
IDNtuplizer::IDNtuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
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
    genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_(),
    packedCands_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCands"))),
    packedCandsH_(),
    lostTracks_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("lostTracks"))),
    lostTracksH_(),
    pfToPackedCands_(consumes<edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("packedCands"))),
    pfToPackedCandsH_(),
    ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
    ebRecHitsH_(),
    eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
    eeRecHitsH_(),
    //ecalTools_(),
    barrelSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("barrelSuperClusters"))),
    barrelSCsH_(),
    endcapSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("endcapSuperClusters"))),
    endcapSCsH_(),
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
    mvaValueLowPtDepth10_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPtDepth10"))),
    mvaValueLowPtDepth10H_(),
    mvaValueLowPtDepth15_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPtDepth15"))),
    mvaValueLowPtDepth15H_(),
    // EGamma collections
    eleSeedsEGamma_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeedsEGamma"))),
    eleSeedsEGammaH_(),
    gsfTracksEGamma_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma"))),
    gsfTracksEGamma_MAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma_MAOD"))),
    gsfTracksEGammaH_(),
    gsfElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma"))),
    patElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectronsEGamma"))),
    gsfElectronsEGammaH_(),
    mvaValueEGamma_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGamma"))),
    mvaValueEGammaH_(),
    mvaValueEGammaRetrained_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGammaRetrained"))),
    mvaValueEGammaRetrainedH_(),
    // Conversions
    //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
    chains_(),
    tracks_(),
    pdgids_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "[IDNtuplizer::IDNtuplizer] Verbosity level: "<< verbose_ << std::endl;
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
bool IDNtuplizer::filter( edm::Event& event, const edm::EventSetup& setup ) { // analyze(const,const)

  // Update all handles - MUST be called every event!
  readCollections(event,setup);
  
  // Clear ElectronChain vector and populate
  chains_.clear();

  // Identify signal electrons, from (GEN) MC or data control regions
  std::set<reco::CandidatePtr> signal_electrons;
  std::set<reco::CandidatePtr> tag_side_muons;
  signalElectrons(signal_electrons,tag_side_muons);

  // Populate std::vector<reco::TrackPtr> tracks_ (from reco::Tracks, PF candidates, lost tracks)
  // Populate std::map<unsigned long,int> pdgids_ (typedef'ed to PdgIds)
  extractTrackPtrs();

  createChains(signal_electrons,tag_side_muons);
  
  // Populate ElectronChain objects using low-pT electrons
  //lowPtElectrons( signal_electrons, sig2trk, other_trk, sig2gsf, trk2gsf, gsf2pfgsf );

  // Populate ElectronChain objects using PF electrons
  //pfElectrons( signal_electrons, sig2trk, other_trk, sig2gsf, trk2gsf, gsf2pfgsf );

  // Debug info for electron "image"
  //debug_image(event,setup);
  build_image(event,setup);
  
  // Fill ntuple
  fill(event,setup);
  
  // Print ElectronChain objects
  if ( verbose_ > 0 ) {
    for ( auto iter : chains_ ) { 
      if ( iter.is_egamma_ || !iter.is_e_ ) { continue; } // skip PF and fakes
      std::cout << iter << std::endl; 
    }
  }

  // Filter based on ElectronChain objects with a PF GSF match
  for ( auto iter : chains_ ) { 
    if ( iter.is_egamma_ ) {
      if (  iter.is_e_                    && iter.pfgsf_match_ && ((event.id().event()%100)==0) ) { return true; }
      if ( !iter.is_e_ && iter.trk_match_ && iter.pfgsf_match_                                  ) { return true; }
    }
  }
  return false;
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "[IDNtuplizer::readCollections] File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(patElectrons_,gsfElectronsH_);
      if ( gsfElectronsH_.isValid() ) { 
	isAOD_ = 0;
	std::cout << "[IDNtuplizer::readCollections] File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << "[IDNtuplizer::readCollections] Failed to find a standard AOD or miniAOD particle collection " 
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
	std::cout << "[IDNtuplizer::readCollections] No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
      event.getByToken(prunedGenParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "[IDNtuplizer::readCollections] No GEN info found in MINIAOD data tier!" << std::endl;
      }
    }
  }

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
    event.getByToken(pfToPackedCands_,pfToPackedCandsH_);
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
  else if ( isAOD_ == 0 ) { event.getByToken(patElectronsEGamma_, gsfElectronsEGammaH_); }

  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
  event.getByToken(mvaValueLowPtDepth10_, mvaValueLowPtDepth10H_);
  event.getByToken(mvaValueLowPtDepth15_, mvaValueLowPtDepth15H_);
  event.getByToken(mvaValueEGamma_, mvaValueEGammaH_);
  event.getByToken(mvaValueEGammaRetrained_, mvaValueEGammaRetrainedH_);

  // Conversions
  //event.getByToken(convVtxFitProb_, convVtxFitProbH_);

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::signalElectrons( std::set<reco::CandidatePtr>& signal_electrons,
				   std::set<reco::CandidatePtr>& tag_side_muons ) {
  
  signal_electrons.clear();
  tag_side_muons.clear();

  if ( isMC_ ) { // Identify "signal" (GEN) electrons from B decays
    std::set<reco::GenParticlePtr> electrons_from_B;
    std::set<reco::GenParticlePtr> gen_muons;
    genElectronsFromB(electrons_from_B,gen_muons);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
    for ( auto gen : gen_muons ) { tag_side_muons.insert(gen); }
  } else { // Identify "signal" electrons from data control regions
    //@@ FOR NOW, A HACK ...
    std::set<reco::GenParticlePtr> electrons_from_B;
    std::set<reco::GenParticlePtr> gen_muons;
    genElectronsFromB(electrons_from_B,gen_muons);
    for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
    for ( auto gen : gen_muons ) { tag_side_muons.insert(gen); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// From AOD
void IDNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
				     std::set<reco::GenParticlePtr>& gen_muons,
				     float tag_muon_pt_threshold, 
				     float tag_muon_eta_threshold ) {
  
  electrons_from_B.clear();
  gen_muons.clear();
  
  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
    
    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "[IDNtuplizer::genElectronsFromB] ERROR! GenParticlePtr:"
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
    bool tag_muon = ( std::abs(gen->pdgId()) == 13 && gen->isLastCopy() 
		      && 
		      ( ( gen->numberOfMothers() >= 1 &&
			  gen->mother() &&
			  std::abs(gen->mother()->pdgId()) > 510 &&
			  std::abs(gen->mother()->pdgId()) < 546 && 
			  gen->mother()->pt() > tag_muon_pt_threshold && 
			  std::abs(gen->mother()->eta()) < tag_muon_eta_threshold ) 
			||
			( gen->numberOfMothers() >= 1 && 
			  gen->mother() &&
			  gen->mother()->numberOfMothers() >= 1 && 
			  gen->mother()->mother() &&
			  std::abs(gen->mother()->mother()->pdgId()) > 510 &&
			  std::abs(gen->mother()->mother()->pdgId()) < 546 && 
			  gen->mother()->mother()->pt() > tag_muon_pt_threshold &&
			  std::abs(gen->mother()->mother()->eta()) < tag_muon_eta_threshold ) ) );
    if ( tag_muon ) { gen_muons.insert(gen); }
    
    // is coming from a B
    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
      electrons_from_B.insert(gen);
      if ( verbose_ > 1 ) {
	std::cout << "[IDNtuplizer::genElectronsFromB] "
		  << " #signal_electrons: " << electrons_from_B.size()
		  << " resonant? " << resonant
		  << " non resonant? " << non_resonant
		  << " tag-side muon? " << !gen_muons.empty()
		  << std::endl;
      }
    }
    
  } // genParticles loop

  if ( gen_muons.empty() ) { electrons_from_B.clear(); }

}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::createChains( std::set<reco::CandidatePtr>& signal_electrons,
				std::set<reco::CandidatePtr>& tag_side_muons ) {
  
  // Match "signal" electrons to generalTracks
  std::vector<SigToTrkDR2> sig2trk;
  std::vector<SigToTrkDR2> other_trk;
  sigToCandLinks<reco::Track>( signal_electrons, 
			       tracks_, 
			       sig2trk, 
			       other_trk );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::createChains] sigToCandLinks<reco::Track>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " tracks_.size(): " << tracks_.size() << std::endl
	      << " sig2trk.size(): " << sig2trk.size() << std::endl
	      << " other_trk.size(): " << other_trk.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : sig2trk ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }
  
  // Match "signal" electrons to low-pT GsfTracks
  std::vector<SigToGsfDR2> sig2gsf;
  std::vector<SigToGsfDR2> other_gsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksH_, 
				  sig2gsf, 
				  other_gsf );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::createChains] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " sig2gsf.size(): " << sig2gsf.size() << std::endl
	      << " other_gsf.size(): " << other_gsf.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : sig2gsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match "signal" electrons to PF GsfTracks
  std::vector<SigToGsfDR2> sig2pfgsf;
  std::vector<SigToGsfDR2> other_pfgsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksEGammaH_,
				  sig2pfgsf, 
				  other_pfgsf );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::createChains] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " sig2pfgsf.size(): " << sig2pfgsf.size() << std::endl
	      << " other_pfgsf.size(): " << other_pfgsf.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : sig2pfgsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match Tracks (including surrogates) to low-pT GsfTracks
  std::vector<TrkToGsfDR2> trk2gsf;
  trkToGsfLinks( tracks_, 
		 gsfTracksH_, 
		 trk2gsf );
  if ( verbose_ > 1 ) {
    int good = 0, surrogate = 0, broken = 0, no_gsf = 0, no_trk = 0, empty = 0;
    for ( auto iter : trk2gsf ) { 
      if ( validPtr(iter.obj1_) && validPtr(iter.obj2_) ) {
	if ( iter.dr2_ >= 0. )                       { good++; } 
	else if ( std::abs(iter.dr2_-IDNtuple::NEG_FLOATSQ) > 1.e-6 ) { surrogate++; } 
	else                                         { broken++; } 
      } 
      else if ( validPtr(iter.obj1_) && !validPtr(iter.obj2_) ) { no_gsf++; }
      else if ( !validPtr(iter.obj1_) && validPtr(iter.obj2_) ) { no_trk++; }
      else                                                      { empty++; }
    }
    std::cout << "[IDNtuplizer::createChains] trkToGsfLinks:" << std::endl
	      << " tracks_.size():      " << tracks_.size() << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " trk2gsf.size():      " << trk2gsf.size() << std::endl
	      << " #good:      " << good << std::endl
	      << " #surrogate: " << surrogate << std::endl
	      << " #no_trk:    " << no_trk << std::endl
	      << " #no_gsf:    " << no_gsf << std::endl
	      << " #broken:    " << broken << std::endl
	      << " #empty:     " << empty << std::endl
	      << " #total:     " << good+surrogate+no_trk+no_gsf+broken+empty << std::endl;
  }

  // Match Tracks (including surrogates) to PF GsfTracks
  std::vector<TrkToGsfDR2> trk2pfgsf;
  trkToGsfLinks( tracks_, 
		 gsfTracksEGammaH_, 
		 trk2pfgsf,
		 true ); // is_egamma
  if ( verbose_ > 1 ) {
    int good = 0, surrogate = 0, broken = 0, no_gsf = 0, no_trk = 0, empty = 0;
    for ( auto iter : trk2pfgsf ) { 
      if ( validPtr(iter.obj1_) && validPtr(iter.obj2_) ) {
	if ( iter.dr2_ >= 0. )                       { 
	  good++; 
	  std::cout << "good " 
		    << validPtr(iter.obj1_) << " " 
		    << validPtr(iter.obj2_) << " " 
		    << iter.dr2_ << std::endl;
	} 
	else if ( std::abs(iter.dr2_-IDNtuple::NEG_FLOATSQ) > 1.e-6 ) { 
	  surrogate++; 
	  std::cout << "surrogate " 
		    << validPtr(iter.obj1_) << " " 
		    << validPtr(iter.obj2_) << " " 
		    << iter.dr2_ << std::endl;
	} 
	else                                         { broken++; } 
      } 
      else if ( validPtr(iter.obj1_) && !validPtr(iter.obj2_) ) { no_gsf++; }
      else if ( !validPtr(iter.obj1_) && validPtr(iter.obj2_) ) { no_trk++; }
      else                                                      { empty++; }
    }
    std::cout << "[IDNtuplizer::createChains] trkToPFGsfLinks:" << std::endl
	      << " tracks_.size():      " << tracks_.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " trk2pfgsf.size():      " << trk2pfgsf.size() << std::endl
	      << " #good:      " << good << std::endl
	      << " #surrogate: " << surrogate << std::endl
	      << " #no_trk:    " << no_trk << std::endl
	      << " #no_gsf:    " << no_gsf << std::endl
	      << " #broken:    " << broken << std::endl
	      << " #empty:     " << empty << std::endl
	      << " #total:     " << good+surrogate+no_trk+no_gsf+broken+empty << std::endl;
  }

  // Hack: PF GsfTracks --> low-pT GsfTracks map
  std::vector<GsfToGsfDR2> gsf2pfgsf;
  gsfToPfGsfLinks( gsfTracksH_, //@@ low-pT GSF tracks!
		   gsfTracksEGammaH_, 
		   gsf2pfgsf );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::chains] gsfToPfGsfLinks:" << std::endl
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " gsf2pfgsf.size(): " << gsf2pfgsf.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : gsf2pfgsf ) { 
	if ( iter.dr2_ >= 0. || validPtr(iter.obj2_) ) { std::cout << iter << std::endl; } 
      }
      std::cout << std::endl;
    }
  }

  // Match low-pT GsfTracks to low-pT GsfElectrons
  std::vector<GsfToEleDR2> gsf2ele;
  gsfToEleLinks( gsfTracksH_, 
		 gsfElectronsH_, 
		 gsf2ele );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::createChains] gsfToEleLinks:" << std::endl 
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " gsfElectronsH_->size(): " << gsfElectronsH_->size() << std::endl
	      << " gsf2ele.size(): " << gsf2ele.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : gsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match PF GsfTracks to PF GsfElectrons
  std::vector<GsfToEleDR2> pfgsf2ele;
  gsfToEleLinks( gsfTracksEGammaH_, 
		 gsfElectronsEGammaH_, 
		 pfgsf2ele );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::createChains] gsfToEleLinks:" << std::endl 
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl
	      << " pfgsf2ele.size(): " << pfgsf2ele.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : pfgsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  signal( signal_electrons, 
	  tag_side_muons, 
	  sig2trk, sig2gsf, sig2pfgsf,
	  other_trk, trk2gsf, trk2pfgsf,
	  gsf2pfgsf, 
	  gsf2ele, pfgsf2ele );

  bkgd( signal_electrons, 
	tag_side_muons,
	sig2trk, sig2gsf, sig2pfgsf,
	other_trk, trk2gsf, trk2pfgsf,
	gsf2pfgsf, 
	gsf2ele, pfgsf2ele );
  
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "signal electrons" 
void IDNtuplizer::signal( std::set<reco::CandidatePtr>& signal_electrons,
			  std::set<reco::CandidatePtr>& tag_side_muons,
			  std::vector<SigToTrkDR2>& sig2trk,
			  std::vector<SigToGsfDR2>& sig2gsf,
			  std::vector<SigToGsfDR2>& sig2pfgsf,
			  std::vector<SigToTrkDR2>& other_trk,
			  std::vector<TrkToGsfDR2>& trk2gsf,
			  std::vector<TrkToGsfDR2>& trk2pfgsf,
			  std::vector<GsfToGsfDR2>& gsf2pfgsf,
			  std::vector<GsfToEleDR2>& gsf2ele,
			  std::vector<GsfToEleDR2>& pfgsf2ele ) {

  // Find highest-pT tag-side muon in event and store
  float tag_pt = IDNtuple::NEG_FLOAT;
  float tag_eta = IDNtuple::NEG_FLOAT;
  for ( auto tag : tag_side_muons ) {
    if ( tag->pt() > tag_pt ) {
      tag_pt = tag->pt();
      tag_eta = tag->eta();
    }
  }
  
  // Iterate through signal electrons
  for ( auto sig : signal_electrons ) {
    
    // Repeat for two options: low pT and PF EGamma reconstruction
    for ( auto is_egamma : std::vector<bool>{ false, true } ) { 
      
      // SIG: Initialise ElectronChain object
      chains_.push_back(ElectronChain());
      ElectronChain& chain = chains_.back();
      chain.is_mc_ = isMC_;
      chain.is_aod_ = isAOD_;
      chain.is_e_ = true;
      chain.is_egamma_ = is_egamma; // Set here if low pT or EGamma!
      chain.sig_ = sig;
      chain.tag_pt_ = tag_pt;
      chain.tag_eta_ = tag_eta;

      // TRK: Store matches between "signal electron" and KF tracks
      match<reco::Track>(sig,
			 sig2trk,
			 chain.trk_, // by ref
			 chain.trk_dr_,
			 chain.trk_match_ );
      
      // GSF: Store matches between "signal electron" and GSF tracks
      match<reco::GsfTrack>(sig,
			    sig2gsf,
			    chain.gsf_, // by ref
			    chain.gsf_dr_,
			    chain.gsf_match_ );
      
      // PFGSF: Store matches between "signal electron" and PF GSF tracks
      match<reco::GsfTrack>(sig, 
			    sig2pfgsf,
			    chain.pfgsf_, // by ref 
			    chain.pfgsf_dr_,
			    chain.pfgsf_match_ );
      
      // Only for matched low-pT GSF tracks ...
      if ( is_egamma == false && chain.gsf_match_ ) { 
	
	// TRK: Update Track info
	reco::TrackPtr trk; 
	if ( gsfToTrk(chain.gsf_,trk) ) {
	  chain.trk_ = trk; 
	  chain.trk_match_ = true;
	  chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	}
	
	// SEED: Store Seed information
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.gsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	    // Store PreId info
	    //chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	    //chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
	  } else {
	    // If no TrackExtra info or maps, then set to false
	    chain.seed_tracker_driven_ = false;
	    chain.seed_ecal_driven_ = false;
	  }
	}
	
	// SEED: Store BDT discrimator outputs for low-pT ElectronSeeds
	chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
	chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

	// GSF: Info is stored if match made above

	// PFGSF: Info is stored if match made above
	
	// ELE: Store GSF electron info if match found with GSF track
	auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
					      gsf2ele.end(), 
					      [chain](const GsfToEleDR2& dr2) { 
						return chain.gsf_ == dr2.obj1_;
					      }
					      );
	if ( match_gsf_to_ele != gsf2ele.end() && 
	     validPtr(match_gsf_to_ele->obj2_) ) { 
	  chain.ele_ = match_gsf_to_ele->obj2_; 
	  chain.ele_match_ = true;
	  chain.ele_dr_ = sqrt(deltaR2(chain.ele_,chain.sig_));
	}
	
      } // !is_egamma and matched
      
      // Only for matched PF GSF tracks ...
      if ( is_egamma == true && chain.pfgsf_match_ ) { 
      
	// TRK: Update Track info (with surrogate track if necessary)
	auto match_pfgsf_to_trk = std::find_if( trk2pfgsf.begin(), 
						trk2pfgsf.end(), 
						[chain](const TrkToGsfDR2& dr2) { 
						  return chain.pfgsf_ == dr2.obj2_; 
						}
						);
	if ( match_pfgsf_to_trk == trk2pfgsf.end() &&
	     validPtr(match_pfgsf_to_trk->obj1_) ) {
	  chain.trk_ = match_pfgsf_to_trk->obj1_;
	  chain.trk_match_ = true;
	  chain.trk_dr_ = sqrt(deltaR2(chain.trk_,chain.sig_)); 
	  if ( match_pfgsf_to_trk->dr2_ > dr_threshold_*dr_threshold_ ) { chain.trk_dr_ *= -1; } //@@ -ve for surrogates!
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	}
	
	// SEED: Store Seed information
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	  } else {
	    // If no TrackExtra info in AOD, then set to false
	    chain.seed_tracker_driven_ = false;
	    chain.seed_ecal_driven_ = false;
	  }
	}

	// SEED: Store BDT discrimator outputs for low-pT ElectronSeeds
	if ( chain.gsf_match_ ) { 
	  chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
	  chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
//	} else {
//	  chain.unbiased_ = 10.;
//	  chain.ptbiased_ = 10.;
	}

	// GSF: Info is stored if match made above
	
	// PFGSF: Info is stored if match made above
	
	// ELE: Store GSF electron info if match found with PF GSF track
	auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
						pfgsf2ele.end(), 
						[chain](const GsfToEleDR2& dr2) { 
						  return chain.pfgsf_ == dr2.obj1_;
						}
						);
	if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
	     validPtr(match_pfgsf_to_ele->obj2_) ) { 
	  chain.ele_ = match_pfgsf_to_ele->obj2_; 
	  chain.ele_match_ = true;
	  chain.ele_dr_ = sqrt(deltaR2(chain.ele_,chain.sig_));
	}
	
      } // is_egamma and matched
      
    } // for ( auto is_egamma : std::vector<bool>{ false, true } ) 
  } // for ( auto sig : signal_electrons )
  
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through fake candidates
void IDNtuplizer::bkgd( std::set<reco::CandidatePtr>& signal_electrons,
			std::set<reco::CandidatePtr>& tag_side_muons,
			std::vector<SigToTrkDR2>& sig2trk,
			std::vector<SigToGsfDR2>& sig2gsf,
			std::vector<SigToGsfDR2>& sig2pfgsf,
			std::vector<SigToTrkDR2>& other_trk,
			std::vector<TrkToGsfDR2>& trk2gsf,
			std::vector<TrkToGsfDR2>& trk2pfgsf,
			std::vector<GsfToGsfDR2>& gsf2pfgsf,
			std::vector<GsfToEleDR2>& gsf2ele,
			std::vector<GsfToEleDR2>& pfgsf2ele ) {

  // Find highest-pT tag-side muon in event and store
  float tag_pt = IDNtuple::NEG_FLOAT;
  float tag_eta = IDNtuple::NEG_FLOAT;
  for ( auto tag : tag_side_muons ) {
    if ( tag->pt() > tag_pt ) {
      tag_pt = tag->pt();
      tag_eta = tag->eta();
    }
  }
  
  // Iterate through tracks
  for ( auto iter : other_trk ) {
    
    // Repeat for two options: low pT and PF EGamma reconstruction
    for ( auto is_egamma : std::vector<bool>{ true, false } ) { 

      // SIG: Initialise ElectronChain object
      chains_.push_back(ElectronChain());
      ElectronChain& chain = chains_.back();
      chain.is_mc_ = isMC_;
      chain.is_aod_ = isAOD_;
      chain.is_e_ = false;
      chain.is_egamma_ = is_egamma; // Set here if low pT or EGamma!
      chain.tag_pt_ = tag_pt;
      chain.tag_eta_ = tag_eta;
      
      // TRK: Store Track info
      chain.trk_ = iter.obj2_; 
      chain.trk_match_ = true;
      chain.trk_dr_ = 0.; // w.r.t. trk
      
      // GSF: Store GSF track info
      auto match_trk_to_gsf = std::find_if( trk2gsf.begin(), 
					    trk2gsf.end(), 
					    [chain](const TrkToGsfDR2& dr2) { 
					      return chain.trk_ == dr2.obj1_;
					    }
					    );
      if ( match_trk_to_gsf != trk2gsf.end() && 
	   validPtr(match_trk_to_gsf->obj2_) ) {
	chain.gsf_ = match_trk_to_gsf->obj2_;
	chain.gsf_match_ = true;
	if ( match_trk_to_gsf->dr2_ >= 0. ) { 
	  chain.gsf_dr_ = sqrt(match_trk_to_gsf->dr2_); 
	} else if ( match_trk_to_gsf->dr2_ > IDNtuple::NEG_FLOATSQ ) { 
	  chain.gsf_dr_ = -1.*sqrt(match_trk_to_gsf->dr2_*-1.); // -ve value for surrogate
	} else { 
	  chain.gsf_dr_ = IDNtuple::NEG_FLOAT;
	}
	//chain.gsf_dr_ = sqrt(deltaR2(chain.gsf_,chain.trk_));
      }
      
      // PFGSF: Store PF GSF track info
      auto match_trk_to_pfgsf = std::find_if( trk2pfgsf.begin(), 
					      trk2pfgsf.end(), 
					      [chain](const TrkToGsfDR2& dr2) { 
						return chain.trk_ == dr2.obj1_; 
					      }
					      );
      if ( match_trk_to_pfgsf != trk2pfgsf.end() ) {
	if ( validPtr(match_trk_to_pfgsf->obj2_) ) {
	  chain.pfgsf_ = match_trk_to_pfgsf->obj2_;
	  chain.pfgsf_match_ = true;
	  if ( match_trk_to_pfgsf->dr2_ >= 0. ) { 
	    chain.pfgsf_dr_ = sqrt(match_trk_to_pfgsf->dr2_); 
	  } else if ( match_trk_to_pfgsf->dr2_ > IDNtuple::NEG_FLOATSQ ) { 
	    chain.pfgsf_dr_ = -1.*sqrt(match_trk_to_pfgsf->dr2_*-1.); // -ve value for surrogate
	  } else { 
	    chain.pfgsf_dr_ = IDNtuple::NEG_FLOAT;
	  }
	  //chain.pfgsf_dr_ = sqrt(match_trk_to_pfgsf->dr2_); // +ve or -ve values for if real or surrogate
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	}
      }
      
      // Only for matched low-pT GSF tracks ...
      if ( is_egamma == false && chain.gsf_match_ ) { 
	  
	// TRK: Info is stored already

	// SEED: Store Seed information
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.gsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	    // Store PreId info
	    //chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	    //chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
	  } else {
	    // If no TrackExtra info or maps, then set to false
	    chain.seed_tracker_driven_ = false;
	    chain.seed_ecal_driven_ = false;
	  }
	}
	
	// SEED: Store BDT discrimator outputs for low-pT ElectronSeeds
	chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
	chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
	
	// GSF: Info is stored already

	// PFGSF: Info is stored already
	
	// ELE: Store GSF electron info if match found with GSF track
	auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
					      gsf2ele.end(), 
					      [chain](const GsfToEleDR2& dr2) { 
						return chain.gsf_ == dr2.obj1_;
					      }
					      );
	if ( match_gsf_to_ele != gsf2ele.end() && 
	     validPtr(match_gsf_to_ele->obj2_) ) { 
	  chain.ele_ = match_gsf_to_ele->obj2_; 
	  chain.ele_match_ = true;
	  chain.ele_dr_ = sqrt(match_gsf_to_ele->dr2_);
	}
	
      } // !is_egamma and matched

      // Only for matched PF GSF tracks ...
      if ( is_egamma == true && chain.pfgsf_match_ ) { 
      
	// TRK: Info is stored already

	// SEED: Store Seed information
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	  } else {
	    // If no TrackExtra info in AOD, then set to false
	    chain.seed_tracker_driven_ = false;
	    chain.seed_ecal_driven_ = false;
	  }
	}
	
	// SEED: Store BDT discrimator outputs for low-pT ElectronSeeds
	if ( chain.gsf_match_ ) { 
	  chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
	  chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
//	} else {
//	  chain.unbiased_ = 10.;
//	  chain.ptbiased_ = 10.;
	}

	// GSF: Info is stored already
	
	// PFGSF: Info is stored already
	
	// ELE: Store GSF electron info if match found with PF GSF track
	auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
						pfgsf2ele.end(), 
						[chain](const GsfToEleDR2& dr2) { 
						  return chain.pfgsf_ == dr2.obj1_;
						}
						);
	if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
	     validPtr(match_pfgsf_to_ele->obj2_) ) { 
	  chain.ele_ = match_pfgsf_to_ele->obj2_; 
	  chain.ele_match_ = true;
	  chain.ele_dr_ = sqrt(match_pfgsf_to_ele->dr2_);
	}
	
      } // is_egamma and matched
      
    } // for ( auto is_egamma : std::vector<bool>{ false, true } ) 
  } // for ( auto trk : other_trk )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::fill( const edm::Event& event,
			const edm::EventSetup& setup ) {

  using value_type = unsigned int;

  // Random seed generator
  static value_type seed = 1; //std::chrono::system_clock::now().time_since_epoch().count(); //@@ fixed seed !!
  static std::default_random_engine generator(seed);

  // Randomise the order of the EletronChain objects
  std::shuffle( chains_.begin(), chains_.end(), generator );

   // Poisson with mean = prescale_ (used when prescale_ is -ve)
  std::poisson_distribution<value_type> poisson(-1.*prescale_);

  std::vector<value_type> cands;
  std::vector<value_type> cands2;
  cands.reserve(chains_.size()); 
  cands2.reserve(chains_.size()); 

  if ( prescale_ < -1.e-9 ) { 

    // If prescale_ is -ve, select cands based on Poisson distr (prescale_ is mean)
    // 1) Store indices of fakes (separately for low pT and EGamma)
    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) { 
      if ( chains_[idx].is_e_ == false ) { 
	if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
	if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
      }
    }
    // 2) Shuffle then truncate to give 'prescale' elements
    unsigned int size = poisson(generator);
    cands.resize(size);
    cands2.resize(size);

  } else if ( prescale_ > 1.e-9 ) { 

    // If prescale_ is +ve, select cands based on flat probability prior (1./prescale_ is threshold)
    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) {
      if ( gRandom->Rndm() < (1./prescale_) ) {
	if ( chains_[idx].is_e_ == false ) { 
	  if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
	  if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
	}
      }
    }
    
  } else {

    // If prescale is 0, select all cands
    for ( unsigned int idx = 0; idx < chains_.size(); ++idx ) {
      if ( chains_[idx].is_e_ == false ) {
	if ( chains_[idx].is_egamma_ == false ) { cands.push_back(idx); }
	if ( chains_[idx].is_egamma_ == true ) { cands2.push_back(idx); }
      }
    }
    
  }

  // Fill ntuple
  for ( size_t idx = 0; idx < chains_.size(); ++idx ) {
    ElectronChain chain = chains_[idx];
    
    // Init tree here
    ntuple_.reset();

    // Apply prescale to fake candidates
    if ( chain.is_e_ == false && 
	 std::find( cands.begin(), cands.end(), idx ) == cands.end() &&
	 std::find( cands2.begin(), cands2.end(), idx ) == cands2.end() ) { continue; }
    
    // Data sample
    ntuple_.is_mc( chain.is_mc_ );
    ntuple_.is_aod( chain.is_aod_ );
    
    // Tag-side muon
    ntuple_.tag_pt( chain.tag_pt_ );
    ntuple_.tag_eta( chain.tag_eta_ );

    // Truth label
    ntuple_.is_e( chain.is_e_ );
    ntuple_.is_other( !chain.is_e_ );
    
    // Event stuff
    ntuple_.fill_evt( event.id() );
    ntuple_.set_rho( *rhoH_ );
    // Number of tracks or charged/neutral PF candidates

    // EGamma of low-pT
    ntuple_.is_egamma( chain.is_egamma_ );
    
    // Set prescale for BACKGROUND
    if ( !chain.is_e_ ) { ntuple_.set_prescale( std::abs(prescale_) > 1.e-6 ? prescale_ : 1. ); }

    //@@ Set weight??

    // "signal" info
    if ( validPtr(chain.sig_) ) {
      ntuple_.fill_gen( chain.sig_ ); // Add sig info
    }

    // Track info
    if ( validPtr(chain.trk_) ) {
      ntuple_.has_trk( chain.trk_match_ );
      ntuple_.fill_trk( chain.trk_, *beamspotH_ );
      ntuple_.trk_dr( chain.trk_dr_ );
      ntuple_.pdg_id( chain.pdg_id_ );
    }

    // ElectronSeed
    if ( validPtr(chain.seed_) ) { ntuple_.has_seed(true); }

    // PreId
//    noZS::EcalClusterLazyTools ecal_tools(event, setup, ebRecHits_, eeRecHits_);
//    ntuple_.fill_preid( *chain.preid_ecal_,
//			*chain.preid_hcal_,
//			*beamspotH_,
//			*rhoH_,
//			ecal_tools );

    // Tracker- or ECAL driven? 
    ntuple_.fill_seed( chain.seed_tracker_driven_, 
		       chain.seed_ecal_driven_ );
    if ( chain.seed_tracker_driven_ == true ||
	 chain.seed_ecal_driven_ == true ) {
    }

    // GsfTrack info
    if ( validPtr(chain.gsf_) ) {
      //if ( !chain.is_egamma_ ) { ntuple_.has_trk( chain.gsf_match_ ); }
      ntuple_.has_gsf( chain.gsf_match_ );
      ntuple_.fill_gsf( chain.gsf_, *beamspotH_ );
      ntuple_.gsf_dr( chain.gsf_dr_ );
      ntuple_.fill_bdt( chain.unbiased_, chain.ptbiased_ );
    }
    
    // PF GsfTrack info
    if ( validPtr(chain.pfgsf_) ) {
      //if ( chain.is_egamma_ ) { ntuple_.has_trk( chain.pfgsf_match_ ); }
      ntuple_.has_pfgsf( chain.pfgsf_match_ );
      ntuple_.fill_pfgsf( chain.pfgsf_, *beamspotH_ );
      ntuple_.pfgsf_dr( chain.pfgsf_dr_ );
    }

    // GsfElectron info
    if ( validPtr(chain.ele_) ) {

      //ntuple_.has_trk( chain.ele_match_ );
      //if ( chain.is_egamma_ ) { ntuple_.has_pfgsf( chain.ele_match_ ); }
      //else { ntuple_.has_gsf( chain.ele_match_ ); }
      ntuple_.has_ele( chain.ele_match_ );
      ntuple_.ele_dr( chain.ele_dr_ );

      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float mva_value = -999.;
      float mva_value_retrained = -999;
      float mva_value_depth10 = -999.;
      float mva_value_depth15 = -999.;
      if ( !chain.is_egamma_ ) {
	if ( mvaValueLowPtH_.isValid() && 
	     mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
	  mva_value = mvaValueLowPtH_->get( chain.ele_.key() );
	  chain.id_ = mva_value; 
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA output to low-pT GsfElectrons!" << std::endl;
	}
	if ( mvaValueLowPtDepth10H_.isValid() && 
	     mvaValueLowPtDepth10H_->size() == gsfElectronsH_->size() ) {
	  mva_value_depth10 = mvaValueLowPtDepth10H_->get( chain.ele_.key() );
	  //} else {
	  //std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA DEPTH10 output to GsfElectrons!" << std::endl;
	}
	if ( mvaValueLowPtDepth15H_.isValid() && 
	     mvaValueLowPtDepth15H_->size() == gsfElectronsH_->size() ) {
	  mva_value_depth15 = mvaValueLowPtDepth15H_->get( chain.ele_.key() );
	  //} else {
	  //std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA DEPTH15 output to GsfElectrons!" << std::endl;
	}
      } else {
	if ( mvaValueEGammaH_.isValid() && 
	     mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_value = mvaValueEGammaH_->get( chain.ele_.key() );
	  chain.id_ = mva_value; 
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching MVA output to PF GsfElectrons!" << std::endl;
	}
	if ( mvaValueEGammaRetrainedH_.isValid() && 
	     mvaValueEGammaRetrainedH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_value_retrained = mvaValueEGammaRetrainedH_->get( chain.ele_.key() );
	} else {
	  std::cout << "[IDNtuplizer::fill] ERROR! Issue matching retrained MVA to PF GsfElectrons!" << std::endl;
	}
      }
      
      //@@ dirty hack as is not in Event nor embedded in pat::Electron
      float conv_vtx_fit_prob = -999.;
      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
      //  conv_vtx_fit_prob = convVtxFitProb->get( chain.ele_.key() );
      //}
      
      ntuple_.fill_ele( chain.ele_, 
			mva_value, mva_value_retrained, mva_value_depth10, mva_value_depth15,
			conv_vtx_fit_prob, *rhoH_, chain.is_egamma_, 
			chain.unbiased_ // Rome ...
			);
      
      //ntuple_.fill_supercluster(chain.ele_);
      
    }
    
    if ( validPtr(chain.gsf_) ) {
      ntuple_.fill_image( chain.gsf_ref_eta_, chain.gsf_ref_phi_, chain.gsf_ref_R_, // Ref 
			  chain.gsf_ref_p_, chain.gsf_ref_pt_,
			  chain.gen_inner_eta_, chain.gen_inner_phi_, chain.gen_inner_R_, // GEN
			  chain.gen_inner_p_, chain.gen_inner_pt_,
			  chain.gen_proj_eta_, chain.gen_proj_phi_, chain.gen_proj_R_,
			  chain.gsf_inner_eta_, chain.gsf_inner_phi_, chain.gsf_inner_R_, // GSF
			  chain.gsf_inner_p_, chain.gsf_inner_pt_, chain.gsf_->charge(),
			  chain.gsf_proj_eta_, chain.gsf_proj_phi_, chain.gsf_proj_R_, 
			  chain.gsf_proj_p_,
			  chain.gsf_atcalo_eta_, chain.gsf_atcalo_phi_, chain.gsf_atcalo_R_, 
			  chain.gsf_atcalo_p_,
			  chain.clu_eta_, chain.clu_phi_, chain.clu_e_, chain.clu_nhit_, // Cluster
			  chain.pf_eta_, chain.pf_phi_, chain.pf_p_, // PFCands
			  chain.pf_pdgid_, chain.pf_matched_, chain.pf_lost_ // PFCands
			  );
    }

    tree_->Fill(); 
    
  }
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Utility methods ... /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::extractTrackPtrs() {
  tracks_.clear();
  pdgids_.clear();
  if ( isAOD_ == 1 ) {
    for ( size_t itrk = 0; itrk < ctfTracksH_->size(); ++itrk ) { 
      reco::TrackPtr ptr(ctfTracksH_,itrk);
      if ( !filterCand<reco::Track>(ptr) ) { continue; }
      tracks_.push_back(ptr); 
      pdgids_.insert(PdgIds::value_type(ptr.key(),0));
    }
  } else if ( isAOD_ == 0 ) {
    size_t iptr = 0;
    for ( const auto& ptr : *packedCandsH_ ) {
      if ( ptr.bestTrack() == nullptr ) { continue; }
      reco::TrackPtr trk(ptr.bestTrack(),iptr);
      if ( !filterCand<reco::Track>(trk) ) { continue; }
      tracks_.push_back(trk);
      pdgids_.insert(PdgIds::value_type(trk.key(),ptr.pdgId()));
      ++iptr;
    }
    for ( const auto& ptr : *lostTracksH_ ) { 
      if ( ptr.bestTrack() == nullptr ) { continue; }
      reco::TrackPtr trk(ptr.bestTrack(),iptr);
      if ( !filterCand<reco::Track>(trk) ) { continue; }
      tracks_.push_back(trk); 
      pdgids_.insert(PdgIds::value_type(trk.key(),ptr.pdgId()));
      ++iptr;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToTrk( reco::GsfTrackPtr& gsf, 
			    reco::TrackPtr& trk, 
			    bool is_egamma ) {

  // Shouldn't happen (but "link maps" may contain a invalid ptr)
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 1 ) {
      std::cout << "[IDNtuplizer::gsfToTrk] ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  
  // Attempt to navigate via Seed (and TrackExtra) to Track
  reco::ElectronSeedPtr seed;
  if ( gsfToSeed(gsf,seed) && seedToTrk(seed,trk) ) { return true; }
  //@@ else { return false; }

  // In the case of low-pT electrons ...
  if ( !is_egamma ) {
    // ... if above fails (e.g. TrackExtra missing), attempt to use following maps:
    //@@ Warning: linked tracks can be missing and/or duplicated!
    
    reco::GsfTrackRef gsf_ref(gsf.id(),gsf.get(),gsf.key());
    if ( gsf_ref.isNull() || !gsf_ref.isAvailable() ) { return false; }

    if ( isAOD_ == 1 ) {
      // ... gsf->track maps in AOD
      reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
      trk = edm::refToPtr(trk_ref);
      if ( validPtr(trk) ) { return true; }
    } else if ( isAOD_ == 0 ) {
      // ... gsf->packed maps in MINIAOD
      pat::PackedCandidateRef packed_ref = (*packedCandLinksH_)[gsf_ref];
      if ( packed_ref.isAvailable() && 
	   packed_ref.isNonnull() && 
	   packed_ref->bestTrack() != nullptr ) { 
	trk = reco::TrackPtr(packed_ref->bestTrack(),packed_ref.key());
	if ( validPtr(trk) ) { return true; }
      }
      // ... gsf->lost maps in MINIAOD
      pat::PackedCandidateRef lost_ref = (*lostTrackLinksH_)[gsf_ref];
      if ( lost_ref.isAvailable() && 
	   lost_ref.isNonnull() && 
	   lost_ref->bestTrack() != nullptr ) { 
	trk = reco::TrackPtr(lost_ref->bestTrack(),lost_ref.key());
	if ( validPtr(trk) ) { return true; }
      }
    }
  }

  return false;

}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::eleToGsf( reco::GsfElectronPtr& ele, reco::GsfTrackPtr& gsf ) {
  if ( !validPtr(ele) ) {
    std::cout << "[IDNtuplizer::eleToGsf] ERROR! GsfElectronPtr:"
	      << " ele.isNull(): " << ele.isNull()
	      << " ele.isAvailable(): " << ele.isAvailable()
	      << std::endl;
    return false;
  }
  gsf = edm::refToPtr(ele->gsfTrack());
  if ( !validPtr(gsf) ) {
    std::cout << "[IDNtuplizer::eleToGsf] ERROR! GsfTrackPtr:"
	      << " gsf.isNull(): " << gsf.isNull()
	      << " gsf.isAvailable(): " << gsf.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 1 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR! GsfTrackPtr:"
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
    if ( verbose_ > 3 ) { // TrackExtra are not stored by default in MINIAOD
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR: TrackExtra:" 
		<< " gsf->extra().isNull(): " << gsf->extra().isNull()
		<< " gsf->extra().isAvailable(): " << gsf->extra().isAvailable()
		<< std::endl; 
    }
    return false;
  }
  if ( traj.isNull() || !traj.isAvailable() ) { 
    if ( verbose_ > 1 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR: TrajectorySeedRef:" 
		<< " traj.isNull(): " << traj.isNull()
		<< " traj.isAvailable(): " << traj.isAvailable()
		<< std::endl; 
    }
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 1 ) {
      std::cout << "[IDNtuplizer::gsfToSeed] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 1 ) {
      std::cout << "[IDNtuplizer::seedToTrk] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  if ( !seed->isTrackerDriven() ) {
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToTrk] INFO! ElectronSeedPtr:"
		<< " seed->isTrackerDriven(): " << seed->isTrackerDriven()
		<< std::endl;
    }
  }
  trk = edm::refToPtr(seed->ctfTrack());
  if ( !validPtr(trk) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToTrk] INFO! TrackPtr:"
		<< " trk.isNull(): " << trk.isNull()
		<< " trk.isAvailable(): " << trk.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool IDNtuplizer::seedToCalo( reco::ElectronSeedPtr& seed, reco::CaloClusterPtr& calo ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToCalo] ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }

  edm::RefToBase<reco::CaloCluster> base = seed->caloCluster();
  if ( base.isNull() || !base.isAvailable() ) { 
    if ( verbose_ > 3 ) {
      std::cout << "[IDNtuplizer::seedToCalo] INFO! edm::RefToBase<reco::CaloCluster>:"
		<< " base.isNull(): " << base.isNull()
		<< " base.isAvailable(): " << base.isAvailable()
		<< std::endl;
    }
    return false;
  }
  calo = edm::refToPtr(seed->caloCluster().castTo<reco::CaloClusterRef>());
  if ( !validPtr(calo) ) { 
    std::cout << "[IDNtuplizer::seedToCalo] ERROR! CaloClusterPtr:"
	      << " calo.isNull(): " << calo.isNull()
	      << " calo.isAvailable(): " << calo.isAvailable()
	      << std::endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::trkToGsfLinks( std::vector<reco::TrackPtr>& ctfTracks,
				 edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 std::vector<TrkToGsfDR2>& trk2gsf,
				 bool is_egamma ) {

  trk2gsf.clear();

  // Record matched tracks via their keys
  std::vector<size_t> keys(ctfTracks.size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );

  // Shuffle keys to allow unbiased selection of surrogate tracks
  static unsigned seed = 1; //std::chrono::system_clock::now().time_since_epoch().count(); //@@ fixed seed !!
  std::shuffle( keys.begin(), keys.end(), std::default_random_engine(seed) );

  //////////
  // 1) Add GsfTracks to trk2gsf map if link (via provenence information) can be made
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackPtr gsf(gsfTracks, idx); 

    // Check validity of GsfTrackPtr (shouldn't ever fail?!)
    if ( !validPtr(gsf) ) { continue; }

    // If gsf not linked to trk (via ElectronSeed or "link maps"), then continue
    reco::TrackPtr trk;
    if ( !gsfToTrk(gsf,trk,is_egamma) ) { continue; }

    // Check if gsf is NOT already stored in map and add entry
    // Warning: multiple tracks can seed same GSF! (SHOULD ADD ALL ENTRIES? ONLY MATTERS FOR FR?)
    auto match_gsf_to_trk = std::find_if( trk2gsf.begin(), 
					  trk2gsf.end(), 
					  [gsf](const TrkToGsfDR2& dr2) { 
					    return gsf == dr2.obj2_; 
					  }
					  );
    if ( match_gsf_to_trk == trk2gsf.end() ) {
      trk2gsf.emplace_back( trk, gsf, deltaR2(trk,gsf) );
      keys.erase( std::remove( keys.begin(), 
			       keys.end(), 
			       trk.key() ), 
		  keys.end() ); // move TrackPtr to end and erase
    } else { 
      // This can happen on occasion with the "link maps"
      std::cout << "[IDNtuplizer::trkToGsfLinks] ERROR! GsfTrackPtr is already in the map!" << std::endl; 
    }
    
  }

  //////////
  // 2) Add entries with GsfTracks linked to (random) surrogate tracks if link cannot be made
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackPtr gsf(gsfTracks, idx); 

    if ( !validPtr(gsf) ) { continue; } //@@ shouldn't ever happen?!
    
    // If gsf is *linked* to trk (i.e. opposite to above), then continue
    reco::TrackPtr trk;
    if ( gsfToTrk(gsf,trk,is_egamma) ) { continue; }
    
    // Check if gsf is NOT already stored in map and add entry with surrogate track
    auto match_gsf_to_trk = std::find_if( trk2gsf.begin(), 
					  trk2gsf.end(), 
					  [gsf](const TrkToGsfDR2& dr2) { 
					    return gsf == dr2.obj2_; 
					  }
					  );
    if ( match_gsf_to_trk == trk2gsf.end() ) {
      
      // Identify (by key) "best" surrogate track to be the closest in pT to the PF GSF track
      auto best = std::min_element(keys.begin(),
				   keys.end(),
				   [gsf,ctfTracks]( const size_t& key1, 
						    const size_t& key2 ) {
				     return
				     std::abs(ctfTracks[key1]->pt()-gsf->ptMode())/gsf->ptMode()
				     <
				     std::abs(ctfTracks[key2]->pt()-gsf->ptMode())/gsf->ptMode();
				   }
				   );
      if ( best != keys.end() && validPtr(ctfTracks[*best]) ) {
	reco::TrackPtr trk = ctfTracks[*best];
	float dr2 = deltaR2(gsf,trk);
	if ( dr2 < dr_threshold_*dr_threshold_ ) {
	  trk2gsf.emplace_back( trk, gsf, deltaR2(gsf,trk) ); // Label as correct match
	} else {
	  trk2gsf.emplace_back( trk, gsf, -1.*deltaR2(gsf,trk) ); // Label as surrogate match
	}
	keys.erase( std::remove( keys.begin(), 
				 keys.end(), 
				 trk.key() ),
		    keys.end() ); // move TrackPtr to end and erase
      } else {
	std::cout << "[IDNtuplizer::trkToGsfLinks] " 
		  << "ERROR: Couldn't find a valid 'best' surrogate track matched in pT to the PF GSF track!";
	trk2gsf.emplace_back( reco::TrackPtr(), gsf, IDNtuple::NEG_FLOATSQ ); // shouldn't ever be called
      }

    } else { std::cout << "[IDNtuplizer::trkToGsfLinks] ERROR! GsfTrackPtr is already in the map!!!" << std::endl; }
  
  } // GsfTracks loop
  
  //////////
  // 3) Add (null) entries for all Tracks without a match to a GsfTrack
  for ( auto trk : ctfTracks ) {
    
    // Check if trk already stored in map and, if not, add "empty" entry
    auto match_trk = std::find_if( trk2gsf.begin(), 
				   trk2gsf.end(), 
				   [trk](const TrkToGsfDR2& dr2) { 
				     return dr2.obj1_ == trk; 
				   }
				   );
    if ( match_trk == trk2gsf.end() ) {
      reco::GsfTrackPtr gsf;
      trk2gsf.emplace_back( trk, gsf, IDNtuple::NEG_FLOATSQ );
    }

  }

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDNtuplizer::gsfToEleLinks( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
				 std::vector<GsfToEleDR2>& gsf2ele ) {
  
  gsf2ele.clear();

  // Record matched GsfTracks via their keys
  std::vector<size_t> keys(gsfTracks->size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );
  
  //////////
  // 1) Add GsfElectrons to gsf2ele map if link (via provenence information) can be made
  for ( size_t idx = 0; idx < gsfElectrons->size(); ++idx ) {
    reco::GsfElectronPtr ele(gsfElectrons, idx);

    // Check validity of GsfElectronPtr (shouldn't ever fail?!)
    if ( !validPtr(ele) ) { continue; }
    
    // If ele not linked to GSF track, then continue
    reco::GsfTrackPtr gsf;
    if ( !eleToGsf(ele,gsf) ) { continue; }
    
    // Check if ele is already stored in map and, if not, add entry
    auto match_ele_to_gsf = std::find_if( gsf2ele.begin(), 
					  gsf2ele.end(), 
					  [ele](const GsfToEleDR2& dr2) { 
					    return ele == dr2.obj2_; 
					  }
					  );
    if ( match_ele_to_gsf == gsf2ele.end() ) {
      gsf2ele.emplace_back( gsf, ele, deltaR2(gsf,ele) );
      keys.erase( std::remove( keys.begin(), 
			       keys.end(), 
			       gsf.key() ), 
		  keys.end() ); // move key to end and erase
    } else { std::cout << "[IDNtuplizer::gsfToEleLinks]" 
		       << " ERROR! GsfElectronPtr is already in the map!" << std::endl; }

  }

  //////////
  // 2) Add (null) entries for all GsfTracks without a match to a GsfElectron
  for ( auto key : keys ) {
    reco::GsfTrackPtr gsf(gsfTracks, key);

    // Check if gsf already stored in map and, if not, add "empty" entry
    auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
					  gsf2ele.end(), 
					  [gsf](const GsfToEleDR2& dr2) { 
					    return gsf == dr2.obj1_;
					  }
					  );
    if ( match_gsf_to_ele == gsf2ele.end() ) {
      gsf2ele.emplace_back( gsf, reco::GsfElectronPtr(), IDNtuple::NEG_FLOATSQ ); // null GSF ele
    }
    
  } 
}

////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::gsfToPfGsfLinks( edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				   edm::Handle< std::vector<reco::GsfTrack> >& gsfTracksEGamma,
				   std::vector<GsfToGsfDR2>& gsf2pfgsf ) {

  if ( gsfTracks->empty() && gsfTracksEGamma->empty() ) { return; }

  // Record matched GsfTracks via their keys
  std::vector<size_t> keys(gsfTracks->size());
  std::generate( keys.begin(), keys.end(), [n=0] () mutable { return n++; } );
  
  //////////
  // 1) Determine deltaR(pfgsg,gsf) for all combinations  
  std::vector<GsfToGsfDR2> gsf2pfgsf_all;
  gsf2pfgsf_all.reserve( gsfTracks->size()*gsfTracksEGamma->size() );
  for ( size_t igsf = 0; igsf < gsfTracks->size(); ++igsf ) {
    reco::GsfTrackPtr gsf(gsfTracks, igsf);
    for ( size_t ipfgsf = 0; ipfgsf < gsfTracksEGamma->size(); ++ipfgsf ) {
      reco::GsfTrackPtr pfgsf(gsfTracksEGamma, ipfgsf);
      gsf2pfgsf_all.emplace_back( gsf, pfgsf, deltaR2(gsf,pfgsf) );
    }
  }

  //////////
  // 2) For each pfgsf, select best match according to DeltaR2 metric
  gsf2pfgsf.clear();
  std::sort( gsf2pfgsf_all.begin(), 
	     gsf2pfgsf_all.end(), 
	     GsfToGsfDR2::compare_by_dr2 );

  for ( size_t ipfgsf = 0; ipfgsf < gsfTracksEGamma->size(); ++ipfgsf ) {
    reco::GsfTrackPtr pfgsf(gsfTracksEGamma, ipfgsf);
    auto match_pfgsf_to_gsf = std::find_if( gsf2pfgsf_all.begin(),
					    gsf2pfgsf_all.end(),
					    [pfgsf](const GsfToGsfDR2& dr2) {
					      return pfgsf == dr2.obj2_;
					    }
					    );
    if ( match_pfgsf_to_gsf != gsf2pfgsf.end() ) { 
      if ( validPtr(match_pfgsf_to_gsf->obj1_) && 
	   match_pfgsf_to_gsf->dr2_ < dr_threshold_*dr_threshold_ ) { 
	gsf2pfgsf.emplace_back(*match_pfgsf_to_gsf); 
	keys.erase( std::remove( keys.begin(), 
				 keys.end(), 
				 match_pfgsf_to_gsf->obj1_.key() ), 
		    keys.end() ); // Erase GSF key (https://en.wikipedia.org/wiki/Erase-remove_idiom)
      } else {
	// If no match within dr_threshold_, then store null GSF track
	gsf2pfgsf.emplace_back( reco::GsfTrackPtr(), match_pfgsf_to_gsf->obj2_, IDNtuple::NEG_FLOATSQ ); // null GSF
      }
    } else {
      std::cout << "[IDNtuplizer::gsfToPfGsfLinks]"
		<< " ERROR! PF GsfTrackPtr is not found in the map!" << std::endl;
    }
    if ( gsf2pfgsf.size() >= gsfTracksEGamma->size() ) { break; }  // found unique match for all PF GsfTracks
  }

  //////////
  // 3) Add (null) entries for all low-pT GsfTracks without a match to a PF GsfTrack
  for ( auto key : keys ) {
    reco::GsfTrackPtr gsf(gsfTracks,key);
    // Check if trk already stored in map and, if not, add "empty" entry
    auto match_gsf = std::find_if( gsf2pfgsf.begin(), 
				   gsf2pfgsf.end(), 
				   [gsf](const GsfToGsfDR2& dr2) { 
				     return dr2.obj1_ == gsf; 
				   }
				   );
    if ( match_gsf == gsf2pfgsf.end() ) {
      gsf2pfgsf.emplace_back( gsf, reco::GsfTrackPtr(), IDNtuple::NEG_FLOATSQ ); // null PF GSF
    }
  }    

}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  edm::Handle< edm::View<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  std::vector< edm::Ptr<T> > cands;
  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
    edm::Ptr<T> cand(candidates,icand);
    if ( validPtr(cand) ) { cands.push_back(cand); }
    else {
      std::cout << "[IDNtuplizer::sigToCandLinks] ERROR! CandidatePtr:"
		<< " cand.isNull(): " << cand.isNull()
		<< " cand.isAvailable(): " << cand.isAvailable()
		<< std::endl;
    }
  }
  sigToCandLinks<T>( signal_electrons, cands, sig2cand, other_cand, append );
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  edm::Handle< std::vector<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  std::vector< edm::Ptr<T> > cands;
  for ( size_t icand = 0; icand < candidates->size(); ++icand ) {
    edm::Ptr<T> cand(candidates,icand);
    if ( validPtr(cand) ) { cands.push_back(cand); }
    else {
      std::cout << "[IDNtuplizer::sigToCandLinks] ERROR! CandidatePtr:"
		<< " cand.isNull(): " << cand.isNull()
		<< " cand.isAvailable(): " << cand.isAvailable()
		<< std::endl;
    }
  }
  sigToCandLinks<T>( signal_electrons, cands, sig2cand, other_cand, append );
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::sigToCandLinks( std::set<reco::CandidatePtr>& signal_electrons,
				  std::vector< edm::Ptr<T> >& candidates,
				  std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
				  std::vector< DeltaR2<reco::Candidate,T> >& other_cand,
				  bool append ) {
  
  if ( !append ) { sig2cand.clear(); }
  if ( !append ) { other_cand.clear(); }
  
  // DeltaR2 matches for all combinations of signal electrons and reco::Candidate
  std::vector< DeltaR2<reco::Candidate,T> > sig2cand_all;
  sig2cand_all.reserve( sig2cand.size() + signal_electrons.size()*candidates.size() );
  unsigned int icand = 0;
  for ( auto cand : candidates ) {
    ++icand;
    if ( !validPtr(cand) ) { continue; }
    if ( !filterCand<T>(cand) ) { continue; }
    for ( auto sig : signal_electrons ) {
      sig2cand_all.emplace_back( sig,
				 cand,
				 deltaR2<reco::Candidate,T>(sig,cand) );
    }
    other_cand.emplace_back( reco::CandidatePtr(), cand, IDNtuple::NEG_FLOATSQ ); // Null signal
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
			       return dr2.obj1_ == iter.obj1_; 
			     }
			     );
    auto cand = std::find_if( sig2cand.begin(), 
			      sig2cand.end(), 
			      [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				return dr2.obj2_ == iter.obj2_; 
			      }
			      );
    if ( sig == sig2cand.end() && cand == sig2cand.end() ) {
      if ( iter.dr2_ < dr_max_*dr_max_ ) { sig2cand.push_back(iter); }
      else { 
	sig2cand.emplace_back( iter.obj1_, edm::Ptr<T>(), IDNtuple::NEG_FLOATSQ ); // Null cand
      }
    }
    if ( sig2cand.size() >= signal_electrons.size() ) { break; }
  }
  
  // Remove matched candidates
  for ( auto iter : sig2cand ) {
    auto match = std::find_if( other_cand.begin(), 
			       other_cand.end(), 
			       [iter](const DeltaR2<reco::Candidate,T>& dr2) { 
				 return dr2.obj2_ == iter.obj2_; 
			       }
			       );
    if ( match != other_cand.end() ) { other_cand.erase(match); }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
void IDNtuplizer::match( reco::CandidatePtr& sig,
			 std::vector< DeltaR2<reco::Candidate,T> >& sig2cand,
			 edm::Ptr<T>& cand_ptr, 
			 float& cand_dr, 
			 bool& cand_match ) { 
  
  // Identify GEN electron matches to candidate
  auto match = std::find_if( sig2cand.begin(), 
			     sig2cand.end(), 
			     [sig](const DeltaR2<reco::Candidate,T>& dr2) { 
			       return dr2.obj1_ == sig; 
			     }
			     );
  if ( match != sig2cand.end() ) {
    if ( validPtr(match->obj2_) ) {
      cand_ptr = match->obj2_; // pass by ref
      if ( deltaR2(sig,cand_ptr) >= 0. ) cand_dr = sqrt(deltaR2(sig,cand_ptr)); // pass by ref
      cand_match = ( cand_dr >= 0. ) && ( cand_dr < dr_threshold_ ); // pass by ref
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T1, typename T2>
float IDNtuplizer::deltaR2( edm::Ptr<T1>& cand1, edm::Ptr<T2>& cand2 ) {
  return reco::deltaR2( cand1->eta(), 
			cand1->phi(),
			cand2->eta(),
			cand2->phi() ); 
}

float IDNtuplizer::deltaR2( reco::TrackPtr& trk,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( trk->eta(), 
			trk->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

float IDNtuplizer::deltaR2( reco::GsfTrackPtr& pfgsf,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( pfgsf->etaMode(), 
			pfgsf->phiMode(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

float IDNtuplizer::deltaR2( reco::CandidatePtr& sig,
			    reco::GsfTrackPtr& gsf ) {
  return reco::deltaR2( sig->eta(), 
			sig->phi(),
			gsf->etaMode(),
			gsf->phiMode() ); 
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T>
bool IDNtuplizer::filterCand( edm::Ptr<T>& cand ) { 
  if ( cand->pt() < minTrackPt_ ) { return false; }
  //if ( std::abs(cand->eta()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  return true;
}

bool IDNtuplizer::filterCand( edm::Ptr<reco::Track>& trk ) { 
  if ( trk->pt() < minTrackPt_ ) { return false; }
  //if ( std::abs(trk->eta()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  if ( !(trk->quality(reco::TrackBase::qualityByName("highPurity"))) ) { return false; }
  return true; 
}

bool IDNtuplizer::filterCand( edm::Ptr<reco::GsfTrack>& gsf ) { 
  if ( gsf->ptMode() < minTrackPt_ ) { return false; }
  //if ( std::abs(gsf->etaMode()) > 2.5 ) { return false; } //@@ Needed for PFCandidates?
  return true; 
}

////////////////////////////////////////////////////////////////////////////////
//
template <typename T> 
bool IDNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

////////////////////////////////////////////////////////////////////////////////
//
BaseParticlePropagator IDNtuplizer::extrapolate_track( const Vector& mom, const Point& pos, int charge,
						       int& reach_ECAL, GlobalPoint& pos_ECAL,
						       int& reach_HCAL, GlobalPoint& pos_HCAL,
						       int& reach_EXIT, GlobalPoint& pos_EXIT ) {
  
  // Propagate 'electron' to ECAL surface
  double energy = mom.R() + 0.000511*0.000511; // electron mass
  RawParticle particle( XYZTLorentzVector( mom.x(), mom.y(), mom.z(), energy ),
			XYZTLorentzVector( pos.x(), pos.y(), pos.z(), 0. ) );
  BaseParticlePropagator propagator( particle, 0., 0., 3.8 ); // fixed Z value for B-field
  particle.setCharge(charge);
  
  // ECAL: true = first half loop; 0 = does not reach ECAL; 1 = yes, barrel; 2 = yes, endcaps
  propagator.propagateToEcalEntrance(true);
  reach_ECAL = propagator.getSuccess();
  pos_ECAL = GlobalPoint( particle.x(),  particle.y(), particle.z() );

  // ECAL: true = first half loop; 0 = does not reach ECAL; 1 = yes, barrel; 2 = yes, endcaps 
  propagator.propagateToHcalEntrance(true);
  reach_HCAL = propagator.getSuccess(); 
  pos_HCAL = GlobalPoint( particle.x(),  particle.y(), particle.z() );

  // ECAL: true = first half loop; 0 = does not reach ECAL; 1 = yes, barrel; 2 = yes, endcaps 
  propagator.propagateToHcalExit(true);
  reach_EXIT = propagator.getSuccess(); 
  pos_EXIT = GlobalPoint( particle.x(),  particle.y(), particle.z() );

  return propagator;

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::build_image( const edm::Event& event,
			       const edm::EventSetup& setup ) {

  if ( isAOD_ == 0 ) { 

    bool first_fake = true;
    for ( auto& chain : chains_ ) {

      // LOW PT ELECTRONS ONLY!!!
      if ( chain.is_egamma_ ) { continue; } 

      //@@@@ TEMP: SIGNAL ONLY !!!! 
      //if ( !chain.is_e_ ) { continue; } 

      // Ensure GSF track is found!
      if ( !validPtr(chain.gsf_) ) { continue; }

      //@@ BARREL ONLY!!! (FOR NOW ...)
      if ( fabs( chain.gsf_->momentumMode().eta() ) > 1.48 ) { continue; } 

      // Print RECO chain 
      std::stringstream ss;
      ss << std::endl
	 << "VALID:      "
	 << std::fixed << std::setprecision(0) 
	 << " GEN: " << std::setw(1) << validPtr(chain.sig_)
	 << " TRK: " << std::setw(1) << validPtr(chain.trk_)
	 << " GSF: " << std::setw(1) << validPtr(chain.gsf_)
	 << " ELE: " << std::setw(1) << validPtr(chain.ele_);
      
      ////////////////////////////////////////
      // REFERENCE coords in eta-phi space (using Inner GSF track P4)
      chain.gsf_ref_eta_ = chain.gsf_->momentumMode().eta();
      chain.gsf_ref_phi_ = chain.gsf_->momentumMode().phi();
      if ( validPtr(chain.ele_) ) { chain.gsf_ref_R_ = chain.ele_->trackPositionAtVtx().R(); }
      chain.gsf_ref_p_ = sqrt(chain.gsf_->momentumMode().Mag2());
      chain.gsf_ref_pt_ = sqrt(chain.gsf_->momentumMode().Perp2());
      
      ////////////////////////////////////////
      // GEN channel (only useful when drawing images!)
      if ( chain.is_mc_ && validPtr(chain.sig_) ) {
	int reach_ECAL = 0; GlobalPoint pos_ECAL;
	int reach_HCAL = 0; GlobalPoint pos_HCAL;
	int reach_EXIT = 0; GlobalPoint pos_EXIT;
	extrapolate_track( chain.sig_->momentum(),
			   chain.sig_->vertex(),
			   chain.sig_->charge(),
			   reach_ECAL, pos_ECAL,
			   reach_HCAL, pos_HCAL,
			   reach_EXIT, pos_EXIT );
	// 1st point: inner P4 (i.e. as if extrapolated to ECAL surface as a neutral)
	chain.gen_inner_eta_ = adj_eta(chain.sig_->momentum().eta(),chain.gsf_ref_eta_,chain.sig_->charge());
	chain.gen_inner_phi_ = adj_phi(chain.sig_->momentum().phi(),chain.gsf_ref_phi_,chain.sig_->charge());
	chain.gen_inner_R_ = chain.sig_->vertex().R();
	chain.gen_inner_p_ = sqrt(chain.sig_->momentum().Mag2());
	chain.gen_inner_pt_ = sqrt(chain.sig_->momentum().Perp2());
	// 2nd point: inner P4 extrapolated to ECAL surface (i.e. as if charged)
	if ( reach_ECAL ) {
	  chain.gen_proj_eta_ = adj_eta(pos_ECAL.eta(),chain.gsf_ref_eta_,chain.sig_->charge());
	  chain.gen_proj_phi_ = adj_phi(pos_ECAL.phi(),chain.gsf_ref_phi_,chain.sig_->charge());
	  chain.gen_proj_R_ = sqrt(pos_ECAL.perp2());
	}
      }
      
      ////////////////////////////////////////
      // WINDOW: GSF-based eta-phi window
      int charge = chain.gsf_->charge();
      float window_eta_min = -5.;
      float window_eta_max = -5.;
      float window_eta_ext =  0.;//05;
      float window_phi_min = -5.;
      float window_phi_max = -5.;
      float window_phi_ext =  0.;//05;

      int reach_ECAL = 0; GlobalPoint pos_ECAL;
      int reach_HCAL = 0; GlobalPoint pos_HCAL;
      int reach_EXIT = 0; GlobalPoint pos_EXIT;
      BaseParticlePropagator propagator = extrapolate_track( chain.gsf_->momentumMode(),
							     chain.gsf_->referencePoint(), // same as deprecated vertex()
							     chain.gsf_->charge(),
							     reach_ECAL, pos_ECAL,
							     reach_HCAL, pos_HCAL,
							     reach_EXIT, pos_EXIT );
      
      // Min/max eta/phi values (w.r.t. ref) based on inner P4
      window_eta_min = adj_eta(chain.gsf_->momentumMode().eta(),chain.gsf_ref_eta_,charge);
      window_phi_min = adj_phi(chain.gsf_->momentumMode().phi(),chain.gsf_ref_phi_,charge);
      window_eta_max = adj_eta(chain.gsf_->momentumMode().eta(),chain.gsf_ref_eta_,charge);
      window_phi_max = adj_phi(chain.gsf_->momentumMode().phi(),chain.gsf_ref_phi_,charge);
      
      // "Maximum" eta/phi values (w.r.t. ref) based on extrapolated P4
      if ( reach_ECAL ) {
	float tmp_eta = adj_eta(pos_ECAL.eta(),chain.gsf_ref_eta_,charge);
	if ( tmp_eta > window_eta_max ) { window_eta_max = tmp_eta; }
	if ( tmp_eta < window_eta_min ) { window_eta_min = tmp_eta; }
	float tmp_phi = adj_phi(pos_ECAL.phi(),chain.gsf_ref_phi_,charge);
	if ( tmp_phi > window_phi_max ) { window_phi_max = tmp_phi; }
	if ( tmp_phi < window_phi_min ) { window_phi_min = tmp_phi; }
      }
      
      if ( validPtr(chain.ele_) ) {
	// Check "maximum" eta/phi values (w.r.t. ref) at the ECAL surface
	float tmp_eta = adj_eta(chain.ele_->trackPositionAtCalo().eta(),chain.gsf_ref_eta_,charge);
	if ( tmp_eta > window_eta_max ) { window_eta_max = tmp_eta; }
	if ( tmp_eta < window_eta_min ) { window_eta_min = tmp_eta; }
	float tmp_phi = adj_phi(chain.ele_->trackPositionAtCalo().phi(),chain.gsf_ref_phi_,charge);
	if ( tmp_phi > window_phi_max ) { window_phi_max = tmp_phi; }
	if ( tmp_phi < window_phi_min ) { window_phi_min = tmp_phi; }
      } 
      
      // // Open up window further by the "_ext" values (0.1?)
      // window_eta_min -= window_eta_ext;
      // window_eta_max += window_eta_ext;
      // window_phi_min -= window_phi_ext;
      // window_phi_max += window_phi_ext;
    
      ////////////////////////////////////////
      // GSF channel
      
      // 1st point: inner P4 (i.e. as if extrapolated to ECAL surface as a neutral)
      chain.gsf_inner_eta_ = adj_eta(chain.gsf_->momentumMode().eta(),chain.gsf_ref_eta_,charge);
      chain.gsf_inner_phi_ = adj_phi(chain.gsf_->momentumMode().phi(),chain.gsf_ref_phi_,charge);
      if ( validPtr(chain.ele_) ) {chain.gsf_inner_R_ = chain.ele_->trackPositionAtVtx().R(); }
      chain.gsf_inner_p_ = sqrt(chain.gsf_->momentumMode().Mag2());
      chain.gsf_inner_pt_ = sqrt(chain.gsf_->momentumMode().Perp2());
      
      // 2nd point: inner P4 extrapolated to ECAL surface (i.e. as if charged)
      chain.gsf_proj_eta_ = adj_eta(pos_ECAL.eta(),chain.gsf_ref_eta_,charge);
      chain.gsf_proj_phi_ = adj_phi(pos_ECAL.phi(),chain.gsf_ref_phi_,charge);
      chain.gsf_proj_R_ = sqrt(pos_ECAL.perp2());
      chain.gsf_proj_p_ = sqrt(chain.gsf_->momentumMode().Mag2());
      
      // 3rd point: outer P4 at ECAL surface (actual entry point)
      if ( validPtr(chain.ele_) ) {
	chain.gsf_atcalo_eta_ = adj_eta(chain.ele_->trackPositionAtCalo().eta(),chain.gsf_ref_eta_,charge);
	chain.gsf_atcalo_phi_ = adj_phi(chain.ele_->trackPositionAtCalo().phi(),chain.gsf_ref_phi_,charge);
	chain.gsf_atcalo_R_ = chain.ele_->trackPositionAtCalo().R();
	chain.gsf_atcalo_p_ = sqrt(chain.ele_->trackMomentumOut().Mag2());
      }
      
      ss << std::endl
	 << "REFERENCE:  "
	 << std::fixed << std::setprecision(2) 
	 << " eta: " << std::setw(5) << chain.gsf_ref_eta_
	 << " phi: " << std::setw(5) << chain.gsf_ref_phi_
	 << "  pt: " << std::setw(5) << chain.gsf_ref_pt_
	 << "   p: " << std::setw(5) << chain.gsf_ref_p_;
      
      if ( validPtr(chain.ele_) ) {
	ss << std::endl
	   << "MomOut:     "
	   << std::fixed << std::setprecision(2) 
	   << " eta: " << std::setw(5) << chain.ele_->trackMomentumOut().eta()
	   << " phi: " << std::setw(5) << chain.ele_->trackMomentumOut().phi()
	   << "  pt: " << std::setw(5) << sqrt(chain.ele_->trackMomentumOut().perp2())
	   << " (atcalo)"
	   << std::endl
	   << "MomAtCalo:  "
	   << std::fixed << std::setprecision(2) 
	   << " eta: " << std::setw(5) << chain.ele_->trackMomentumAtCalo().eta()
	   << " phi: " << std::setw(5) << chain.ele_->trackMomentumAtCalo().phi()
	   << "  pt: " << std::setw(5) << sqrt(chain.ele_->trackMomentumAtCalo().perp2())
	   << " (not used)";
      }

      ss << std::endl
	 << "MomMode:    "
	 << std::fixed << std::setprecision(2) 
	 << " eta: " << std::setw(5) << chain.gsf_->momentumMode().eta()
	 << " phi: " << std::setw(5) << chain.gsf_->momentumMode().phi()
	 << "  pt: " << std::setw(5) << sqrt(chain.gsf_->momentumMode().perp2())
	 << " (#1 inner)"
	 << std::endl
	 << "posECAL:    "
	 << std::fixed << std::setprecision(2) 
	 << " eta: " << std::setw(5) << pos_ECAL.eta()
	 << " phi: " << std::setw(5) << pos_ECAL.phi()
	 << "   R: " << std::setw(6) << sqrt(pos_ECAL.perp2())
	 << "   z: " << std::setw(7) << pos_ECAL.z()
	 << " Reach? " << std::setw(1) << int(reach_ECAL)
	 << " (#2 proj)";
      
      if ( validPtr(chain.ele_) ) {
	ss << std::endl
	   << "PosAtCalo:  "
	   << std::fixed << std::setprecision(2) 
	   << " eta: " << std::setw(5) << chain.ele_->trackPositionAtCalo().eta()
	   << " phi: " << std::setw(5) << chain.ele_->trackPositionAtCalo().phi()
	   << "   R: " << std::setw(6) << sqrt(chain.ele_->trackPositionAtCalo().perp2())
	   << "   z: " << std::setw(7) << chain.ele_->trackPositionAtCalo().z()
	   << " (#3 atcalo)";
      }
      
      ss << std::endl
	 << "posHCAL:    "
	 << std::fixed << std::setprecision(2) 
	 << " eta: " << std::setw(5) << pos_HCAL.eta()
	 << " phi: " << std::setw(5) << pos_HCAL.phi()
	 << "   R: " << std::setw(6) << sqrt(pos_HCAL.perp2())
	 << "   z: " << std::setw(7) << pos_HCAL.z()
	 << " Reach? " << std::setw(1) << int(reach_HCAL)
	 << std::endl
	 << "posEXIT:    "
	 << std::fixed << std::setprecision(2) 
	 << " eta: " << std::setw(5) << pos_EXIT.eta()
	 << " phi: " << std::setw(5) << pos_EXIT.phi()
	 << "   R: " << std::setw(6) << sqrt(pos_EXIT.perp2())
	 << "   z: " << std::setw(7) << pos_EXIT.z()
	 << " Reach? " << std::setw(1) << int(reach_EXIT);

      if ( validPtr(chain.ele_) ) {
	ss << std::endl
	   << "fbrem:      "
	   << std::fixed << std::setprecision(2) 
	   << std::setw(5) << chain.ele_->fbrem()
	   << ", fbrem(SC): "
	   << std::setw(5) << chain.ele_->superClusterFbrem()
	   << ", IPxy: "
	   << std::setw(5) << propagator.xyImpactParameter()
	   << ", helixRadius: "
	   << std::setw(5) << propagator.helixRadius()
	   << ", onBarrel: "
	   << std::setw(5) << propagator.onBarrel()
	   << ", onEndcap: "
	   << std::setw(5) << propagator.onEndcap();
      }

      // Print eta-phi window
      ss << std::endl
	 << "WINDOW eta: "
	 << std::fixed << std::setprecision(2) 
	 << " min: " << std::setw(5) << window_eta_min
	 << " max: " << std::setw(5) << window_eta_max
	 << " del: " << std::setw(5) << window_eta_max - window_eta_min
	 << " ext: " << std::setw(5) << window_eta_ext
	 << std::fixed << std::setprecision(0) 
	 << " charge: " << std::setw(2) << charge << std::endl
	 << "WINDOW phi: "
	 << std::fixed << std::setprecision(2) 
	 << " min: " << std::setw(5) << window_phi_min
	 << " max: " << std::setw(5) << window_phi_max
	 << " del: " << std::setw(5) << reco::deltaPhi(window_phi_max,window_phi_min)
	 << " ext: " << std::setw(5) << window_phi_ext;
      
      ////////////////////////////////////////
      // TRK channel???
      
      ////////////////////////////////////////
      // Clusters channel
      if ( validPtr(chain.ele_) ) {
	const reco::SuperClusterRef& sc = chain.ele_->superCluster();
	if ( sc.isNonnull() ) { 
	  int charge = chain.gsf_->charge();
	  for ( auto& cluster : sc->clusters() ) {
	    chain.clu_eta_.push_back(adj_eta(cluster->eta(),chain.gsf_ref_eta_,charge));
	    chain.clu_phi_.push_back(adj_phi(cluster->phi(),chain.gsf_ref_phi_,charge));
	    chain.clu_e_.push_back(cluster->correctedEnergy());
	    chain.clu_nhit_.push_back(cluster->hitsAndFractions().size());
	    ss << std::endl
	       << " CLUSTER:   "
	       << std::fixed << std::setprecision(2) 
	       << " eta: " << std::setw(5) << adj_eta(cluster->eta(),chain.gsf_ref_eta_,charge)
	       << " phi: " << std::setw(5) << adj_phi(cluster->phi(),chain.gsf_ref_phi_,charge)
	       << "   e: " << std::setw(5) << cluster->correctedEnergy()
	       << std::fixed << std::setprecision(0) 
	       << " nhit: " << std::setw(2) << cluster->hitsAndFractions().size();
	  }
	}
      }
      
      ////////////////////////////////////////
      // PF candidates channel
      edm::Ptr<pat::PackedCandidate> ele = edm::refToPtr((*packedCandLinksH_)[chain.gsf_]);
      if ( ele.isNonnull() ) { 
	size_t total_size = packedCandsH_->size() + lostTracksH_->size();
	for ( size_t idx = 0; idx < total_size; ++idx ) {
	  edm::Ptr<pat::PackedCandidate> cand;
	  if ( idx < packedCandsH_->size() ) { 
	    cand = edm::Ptr<pat::PackedCandidate>(packedCandsH_,idx);
	  } else {
	    cand = edm::Ptr<pat::PackedCandidate>(lostTracksH_,idx-packedCandsH_->size());
	  }
	  if ( cand.isNull() ) { continue; }
	  //if ( cand == ele ) { continue; }
	  int matched = (cand==ele) ? 1 : 0;
	  int charge = cand->charge();
	  float adjusted_eta = adj_eta(cand->eta(),chain.gsf_ref_eta_,charge);
	  float adjusted_phi = adj_phi(cand->phi(),chain.gsf_ref_phi_,charge);
	  if ( adjusted_eta < window_eta_min || adjusted_eta > window_eta_max || 
	       adjusted_phi < window_phi_min || adjusted_phi > window_phi_max ) { continue; }
	  chain.pf_eta_.push_back(adjusted_eta);
	  chain.pf_phi_.push_back(adjusted_phi);
	  chain.pf_p_.push_back(cand->p());
	  chain.pf_pdgid_.push_back(cand->pdgId());
	  chain.pf_matched_.push_back(matched);
	  chain.pf_lost_.push_back( idx < packedCandsH_->size() ? 0 : 1 );
	  //@@
	  if ( std::find( pf_pdgids_.begin(),
			  pf_pdgids_.end(),
			  cand->pdgId() ) == pf_pdgids_.end() ) {
	    pf_pdgids_.push_back(cand->pdgId());
	  }
	  
	  ss << std::endl
	     << " PFCAND:    "
	     << std::fixed << std::setprecision(2) 
	     << " eta: " << std::setw(5) << adjusted_eta
	     << " phi: " << std::setw(5) << adjusted_phi
	     << "  pt: " << std::setw(5) << cand->pt()
	     << "   p: " << std::setw(5) << cand->p()
	     << std::fixed << std::setprecision(0) 
	     << " pdg " << std::setw(3) << cand->pdgId()
	     << " Q: " << std::setw(3) << cand->charge()
	     << " idx: " << std::setw(2) << idx
	     << " match: " << std::setw(1) << matched
	     << " lost: " << std::setw(1) << ( idx < packedCandsH_->size() ? 0 : 1 );
	}
      }

      // Print out image details
      if ( chain.is_e_ || first_fake ) { 
	if ( !chain.is_e_ ) { first_fake = false; } // Only print first fake candidate
	//std::cout << ss.str() << std::endl; 
      }
      
    } // if ( isAOD_ == 0 )
    
  } // for ( auto& chain : chains_ )
  
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDNtuplizer);

////////////////////
////////////////////
////////////////////
// OBSOLETE METHODS
////////////////////
////////////////////
////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::lowPtElectrons( std::set<reco::CandidatePtr>& signal_electrons,
				  std::vector<SigToTrkDR2>& sig2trk,
				  std::vector<SigToTrkDR2>& other_trk,
				  std::vector<SigToGsfDR2>& sig2gsf,
				  std::vector<TrkToGsfDR2>& trk2gsf,
				  std::vector<GsfToGsfDR2>& gsf2pfgsf ) {

  // Match low-pT GsfTracks to low-pT GsfElectrons
  std::vector<GsfToEleDR2> gsf2ele;
  gsfToEleLinks( gsfTracksH_, 
		 gsfElectronsH_, 
		 gsf2ele );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::lowPtElectrons] gsfToEleLinks:" << std::endl 
	      << " gsfTracksH_->size(): " << gsfTracksH_->size() << std::endl
	      << " gsfElectronsH_->size(): " << gsfElectronsH_->size() << std::endl
	      << " gsf2ele.size(): " << gsf2ele.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : gsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  lowPtElectrons_signal( signal_electrons, sig2trk, other_trk, sig2gsf, gsf2pfgsf, gsf2ele );
  lowPtElectrons_fakes( other_trk, trk2gsf, gsf2pfgsf, gsf2ele );
    
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "signal electrons" 
void IDNtuplizer::lowPtElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
					 std::vector<SigToTrkDR2>& sig2trk,
					 std::vector<SigToTrkDR2>& other_trk,
					 std::vector<SigToGsfDR2>& sig2gsf,
					 std::vector<GsfToGsfDR2>& gsf2pfgsf,
					 std::vector<GsfToEleDR2>& gsf2ele ) {

  for ( auto sig : signal_electrons ) {

    // SIG: Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = true;
    chain.is_egamma_ = false;
    chain.sig_ = sig;

    // TRK: Store matches between "signal electron" and KF tracks
    match<reco::Track>(sig,
		       sig2trk,
		       chain.trk_, // by ref
		       chain.trk_dr_,
		       chain.trk_match_ );
    
    // GSF: Store matches between "signal electron" and GSF tracks
    match<reco::GsfTrack>(sig,
			  sig2gsf,
			  chain.gsf_, // by ref
			  chain.gsf_dr_,
			  chain.gsf_match_ );

    // ELE: No matching between signal electron and GSF electron here, done later

    // GSF: If not matched to GsfTrack, then move onto next "signal electron"
    if ( !chain.gsf_match_ ) { continue; }

    // TRK: Update Track info
    reco::TrackPtr trk; 
    if ( gsfToTrk(chain.gsf_,trk) ) {
      chain.trk_ = trk; 
      chain.trk_match_ = true;
      chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
      PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
      if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
    }

    // SEED: Store Seed information
    if ( isAOD_ == 0 ) {
      // If no TrackExtra info, then assume tracker-driven
      chain.seed_tracker_driven_ = true;
      chain.seed_ecal_driven_ = false;
    } else if ( isAOD_ == 1 ) {
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) { // Store ElectronSeed info
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
	reco::CaloClusterPtr calo;
	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	// Store PreId info
	//chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	//chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
      } else {
	// If no TrackExtra info, then assume tracker-driven
	chain.seed_tracker_driven_ = true;
	chain.seed_ecal_driven_ = false;
      }
    }
    
    // SEED: Store ElectronSeed BDT discrimator outputs
    chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
    chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

    // PFGSF: Store PF GSF track info if match found with GSF track
    auto match_gsf_to_pfgsf = std::find_if( gsf2pfgsf.begin(),
					    gsf2pfgsf.end(),
					    [chain](const GsfToGsfDR2& dr2) {
					      return chain.gsf_ == dr2.obj1_;
					    }
					    );
    if ( match_gsf_to_pfgsf != gsf2pfgsf.end() && 
	 validPtr(match_gsf_to_pfgsf->obj2_) ) {
      chain.pfgsf_ = match_gsf_to_pfgsf->obj2_;
      chain.pfgsf_match_ = true;
      chain.pfgsf_dr_ = sqrt(deltaR2(chain.sig_,chain.pfgsf_));
    }
    
    // PFGSF: No check here on if PF GSF is found, as this is for information only!

    // ELE: Store GSF electron info if match found with GSF track
    auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
					  gsf2ele.end(), 
					  [chain](const GsfToEleDR2& dr2) { 
					    return chain.gsf_ == dr2.obj1_;
					  }
					  );
    if ( match_gsf_to_ele != gsf2ele.end() && 
	 validPtr(match_gsf_to_ele->obj2_) ) { 
      chain.ele_ = match_gsf_to_ele->obj2_; 
      chain.ele_match_ = true;
      chain.ele_dr_ = sqrt(deltaR2(chain.sig_,chain.ele_));
    }
    
  } // for ( auto sig : signal_electrons )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "fake tracks" 
void IDNtuplizer::lowPtElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
					std::vector<TrkToGsfDR2>& trk2gsf,
					std::vector<GsfToGsfDR2>& gsf2pfgsf,
					std::vector<GsfToEleDR2>& gsf2ele ) {

  // Iterate through tracks
  for ( auto iter : other_trk ) {

    if ( !validPtr(iter.obj2_) ) { continue; } // Shouldn't happen

    // SIG: Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = false;
    chain.is_egamma_ = false;
    
    // TRK: Store Track info
    chain.trk_ = iter.obj2_;
    chain.trk_match_ = true;
    chain.trk_dr_ = 0.; // w.r.t. trk
     
    // GSF: Store GSF track info if match found with KF track (either via ElectronSeed or just deltaR)
    auto match_trk_to_gsf = std::find_if( trk2gsf.begin(), 
					  trk2gsf.end(), 
					  [chain](const TrkToGsfDR2& dr2) { 
					    return chain.trk_ == dr2.obj1_;
					  }
					  );
    if ( match_trk_to_gsf != trk2gsf.end() && 
	 validPtr(match_trk_to_gsf->obj2_) ) {
      chain.gsf_ = match_trk_to_gsf->obj2_;
      chain.gsf_match_ = true;
      chain.gsf_dr_ = sqrt(deltaR2(chain.trk_,chain.gsf_));
    }

    // GSF: If not matched to GsfTrack, then move on to next "fake candidate"
    if ( !chain.gsf_match_ ) { continue; }
    // Otherwise ... 

    // SEED: Store Seed information
    if ( isAOD_ == 0 ) {
      // If no TrackExtra info, then assume tracker-driven
      chain.seed_tracker_driven_ = true;
      chain.seed_ecal_driven_ = false;
    } else if ( isAOD_ == 1 ) {
      reco::ElectronSeedPtr seed;
      if ( gsfToSeed(chain.gsf_,seed) ) { // Store ElectronSeed info
	chain.seed_ = seed;
	chain.seed_tracker_driven_ = seed->isTrackerDriven();
	chain.seed_ecal_driven_ = seed->isEcalDriven();
	reco::CaloClusterPtr calo;
	if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	// Store PreId info
	//chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
	//chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
      } else {
	// If no TrackExtra info in AOD, then assume tracker-driven
	chain.seed_tracker_driven_ = true;
	chain.seed_ecal_driven_ = false;
      }
    }

    // GSF: Store ElectronSeed BDT discrimator outputs
    chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
    chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];

    // PFGSF: Store PF GSF track info if match found with GSF track
    auto match_gsf_to_pfgsf = std::find_if( gsf2pfgsf.begin(),
					    gsf2pfgsf.end(),
					    [chain](const GsfToGsfDR2& dr2) {
					      return chain.gsf_ == dr2.obj1_;
					    }
					    );
    if ( match_gsf_to_pfgsf != gsf2pfgsf.end() && 
	 validPtr(match_gsf_to_pfgsf->obj2_) ) {
      chain.pfgsf_ = match_gsf_to_pfgsf->obj2_;
      chain.pfgsf_match_ = true;
      chain.pfgsf_dr_ = sqrt(deltaR2(chain.trk_,chain.pfgsf_));
    }
    
    // ELE: Store GSF electron info if match found with GSF track
    auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
					  gsf2ele.end(), 
					  [chain](const GsfToEleDR2& dr2) { 
					    return chain.gsf_ == dr2.obj1_;
					  }
					  );
    if ( match_gsf_to_ele != gsf2ele.end() && 
	 validPtr(match_gsf_to_ele->obj2_) ) { 
      chain.ele_ = match_gsf_to_ele->obj2_; 
      chain.ele_match_ = true;
      chain.ele_dr_ = sqrt(deltaR2(chain.trk_,chain.ele_));
    }
    
  } // for ( auto iter : other_trk )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
void IDNtuplizer::pfElectrons( std::set<reco::CandidatePtr>& signal_electrons,
			       std::vector<SigToTrkDR2>& sig2trk,
			       std::vector<SigToTrkDR2>& other_trk,
			       std::vector<SigToGsfDR2>& sig2gsf,
			       std::vector<TrkToGsfDR2>& trk2gsf,
			       std::vector<GsfToGsfDR2>& gsf2pfgsf ) {

  // Match "signal" electrons to PF GsfTracks
  std::vector<SigToGsfDR2> sig2pfgsf;
  std::vector<SigToGsfDR2> other_pfgsf;
  sigToCandLinks<reco::GsfTrack>( signal_electrons, 
				  gsfTracksEGammaH_,
				  sig2pfgsf, 
				  other_pfgsf );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::pfElectrons] sigToCandLinks<reco::GsfTrack>:" << std::endl
	      << " signal_electrons.size(): " << signal_electrons.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " sig2pfgsf.size(): " << sig2pfgsf.size() << std::endl
	      << " other_pfgsf.size(): " << other_pfgsf.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : sig2pfgsf ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }

  // Match Tracks (including surrogates) to PF GsfTracks
  std::vector<TrkToGsfDR2> trk2pfgsf;
  trkToGsfLinks( tracks_, 
		 gsfTracksEGammaH_, 
		 trk2pfgsf,
		 true ); // is_egamma
  if ( verbose_ > 1 ) {
    int good = 0, surrogate = 0, broken = 0, no_gsf = 0, no_trk = 0, empty = 0;
    for ( auto iter : trk2pfgsf ) { 
      if ( validPtr(iter.obj1_) && validPtr(iter.obj2_) ) {
	if ( iter.dr2_ >= 0. )                       { good++; } 
	else if ( std::abs(iter.dr2_-IDNtuple::NEG_FLOATSQ) > 1.e-6 ) { surrogate++; } 
	else                                         { broken++; } 
      } 
      else if ( validPtr(iter.obj1_) && !validPtr(iter.obj2_) ) { no_gsf++; }
      else if ( !validPtr(iter.obj1_) && validPtr(iter.obj2_) ) { no_trk++; }
      else                                                      { empty++; }
    }
    std::cout << "[IDNtuplizer::pfElectrons] trkToPFGsfLinks:" << std::endl
	      << " tracks_.size():      " << tracks_.size() << std::endl
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " trk2pfgsf.size():      " << trk2pfgsf.size() << std::endl
	      << " #good:      " << good << std::endl
	      << " #surrogate: " << surrogate << std::endl
	      << " #no_trk:    " << no_trk << std::endl
	      << " #no_gsf:    " << no_gsf << std::endl
	      << " #broken:    " << broken << std::endl
	      << " #empty:     " << empty << std::endl
	      << " #total:     " << good+surrogate+no_trk+no_gsf+broken+empty << std::endl;
  }

  // Match PF GsfTracks to PF GsfElectrons
  std::vector<GsfToEleDR2> pfgsf2ele;
  gsfToEleLinks( gsfTracksEGammaH_, 
		 gsfElectronsEGammaH_, 
		 pfgsf2ele );
  if ( verbose_ > 1 ) {
    std::cout << "[IDNtuplizer::pfElectrons] gsfToEleLinks:" << std::endl 
	      << " gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << " gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl
	      << " pfgsf2ele.size(): " << pfgsf2ele.size() << std::endl;
    if ( verbose_ > 2 ) {
      for ( auto iter : pfgsf2ele ) { if ( iter.dr2_ >= 0. ) { std::cout << iter << std::endl; } }
    }
    std::cout << std::endl;
  }
  
  pfElectrons_signal( signal_electrons, sig2trk, other_trk, sig2gsf, sig2pfgsf, trk2pfgsf, gsf2pfgsf, pfgsf2ele );
  pfElectrons_fakes( other_trk, trk2gsf, trk2pfgsf, gsf2pfgsf, pfgsf2ele );
  
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "signal electrons" 
void IDNtuplizer::pfElectrons_signal( std::set<reco::CandidatePtr>& signal_electrons,
				      std::vector<SigToTrkDR2>& sig2trk,
				      std::vector<SigToTrkDR2>& other_trk,
				      std::vector<SigToGsfDR2>& sig2gsf,
				      std::vector<SigToGsfDR2>& sig2pfgsf,
			       	      std::vector<TrkToGsfDR2>& trk2pfgsf,
			       	      std::vector<GsfToGsfDR2>& gsf2pfgsf,
				      std::vector<GsfToEleDR2>& pfgsf2ele ) {

  for ( auto sig : signal_electrons ) {

    // SIG: Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = true;
    chain.is_egamma_ = true;
    chain.sig_ = sig;

    // TRK: Store matches between "signal electron" and KF tracks
    match<reco::Track>(sig,
		       sig2trk,
		       chain.trk_, // by ref
		       chain.trk_dr_,
		       chain.trk_match_ );
    
    // GSF: Store matches between "signal electron" and (low-pT!) GSF tracks
    match<reco::GsfTrack>(sig,
			  sig2gsf,
			  chain.gsf_, // by ref
			  chain.gsf_dr_,
			  chain.gsf_match_ );
	
    // PFGSF: Store matches between "signal electron" and PF GSF tracks
    match<reco::GsfTrack>(sig, 
			  sig2pfgsf,
			  chain.pfgsf_, // by ref 
			  chain.pfgsf_dr_,
			  chain.pfgsf_match_ );

    // ELE: No matching between signal electron and GSF electron here, done later

    // GSF: Check if matched to a (low-pT!) GsfTrack
    if ( chain.gsf_match_ ) {

      // SEED: Store ElectronSeed BDT discrimator outputs
      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
      
      // PFGSF: Store PF GSF track info if match found with GSF track
      auto match_gsf_to_pfgsf = std::find_if( gsf2pfgsf.begin(),
					      gsf2pfgsf.end(),
					      [chain](const GsfToGsfDR2& dr2) {
						return chain.gsf_ == dr2.obj1_;
					      }
					      );
      if ( match_gsf_to_pfgsf != gsf2pfgsf.end() && 
	   validPtr(match_gsf_to_pfgsf->obj2_) ) {

	// TRK: Update Track info (i.e. calc deltaR with seed track, probably identical to original trk...)
	reco::TrackPtr trk; 
	if ( gsfToTrk(chain.pfgsf_,trk,true) ) {
	  chain.trk_ = trk; 
	  chain.trk_match_ = true;
	  chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	}

	// SEED: Store Seed information associated with PF GSF track
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	  } else {
	    // If no TrackExtra info in AOD, then assume tracker-driven
	    chain.seed_tracker_driven_ = true;
	    chain.seed_ecal_driven_ = false;
	  }
	}

	// PFGSF: Store PF GSF track
	chain.pfgsf_ = match_gsf_to_pfgsf->obj2_;
	chain.pfgsf_match_ = true;
	chain.pfgsf_dr_ = sqrt(deltaR2(chain.sig_,chain.pfgsf_));

      } else {
	// PFGSF: If do not match GSF to PFGSF, then deal with this below
      }
      
    } else if  ( chain.pfgsf_match_ ) {

      // PFGSF: Check if PFGSF track is NOT matched to GSF track
      auto match_pfgsf_to_gsf = std::find_if( gsf2pfgsf.begin(),
					      gsf2pfgsf.end(),
					      [chain](const GsfToGsfDR2& dr2) {
						return chain.pfgsf_ == dr2.obj2_;
					      }
					      );
      if ( match_pfgsf_to_gsf != gsf2pfgsf.end() && // Find PFGSF entry
	   validPtr(match_pfgsf_to_gsf->obj2_) &&   // Ensure nonnull PFGSF
	   !validPtr(match_pfgsf_to_gsf->obj1_) ) { // Ensure null GSF
      
	// TRK: Check if PF GSF track is matched to a surrogate PF track (closest in pT)
	reco::GsfTrackPtr pfgsf = match_pfgsf_to_gsf->obj2_;
	auto match_pfgsf_to_trk = std::find_if( trk2pfgsf.begin(), 
						trk2pfgsf.end(), 
						[pfgsf](const TrkToGsfDR2& dr2) { 
						  return pfgsf == dr2.obj2_; 
						}
						);
	if ( match_pfgsf_to_trk == trk2pfgsf.end() &&
	     validPtr(match_pfgsf_to_trk->obj1_) ) {
	  chain.trk_ = match_pfgsf_to_trk->obj1_; // Store surrogate track
	  chain.trk_match_ = true;
	  chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
	  PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
	  if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
	}

	// SEED: Store Seed information
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume ECAL-driven
	  chain.seed_tracker_driven_ = false;
	  chain.seed_ecal_driven_ = true;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	  } else {
	    // If no TrackExtra info in AOD, then assume ECAL-driven!
	    chain.seed_tracker_driven_ = false;
	    chain.seed_ecal_driven_ = true;
	  }
	}

	// PFGSF: Information is already stored

      } else {
	std::cout << "[IDNtuplizer::pfElectrons_signal] " 
		  << "ERROR: Couldn't find PF GSF track in gsf2pfgsf map!";
      }
      
    }
    
    // GSF: No check here on if GSF track is found, as this we care only about PF GSF!

    // PFGSF: If not matched to PF GsfTrack, then move onto next "signal electron"
    if ( !chain.pfgsf_match_ ) { continue; }
    // Otherwise ... 

    // ELE: Store GSF electron info if match found with PF GSF track
    auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
					    pfgsf2ele.end(), 
					    [chain](const GsfToEleDR2& dr2) { 
					      return chain.pfgsf_ == dr2.obj1_;
					    }
					    );
    if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
	 validPtr(match_pfgsf_to_ele->obj2_) ) { 
      chain.ele_ = match_pfgsf_to_ele->obj2_; 
      chain.ele_match_ = true;
      chain.ele_dr_ = sqrt(deltaR2(chain.sig_,chain.ele_));
    }
    
  } // for ( auto sig : signal_electrons )

}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Iterate through "fake tracks" 
void IDNtuplizer::pfElectrons_fakes( std::vector<SigToTrkDR2>& other_trk,
				     std::vector<TrkToGsfDR2>& trk2gsf,
				     std::vector<TrkToGsfDR2>& trk2pfgsf,
				     std::vector<GsfToGsfDR2>& gsf2pfgsf,
				     std::vector<GsfToEleDR2>& pfgsf2ele ) {
  
  // Iterate through tracks
  for ( auto iter : other_trk ) {
    
    if ( !validPtr(iter.obj2_) ) { continue; } // Shouldn't happen
    
    // SIG: Initialise ElectronChain object
    chains_.push_back(ElectronChain());
    ElectronChain& chain = chains_.back();
    chain.is_mc_ = isMC_;
    chain.is_aod_ = isAOD_;
    chain.is_e_ = false;
    chain.is_egamma_ = true;
    
    // TRK: Store Track info
    chain.trk_ = iter.obj2_;
    chain.trk_match_ = true;
    chain.trk_dr_ = 0.; // w.r.t. trk
    
    // GSF: Store GSF track info if match found with KF track (either via ElectronSeed or just deltaR)
    auto match_trk_to_gsf = std::find_if( trk2gsf.begin(), 
					  trk2gsf.end(), 
					  [chain](const TrkToGsfDR2& dr2) { 
					    return dr2.obj1_ == chain.trk_;
					  }
					  );
    if ( match_trk_to_gsf != trk2gsf.end() && 
	 validPtr(match_trk_to_gsf->obj2_) ) {

      // GSF: Store GSF track info
      chain.gsf_ = match_trk_to_gsf->obj2_;
      chain.gsf_match_ = true;
      chain.gsf_dr_ = sqrt(deltaR2(chain.trk_,chain.gsf_));
      
      // SEED: Store ElectronSeed BDT discrimator outputs
      chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
      chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
      
      // PFGSF: Store PF GSF track info if match found with GSF track
      auto match_gsf_to_pfgsf = std::find_if( gsf2pfgsf.begin(),
					      gsf2pfgsf.end(),
					      [chain](const GsfToGsfDR2& dr2) {
						return chain.gsf_ == dr2.obj1_;
					      }
					      );
      if ( match_gsf_to_pfgsf != gsf2pfgsf.end() && 
	   validPtr(match_gsf_to_pfgsf->obj2_) ) {

	// TRK: KF track info already stored

	// PFGSF: Store PF GSF track
	chain.pfgsf_ = match_gsf_to_pfgsf->obj2_;
	chain.pfgsf_match_ = true;
	chain.pfgsf_dr_ = sqrt(deltaR2(chain.pfgsf_,chain.trk_));

	// SEED: Store Seed information associated with PF GSF track
	if ( isAOD_ == 0 ) {
	  // If no TrackExtra info, then assume tracker-driven
	  chain.seed_tracker_driven_ = true;
	  chain.seed_ecal_driven_ = false;
	} else if ( isAOD_ == 1 ) {
	  reco::ElectronSeedPtr seed;
	  if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	    chain.seed_ = seed;
	    chain.seed_tracker_driven_ = seed->isTrackerDriven();
	    chain.seed_ecal_driven_ = seed->isEcalDriven();
	    reco::CaloClusterPtr calo;
	    if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	  } else {
	    // If no TrackExtra info in AOD, then assume tracker-driven
	    chain.seed_tracker_driven_ = true;
	    chain.seed_ecal_driven_ = false;
	  }
	}

	// ELE: Store GSF electron info if match found with PF GSF track
	auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
						pfgsf2ele.end(), 
						[chain](const GsfToEleDR2& dr2) { 
						  return chain.pfgsf_ == dr2.obj1_;
						}
						);
	if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
	     validPtr(match_pfgsf_to_ele->obj2_) ) { 
	  chain.ele_ = match_pfgsf_to_ele->obj2_; 
	  chain.ele_match_ = true;
	  chain.ele_dr_ = sqrt(deltaR2(chain.trk_,chain.ele_));
	}

      } // if PF GSF matched to GSF then ...
    } // if GSF matched to TRK then ...
  } // for ( auto iter : other_trk )

  // Loop through PF GSF tracks with no match to a GSF track and a surrogate KF track
  for ( auto match_trk_to_pfgsf : trk2pfgsf ) {
    reco::TrackPtr trk = match_trk_to_pfgsf.obj1_;
    reco::GsfTrackPtr pfgsf = match_trk_to_pfgsf.obj2_;
    
    // PFGSF: Check if PF GSF track is NOT matched to GSF track
    auto match_pfgsf_to_gsf = std::find_if( gsf2pfgsf.begin(),
					    gsf2pfgsf.end(),
					    [pfgsf](const GsfToGsfDR2& dr2) {
					      return pfgsf == dr2.obj2_;
					    }
					    );
    if ( match_pfgsf_to_gsf != gsf2pfgsf.end() &&
	 validPtr(match_pfgsf_to_gsf->obj2_) &&
	 validPtr(match_pfgsf_to_gsf->obj1_) ) { continue; }
      
    // Check if PF GSF match to surrogate TRK is found 
    if ( validPtr(trk) && // Valid TRK Ptr
	 validPtr(pfgsf) && // Valid PF GSF Ptr
	 std::abs( match_trk_to_pfgsf.dr2_-IDNtuple::NEG_FLOATSQ) > 1.e-6 && // Implies valid match
	 match_trk_to_pfgsf.dr2_ < 0. ) { // Implies surrogate trk
	
      // SIG: Initialise ElectronChain object
      chains_.push_back(ElectronChain());
      ElectronChain& chain = chains_.back();
      chain.is_mc_ = isMC_;
      chain.is_aod_ = isAOD_;
      chain.is_e_ = false;
      chain.is_egamma_ = true;
      
      // Store (surrogate) TRK info
      chain.trk_ = trk; // Store surrogate track
      chain.trk_match_ = true;
      chain.trk_dr_ = 0.; // w.r.t. trk

      // GSF: No match by definition, so not set

      // PFGSF: Store PF GSF track
      chain.pfgsf_ = pfgsf;
      chain.pfgsf_match_ = true;
      chain.pfgsf_dr_ = sqrt(deltaR2(chain.trk_,chain.pfgsf_));
      
      // SEED: Store Seed information associated with PF GSF track
      if ( isAOD_ == 0 ) {
	// If no TrackExtra info, then assume ECAL-driven
	chain.seed_tracker_driven_ = false;
	chain.seed_ecal_driven_ = true;
      } else if ( isAOD_ == 1 ) {
	reco::ElectronSeedPtr seed;
	if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
	  chain.seed_ = seed;
	  chain.seed_tracker_driven_ = seed->isTrackerDriven();
	  chain.seed_ecal_driven_ = seed->isEcalDriven();
	  reco::CaloClusterPtr calo;
	  if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
	} else {
	  // If no TrackExtra info in AOD, then assume ECAL-driven!
	  chain.seed_tracker_driven_ = false;
	  chain.seed_ecal_driven_ = true;
	}
      }

      // ELE: Store GSF electron info if match found with PF GSF track
      auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
					      pfgsf2ele.end(), 
					      [chain](const GsfToEleDR2& dr2) { 
						return chain.pfgsf_ == dr2.obj1_;
					      }
					      );
      if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
	   validPtr(match_pfgsf_to_ele->obj2_) ) { 
	chain.ele_ = match_pfgsf_to_ele->obj2_; 
	chain.ele_match_ = true;
	chain.ele_dr_ = sqrt(deltaR2(chain.trk_,chain.ele_));
      }

    } // if PF GSF matched to (surrogate) TRK then ...

  } // trk2pfgsf loop

////////////////////
////////////////////
////////////////////
//    } else {
//
//      std::cout << "!!!!! NO GSF MATCH !!!!!" << std::endl;
//
//      // TRK: Check if PF GSF track is matched to a surrogate PF track (closest in pT)
//      reco::GsfTrackPtr pfgsf = match_pfgsf_to_gsf->obj2_;
//      auto match_gsf_to_trk = std::find_if( trk2pfgsf.begin(), 
//					    trk2pfgsf.end(), 
//					    [pfgsf](const TrkToGsfDR2& dr2) { 
//					      return pfgsf == dr2.obj2_; 
//					    }
//					    );
//      if ( match_gsf_to_trk == trk2gsf.end() &&
//	   validPtr(match_gsf_to_trk->obj1_) ) {
//	chain.trk_ = match_gsf_to_trk->obj1_; // Store surrogate track
//	chain.trk_match_ = true;
//	chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
//	PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//      }
//
////      // PFGSF: Check if PFGSF track IS matched to KF track
////      auto match_trk_to_pfgsf = std::find_if( trk2pfgsf.begin(), 
////					      trk2pfgsf.end(), 
////					      [chain](const TrkToGsfDR2& dr2) { 
////						return dr2.obj1_ == chain.trk_;
////					      }
////					      );
////      if ( match_trk_to_pfgsf != trk2pfgsf.end() &&
////	   validPtr(match_trk_to_pfgsf->obj2_) ) { 
//
//	std::cout << "!!!!! TRK MATCHED TO PFGSF !!!!!" << std::endl;
//      
//	// PFGSF: Check if PFGSF track is NOT matched to GSF track
//	auto match_pfgsf_to_gsf = std::find_if( gsf2pfgsf.begin(),
//						gsf2pfgsf.end(),
//						[chain](const GsfToGsfDR2& dr2) {
//						  return chain.pfgsf_ == dr2.obj2_;
//						}
//						);
//	if ( match_pfgsf_to_gsf != gsf2pfgsf.end() && // Find PFGSF entry
//	     validPtr(match_pfgsf_to_gsf->obj2_) &&   // Ensure nonnull PFGSF
//	     !validPtr(match_pfgsf_to_gsf->obj1_) ) { // Ensure null GSF
//
//	  std::cout << "!!!!! PFGSF NOT MATCHED TO GSF !!!!!" << std::endl;
//	  
//	  // TRK: KF track info already stored
//	  
//	  // PFGSF: Store PF GSF track
//	  chain.pfgsf_ = match_pfgsf_to_gsf->obj2_;
//	  chain.pfgsf_match_ = true;
//	  chain.pfgsf_dr_ = sqrt(deltaR2(chain.trk_,chain.pfgsf_));
//	  
//	  // SEED: Store Seed information associated with PF GSF track
//	  if ( isAOD_ == 0 ) {
//	    // If no TrackExtra info, then assume ECAL-driven
//	    chain.seed_tracker_driven_ = false;
//	    chain.seed_ecal_driven_ = true;
//	  } else if ( isAOD_ == 1 ) {
//	    reco::ElectronSeedPtr seed;
//	    if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
//	      chain.seed_ = seed;
//	      chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	      chain.seed_ecal_driven_ = seed->isEcalDriven();
//	      reco::CaloClusterPtr calo;
//	      if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
//	    } else {
//	      // If no TrackExtra info in AOD, then assume ECAL-driven!
//	      chain.seed_tracker_driven_ = false;
//	      chain.seed_ecal_driven_ = true;
//	    }
//	  }
//	  
//	} else {
//	  std::cout << "!!!!! PFGSF is matched to GSF !!!!!" << std::endl;
//	}
//      } else {
//	std::cout << "!!!!! PFGSF is not matched to TRK !!!!!" << std::endl;
//      }
//      
//    } // if GSF matched to TRK then ...
//
//    // PFGSF: If not matched to KF track, then move onto next "fake candidate"
//    if ( !chain.pfgsf_match_ ) { continue; }
//
//    // ELE: Store GSF electron info if match found with PF GSF track
//    auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
//					    pfgsf2ele.end(), 
//					    [chain](const GsfToEleDR2& dr2) { 
//					      return chain.pfgsf_ == dr2.obj1_;
//					    }
//					    );
//    if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
//	 validPtr(match_pfgsf_to_ele->obj2_) ) { 
//      chain.ele_ = match_pfgsf_to_ele->obj2_; 
//      chain.ele_match_ = true;
//      chain.ele_dr_ = sqrt(deltaR2(chain.trk_,chain.ele_));
//    }
//
//  } // for ( auto iter : other_trk )

}


////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//// Iterate through "signal electrons" 
//void IDNtuplizer::signal( std::set<reco::CandidatePtr>& signal_electrons,
//			  std::vector<SigToTrkDR2>& sig2trk,
//			  std::vector<SigToGsfDR2>& sig2gsf,
//			  std::vector<SigToGsfDR2>& sig2pfgsf,
//			  std::vector<SigToTrkDR2>& other_trk,
//			  std::vector<TrkToGsfDR2>& trk2gsf,
//			  std::vector<TrkToGsfDR2>& trk2pfgsf,
//			  std::vector<GsfToGsfDR2>& gsf2pfgsf,
//			  std::vector<GsfToEleDR2>& gsf2ele,
//			  std::vector<GsfToEleDR2>& pfgsf2ele ) {
//
//  // Iterate through signal electrons
//  for ( auto sig : signal_electrons ) {
//    
//    // Repeat for two options: low pT and PF EGamma reconstruction
//    for ( auto is_egamma : std::vector<bool>{ false, true } ) { 
//      
//      // SIG: Initialise ElectronChain object
//      chains_.push_back(ElectronChain());
//      ElectronChain& chain = chains_.back();
//      chain.is_mc_ = isMC_;
//      chain.is_aod_ = isAOD_;
//      chain.is_e_ = true;
//      chain.is_egamma_ = is_egamma; // Set here if low pT or EGamma!
//      chain.sig_ = sig;
//      
//      // TRK: Store matches between "signal electron" and KF tracks
//      match<reco::Track>(sig,
//			 sig2trk,
//			 chain.trk_, // by ref
//			 chain.trk_dr_,
//			 chain.trk_match_ );
//      
//      // GSF: Store matches between "signal electron" and GSF tracks
//      match<reco::GsfTrack>(sig,
//			    sig2gsf,
//			    chain.gsf_, // by ref
//			    chain.gsf_dr_,
//			    chain.gsf_match_ );
//	
//      // PFGSF: Store matches between "signal electron" and PF GSF tracks
//      match<reco::GsfTrack>(sig, 
//			    sig2pfgsf,
//			    chain.pfgsf_, // by ref 
//			    chain.pfgsf_dr_,
//			    chain.pfgsf_match_ );
//      
//      // ELE: No matching between signal electron and GSF electron
//
//      // GSF: Check if signal electron is matched to GsfTrack
//      if ( chain.gsf_match_ ) { 
//
//	// Only for low-pT GSF tracks ...
//	if ( is_egamma == false ) {
//	  
//	  // TRK: *Update* Track info
//	  reco::TrackPtr trk; 
//	  if ( gsfToTrk(chain.gsf_,trk) ) {
//	    chain.trk_ = trk; 
//	    chain.trk_match_ = true;
//	    chain.trk_dr_ = sqrt(deltaR2(chain.sig_,chain.trk_));
//	    PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	    if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//	  }
//	  
//	  // SEED: Store Seed information
//	  if ( isAOD_ == 0 ) {
//	    // If no TrackExtra info, then assume tracker-driven
//	    chain.seed_tracker_driven_ = true;
//	    chain.seed_ecal_driven_ = false;
//	  } else if ( isAOD_ == 1 ) {
//	    reco::ElectronSeedPtr seed;
//	    if ( gsfToSeed(chain.gsf_,seed) ) { // Store ElectronSeed info
//	      chain.seed_ = seed;
//	      chain.seed_tracker_driven_ = seed->isTrackerDriven();
//	      chain.seed_ecal_driven_ = seed->isEcalDriven();
//	      reco::CaloClusterPtr calo;
//	      if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
//	      // Store PreId info
//	      //chain.preid_ecal_ = edm::refToPtr((*preIdRefsH_)[seed->ctfTrack()]);
//	      //chain.preid_hcal_ = reco::PreIdPtr( preIdsHcalH_, chain.preid_ecal_.key() );
//	    } else {
//	      // If no TrackExtra info, then assume tracker-driven
//	      chain.seed_tracker_driven_ = true;
//	      chain.seed_ecal_driven_ = false;
//	    }
//	  }
//
//	} // !is_egamma
//      
//	// SEED: Store BDT discrimator outputs for low-pT ElectronSeeds
//	chain.unbiased_ = (*mvaUnbiasedH_)[chain.gsf_];
//	chain.ptbiased_ = (*mvaPtbiasedH_)[chain.gsf_];
//
//	// PFGSF: Store PF GSF track info if match found with GSF track
//	auto match_gsf_to_pfgsf = std::find_if( gsf2pfgsf.begin(),
//						gsf2pfgsf.end(),
//						[chain](const GsfToGsfDR2& dr2) {
//						  return chain.gsf_ == dr2.obj1_;
//						}
//						);
//	if ( match_gsf_to_pfgsf != gsf2pfgsf.end() && 
//	     validPtr(match_gsf_to_pfgsf->obj2_) ) {
//
//	  // PFGSF: Store PF GSF track info
//	  chain.pfgsf_ = match_gsf_to_pfgsf->obj2_;
//	  chain.pfgsf_match_ = true;
//	  chain.pfgsf_dr_ = sqrt(deltaR2(chain.sig_,chain.pfgsf_));
//
//	  // Only for PF GSF tracks ...
//	  if ( is_egamma == true ) {
//
//	    // TRK: Update Track info (with surrogate track if necessary)
//	    auto match_pfgsf_to_trk = std::find_if( trk2pfgsf.begin(), 
//						    trk2pfgsf.end(), 
//						    [chain](const TrkToGsfDR2& dr2) { 
//						      return chain.pfgsf_ == dr2.obj2_; 
//						    }
//						    );
//	    if ( match_pfgsf_to_trk == trk2pfgsf.end() &&
//		 validPtr(match_pfgsf_to_trk->obj1_) ) {
//	      chain.trk_ = match_pfgsf_to_trk->obj1_; // Store surrogate track
//	      chain.trk_match_ = true;
//	      chain.trk_dr_ = sqrt(deltaR2(chain.trk_,chain.sig_));
//	      PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	      if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//	    }
//
//	    // SEED: Store Seed information
//	    if ( isAOD_ == 0 ) {
//	      // If no TrackExtra info, then assume tracker-driven
//	      chain.seed_tracker_driven_ = true;
//	      chain.seed_ecal_driven_ = false;
//	    } else if ( isAOD_ == 1 ) {
//	      reco::ElectronSeedPtr seed;
//	      if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
//		chain.seed_ = seed;
//		chain.seed_tracker_driven_ = seed->isTrackerDriven();
//		chain.seed_ecal_driven_ = seed->isEcalDriven();
//		reco::CaloClusterPtr calo;
//		if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
//	      } else {
//		// If no TrackExtra info in AOD, then check if surrogate track
//		if ( chain.trk_dr_ >= 0. ) { 
//		  chain.seed_tracker_driven_ = true; 
//		  chain.seed_ecal_driven_ = false; 
//		} else if ( chain.trk_dr_ > IDNtuple::NEG_FLOATSQ ) { 
//		  chain.seed_tracker_driven_ = false; 
//		  chain.seed_ecal_driven_ = true; 
//		}
//	      }
//	    }
//	    
//	    // ELE: Store GSF electron info if match found with PF GSF track
//	    auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
//						    pfgsf2ele.end(), 
//						    [chain](const GsfToEleDR2& dr2) { 
//						      return chain.pfgsf_ == dr2.obj1_;
//						    }
//						    );
//	    if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
//		 validPtr(match_pfgsf_to_ele->obj2_) ) { 
//	      chain.ele_ = match_pfgsf_to_ele->obj2_; 
//	      chain.ele_match_ = true;
//	      chain.ele_dr_ = sqrt(deltaR2(chain.ele_,chain.sig_));
//	    }
//
//	  } // is_egamma
//	} // If PF GSF
//    
//	// ELE: Store GSF electron info if match found with GSF track
//	if ( is_egamma == false ) { 
//	  auto match_gsf_to_ele = std::find_if( gsf2ele.begin(), 
//						gsf2ele.end(), 
//						[chain](const GsfToEleDR2& dr2) { 
//						  return chain.gsf_ == dr2.obj1_;
//						}
//						);
//	  if ( match_gsf_to_ele != gsf2ele.end() && 
//	       validPtr(match_gsf_to_ele->obj2_) ) { 
//	    chain.ele_ = match_gsf_to_ele->obj2_; 
//	    chain.ele_match_ = true;
//	    chain.ele_dr_ = sqrt(deltaR2(chain.ele_,chain.sig_));
//	  }
//	}
//
//      } else if  ( chain.pfgsf_match_ ) { // i.e. if chain.gsf_match_ == false
//
//	// Only for PF GSF tracks ...
//	if ( is_egamma == true ) {
//	
//	  // PFGSF: Check if PFGSF track is NOT matched to GSF track
//	  auto match_pfgsf_to_gsf = std::find_if( gsf2pfgsf.begin(),
//						  gsf2pfgsf.end(),
//						  [chain](const GsfToGsfDR2& dr2) {
//						    return chain.pfgsf_ == dr2.obj2_;
//						  }
//						  );
//	  if ( match_pfgsf_to_gsf != gsf2pfgsf.end() && // Find PFGSF entry
//	       validPtr(match_pfgsf_to_gsf->obj2_) &&   // Ensure nonnull PFGSF
//	       !validPtr(match_pfgsf_to_gsf->obj1_) ) { // Ensure null GSF
//	    
//	    // TRK: Check if PF GSF track is matched to a surrogate PF track (closest in pT)
//	    auto match_pfgsf_to_trk = std::find_if( trk2pfgsf.begin(), 
//						    trk2pfgsf.end(), 
//						    [chain](const TrkToGsfDR2& dr2) { 
//						      return chain.pfgsf_ == dr2.obj2_; 
//						    }
//						    );
//	    if ( match_pfgsf_to_trk == trk2pfgsf.end() &&
//		 validPtr(match_pfgsf_to_trk->obj1_) ) {
//	      chain.trk_ = match_pfgsf_to_trk->obj1_; // Store surrogate track
//	      chain.trk_match_ = true;
//	      chain.trk_dr_ = sqrt(deltaR2(chain.trk_,chain.sig_));
//	      PdgIds::const_iterator pos = pdgids_.find(chain.trk_.key());
//	      if ( pos != pdgids_.end() ) { chain.pdg_id_ = pos->second; }
//	    }
//	    
//	    // SEED: Store Seed information
//	    if ( isAOD_ == 0 ) {
//	      // If no TrackExtra info, then assume ECAL-driven
//	      chain.seed_tracker_driven_ = false;
//	      chain.seed_ecal_driven_ = true;
//	    } else if ( isAOD_ == 1 ) {
//	      reco::ElectronSeedPtr seed;
//	      if ( gsfToSeed(chain.pfgsf_,seed) ) { // Store ElectronSeed info
//		chain.seed_ = seed;
//		chain.seed_tracker_driven_ = seed->isTrackerDriven();
//		chain.seed_ecal_driven_ = seed->isEcalDriven();
//		reco::CaloClusterPtr calo;
//		if ( seedToCalo(seed,calo) ) { chain.calo_ = calo; } // Store CaloCluster info
//	      } else {
//		// If no TrackExtra info in AOD, then check if surrogate track
//		if ( chain.trk_dr_ >= 0. ) { 
//		  chain.seed_tracker_driven_ = true; 
//		  chain.seed_ecal_driven_ = false; 
//		} else if ( chain.trk_dr_ > IDNtuple::NEG_FLOATSQ ) { 
//		  chain.seed_tracker_driven_ = false; 
//		  chain.seed_ecal_driven_ = true; 
//		}
//	      }
//	    }
//	    
//	    // ELE: Store GSF electron info if match found with PF GSF track
//	    auto match_pfgsf_to_ele = std::find_if( pfgsf2ele.begin(), 
//						    pfgsf2ele.end(), 
//						    [chain](const GsfToEleDR2& dr2) { 
//						      return chain.pfgsf_ == dr2.obj1_;
//						    }
//						    );
//	    if ( match_pfgsf_to_ele != pfgsf2ele.end() && 
//		 validPtr(match_pfgsf_to_ele->obj2_) ) { 
//	      chain.ele_ = match_pfgsf_to_ele->obj2_; 
//	      chain.ele_match_ = true;
//	      chain.ele_dr_ = sqrt(deltaR2(chain.ele_,chain.sig_));
//	    }
//	    
//	  } // PF GSF track not matched to GSF track
//	} // is_egamma
//      } // if ( chain.pfgsf_match_ ) ...
//      
//    } // for ( auto is_egamma : std::vector<bool>{ false, true } ) 
//  } // for ( auto sig : signal_electrons )
//
//}
