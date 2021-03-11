#ifndef LowPtElectrons_LowPtElectrons_IDSlimNtuple
#define LowPtElectrons_LowPtElectrons_IDSlimNtuple

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

// Small class to provide fillers and hide tree I/O
class IDSlimNtuple {

 public:

  static constexpr size_t NHITS_MAX = 30;
  static constexpr int NCLUS_MAX = 50;
  static constexpr int NEG_INT = -999;
  static constexpr float NEG_FLOAT = -999.;
  
  IDSlimNtuple() {}
  
  void reset() {
    IDSlimNtuple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );

  void fill_evt( const edm::EventID& id );
  void fill_gen( const pat::PackedGenParticleRef );   
  void fill_gen( const reco::GenParticlePtr );
  void fill_gen_default( );

  //void fill_trk( const reco::TrackPtr& trk,
  //const reco::BeamSpot& spot );
  void fill_trk( const reco::TrackPtr& trk );

  void fill_bdt( double seed_unbiased,
		 double seed_ptbiased );

  void fill_gsf( const reco::GsfTrackPtr trk );
  //void fill_gsf( const reco::GsfTrackPtr trk,
  //const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value, int mva_id,
		 const double rho, float seed_unbiased, float field_z);

  void fill_supercluster_miniAOD(const reco::GsfElectronPtr ele); 

  template < typename T> bool validPtr( edm::Ptr<T>& ptr);

 public:

  // Event
  unsigned long long evt_ = 0;
  float weight_ = 1.;   

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;  
  bool is_other_ = false;             
  bool is_egamma_ = false;
  bool has_trk_  = true;
  bool has_seed_ = true;
  bool has_gsf_  = true;
  bool has_ele_  = true;

  // GEN electrons
  float gen_dR_  = IDSlimNtuple::NEG_FLOAT;
  float genOther_dR_  = IDSlimNtuple::NEG_FLOAT;
  float gen_pt_  = IDSlimNtuple::NEG_FLOAT;
  float gen_eta_ = IDSlimNtuple::NEG_FLOAT;
  float gen_phi_ = IDSlimNtuple::NEG_FLOAT;
  float gen_p_   = IDSlimNtuple::NEG_FLOAT;
  int gen_tag_side_ = IDSlimNtuple::NEG_INT;  

  // GSF tracks: kine
  float gsf_pt_   = IDSlimNtuple::NEG_FLOAT;
  float gsf_eta_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_phi_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_p_    = IDSlimNtuple::NEG_FLOAT;

  // GSF tracks: kine (mode)
  float gsf_mode_pt_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_mode_eta_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_mode_phi_ = IDSlimNtuple::NEG_FLOAT;

  // GSF tracks: quality
  int gsf_missing_inner_hits_ = IDSlimNtuple::NEG_INT;

  // GSF tracks: tangents
  int gsf_ntangents_ = 0; 

  // Seed BDT discriminator values at GsfTrack level
  float seed_unbiased_ = IDSlimNtuple::NEG_FLOAT;
  float seed_ptbiased_ = IDSlimNtuple::NEG_FLOAT;

  // KF tracks: kine
  float trk_pt_  = IDSlimNtuple::NEG_FLOAT;
  float trk_eta_  = IDSlimNtuple::NEG_FLOAT;
  float trk_phi_  = IDSlimNtuple::NEG_FLOAT;
  float trk_p_    = IDSlimNtuple::NEG_INT;

  // GSF electrons: kinematics
  float ele_p_   = IDSlimNtuple::NEG_FLOAT;
  float ele_pt_  = IDSlimNtuple::NEG_FLOAT;
  float ele_eta_ = IDSlimNtuple::NEG_FLOAT;
  float ele_phi_ = IDSlimNtuple::NEG_FLOAT;
  int fiducial_isEB_      = IDSlimNtuple::NEG_INT;
  int fiducial_isEE_      = IDSlimNtuple::NEG_INT; 
  int fiducial_isEBEEGap_ = IDSlimNtuple::NEG_INT;
    
  // Electrons: IDs
  float ele_mva_value_ = IDSlimNtuple::NEG_FLOAT;
  int ele_mva_id_      = IDSlimNtuple::NEG_INT;

  float eid_rho_    = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_eta_ = IDSlimNtuple::NEG_FLOAT;
  float eid_shape_full5x5_r9_            = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_etaWidth_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_phiWidth_ = IDSlimNtuple::NEG_FLOAT;
  float eid_shape_full5x5_HoverE_ = IDSlimNtuple::NEG_FLOAT;
  float eid_trk_nhits_   = IDSlimNtuple::NEG_FLOAT;
  float eid_trk_chi2red_ = IDSlimNtuple::NEG_FLOAT;
  float eid_gsf_chi2red_ = IDSlimNtuple::NEG_FLOAT;
  float eid_brem_frac_   = IDSlimNtuple::NEG_FLOAT;
  float eid_gsf_nhits_   = IDSlimNtuple::NEG_FLOAT;
  float eid_match_SC_EoverP_   = IDSlimNtuple::NEG_FLOAT;
  float eid_match_eclu_EoverP_ = IDSlimNtuple::NEG_FLOAT;
  float eid_match_SC_dEta_   = IDSlimNtuple::NEG_FLOAT;
  float eid_match_SC_dPhi_   = IDSlimNtuple::NEG_FLOAT;
  float eid_match_seed_dEta_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_E_  = IDSlimNtuple::NEG_FLOAT;
  float eid_trk_p_ = IDSlimNtuple::NEG_FLOAT;
  float eid_gsf_mode_p_ = IDSlimNtuple::NEG_FLOAT;  
  float eid_core_shFracHits_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_gsf_bdtout1_ = IDSlimNtuple::NEG_FLOAT;  
  float eid_gsf_dr_ = IDSlimNtuple::NEG_FLOAT;    
  float eid_trk_dr_ = IDSlimNtuple::NEG_FLOAT;    
  float eid_sc_Nclus_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus1_nxtal_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_clus1_dphi_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus2_dphi_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus1_deta_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus2_deta_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus1_E_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_clus2_E_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_clus1_E_ov_p_ = IDSlimNtuple::NEG_FLOAT; 
  float eid_sc_clus2_E_ov_p_ = IDSlimNtuple::NEG_FLOAT;   

  // Electron, brem fractions
  float brem_fracTrk_ = IDSlimNtuple::NEG_FLOAT;
  float brem_fracSC_  = IDSlimNtuple::NEG_FLOAT;
  int brem_N_         = IDSlimNtuple::NEG_INT;

  // Energy regression 
  float pre_ecal_     = IDSlimNtuple::NEG_FLOAT;
  float pre_ecaltrk_  = IDSlimNtuple::NEG_FLOAT;
  float post_ecal_    = IDSlimNtuple::NEG_FLOAT;
  float post_ecaltrk_ = IDSlimNtuple::NEG_FLOAT;
  float sc_raw_energy_= IDSlimNtuple::NEG_FLOAT;
  float sc_energy_    = IDSlimNtuple::NEG_FLOAT;

  // SuperClusters 
  int sc_Nclus_deta01_ = IDSlimNtuple::NEG_INT;
  int sc_Nclus_deta02_ = IDSlimNtuple::NEG_INT;
  int sc_Nclus_deta03_ = IDSlimNtuple::NEG_INT;  
  bool sc_goodSeed_ = false;

  // When running on miniaod, only 3 highest ET clusters only
  float sc_clus1_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus1_nxtal_    = IDSlimNtuple::NEG_INT;
  //
  float sc_clus2_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus2_nxtal_    = IDSlimNtuple::NEG_INT;
  //
  float sc_clus3_E_      = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_E_ov_p_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus3_nxtal_    = IDSlimNtuple::NEG_INT;
  float sc_clus3_dphi_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_deta_   = IDSlimNtuple::NEG_FLOAT;

  // Distance wrt muons
  float minDrWithMu_ = IDSlimNtuple::NEG_FLOAT;
};

#endif // LowPtElectrons_LowPtElectrons_IDSlimNtuple
