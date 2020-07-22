#ifndef LowPtElectrons_LowPtElectrons_IDNtuple
#define LowPtElectrons_LowPtElectrons_IDNtuple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

constexpr size_t ARRAY_SIZE = 20;

// Small class to provide fillers and hide tree I/O
class IDNtuple {

 public:

  static constexpr size_t NHITS_MAX = 30;
  static constexpr int NEG_INT = -10;
  static constexpr float NEG_FLOAT = -10.;
  static constexpr float NEG_FLOATSQ = -1.*NEG_FLOAT*NEG_FLOAT;
  
  IDNtuple() {}

  void reset() {
    IDNtuple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );
  
  void set_weight( float w ) { weight_ = w; }
  void set_prescale( float p ) { prescale_ = p; }
  void set_rho( float r ) { rho_ = r; }

  void is_aod( int aod ) { is_aod_ = aod; }
  void is_mc( int mc ) { is_mc_ = mc; }

  void tag_pt( float x ) { tag_pt_ = x; }
  void tag_eta( float x ) { tag_eta_ = x; }

  void is_e( bool t = true ) { is_e_ = t; }
  void is_e_not_matched( bool t = true ) { is_e_not_matched_ = t; }
  void is_other( bool t = true ) { is_other_ = t; }
  void is_egamma( bool t = true ) { is_egamma_ = t; }

  void has_trk( bool f = false ) { has_trk_ = f; }
  void has_seed( bool f = false ) { has_seed_ = f; }
  void has_gsf( bool f = false ) { has_gsf_ = f; }
  void has_pfgsf( bool f = false ) { has_pfgsf_ = f; }
  void has_ele( bool f = false ) { has_ele_ = f; }

  void trk_dr( float dr ) { trk_dr_ = dr; }
  void gsf_dr( float dr ) { gsf_dr_ = dr; }
  void pfgsf_dr( float dr ) { pfgsf_dr_ = dr; }
  void ele_dr( float dr ) { ele_dr_ = dr; }

  void fill_evt( const edm::EventID& id );

  void fill_gen( const pat::PackedGenParticleRef );
  void fill_gen( const reco::GenParticlePtr ); //@@ AOD
  void fill_gen( const reco::CandidatePtr );

  void fill_trk( const reco::TrackPtr& trk,
		 const reco::BeamSpot& spot );

  void pdg_id( int id ) { pdg_id_ = id; }

  void fill_seed( bool seed_trk_driven, 
		  bool seed_ecal_driven );

   // to be deprecated
  void fill_seed( double seed_unbiased,
		  double seed_ptbiased );

  void fill_bdt( double seed_unbiased,
		 double seed_ptbiased );
  
  void fill_preid( const reco::PreId& preid_ecal,
		   const reco::PreId& preid_hcal,
		   const reco::BeamSpot& spot, 
		   const double rho, 
		   noZS::EcalClusterLazyTools& ecalTools );

  void fill_gsf( const reco::GsfTrackPtr gsf,
		 const reco::BeamSpot& spot );

  void fill_pfgsf( const reco::GsfTrackPtr pfgsf,
		 const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value,
		 float mva_value_retrained,
		 float mva_value_depth10,
		 float mva_value_depth15,
		 float ele_conv_vtx_fit_prob,
		 const double rho,
		 bool is_egamma = false,
		 float unbiased = 0. );
  
  void fill_image( const float gsf_ref_eta, const float gsf_ref_phi, const float gsf_ref_R,
		   const float gsf_ref_p, const float gsf_ref_pt,
		   const float gen_inner_eta, const float gen_inner_phi, const float gen_inner_R,
		   const float gen_inner_p, const float gen_inner_pt,
		   const float gen_proj_eta, const float gen_proj_phi, const float gen_proj_R,
		   const float gsf_inner_eta, const float gsf_inner_phi, const float gsf_inner_R,
		   const float gsf_inner_p, const float gsf_inner_pt, const int gsf_charge,
		   const float gsf_proj_eta, const float gsf_proj_phi, const float gsf_proj_R,
		   const float gsf_proj_p,
		   const float gsf_atcalo_eta, const float gsf_atcalo_phi, const float gsf_atcalo_R,
		   const float gsf_atcalo_p,
		   const std::vector<float>& clu_eta,
		   const std::vector<float>& clu_phi,
		   const std::vector<float>& clu_e,
		   const std::vector<int>& clu_nhit,
		   const std::vector<float>& pf_eta,
		   const std::vector<float>& pf_phi,
		   const std::vector<float>& pf_p,
		   const std::vector<int>& pf_pdgid,
		   const std::vector<int>& pf_matched,
		   const std::vector<int>& pf_lost );
  
 public:

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float prescale_ = 0.;
  float weight_ = 1.;
  float rho_ = IDNtuple::NEG_FLOAT;

  // Data sample
  int is_aod_ = -1;
  int is_mc_ = -1;

  // Tag-side muon
  float tag_pt_ = IDNtuple::NEG_FLOAT;
  float tag_eta_ = IDNtuple::NEG_FLOAT;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;
  bool is_other_ = false;
  bool is_egamma_ = false;

  // RECO steps
  bool has_trk_ = false;
  bool has_seed_ = false;
  bool has_gsf_ = false;
  bool has_pfgsf_ = false;
  bool has_ele_ = false;

  float trk_dr_ = IDNtuple::NEG_FLOAT;
  float gsf_dr_ = IDNtuple::NEG_FLOAT;
  float pfgsf_dr_ = IDNtuple::NEG_FLOAT;
  float ele_dr_ = IDNtuple::NEG_FLOAT;

  // GEN electrons
  float gen_pt_ = IDNtuple::NEG_FLOAT;
  float gen_eta_ = IDNtuple::NEG_FLOAT;
  float gen_phi_ = IDNtuple::NEG_FLOAT;
  float gen_e_ = IDNtuple::NEG_FLOAT;
  float gen_p_ = IDNtuple::NEG_FLOAT;
  int gen_charge_ = IDNtuple::NEG_INT;
  int gen_pdgid_ = 0;
  int gen_mom_pdgid_ = 0;
  int gen_gran_pdgid_ = 0;

  // KF tracks: kine
  float trk_pt_ = IDNtuple::NEG_FLOAT;
  float trk_eta_ = IDNtuple::NEG_FLOAT;
  float trk_phi_ = IDNtuple::NEG_FLOAT;
  float trk_p_ = IDNtuple::NEG_INT;
  int trk_charge_ = 0; //@@ IDNtuple::NEG_INT;
  float trk_inp_ = IDNtuple::NEG_FLOAT;
  float trk_outp_ = IDNtuple::NEG_FLOAT;
  float trk_dpt_ = IDNtuple::NEG_FLOAT;

  int pdg_id_ = 0;

  // KF tracks: quality
  int trk_nhits_ = IDNtuple::NEG_INT;
  int trk_missing_inner_hits_ = IDNtuple::NEG_INT;
  int trk_high_purity_ = 0; //@@ IDNtuple::NEG_INT;
  float trk_chi2red_ = IDNtuple::NEG_FLOAT;

  // KF tracks: displ
  float trk_dxy_ = IDNtuple::NEG_FLOAT;
  float trk_dxy_err_ = IDNtuple::NEG_FLOAT;
  float trk_dz_ = IDNtuple::NEG_FLOAT;
  float trk_dz_err_ = IDNtuple::NEG_FLOAT;
  
  // Seed BDT discriminator values
  float preid_unbiased_ = IDNtuple::NEG_FLOAT;
  float preid_ptbiased_ = IDNtuple::NEG_FLOAT;

  // Seed BDT discriminator values at GsfTrack level
  float seed_unbiased_ = IDNtuple::NEG_FLOAT;
  float seed_ptbiased_ = IDNtuple::NEG_FLOAT;

  bool seed_trk_driven_ = false;
  bool seed_ecal_driven_ = false;

  // GSF tracks: kine
  float gsf_pt_ = IDNtuple::NEG_FLOAT;
  float gsf_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_p_ = IDNtuple::NEG_FLOAT;
  int gsf_charge_ = 0; //@@ IDNtuple::NEG_INT;
  float gsf_inp_ = IDNtuple::NEG_FLOAT;
  float gsf_outp_ = IDNtuple::NEG_FLOAT;
  float gsf_dpt_ = IDNtuple::NEG_FLOAT;

  // GSF tracks: kine (mode)
  float gsf_mode_pt_ = IDNtuple::NEG_FLOAT;
  float gsf_mode_eta_ = IDNtuple::NEG_FLOAT;
  float gsf_mode_phi_ = IDNtuple::NEG_FLOAT;
  float gsf_mode_p_ = IDNtuple::NEG_FLOAT;

  // GSF tracks: quality
  int gsf_nhits_ = IDNtuple::NEG_INT;
  int gsf_missing_inner_hits_ = IDNtuple::NEG_INT;
  float gsf_chi2red_ = IDNtuple::NEG_FLOAT;

  // GSF tracks: displacement
  float gsf_dxy_ = IDNtuple::NEG_FLOAT;
  float gsf_dxy_err_ = IDNtuple::NEG_FLOAT;
  float gsf_dz_ = IDNtuple::NEG_FLOAT;
  float gsf_dz_err_ = IDNtuple::NEG_FLOAT;

  // GSF tracks: tangents
  int gsf_ntangents_ = 0; //@@ IDNtuple::NEG_INT;
  float gsf_hit_dpt_[NHITS_MAX] = {0}; //@@ {IDNtuple::NEG_FLOAT};
  float gsf_hit_dpt_unc_[NHITS_MAX] = {0}; //@@ {IDNtuple::NEG_FLOAT};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // PF GSF tracks: kine
  float pfgsf_pt_ = IDNtuple::NEG_FLOAT;
  float pfgsf_eta_ = IDNtuple::NEG_FLOAT;
  float pfgsf_phi_ = IDNtuple::NEG_FLOAT;
  float pfgsf_p_ = IDNtuple::NEG_FLOAT;
  int pfgsf_charge_ = 0; //@@ IDNtuple::NEG_INT;
  float pfgsf_inp_ = IDNtuple::NEG_FLOAT;
  float pfgsf_outp_ = IDNtuple::NEG_FLOAT;
  float pfgsf_dpt_ = IDNtuple::NEG_FLOAT;

  // PF GSF tracks: kine (mode)
  float pfgsf_mode_pt_ = IDNtuple::NEG_FLOAT;
  float pfgsf_mode_eta_ = IDNtuple::NEG_FLOAT;
  float pfgsf_mode_phi_ = IDNtuple::NEG_FLOAT;
  float pfgsf_mode_p_ = IDNtuple::NEG_FLOAT;

  // PF GSF tracks: quality
  int pfgsf_nhits_ = IDNtuple::NEG_INT;
  int pfgsf_missing_inner_hits_ = IDNtuple::NEG_INT;
  float pfgsf_chi2red_ = IDNtuple::NEG_FLOAT;

  // PF GSF tracks: displacement
  float pfgsf_dxy_ = IDNtuple::NEG_FLOAT;
  float pfgsf_dxy_err_ = IDNtuple::NEG_FLOAT;
  float pfgsf_dz_ = IDNtuple::NEG_FLOAT;
  float pfgsf_dz_err_ = IDNtuple::NEG_FLOAT;

  // PF GSF tracks: tangents
  int pfgsf_ntangents_ = 0; //@@ IDNtuple::NEG_INT;
  float pfgsf_hit_dpt_[NHITS_MAX] = {0}; //@@ {IDNtuple::NEG_FLOAT};
  float pfgsf_hit_dpt_unc_[NHITS_MAX] = {0}; //@@ {IDNtuple::NEG_FLOAT};
  //std::vector<float> pfgsf_extapolated_eta_;
  //std::vector<float> pfgsf_extapolated_phi_;

  // GSF electrons: kinematics
  float ele_pt_ = IDNtuple::NEG_FLOAT;
  float ele_eta_ = IDNtuple::NEG_FLOAT;
  float ele_phi_ = IDNtuple::NEG_FLOAT;
  float ele_p_ = IDNtuple::NEG_FLOAT;

  // Electrons: IDs
  float ele_mva_value_ = -999.; //@ IDNtuple::NEG_FLOAT;
  float ele_mva_value_retrained_ = -999.;
  float ele_conv_vtx_fit_prob_ = IDNtuple::NEG_FLOAT;
  float ele_mva_value_depth10_ = -999.; //@ IDNtuple::NEG_FLOAT;
  float ele_mva_value_depth15_ = -999.; //@ IDNtuple::NEG_FLOAT;

  // Electrons: MVA variables
  float eid_rho_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_ele_pt_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_trk_p_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_trk_nhits_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_trk_chi2red_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_gsf_nhits_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_gsf_chi2red_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_sc_E_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_sc_eta_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_sc_etaWidth_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_sc_phiWidth_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_match_seed_dEta_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_match_eclu_EoverP_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_match_SC_EoverP_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_match_SC_dEta_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_match_SC_dPhi_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_shape_full5x5_sigmaIetaIeta_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_shape_full5x5_sigmaIphiIphi_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_shape_full5x5_HoverE_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_shape_full5x5_r9_ = -666; //@@ IDNtuple::NEG_FLOAT;
  float eid_shape_full5x5_circularity_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float eid_brem_frac_ = -666; //@@ IDNtuple::NEG_FLOAT;

  float image_gsf_ref_eta_ = IDNtuple::NEG_FLOAT;
  float image_gsf_ref_phi_ = IDNtuple::NEG_FLOAT;
  float image_gsf_ref_R_ = IDNtuple::NEG_FLOAT;
  float image_gsf_ref_p_ = IDNtuple::NEG_FLOAT;
  float image_gsf_ref_pt_ = IDNtuple::NEG_FLOAT;

  float image_gen_inner_eta_ = IDNtuple::NEG_FLOAT;
  float image_gen_inner_phi_ = IDNtuple::NEG_FLOAT;
  float image_gen_inner_R_ = IDNtuple::NEG_FLOAT;
  float image_gen_inner_p_ = IDNtuple::NEG_FLOAT;
  float image_gen_inner_pt_ = IDNtuple::NEG_FLOAT;

  float image_gsf_inner_eta_ = IDNtuple::NEG_FLOAT;
  float image_gsf_inner_phi_ = IDNtuple::NEG_FLOAT;
  float image_gsf_inner_R_ = IDNtuple::NEG_FLOAT;
  float image_gsf_inner_p_ = IDNtuple::NEG_FLOAT;
  float image_gsf_inner_pt_ = IDNtuple::NEG_FLOAT;
  int image_gsf_charge_ = IDNtuple::NEG_INT*10;

  float image_gsf_proj_eta_ = IDNtuple::NEG_FLOAT;
  float image_gsf_proj_phi_ = IDNtuple::NEG_FLOAT;
  float image_gsf_proj_R_ = IDNtuple::NEG_FLOAT;
  float image_gsf_proj_p_ = IDNtuple::NEG_FLOAT;

  float image_gsf_atcalo_eta_ = IDNtuple::NEG_FLOAT;
  float image_gsf_atcalo_phi_ = IDNtuple::NEG_FLOAT;
  float image_gsf_atcalo_R_ = IDNtuple::NEG_FLOAT;
  float image_gsf_atcalo_p_ = IDNtuple::NEG_FLOAT;

  unsigned int image_clu_n_ = 0;
  //std::vector<float> image_clu_eta_ = {};
  //std::vector<float> image_clu_phi_ = {};
  //std::vector<float> image_clu_e_ = {};
  float image_clu_eta_[ARRAY_SIZE] = {};
  float image_clu_phi_[ARRAY_SIZE] = {};
  float image_clu_e_[ARRAY_SIZE] = {};
  int image_clu_nhit_[ARRAY_SIZE] = {};

  unsigned int image_pf_n_ = 0;
  //std::vector<float> image_pf_eta_ = {};
  //std::vector<float> image_pf_phi_ = {};
  //std::vector<float> image_pf_p_ = {};
  //std::vector<int> image_pf_pdgid_ = {};
  float image_pf_eta_[ARRAY_SIZE] = {};
  float image_pf_phi_[ARRAY_SIZE] = {};
  float image_pf_p_[ARRAY_SIZE] = {};
  int image_pf_pdgid_[ARRAY_SIZE] = {};
  int image_pf_matched_[ARRAY_SIZE] = {};
  int image_pf_lost_[ARRAY_SIZE] = {};

};

#endif // LowPtElectrons_LowPtElectrons_IDNtuple
