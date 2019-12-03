#ifndef LowPtElectrons_LowPtElectrons_IDSlimNtuple
#define LowPtElectrons_LowPtElectrons_IDSlimNtuple

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
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
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

  void set_rho( float r ) { rho_ = r; }

  void fill_evt( const edm::EventID& id );
  void fill_gen( const pat::PackedGenParticleRef );   
  void fill_gen( const reco::GenParticlePtr );
  void fill_gen_default( );

  void fill_trk( const reco::TrackPtr& trk,
		 const reco::BeamSpot& spot );

  void fill_trk_dEdx( const reco::TrackPtr& trk,
                      std::vector< const edm::ValueMap< reco::DeDxData > * > & x );
  void fill_trk_dEdx_default();

  void fill_bdt( double seed_unbiased,
		 double seed_ptbiased );

  void fill_gsf( const reco::GsfTrackPtr trk,
		 const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value,
		 int mva_id,
		 float ele_conv_vtx_fit_prob,
		 const double rho );

  void fill_supercluster(const reco::GsfElectronPtr ele, noZS::EcalClusterLazyTools *ecalTools_);
  void fill_supercluster_miniAOD(const reco::GsfElectronPtr ele); 

  template < typename T> bool validPtr( edm::Ptr<T>& ptr);

 public:

  // To reduce ntuples size
  bool largeNtuple = false;

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float weight_ = 1.;   
  float rho_ = IDSlimNtuple::NEG_FLOAT;

  // Data sample
  int is_aod_ = -1;
  int is_mc_ = -1;

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
  float gen_pt_  = IDSlimNtuple::NEG_FLOAT;
  float gen_eta_ = IDSlimNtuple::NEG_FLOAT;
  float gen_phi_ = IDSlimNtuple::NEG_FLOAT;
  float gen_p_   = IDSlimNtuple::NEG_FLOAT;
  int gen_charge_ = IDSlimNtuple::NEG_INT;
  int gen_pdgid_  = IDSlimNtuple::NEG_INT;
  int gen_mom_pdgid_ = IDSlimNtuple::NEG_INT;
  int gen_gran_pdgid_ = IDSlimNtuple::NEG_INT;
  int gen_tag_side_ = IDSlimNtuple::NEG_INT;  
  float gen_trk_dr_  = IDSlimNtuple::NEG_FLOAT;
  float gen_gsf_dr_  = IDSlimNtuple::NEG_FLOAT;

  // RECO steps
  float gsf_dr_ = IDSlimNtuple::NEG_FLOAT;
  float trk_dr_ = IDSlimNtuple::NEG_FLOAT;

  // GSF tracks: kine
  float gsf_pt_   = IDSlimNtuple::NEG_FLOAT;
  float gsf_eta_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_phi_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_p_    = IDSlimNtuple::NEG_FLOAT;
  int gsf_charge_ = IDSlimNtuple::NEG_INT;
  float gsf_inp_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_outp_ = IDSlimNtuple::NEG_FLOAT;

  // GSF tracks: kine (mode)
  float gsf_mode_pt_  = IDSlimNtuple::NEG_FLOAT;
  float gsf_mode_eta_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_mode_phi_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_mode_p_   = IDSlimNtuple::NEG_FLOAT;

  // GSF tracks: quality
  int gsf_missing_inner_hits_ = IDSlimNtuple::NEG_INT;

  // GSF tracks: displacement
  float gsf_dxy_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_dxy_err_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_dz_ = IDSlimNtuple::NEG_FLOAT;
  float gsf_dz_err_ = IDSlimNtuple::NEG_FLOAT;

  // GSF pos at ECAL
  float gsf_x_ = IDSlimNtuple::NEG_FLOAT;     
  float gsf_y_ = IDSlimNtuple::NEG_FLOAT;      
  float gsf_z_ = IDSlimNtuple::NEG_FLOAT;     

  // GSF tracks: tangents
  int gsf_ntangents_ = 0; //@@ IDSlimNtuple::NEG_INT;
  //float gsf_hit_dpt_[NHITS_MAX] = {0}; //@@ {IDSlimNtuple::NEG_FLOAT};
  //float gsf_hit_dpt_unc_[NHITS_MAX] = {0}; //@@ {IDSlimNtuple::NEG_FLOAT};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // Seed BDT discriminator values at GsfTrack level
  float seed_unbiased_ = IDSlimNtuple::NEG_FLOAT;
  float seed_ptbiased_ = IDSlimNtuple::NEG_FLOAT;

  // KF tracks: kine
  float trk_pt_  = IDSlimNtuple::NEG_FLOAT;
  float trk_eta_  = IDSlimNtuple::NEG_FLOAT;
  float trk_phi_  = IDSlimNtuple::NEG_FLOAT;
  float trk_p_    = IDSlimNtuple::NEG_INT;
  int trk_charge_ = IDSlimNtuple::NEG_INT;
  float trk_inp_  = IDSlimNtuple::NEG_FLOAT;
  float trk_outp_ = IDSlimNtuple::NEG_FLOAT;
  int pdg_id_     = IDSlimNtuple::NEG_INT;

  // KF tracks: quality
  int trk_nhits_ = IDSlimNtuple::NEG_INT;
  int trk_missing_inner_hits_ = IDSlimNtuple::NEG_INT;
  float trk_chi2red_   = IDSlimNtuple::NEG_FLOAT;
  int trk_high_purity_ = IDSlimNtuple::NEG_INT;

  // KF tracks: displ
  float trk_dxy_     = IDSlimNtuple::NEG_FLOAT;
  float trk_dxy_err_ = IDSlimNtuple::NEG_FLOAT;
  float trk_dz_      = IDSlimNtuple::NEG_FLOAT;
  float trk_dz_err_  = IDSlimNtuple::NEG_FLOAT;

  // KF tracks: dE/dx
  float trk_dEdx1_   = IDSlimNtuple::NEG_FLOAT;
  int trk_dEdx1_Nm_  = IDSlimNtuple::NEG_INT;
  int trk_dEdx1_NSm_ = IDSlimNtuple::NEG_INT;

  // GSF electrons: kinematics
  float ele_p_   = IDSlimNtuple::NEG_FLOAT;
  float ele_pt_  = IDSlimNtuple::NEG_FLOAT;
  float ele_eta_ = IDSlimNtuple::NEG_FLOAT;
  float ele_phi_ = IDSlimNtuple::NEG_FLOAT;
  int p4kind_    = IDSlimNtuple::NEG_INT;
  float core_shFracHits_ = IDSlimNtuple::NEG_FLOAT;
  float ele_p_atvtx_  = IDSlimNtuple::NEG_FLOAT;
  float ele_p_atcalo_ = IDSlimNtuple::NEG_FLOAT;
  int fiducial_isEB_      = IDSlimNtuple::NEG_INT;
  int fiducial_isEE_      = IDSlimNtuple::NEG_INT; 
  int fiducial_isEBEEGap_ = IDSlimNtuple::NEG_INT;
    
  // Electron: Charge
  int chPix_ = IDSlimNtuple::NEG_INT;
  int chGCP_ = IDSlimNtuple::NEG_INT;
  int chGP_  = IDSlimNtuple::NEG_INT;
  int chGC_  = IDSlimNtuple::NEG_INT;  

  // Electrons: IDs
  float ele_mva_value_ = IDSlimNtuple::NEG_FLOAT;
  int ele_mva_id_      = IDSlimNtuple::NEG_INT;
  float ele_conv_vtx_fit_prob_ = IDSlimNtuple::NEG_FLOAT;

  float eid_rho_    = IDSlimNtuple::NEG_FLOAT;
  float eid_ele_pt_ = IDSlimNtuple::NEG_FLOAT;
  float eid_sc_eta_ = IDSlimNtuple::NEG_FLOAT;
  float eid_shape_full5x5_sigmaIetaIeta_ = IDSlimNtuple::NEG_FLOAT;
  float eid_shape_full5x5_sigmaIphiIphi_ = IDSlimNtuple::NEG_FLOAT;
  float eid_shape_full5x5_circularity_   = IDSlimNtuple::NEG_FLOAT;
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

  // Electron: firther track-Cluster matching
  float match_seed_EoverP_    = IDSlimNtuple::NEG_FLOAT;
  float match_seed_EoverPout_ = IDSlimNtuple::NEG_FLOAT;
  float match_seed_dPhi_      = IDSlimNtuple::NEG_FLOAT;
  float match_seed_dEta_vtx_  = IDSlimNtuple::NEG_FLOAT;
  float match_eclu_EoverPout_ = IDSlimNtuple::NEG_FLOAT;
  float match_eclu_dEta_      = IDSlimNtuple::NEG_FLOAT;
  float match_eclu_dPhi_      = IDSlimNtuple::NEG_FLOAT;

  // Further full 5x5 shower shape 
  float shape_full5x5_e1x5_    = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_e2x5Max_ = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_e5x5_    = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_HoverEBc_     = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth1_   = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth2_   = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth1Bc_ = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth2Bc_ = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_eLeft_   = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_eRight_  = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_eTop_    = IDSlimNtuple::NEG_FLOAT;
  float shape_full5x5_eBottom_ = IDSlimNtuple::NEG_FLOAT;
  
  // Electron, brem fractions
  float brem_fracTrk_ = IDSlimNtuple::NEG_FLOAT;
  float brem_fracSC_  = IDSlimNtuple::NEG_FLOAT;
  int brem_N_         = IDSlimNtuple::NEG_INT;

  // Electron, isolation
  float ele_sumPhotonEt_;
  float ele_sumChargedHadronPt_;
  float ele_sumNeutralHadronEt_;

  // SuperClusters 
  float sc_Et_  = IDSlimNtuple::NEG_FLOAT;
  int sc_Nclus_ = IDSlimNtuple::NEG_INT;
  int sc_Nclus_deta01_ = IDSlimNtuple::NEG_INT;
  int sc_Nclus_deta02_ = IDSlimNtuple::NEG_INT;
  int sc_Nclus_deta03_ = IDSlimNtuple::NEG_INT;  
  bool sc_goodSeed_ = false;
  float sc_E_ps_ = IDSlimNtuple::NEG_FLOAT;
  float sc_E_ps1_= IDSlimNtuple::NEG_FLOAT;
  float sc_E_ps2_= IDSlimNtuple::NEG_FLOAT;

  // Clusters 
  float sc_cluster_et_[NCLUS_MAX]    = {IDSlimNtuple::NEG_FLOAT}; 
  float sc_cluster_E_[NCLUS_MAX]     = {IDSlimNtuple::NEG_FLOAT}; 
  float sc_cluster_eta_[NCLUS_MAX]   = {IDSlimNtuple::NEG_FLOAT}; 
  float sc_cluster_phi_[NCLUS_MAX]   = {IDSlimNtuple::NEG_FLOAT}; 
  int   sc_cluster_nxtal_[NCLUS_MAX] = {IDSlimNtuple::NEG_INT}; 
  float sc_cluster_e1x3_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e1x5_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2x2_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e3x3_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e5x5_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eMax_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2nd_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Right_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Left_[NCLUS_MAX]   = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Top_[NCLUS_MAX]    = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Bottom_[NCLUS_MAX] = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eRight_[NCLUS_MAX]  = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eLeft_[NCLUS_MAX]   = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eTop_[NCLUS_MAX]    = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eBottom_[NCLUS_MAX] = {IDSlimNtuple::NEG_FLOAT};  
  float sc_cluster_eMaxOver2x2_[NCLUS_MAX] = {IDSlimNtuple::NEG_FLOAT};
  float sc_cluster_eMaxOver3x3_[NCLUS_MAX] = {IDSlimNtuple::NEG_FLOAT};
  float sc_cluster_eMaxOver1x3_[NCLUS_MAX] = {IDSlimNtuple::NEG_FLOAT};

  // When running on miniaod, only 3 highest ET clusters only
  float sc_clus1_et_     = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_E_      = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_E_ov_p_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus1_nxtal_    = IDSlimNtuple::NEG_INT;
  float sc_clus1_dphi_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_deta_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus1_ntrk_deta01_ = IDSlimNtuple::NEG_FLOAT;
  //
  float sc_clus2_et_     = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_E_      = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_E_ov_p_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus2_nxtal_    = IDSlimNtuple::NEG_INT;
  float sc_clus2_dphi_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_deta_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus2_ntrk_deta01_ =IDSlimNtuple::NEG_FLOAT;
  //
  float sc_clus3_et_     = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_E_      = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_E_ov_p_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_E_ov_E_ = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_eta_    = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_phi_    = IDSlimNtuple::NEG_FLOAT;
  int sc_clus3_nxtal_    = IDSlimNtuple::NEG_INT;
  float sc_clus3_dphi_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_deta_   = IDSlimNtuple::NEG_FLOAT;
  float sc_clus3_ntrk_deta01_ = IDSlimNtuple::NEG_FLOAT;
};

#endif // LowPtElectrons_LowPtElectrons_IDSlimNtuple
