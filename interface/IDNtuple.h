#ifndef LowPtElectrons_LowPtElectrons_IDNtuple
#define LowPtElectrons_LowPtElectrons_IDNtuple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

constexpr size_t NHITS_MAX = 30;
constexpr float NEG_INT = -10;
constexpr float NEG_FLOAT = -10.;

// Small class to provide fillers and hide tree I/O
class IDNtuple {

 public:

  IDNtuple() {}

  void reset() {
    IDNtuple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }
  
  void link_tree( TTree* tree );
  
  void set_rho( float r ) { rho_ = r; }

  void is_e( bool t = true ) { is_e_ = t; }
  void is_e_not_matched( bool t = true ) { is_e_not_matched_ = t; }
  void is_other( bool t = true ) { is_other_ = t; }
  
  void has_egamma_gsf( bool f = false ) { has_egamma_gsf_ = f; }
  void has_egamma_ele( bool f = false ) { has_egamma_ele_ = f; }
	void is_egamma( bool f = false ) { is_egamma_ = f; }

  void fill_evt( const edm::EventID& id );

  void fill_gen( const pat::PackedGenParticleRef );
  void fill_gen( const reco::GenParticleRef ); //@@ AOD
  void fill_gen( const reco::CandidatePtr );

  void fill_seed( double seed_bdt_unbiased,
		  double seed_bdt_ptbiased );

  void fill_gsf( const reco::GsfTrackRef trk,
		 const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float id_lowpt,
		 float id_v2,
		 float ele_conv_vtx_fit_prob,
		 const double rho );

  void fill_ele( const pat::ElectronRef ele,
		 float id_lowpt,
		 float id_v2,
		 float ele_conv_vtx_fit_prob,
		 const double rho );
  
 private:

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float rho_ = NEG_FLOAT;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;
  bool is_other_ = false;

  // CMS PF reco
  bool has_egamma_gsf_ = false;
  bool has_egamma_ele_ = false;
  bool is_egamma_ = false;
  
  // GEN electrons
  float gen_pt_ = NEG_FLOAT;
  float gen_eta_ = NEG_FLOAT;
  float gen_phi_ = NEG_FLOAT;
  float gen_e_ = NEG_FLOAT;
  float gen_p_ = NEG_FLOAT;
  int gen_charge_ = NEG_INT;
  int gen_pdgid_ = 0;
  int gen_mom_pdgid_ = 0;
  int gen_gran_pdgid_ = 0;

  // Seed BDT discriminator values
  float seed_bdt_unbiased_ = NEG_FLOAT;
  float seed_bdt_ptbiased_ = NEG_FLOAT;

  // GSF tracks: kine
  float gsf_pt_ = NEG_FLOAT;
  float gsf_eta_ = NEG_FLOAT;
  float gsf_phi_ = NEG_FLOAT;
  float gsf_p_ = NEG_FLOAT;
  int gsf_charge_ = 0;
  float gsf_inp_ = NEG_FLOAT;
  float gsf_outp_ = NEG_FLOAT;
  float gsf_dpt_ = NEG_FLOAT;

  // GSF tracks: kine (mode)
  float mode_pt_ = NEG_FLOAT;
  float mode_eta_ = NEG_FLOAT;
  float mode_phi_ = NEG_FLOAT;
  float mode_p_ = NEG_FLOAT;

  // GSF tracks: quality
  int gsf_nhits_ = NEG_INT;
  int gsf_missing_inner_hits_ = NEG_INT;
  float gsf_chi2red_ = NEG_FLOAT;

  // GSF tracks: displacement
  float gsf_dxy_ = NEG_FLOAT;
  float gsf_dxy_err_ = NEG_FLOAT;
  float gsf_dz_ = NEG_FLOAT;
  float gsf_dz_err_ = NEG_FLOAT;

  // GSF tracks: tangents
  int gsf_ntangents_ = 0;
  float gsf_hit_dpt_[NHITS_MAX] = {NEG_FLOAT};
  float gsf_hit_dpt_unc_[NHITS_MAX] = {NEG_FLOAT};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // GSF electrons: kinematics
  float ele_pt_ = NEG_FLOAT;
  float ele_eta_ = NEG_FLOAT;
  float ele_phi_ = NEG_FLOAT;
  float ele_p_ = NEG_FLOAT;

  // Electrons: IDs
  float ele_mvaIdV2_ = NEG_FLOAT;
  float ele_lowPtMva_ = NEG_FLOAT;
  float ele_conv_vtx_fit_prob_ = NEG_FLOAT;

  // Electrons: MVA variables
  float eid_rho_ = NEG_FLOAT;
  float eid_ele_pt_ = NEG_FLOAT;

  float eid_trk_p_ = NEG_FLOAT;
  float eid_trk_nhits_ = NEG_FLOAT;
  float eid_trk_chi2red_ = NEG_FLOAT;

  float eid_gsf_nhits_ = NEG_FLOAT;
  float eid_gsf_chi2red_ = NEG_FLOAT;

  float eid_sc_E_ = NEG_FLOAT;
  float eid_sc_eta_ = NEG_FLOAT;
  float eid_sc_etaWidth_ = NEG_FLOAT;
  float eid_sc_phiWidth_ = NEG_FLOAT;

  float eid_match_seed_dEta_ = NEG_FLOAT;
  float eid_match_eclu_EoverP_ = NEG_FLOAT;
  float eid_match_SC_EoverP_ = NEG_FLOAT;
  float eid_match_SC_dEta_ = NEG_FLOAT;
  float eid_match_SC_dPhi_ = NEG_FLOAT;

  float eid_shape_full5x5_sigmaIetaIeta_ = NEG_FLOAT;
  float eid_shape_full5x5_sigmaIphiIphi_ = NEG_FLOAT;
  float eid_shape_full5x5_HoverE_ = NEG_FLOAT;
  float eid_shape_full5x5_r9_ = NEG_FLOAT;
  float eid_shape_full5x5_circularity_ = NEG_FLOAT;

  float eid_brem_frac_ = NEG_FLOAT;
  
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace eleid {
  
  class Features {
  public:
    // KF track
    float trk_p_ = NEG_FLOAT;
    float trk_nhits_ = NEG_FLOAT;
    float trk_chi2red_ = NEG_FLOAT;
    // GSF track
    float gsf_nhits_ = NEG_FLOAT;
    float gsf_chi2red_ = NEG_FLOAT;
    // SC 
    float sc_E_ = NEG_FLOAT;
    float sc_eta_ = NEG_FLOAT;
    float sc_etaWidth_ = NEG_FLOAT;
    float sc_phiWidth_ = NEG_FLOAT;
    // Track-cluster matching
    float match_seed_dEta_ = NEG_FLOAT;
    float match_eclu_EoverP_ = NEG_FLOAT;
    float match_SC_EoverP_ = NEG_FLOAT;
    float match_SC_dEta_ = NEG_FLOAT;
    float match_SC_dPhi_ = NEG_FLOAT;
    // Shower shape vars
    float shape_full5x5_sigmaIetaIeta_ = NEG_FLOAT;
    float shape_full5x5_sigmaIphiIphi_ = NEG_FLOAT;
    float shape_full5x5_HoverE_ = NEG_FLOAT;
    float shape_full5x5_r9_ = NEG_FLOAT;
    float shape_full5x5_circularity_ = NEG_FLOAT;
    // Misc
    float rho_ = NEG_FLOAT;
    float brem_frac_ = NEG_FLOAT;
    float ele_pt_ = NEG_FLOAT;
  public:
    std::vector<float> get();
    void set( const pat::ElectronRef& ele, double rho );
    void set( const reco::GsfElectronPtr& ele, double rho );
  };

}

#endif // LowPtElectrons_LowPtElectrons_IDNtuple
