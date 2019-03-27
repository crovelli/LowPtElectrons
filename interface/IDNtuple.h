#ifndef LowPtElectrons_LowPtElectrons_IDNtuple
#define LowPtElectrons_LowPtElectrons_IDNtuple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include <vector>
class TTree;

constexpr size_t NHITS_MAX = 30;

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

  void fill_evt( const edm::EventID& id );
  void fill_gen( const pat::PackedGenParticleRef );
  void fill_gen( const reco::GenParticleRef ); //@@ AOD

  void fill_gsf_trk( const reco::GsfTrackRef trk,
		     const reco::BeamSpot& spot );

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
  float rho_ = -1.;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;
  bool is_other_ = false;
  
  // GEN electrons
  float gen_pt_ = -1.;
  float gen_eta_ = -1.;
  float gen_phi_ = -1.;
  float gen_e_ = -1.;
  float gen_p_ = -1.;
  int gen_charge_ = -1;
  int gen_pdgid_ = 0;
  int gen_mom_pdgid_ = 0;
  int gen_gran_pdgid_ = 0;

  // GSF tracks: kine
  float gsf_pt_ = -1.;
  float gsf_eta_ = -1.;
  float gsf_phi_ = -1.;
  float gsf_p_ = -1.;
  int gsf_charge_ = 0;
  float gsf_inp_ = -1.;
  float gsf_outp_ = -1.;
  float gsf_dpt_ = -1.;

  // GSF tracks: quality
  int gsf_nhits_ = -1;
  int gsf_missing_inner_hits_ = -1;
  float gsf_chi2red_ = -1.;

  // GSF tracks: displacement
  float gsf_dxy_ = -1;
  float gsf_dxy_err_ = -1;
  float gsf_dz_ = -1;
  float gsf_dz_err_ = -1;

  // GSF tracks: tangents
  int gsf_ntangents_ = 0;
  float gsf_hit_dpt_[NHITS_MAX] = {0.};
  float gsf_hit_dpt_unc_[NHITS_MAX] = {0.};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // GSF electrons: kinematics
  float ele_pt_ = -1.;
  float ele_eta_ = -1.;
  float ele_phi_ = -1.;
  float ele_p_ = -1.;

  // Electrons: IDs
  float ele_mvaIdV2_ = -2.;
  float ele_lowPtMva_ = -999.;
  float ele_conv_vtx_fit_prob_ = -1.;

  // Electrons: MVA variables
  float eid_rho_ = -666.;
  float eid_ele_pt_ = -666.;

  float eid_trk_p_ = -666.;
  float eid_trk_nhits_ = -666.;
  float eid_trk_chi2red_ = -666.;

  float eid_gsf_nhits_ = -666.;
  float eid_gsf_chi2red_ = -666.;

  float eid_sc_E_ = -666.;
  float eid_sc_eta_ = -666.;
  float eid_sc_etaWidth_ = -666.;
  float eid_sc_phiWidth_ = -666.;

  float eid_match_seed_dEta_ = -666.;
  float eid_match_eclu_EoverP_ = -666.;
  float eid_match_SC_EoverP_ = -666.;
  float eid_match_SC_dEta_ = -666.;
  float eid_match_SC_dPhi_ = -666.;

  float eid_shape_full5x5_sigmaIetaIeta_ = -666.;
  float eid_shape_full5x5_sigmaIphiIphi_ = -666.;
  float eid_shape_full5x5_HoverE_ = -666.;
  float eid_shape_full5x5_r9_ = -666.;
  float eid_shape_full5x5_circularity_ = -666.;

  float eid_brem_frac_ = -666.;
  
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace eleid {
  
  class Features {
  public:
    // KF track
    float trk_p_ = -1.;
    float trk_nhits_ = -1.;
    float trk_chi2red_ = -1.;
    // GSF track
    float gsf_nhits_ = -1.;
    float gsf_chi2red_ = -1.;
    // SC 
    float sc_E_ = -1.;
    float sc_eta_ = -1.;
    float sc_etaWidth_ = -1.;
    float sc_phiWidth_ = -1.;
    // Track-cluster matching
    float match_seed_dEta_ = -1.;
    float match_eclu_EoverP_ = -1.;
    float match_SC_EoverP_ = -1.;
    float match_SC_dEta_ = -1.;
    float match_SC_dPhi_ = -1.;
    // Shower shape vars
    float shape_full5x5_sigmaIetaIeta_ = -1.;
    float shape_full5x5_sigmaIphiIphi_ = -1.;
    float shape_full5x5_HoverE_ = -1.;
    float shape_full5x5_r9_ = -1.;
    float shape_full5x5_circularity_ = -1.;
    // Misc
    float rho_ = -1.;
    float brem_frac_ = -1.;
    float ele_pt_ = -1.;
  public:
    std::vector<float> get();
    void set( const pat::ElectronRef& ele, double rho );
  };

}

#endif // LowPtElectrons_LowPtElectrons_IDNtuple
