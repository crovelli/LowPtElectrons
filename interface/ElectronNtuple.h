#ifndef LowPtElectrons_LowPtElectrons_ElectronNtuple
#define LowPtElectrons_LowPtElectrons_ElectronNtuple

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include <vector>
class TTree;
/*namespace reco {
	class GsfTrackRef;
	class GsfElectronRef;
	class GenParticleRef;
	class PreId;
	class TrackRef;
	class BeamSpot;
}

namespace edm {
	class EventID;
	}*/

//constants
// #define NHITS_MAX 30
// #define BREM_WINDOW_ETA 7
// #define BREM_WINDOW_PHI 12
// #define CLUSTER_WINDOW_ETA 10
// #define CLUSTER_WINDOW_PHI 15

constexpr size_t NHITS_MAX = 30;
constexpr size_t BREM_WINDOW_ETA = 7;
constexpr size_t BREM_WINDOW_PHI = 12;
constexpr size_t CLUSTER_WINDOW_ETA = 10;
constexpr size_t CLUSTER_WINDOW_PHI = 15;

constexpr size_t NECAL_PFCLUSTERS = 30;
constexpr size_t NHCAL_PFCLUSTERS = 30;
constexpr size_t ECAL_CLUSTER_SIZE = 7;


class ElectronNtuple {
	/*Small class to provide fillers and hide tree I/O*/
public:
	ElectronNtuple() {}
	
	void reset() {
		ElectronNtuple dummy; //create a new one
		*this = dummy; //use assignment to reset
	}

	void link_tree(TTree *tree);

	//fillers
	void fill_evt(const edm::EventID &id);
	void fill_gen(const reco::GenParticleRef genp);
	void fill_gsf_trk(const reco::GsfTrackRef trk, const reco::BeamSpot &spot);
	void fill_pfgsf_trk(const reco::GsfPFRecTrackRef pfgsf);
	void fill_preid(const reco::PreId &preid_ecal, const reco::PreId &preid_hcal, 
									const reco::BeamSpot &spot, const double rho, const int num_gsf,
									noZS::EcalClusterLazyTools& ecalTools);
	void fill_ele(const reco::GsfElectronRef ele, float mvaid_v1=-2, float mvaid_v2=-2, 
								float ele_conv_vtx_fit_prob = -1., const std::vector<float>& iso_rings={0., 0., 0., 0.});
	void fill_supercluster(const reco::GsfElectronRef ele);
	void fill_ktf_trk(const reco::TrackRef trk, const reco::BeamSpot &spot);
	void fill_GSF_ECAL_cluster_info(
		const reco::PFClusterRef cluster,
		const reco::PFTrajectoryPoint &gsf,
		noZS::EcalClusterLazyTools& tools
		);
	void fill_GSF_HCAL_cluster_info(
		const reco::PFClusterRef cluster,
		const reco::PFTrajectoryPoint &gsf
		//noZS::EcalClusterLazyTools& tools //Something similar for HCAL?
		);

	void print_val(int i=0) {
		std::cout<< "  Fill ntuple: " << i << " " << preid_ktf_ecal_cluster_e_ << std::endl;
	}

	//TODO: refactor to avoid repetitions
	//TODO: refactor to avoid repetitions
	void fill_KTF_ECAL_cluster_info(
		const reco::PFClusterRef cluster,
		const reco::PFTrajectoryPoint &ktf,
		noZS::EcalClusterLazyTools& tools
		);
	void fill_KTF_HCAL_cluster_info(
		const reco::PFClusterRef cluster,
		const reco::PFTrajectoryPoint &ktf
		//noZS::EcalClusterLazyTools& tools //Something similar for HCAL?
		);

	void is_ECAL_cluster_same(bool t=true){gsf_ktf_same_ecal_ = t;}
	void is_HCAL_cluster_same(bool t=true){gsf_ktf_same_hcal_ = t;}

	void is_e(bool t=true) {is_e_=t;}
	void is_e_not_matched(bool t=true) {is_e_not_matched_=t;}
	void is_other(bool t=true) {is_other_=t;}
	void has_ele_core(bool t=true) {has_ele_core_=t;}
	void has_pfEgamma(bool t=true) {has_pfEgamma_=t;}
	void has_pfBlock_with_SC(bool t=true) {has_pfBlock_with_SC_ = t;}
	void has_pfBlock_with_ECAL(bool t=true) {has_pfBlock_with_ECAL_ = t;}
	void has_pfBlock(bool t=true) {has_pfBlock_ = t;}
	void has_pfBlock_size(float f) {has_pfBlock_size_ = f;}
	void has_pfBlock_dr(float f) {has_pfBlock_dr_ = f;}
	void has_pfGSFTrk(bool t=true) {has_pfGSF_trk_=t;}
	void unpack_pfgsf_flags(int flags);
	void set_rho(float r) {rho_=r;}

private:

	bool get_ith_bit(const int val, const size_t ibit) const {return (val >> ibit) & 1;}

	//only simple types (no arrays allowed, otherwise the reset() method fails;
	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	unsigned long long evt_ = 0;

	bool is_e_ = false;
	bool is_e_not_matched_ = false;
	bool is_other_ = false;

	float rho_ = -1;
	
	// GEN electrons
	float gen_pt_ = -1;
	float gen_eta_ = -1;
	float gen_phi_ = -1;
	float gen_e_ = -1;
	float gen_p_ = -1;
	int   gen_charge_ = 0;

	// KF tracks: kine
	float trk_pt_ = -1.;
	float trk_eta_ = -1.;
	float trk_phi_ = -1.;
	float trk_p_ = -1;
	int   trk_charge_ = 0;
	float trk_inp_ = -1.;
	float trk_outp_ = -1.;
	float trk_dpt_ = -1.;
	// KF tracks: quality
	int trk_nhits_ = -1;
	int trk_missing_inner_hits_ = -1;
	int   trk_high_purity_ = 0;
	float trk_chi2red_ = -1.;
	// KF tracks: displ
	float trk_dxy_ = -1;
	float trk_dxy_err_ = -1;
	float trk_dz_ = -1;
	float trk_dz_err_ = -1;

	// PreId: ECAL/track matching
	float preid_trk_pt_  = -1.;
	float preid_trk_eta_ = -1.;
	float preid_trk_phi_ = -1.;
  float preid_trk_p_   = -1.;
  float preid_trk_nhits_ = -1.;
  float preid_trk_high_quality_ = -1.;
  float preid_trk_chi2red_ = -1.;
  float preid_rho_ = -1.;
  float preid_ktf_ecal_cluster_e_ = -1.;
  float preid_ktf_ecal_cluster_deta_ = -1.;
  float preid_ktf_ecal_cluster_dphi_ = -1.;
  float preid_ktf_ecal_cluster_e3x3_ = -1.;
  float preid_ktf_ecal_cluster_e5x5_ = -1.;
  float preid_ktf_ecal_cluster_covEtaEta_ = -1.;
  float preid_ktf_ecal_cluster_covEtaPhi_ = -1.;
  float preid_ktf_ecal_cluster_covPhiPhi_ = -1.;
  float preid_ktf_ecal_cluster_r9_ = -1.;
  float preid_ktf_ecal_cluster_circularity_ = -1.;
  float preid_ktf_hcal_cluster_e_ = -1.;
  float preid_ktf_hcal_cluster_deta_ = -1.;
  float preid_ktf_hcal_cluster_dphi_ = -1.;
  float preid_gsf_dpt_ = -1.;
  float preid_trk_gsf_chiratio_ = -1.;
  float preid_gsf_chi2red_ = -1.;
  float preid_trk_dxy_sig_ = -1.; // must be last (not used by unbiased model)

	// PreId: MVA output
	float preid_bdtout1_ = -1.;
	float preid_bdtout2_ = -1.;
	bool preid_mva1_pass_ = false;
	bool preid_mva2_pass_ = false;
	
	// GSF tracks: kine
	float gsf_pt_ = -1.;
	float gsf_eta_ = -1.;
	float gsf_phi_ = -1.;
	float gsf_p_ = -1;
	int   gsf_charge_ = 0;
	float gsf_inp_ = -1.;
	float gsf_outp_ = -1.;
	float gsf_dpt_ = -1.;
	// GSF tracks: quality
	int gsf_nhits_ = -1;
	int gsf_missing_inner_hits_ = -1;
	float gsf_chi2red_ = -1.;
	// GSF tracks: displ
	float gsf_dxy_ = -1;
	float gsf_dxy_err_ = -1;
	float gsf_dz_ = -1;
	float gsf_dz_err_ = -1;
	int gsf_ntangents_ = 0;
 	float gsf_hit_dpt_[NHITS_MAX] = {0};
 	float gsf_hit_dpt_unc_[NHITS_MAX] = {0};
	std::vector<float> gsf_extapolated_eta_;
	std::vector<float> gsf_extapolated_phi_;


	//PFGSFTrack internal steps flags
	bool pfgsf_gsf_has_ktf_ = false;
	bool pfgsf_ktf_is_fifthStep_ = false;
	bool pfgsf_gsf_ecalDriven_ = false;
	bool pfgsf_gsf_trackerDriven_ = false;
	bool pfgsf_valid_gsf_brem_ = false;
	bool pfgsf_passes_preselection_ = false;
	bool pfgsf_passes_selection_ = false;

	bool pfgsf_xclean_seedref_ = false;
	bool pfgsf_xclean_ECALDriven_too_few_hits_ = false;
	bool pfgsf_xclean_ECALDriven_bad_EoverP_ = false;
	bool pfgsf_xclean_TrkDriven_too_few_hits_ = false;
	bool pfgsf_xclean_TrkDriven_vs_ECALDriven_ = false;
	bool pfgsf_xclean_BothTrk_bad_EoverP_ = false;
	bool pfgsf_xclean_BothTrk_noECAL_match_ = false;
	bool pfgsf_xclean_AngularGsfCleaning_ = false;
	bool pfgsf_xclean_noECAL_match_AGAIN_ = false;
	bool pfgsf_xclean_FINAL_ = false;

//	bool pfgsf_gsf_cntr_ = false;
//	bool pfgsf_findpfref_seednullptr_ = false;
//	bool pfgsf_findpfref_castnullptr_ = false;
//	bool pfgsf_findpfref_trknullptr_ = false;
//	bool pfgsf_findpfref_nosharedhits_ = false;
//	bool pfgsf_findpfref_nodrmatch_ = false;
//	bool pfgsf_findpfref_blah_ = false;
//	bool pfgsf_findpfref_trkmatch_ = false;
//	bool pfgsf_findpfref_end_ = false;

	//Middle steps
	bool has_ele_core_ = false;
	bool has_pfEgamma_ = false;
	bool has_pfBlock_with_SC_ = false;
	bool has_pfBlock_with_ECAL_ = false;
	bool has_pfBlock_ = false;
	float has_pfBlock_size_ = -1.;
	float has_pfBlock_dr_ = -1.;
	bool has_pfGSF_trk_ = false;

	// GSF electrons
	float ele_pt_ = -1.;
	float ele_eta_ = -1.;
	float ele_phi_ = -1.;
	float ele_p_ = -1.;
	float ele_mvaIdV1_ = -2.;
	float ele_mvaIdV2_ = -2.;
	float ele_conv_vtx_fit_prob_ = -1.;
	float ele_iso01_ = 0.;
	float ele_iso02_ = 0.;
	float ele_iso03_ = 0.;
	float ele_iso04_ = 0.;

	// Bottom up approach	
	float gsf_ecal_cluster_e_ = -1;
	float gsf_ecal_cluster_ecorr_ = -1;
	float gsf_ecal_cluster_eta_ = -1;
	float gsf_ecal_cluster_deta_ = -42;
	float gsf_ecal_cluster_dphi_ = -42;
	float gsf_ecal_cluster_e3x3_ = -1;
	float gsf_ecal_cluster_e5x5_ = -1;
	float gsf_ecal_cluster_covEtaEta_ = -42;
	float gsf_ecal_cluster_covEtaPhi_ = -42;
	float gsf_ecal_cluster_covPhiPhi_ = -42;
	float gsf_ecal_cluster_ematrix_[ECAL_CLUSTER_SIZE][ECAL_CLUSTER_SIZE] = {{0}};

	float gsf_hcal_cluster_e_ = -1;
	float gsf_hcal_cluster_eta_ = -1;
	float gsf_hcal_cluster_deta_ = -42;
	float gsf_hcal_cluster_dphi_ = -42;

	bool gsf_ktf_same_ecal_ = false;
	bool gsf_ktf_same_hcal_ = false;

	float ktf_ecal_cluster_e_ = -1;
	float ktf_ecal_cluster_ecorr_ = -1;
	float ktf_ecal_cluster_eta_ = -1;
	float ktf_ecal_cluster_deta_ = -42;
	float ktf_ecal_cluster_dphi_ = -42;
	float ktf_ecal_cluster_e3x3_ = -1;
	float ktf_ecal_cluster_e5x5_ = -1;
	float ktf_ecal_cluster_covEtaEta_ = -42;
	float ktf_ecal_cluster_covEtaPhi_ = -42;
	float ktf_ecal_cluster_covPhiPhi_ = -42;
	float ktf_ecal_cluster_ematrix_[ECAL_CLUSTER_SIZE][ECAL_CLUSTER_SIZE] = {{0}};
	float ktf_ecal_cluster_r9_ = -0.1;
	float ktf_ecal_cluster_circularity_ = -0.1;

	float ktf_hcal_cluster_e_ = -1;
	float ktf_hcal_cluster_eta_ = -1;
	float ktf_hcal_cluster_deta_ = -42;
	float ktf_hcal_cluster_dphi_ = -42;


// 	float gsf_dEdx_[NHITS_MAX] = {0};
// 	float gsf_hit_costh_impact_[NHITS_MAX] = {0};
// 	float gsf_brem_ecal_map_[NHITS_MAX][BREM_WINDOW_ETA][BREM_WINDOW_PHI] = {{{0}}};
// 	float gsf_brem_hcal_map_[NHITS_MAX][BREM_WINDOW_ETA][BREM_WINDOW_PHI] = {{{0}}};
// 	float gsf_brem_hit_map_ [NHITS_MAX][BREM_WINDOW_ETA][BREM_WINDOW_PHI] = {{{0}}};

// 	float gsf_cluster_ecal_map_[CLUSTER_WINDOW_ETA][CLUSTER_WINDOW_PHI]    = {{0}};
// 	float gsf_cluster_hit_map_[CLUSTER_WINDOW_ETA][CLUSTER_WINDOW_PHI]     = {{0}};
// 	float gsf_ktf_cluster_hit_map_[CLUSTER_WINDOW_ETA][CLUSTER_WINDOW_PHI] = {{0}};
// 	float gsf_cluster_hcal_map_[CLUSTER_WINDOW_ETA][CLUSTER_WINDOW_PHI]    = {{0}};
// 	float gsf_cluster_other_map_[CLUSTER_WINDOW_ETA][CLUSTER_WINDOW_PHI]   = {{0}};

// 	float pfecal_correctedEnergy[NECAL_PFCLUSTERS] = {0};
// 	float pfecal_deta_impact[NECAL_PFCLUSTERS] = {0};

	// 


  float core_shFracHits_ = -1.;

  // Track-Cluster matching //////////

  float match_SC_EoverP_ = -1.;

  float match_SC_dEta_ = -1.;
  float match_SC_dPhi_ = -1.;

  float match_seed_EoverP_ = -1.;
  float match_seed_EoverPout_ = -1.;

  float match_seed_dEta_ = -1.;
  float match_seed_dPhi_ = -1.;
  float match_seed_dEta_vtx_ = -1.;

  float match_eclu_EoverP_ = -1.;
  float match_eclu_EoverPout_ = -1.;

  float match_eclu_dEta_ = -1.;
  float match_eclu_dPhi_ = -1.;

  // Fiducial flags (booleans) //////////

  int fiducial_isEB_ = -1;
  int fiducial_isEE_ = -1;

  int fiducial_isGap_ = -1;
  int fiducial_isEBEEGap_ = -1;
  int fiducial_isEBGap_ = -1;
  int fiducial_isEBEtaGap_ = -1;
  int fiducial_isEBPhiGap_ = -1;
  int fiducial_isEEGap_ = -1;
  int fiducial_isEEDeeGap_ = -1;
  int fiducial_isEERingGap_ = -1;

  // Shower shape //////////
  
  float shape_sigmaEtaEta_ = -1.;
  float shape_sigmaIetaIeta_ = -1.;
  float shape_sigmaIphiIphi_ = -1.;

  float shape_e1x5_ = -1.;
  float shape_e2x5Max_ = -1.;
  float shape_e5x5_ = -1.;

  float shape_r9_ = -1.;

  float shape_HoverE_ = -1.;
  float shape_HoverEBc_ = -1.;

  float shape_hcalDepth1_ = -1.;
  float shape_hcalDepth2_ = -1.;
  float shape_hcalDepth1Bc_ = -1.;
  float shape_hcalDepth2Bc_ = -1.;
  int shape_nHcalTowersBc_ = -1;

  float shape_eLeft_ = -1.;
  float shape_eRight_ = -1.;
  float shape_eTop_ = -1.;
  float shape_eBottom_ = -1.;

  // full 5x5

  float shape_full5x5_sigmaEtaEta_ = -1.;
  float shape_full5x5_sigmaIetaIeta_ = -1.;
  float shape_full5x5_sigmaIphiIphi_ = -1.;
  float shape_full5x5_circularity_ = -1.;

  float shape_full5x5_e1x5_ = -1.;
  float shape_full5x5_e2x5Max_ = -1.;
  float shape_full5x5_e5x5_ = -1.;

  float shape_full5x5_r9_ = -1.;

  float shape_full5x5_HoverE_ = -1.;
  float shape_full5x5_HoverEBc_ = -1.;

  float shape_full5x5_hcalDepth1_ = -1.;
  float shape_full5x5_hcalDepth2_ = -1.;
  float shape_full5x5_hcalDepth1Bc_ = -1.;
  float shape_full5x5_hcalDepth2Bc_ = -1.;

  float shape_full5x5_eLeft_ = -1.;
  float shape_full5x5_eRight_ = -1.;
  float shape_full5x5_eTop_ = -1.;
  float shape_full5x5_eBottom_ = -1.;
  
  // Isolation variables //////////
  // Conversion rejection //////////
  // PFlow info //////////
  // Preselection and ambiguity //////////
  // Corrections //////////
  // ???
  
  // Brem fractions and classification //////////

  float brem_frac_ = -1.;
  float brem_fracTrk_ = -1.;
  float brem_fracSC_ = -1.;
  int brem_N_ = -1;
  int p4kind_ = -1;

  // SuperClusters //////////
	std::vector<float> sc_cluster_eta_;
	std::vector<float> sc_cluster_phi_;

  float sc_etaWidth_ = -1.;
  float sc_phiWidth_ = -1.;

  float sc_ps_EoverEraw_ = -1.;
  float sc_E_ = -1.;
  float sc_Et_ = -1.;

  float sc_eta_ = -1.;
  float sc_phi_ = -1.;

  float sc_RawE_ = -1.;
  int sc_Nclus_ = -1;


};

#endif
