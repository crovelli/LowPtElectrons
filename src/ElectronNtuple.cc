#include "LowPtElectrons/LowPtElectrons/interface/ElectronNtuple.h"
#include "TTree.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <sstream>
#include <cmath>

using namespace reco;
using namespace edm;
void ElectronNtuple::link_tree(TTree *tree) {
	tree->Branch("run",  &run_ , "run/i");
	tree->Branch("lumi", &lumi_, "lumi/i");
	tree->Branch("evt",  &evt_ , "evt/i");

	tree->Branch("is_e", &is_e_, "is_e/O");
	tree->Branch("is_e_not_matched", &is_e_not_matched_, "is_e_not_matched/O");
	tree->Branch("is_other", &is_other_, "is_other/O");

	tree->Branch("rho", &rho_, "rho/f");

	tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
	tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
	tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
	tree->Branch("gen_e", &gen_e_, "gen_e/f");
	tree->Branch("gen_p", &gen_p_, "gen_p/f");
	tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");

	tree->Branch("trk_pt",				 	&trk_pt_				   , "trk_pt/f");
	tree->Branch("trk_eta",		 	  &trk_eta_		       , "trk_eta/f");
	tree->Branch("trk_phi",		 	  &trk_phi_     		 , "trk_phi/f");
	tree->Branch("trk_p", &trk_p_, "trk_p/f");
	tree->Branch("trk_charge", &trk_charge_, "trk_charge/I");
	tree->Branch("trk_nhits",			&trk_nhits_        , "trk_nhits/I");
	tree->Branch("trk_missing_inner_hits", &trk_missing_inner_hits_, "trk_missing_inner_hits/I");

	tree->Branch("trk_high_purity",&trk_high_purity_	 , "trk_high_purity/i");
	tree->Branch("trk_dxy",			  &trk_dxy_		  		 , "trk_dxy/f");
	tree->Branch("trk_dxy_err",		&trk_dxy_err_			 , "trk_dxy_err/f");
	tree->Branch("trk_dz",			  &trk_dz_		  		 , "trk_dz/f");
	tree->Branch("trk_dz_err",		&trk_dz_err_			 , "trk_dz_err/f");
	tree->Branch("trk_inp",  			&trk_inp_  	    	 , "trk_inp/f");
	tree->Branch("trk_outp",	  		&trk_outp_	  	   , "trk_outp/f");
	tree->Branch("trk_chi2red",    &trk_chi2red_      , "trk_chi2red/f"); 

	tree->Branch("preid_bdtout",				 	&preid_bdtout_		  	, "preid_bdtout/f");
	tree->Branch("preid_trk_ecal_Deta",	&preid_trk_ecal_Deta_ , "preid_trk_ecal_Deta/f");
	tree->Branch("preid_trk_ecal_Dphi",	&preid_trk_ecal_Dphi_ , "preid_trk_ecal_Dphi/f");
	tree->Branch("preid_e_over_p",			 	&preid_e_over_p_			, "preid_e_over_p/f");
	//stage 2, with GSF
	tree->Branch("preid_gsf_success"			, &preid_gsf_success_     , "preid_gsf_success/O");
	tree->Branch("preid_gsf_dpt"					, &preid_gsf_dpt_		  		, "preid_gsf_dpt/f");					
	tree->Branch("preid_trk_gsf_chiratio", &preid_trk_gsf_chiratio_, "preid_trk_gsf_chiratio/f");
	tree->Branch("preid_gsf_chi2red"     , &preid_gsf_chi2red_     , "preid_gsf_chi2red/f");     
	tree->Branch("preid_numGSF", &preid_numGSF_, "preid_gsf_chi2red/i");
	//step-wise standard selection
	tree->Branch("preid_trk_ecal_match", &preid_trk_ecal_match_, "preid_trk_ecal_match/O");
	tree->Branch("preid_trkfilter_pass", &preid_trkfilter_pass_, "preid_trkfilter_pass/O");
	tree->Branch("preid_mva_pass", &preid_mva_pass_, "preid_mva_pass/O");
	
	tree->Branch("gsf_pt",				 	&gsf_pt_				   , "gsf_pt/f");
	tree->Branch("gsf_eta",		 	  &gsf_eta_		       , "gsf_eta/f");
	tree->Branch("gsf_phi",		 	  &gsf_phi_     		 , "gsf_phi/f");
	tree->Branch("gsf_p", &gsf_p_, "gsf_p/f");
	tree->Branch("gsf_charge", &gsf_charge_, "gsf_charge/I");
	tree->Branch("gsf_nhits",			&gsf_nhits_        , "gsf_nhits/I");
	tree->Branch("gsf_missing_inner_hits", &gsf_missing_inner_hits_, "gsf_missing_inner_hits/I");

	tree->Branch("gsf_dxy",			  &gsf_dxy_		  		 , "gsf_dxy/f");
	tree->Branch("gsf_dxy_err",		&gsf_dxy_err_			 , "gsf_dxy_err/f");
	tree->Branch("gsf_dz",			  &gsf_dz_		  		 , "gsf_dz/f");
	tree->Branch("gsf_dz_err",		&gsf_dz_err_			 , "gsf_dz_err/f");
	tree->Branch("gsf_inp",  			&gsf_inp_  	    	 , "gsf_inp/f");
	tree->Branch("gsf_outp",  		&gsf_outp_	  	   , "gsf_outp/f");
	tree->Branch("gsf_chi2red",    &gsf_chi2red_      , "gsf_chi2red/f"); 
	tree->Branch("gsf_ntangents", &gsf_ntangents_, "gsf_ntangents/i");
	tree->Branch("gsf_hit_dpt", gsf_hit_dpt_, "gsf_hit_dpt[gsf_ntangents]/f");
	tree->Branch("gsf_hit_dpt_unc", gsf_hit_dpt_unc_, "gsf_hit_dpt_unc[gsf_ntangents]/f");

	//PFGSFTrack internal steps flags
	tree->Branch("pfgsf_gsf_has_ktf", &pfgsf_gsf_has_ktf_, "pfgsf_gsf_has_ktf/O");
	tree->Branch("pfgsf_ktf_is_fifthStep", &pfgsf_ktf_is_fifthStep_, "pfgsf_ktf_is_fifthStep/O");
	tree->Branch("pfgsf_gsf_ecalDriven", &pfgsf_gsf_ecalDriven_, "pfgsf_gsf_ecalDriven/O");
	tree->Branch("pfgsf_gsf_trackerDriven", &pfgsf_gsf_trackerDriven_, "pfgsf_gsf_trackerDriven/O");
	tree->Branch("pfgsf_valid_gsf_brem", &pfgsf_valid_gsf_brem_, "pfgsf_valid_gsf_brem/O");
	tree->Branch("pfgsf_passes_preselection", &pfgsf_passes_preselection_, "pfgsf_passes_preselection/O");
	tree->Branch("pfgsf_passes_selection", &pfgsf_passes_selection_, "pfgsf_passes_selection/O");

	tree->Branch("pfgsf_xclean_seedref", &pfgsf_xclean_seedref_, "pfgsf_xclean_seedref/O");
	tree->Branch("pfgsf_xclean_ECALDriven_too_few_hits", &pfgsf_xclean_ECALDriven_too_few_hits_, "pfgsf_xclean_ECALDriven_too_few_hits/O");
	tree->Branch("pfgsf_xclean_ECALDriven_bad_EoverP", &pfgsf_xclean_ECALDriven_bad_EoverP_, "pfgsf_xclean_ECALDriven_bad_EoverP/O");
	tree->Branch("pfgsf_xclean_TrkDriven_too_few_hits", &pfgsf_xclean_TrkDriven_too_few_hits_, "pfgsf_xclean_TrkDriven_too_few_hits/O");
	tree->Branch("pfgsf_xclean_TrkDriven_vs_ECALDriven", &pfgsf_xclean_TrkDriven_vs_ECALDriven_, "pfgsf_xclean_TrkDriven_vs_ECALDriven/O");
	tree->Branch("pfgsf_xclean_BothTrk_bad_EoverP", &pfgsf_xclean_BothTrk_bad_EoverP_, "pfgsf_xclean_BothTrk_bad_EoverP/O");
	tree->Branch("pfgsf_xclean_BothTrk_noECAL_match", &pfgsf_xclean_BothTrk_noECAL_match_, "pfgsf_xclean_BothTrk_noECAL_match/O");
	tree->Branch("pfgsf_xclean_AngularGsfCleaning", &pfgsf_xclean_AngularGsfCleaning_, "pfgsf_xclean_AngularGsfCleaning/O");
	tree->Branch("pfgsf_xclean_noECAL_match_AGAIN", &pfgsf_xclean_noECAL_match_AGAIN_, "pfgsf_xclean_noECAL_match_AGAIN/O");
	tree->Branch("pfgsf_xclean_FINAL", &pfgsf_xclean_FINAL_, "pfgsf_xclean_FINAL/O");

//	tree->Branch("pfgsf_gsf_cntr", &pfgsf_gsf_cntr_, "pfgsf_gsf_cntr/O");
//	tree->Branch("pfgsf_findpfref_seednullptr", &pfgsf_findpfref_seednullptr_, "pfgsf_findpfref_seednullptr/O");
//	tree->Branch("pfgsf_findpfref_castnullptr", &pfgsf_findpfref_castnullptr_, "pfgsf_findpfref_castnullptr/O");
//	tree->Branch("pfgsf_findpfref_trknullptr", &pfgsf_findpfref_trknullptr_, "pfgsf_findpfref_trknullptr/O");
//	tree->Branch("pfgsf_findpfref_nosharedhits", &pfgsf_findpfref_nosharedhits_, "pfgsf_findpfref_nosharedhits/O");
//	tree->Branch("pfgsf_findpfref_nodrmatch", &pfgsf_findpfref_nodrmatch_, "pfgsf_findpfref_nodrmatch/O");
//	tree->Branch("pfgsf_findpfref_blah", &pfgsf_findpfref_blah_, "pfgsf_findpfref_blah/O");
//	tree->Branch("pfgsf_findpfref_trkmatch", &pfgsf_findpfref_trkmatch_, "pfgsf_findpfref_trkmatch/O");
//	tree->Branch("pfgsf_findpfref_end", &pfgsf_findpfref_end_, "pfgsf_findpfref_end/O");

	//Middle steps
	tree->Branch("has_ele_core", &has_ele_core_, "has_ele_core/O");
	tree->Branch("has_pfEgamma", &has_pfEgamma_,  "has_pfEgamma/O");
	tree->Branch("has_pfGSF_trk", &has_pfGSF_trk_, "has_pfGSF_trk/O");
	tree->Branch("has_pfBlock_with_SC", &has_pfBlock_with_SC_, "has_pfBlock_with_SC/O");
	tree->Branch("has_pfBlock_with_ECAL", &has_pfBlock_with_ECAL_, "has_pfBlock_with_ECAL/O");
	tree->Branch("has_pfBlock", &has_pfBlock_, "has_pfBlock/O");
	tree->Branch("has_pfBlock_size", &has_pfBlock_size_, "has_pfBlock_size/f");
	tree->Branch("has_pfBlock_dr", &has_pfBlock_dr_, "has_pfBlock_dr/f");

	//bool has_pfEgamma_ = false;
	tree->Branch("ele_pt",				 	&ele_pt_				   , "ele_pt/f");
	tree->Branch("ele_p",				 	&ele_p_				   , "ele_p/f");
	tree->Branch("ele_eta",		 	  &ele_eta_		       , "ele_eta/f");
	tree->Branch("ele_phi",		 	  &ele_phi_     		 , "ele_phi/f");
	tree->Branch("ele_mvaIdV1",		 	  &ele_mvaIdV1_     		 , "ele_mvaIdV1/f");
	tree->Branch("ele_mvaIdV2",		 	  &ele_mvaIdV2_     		 , "ele_mvaIdV2/f");
	tree->Branch("ele_conv_vtx_fit_prob", &ele_conv_vtx_fit_prob_, "ele_conv_vtx_fit_prob/f");
	tree->Branch("ele_iso01", &ele_iso01_, "ele_iso01/f");
	tree->Branch("ele_iso02", &ele_iso02_, "ele_iso02/f");
	tree->Branch("ele_iso03", &ele_iso03_, "ele_iso03/f");
	tree->Branch("ele_iso04", &ele_iso04_, "ele_iso04/f");

	//Bottom up approach
	tree->Branch("gsf_ecal_cluster_e", &gsf_ecal_cluster_e_, "gsf_ecal_cluster_e/f");
	tree->Branch("gsf_ecal_cluster_ecorr", &gsf_ecal_cluster_ecorr_, "gsf_ecal_cluster_ecorr/f");
	tree->Branch("gsf_ecal_cluster_eta", &gsf_ecal_cluster_eta_, "gsf_ecal_cluster_eta/f");
	tree->Branch("gsf_ecal_cluster_deta", &gsf_ecal_cluster_deta_, "gsf_ecal_cluster_deta/f");
	tree->Branch("gsf_ecal_cluster_dphi", &gsf_ecal_cluster_dphi_, "gsf_ecal_cluster_dphi/f");
	tree->Branch("gsf_ecal_cluster_e3x3", &gsf_ecal_cluster_e3x3_, "gsf_ecal_cluster_e3x3/f");
	tree->Branch("gsf_ecal_cluster_e5x5", &gsf_ecal_cluster_e5x5_, "gsf_ecal_cluster_e5x5/f");
	tree->Branch("gsf_ecal_cluster_covEtaEta", &gsf_ecal_cluster_covEtaEta_, "gsf_ecal_cluster_covEtaEta/f");
	tree->Branch("gsf_ecal_cluster_covEtaPhi", &gsf_ecal_cluster_covEtaPhi_, "gsf_ecal_cluster_covEtaPhi/f");
	tree->Branch("gsf_ecal_cluster_covPhiPhi", &gsf_ecal_cluster_covPhiPhi_, "gsf_ecal_cluster_covPhiPhi/f");
	std::stringstream buffer;
	buffer << "[" << ECAL_CLUSTER_SIZE << "][" << ECAL_CLUSTER_SIZE << "]/f";
//	tree->Branch("gsf_ecal_cluster_ematrix", &gsf_ecal_cluster_ematrix_, ("gsf_ecal_cluster_ematrix"+buffer.str()).c_str());

	tree->Branch("gsf_hcal_cluster_e",    &gsf_hcal_cluster_e_,    "gsf_hcal_cluster_e/f");
	tree->Branch("gsf_hcal_cluster_eta",  &gsf_hcal_cluster_eta_,  "gsf_hcal_cluster_eta/f");
	tree->Branch("gsf_hcal_cluster_deta", &gsf_hcal_cluster_deta_, "gsf_hcal_cluster_deta/f");
	tree->Branch("gsf_hcal_cluster_dphi", &gsf_hcal_cluster_dphi_, "gsf_hcal_cluster_dphi/f");

	tree->Branch("gsf_ktf_same_ecal", &gsf_ktf_same_ecal_, "gsf_ktf_same_ecal/O");
	tree->Branch("gsf_ktf_same_hcal", &gsf_ktf_same_hcal_, "gsf_ktf_same_hcal/O");

	tree->Branch("ktf_ecal_cluster_e", &ktf_ecal_cluster_e_, "ktf_ecal_cluster_e/f");
	tree->Branch("ktf_ecal_cluster_ecorr", &ktf_ecal_cluster_ecorr_, "ktf_ecal_cluster_ecorr/f");
	tree->Branch("ktf_ecal_cluster_eta", &ktf_ecal_cluster_eta_, "ktf_ecal_cluster_eta/f");
	tree->Branch("ktf_ecal_cluster_deta", &ktf_ecal_cluster_deta_, "ktf_ecal_cluster_deta/f");
	tree->Branch("ktf_ecal_cluster_dphi", &ktf_ecal_cluster_dphi_, "ktf_ecal_cluster_dphi/f");
	tree->Branch("ktf_ecal_cluster_e3x3", &ktf_ecal_cluster_e3x3_, "ktf_ecal_cluster_e3x3/f");
	tree->Branch("ktf_ecal_cluster_e5x5", &ktf_ecal_cluster_e5x5_, "ktf_ecal_cluster_e5x5/f");
	tree->Branch("ktf_ecal_cluster_covEtaEta", &ktf_ecal_cluster_covEtaEta_, "ktf_ecal_cluster_covEtaEta/f");
	tree->Branch("ktf_ecal_cluster_covEtaPhi", &ktf_ecal_cluster_covEtaPhi_, "ktf_ecal_cluster_covEtaPhi/f");
	tree->Branch("ktf_ecal_cluster_covPhiPhi", &ktf_ecal_cluster_covPhiPhi_, "ktf_ecal_cluster_covPhiPhi/f");
//	tree->Branch("ktf_ecal_cluster_ematrix", &ktf_ecal_cluster_ematrix_, ("ktf_ecal_cluster_ematrix"+buffer.str()).c_str());
	tree->Branch("ktf_ecal_cluster_r9", &ktf_ecal_cluster_r9_, "ktf_ecal_cluster_r9/f");
	tree->Branch("ktf_ecal_cluster_circularity_", &ktf_ecal_cluster_circularity_, "ktf_ecal_cluster_circularity/f");

	tree->Branch("ktf_hcal_cluster_e",    &ktf_hcal_cluster_e_,    "ktf_hcal_cluster_e/f");
	tree->Branch("ktf_hcal_cluster_eta",  &ktf_hcal_cluster_eta_,  "ktf_hcal_cluster_eta/f");
	tree->Branch("ktf_hcal_cluster_deta", &ktf_hcal_cluster_deta_, "ktf_hcal_cluster_deta/f");
	tree->Branch("ktf_hcal_cluster_dphi", &ktf_hcal_cluster_dphi_, "ktf_hcal_cluster_dphi/f");

	// tree->Branch("gsf_dEdx", gsf_dEdx_, "gsf_dEdx[gsf_nhits]/f");
	// tree->Branch("gsf_hit_costh_impact", gsf_hit_costh_impact_, "gsf_hit_costh_impact[gsf_nhits]/f");

	// std::stringstream dim;
	// dim << "[" << BREM_WINDOW_ETA << "][" << BREM_WINDOW_PHI << "]/f";
	// tree->Branch("gsf_brem_ecal_map", gsf_brem_ecal_map_, ("gsf_brem_ecal_map[gsf_nhits]"+dim.str()).c_str());
	// tree->Branch("gsf_brem_hcal_map", gsf_brem_hcal_map_, ("gsf_brem_hcal_map[gsf_nhits]"+dim.str()).c_str());
	// tree->Branch("gsf_brem_hit_map", gsf_brem_hit_map_, ("gsf_brem_hit_map[gsf_nhits]"+dim.str()).c_str());

	// std::stringstream dims;
	// dims << "[" << CLUSTER_WINDOW_ETA << "][" << CLUSTER_WINDOW_PHI << "]/f";
	// tree->Branch("gsf_cluster_ecal_map", gsf_cluster_ecal_map_, ("gsf_cluster_ecal_map"+dims.str()).c_str());
	// tree->Branch("gsf_cluster_hit_map", gsf_cluster_hit_map_, ("gsf_cluster_hit_map"+dims.str()).c_str());
	// tree->Branch("gsf_ktf_cluster_hit_map", gsf_ktf_cluster_hit_map_, ("gsf_ktf_cluster_hit_map"+dims.str()).c_str());
	// tree->Branch("gsf_cluster_hcal_map", gsf_cluster_hcal_map_, ("gsf_cluster_hcal_map"+dims.str()).c_str());
	// tree->Branch("gsf_cluster_other_map", gsf_cluster_other_map_, ("gsf_cluster_other_map"+dims.str()).c_str());

  // Core //////////

  tree->Branch("core_shFracHits",&core_shFracHits_);

  // Track-Cluster matching //////////

  tree->Branch("match_SC_EoverP",&match_SC_EoverP_); 

  tree->Branch("match_SC_dEta",&match_SC_dEta_); 
  tree->Branch("match_SC_dPhi",&match_SC_dPhi_); 

  tree->Branch("match_seed_EoverP",&match_seed_EoverP_); 
  tree->Branch("match_seed_EoverPout",&match_seed_EoverPout_); 

  tree->Branch("match_seed_dEta",&match_seed_dEta_); 
  tree->Branch("match_seed_dPhi",&match_seed_dPhi_); 
  tree->Branch("match_seed_dEta_vtx",&match_seed_dEta_vtx_); 

  tree->Branch("match_eclu_EoverP",&match_eclu_EoverP_); 
  tree->Branch("match_eclu_EoverPout",&match_eclu_EoverPout_); 

  tree->Branch("match_eclu_dEta",&match_eclu_dEta_); 
  tree->Branch("match_eclu_dPhi",&match_eclu_dPhi_); 

  // Fiducial flags (booleans) //////////

  tree->Branch("fiducial_isEB",&fiducial_isEB_,"fiducial_isEB/I");
  tree->Branch("fiducial_isEE",&fiducial_isEE_,"fiducial_isEE/I"); 

  tree->Branch("fiducial_isGap",&fiducial_isGap_,"fiducial_isGap/I"); 
  tree->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap_,"fiducial_isEBEEGap/I"); 
  tree->Branch("fiducial_isEBGap",&fiducial_isEBGap_,"fiducial_isEBGap/I"); 
  tree->Branch("fiducial_isEBEtaGap",&fiducial_isEBEtaGap_,"fiducial_isEBEtaGap/I"); 
  tree->Branch("fiducial_isEBPhiGap",&fiducial_isEBPhiGap_,"fiducial_isEBPhiGap/I"); 
  tree->Branch("fiducial_isEEGap",&fiducial_isEEGap_,"fiducial_isEEGap/I"); 
  tree->Branch("fiducial_isEEDeeGap",&fiducial_isEEDeeGap_,"fiducial_isEEDeeGap/I"); 
  tree->Branch("fiducial_isEERingGap",&fiducial_isEERingGap_,"fiducial_isEERingGap/I"); 

  // Shower shape //////////

  tree->Branch("shape_sigmaEtaEta",&shape_sigmaEtaEta_); 
  tree->Branch("shape_sigmaIetaIeta",&shape_sigmaIetaIeta_); 
  tree->Branch("shape_sigmaIphiIphi",&shape_sigmaIphiIphi_); 

  tree->Branch("shape_e1x5",&shape_e1x5_); 
  tree->Branch("shape_e2x5Max",&shape_e2x5Max_); 
  tree->Branch("shape_e5x5",&shape_e5x5_); 

  tree->Branch("shape_r9",&shape_r9_); 

  tree->Branch("shape_HoverE",&shape_HoverE_); 
  tree->Branch("shape_HoverEBc",&shape_HoverEBc_); 

  tree->Branch("shape_hcalDepth1",&shape_hcalDepth1_); 
  tree->Branch("shape_hcalDepth2",&shape_hcalDepth2_); 
  tree->Branch("shape_hcalDepth1Bc",&shape_hcalDepth1Bc_); 
  tree->Branch("shape_hcalDepth2Bc",&shape_hcalDepth2Bc_); 
  tree->Branch("shape_nHcalTowersBc",&shape_nHcalTowersBc_,"shape_nHcalTowersBc/I"); 

  tree->Branch("shape_eLeft",&shape_eLeft_); 
  tree->Branch("shape_eRight",&shape_eRight_); 
  tree->Branch("shape_eTop",&shape_eTop_); 
  tree->Branch("shape_eBottom",&shape_eBottom_); 

  // full 5x5

  tree->Branch("shape_full5x5_sigmaEtaEta",&shape_full5x5_sigmaEtaEta_); 
  tree->Branch("shape_full5x5_sigmaIetaIeta",&shape_full5x5_sigmaIetaIeta_); 
  tree->Branch("shape_full5x5_sigmaIphiIphi",&shape_full5x5_sigmaIphiIphi_); 
  tree->Branch("shape_full5x5_circularity",&shape_full5x5_circularity_); 

  tree->Branch("shape_full5x5_e1x5",&shape_full5x5_e1x5_); 
  tree->Branch("shape_full5x5_e2x5Max",&shape_full5x5_e2x5Max_); 
  tree->Branch("shape_full5x5_e5x5",&shape_full5x5_e5x5_); 

  tree->Branch("shape_full5x5_r9",&shape_full5x5_r9_); 

  tree->Branch("shape_full5x5_HoverE",&shape_full5x5_HoverE_); 
  tree->Branch("shape_full5x5_HoverEBc",&shape_full5x5_HoverEBc_); 

  tree->Branch("shape_full5x5_hcalDepth1",&shape_full5x5_hcalDepth1_); 
  tree->Branch("shape_full5x5_hcalDepth2",&shape_full5x5_hcalDepth2_); 
  tree->Branch("shape_full5x5_hcalDepth1Bc",&shape_full5x5_hcalDepth1Bc_); 
  tree->Branch("shape_full5x5_hcalDepth2Bc",&shape_full5x5_hcalDepth2Bc_); 

  tree->Branch("shape_full5x5_eLeft",&shape_full5x5_eLeft_); 
  tree->Branch("shape_full5x5_eRight",&shape_full5x5_eRight_); 
  tree->Branch("shape_full5x5_eTop",&shape_full5x5_eTop_); 
  tree->Branch("shape_full5x5_eBottom",&shape_full5x5_eBottom_); 
  
  // Brem fractions and classification //////////

  tree->Branch("brem_frac",&brem_frac_); 
  tree->Branch("brem_fracTrk",&brem_fracTrk_); 
  tree->Branch("brem_fracSC",&brem_fracSC_); 
  tree->Branch("brem_N",&brem_N_,"brem_N/I"); 
  tree->Branch("p4kind",&p4kind_,"p4kind/I"); 
  
  // SuperClusters //////////

  tree->Branch("sc_etaWidth",&sc_etaWidth_); 
  tree->Branch("sc_phiWidth",&sc_phiWidth_); 

  tree->Branch("sc_ps_EoverEraw",&sc_ps_EoverEraw_); 
  tree->Branch("sc_E",&sc_E_); 
  tree->Branch("sc_Et",&sc_Et_); 

  tree->Branch("sc_eta",&sc_eta_); 
  tree->Branch("sc_phi",&sc_phi_); 

  tree->Branch("sc_RawE",&sc_RawE_); 
  tree->Branch("sc_Nclus",&sc_Nclus_,"sc_Nclus/I"); 

}

//fillers
void ElectronNtuple::fill_evt(const edm::EventID &id) {
	run_  = id.run();						
  lumi_ = id.luminosityBlock();
  evt_  = id.event();          
}

void ElectronNtuple::fill_gen(const GenParticleRef genp) {
	gen_pt_  = genp->pt();
	gen_eta_ = genp->eta();
	gen_phi_ = genp->phi();
	gen_e_ = genp->energy();
	gen_p_ = genp->p();
	gen_charge_ = genp->charge();
}

void ElectronNtuple::fill_gsf_trk(const GsfTrackRef trk, const reco::BeamSpot &spot) {
  if ( trk.isNonnull() ) {
    // kine
    gsf_pt_ = trk->pt();
    gsf_eta_ = trk->eta();
    gsf_phi_ = trk->phi();
		gsf_p_ = trk->p();
		gsf_charge_ = trk->charge();
    gsf_inp_ = sqrt(trk->innerMomentum().mag2());
    gsf_outp_ = sqrt(trk->outerMomentum().mag2());
    gsf_dpt_ = ( gsf_inp_ > 0. ) ? fabs( gsf_outp_ - gsf_inp_ ) / gsf_inp_ : 0.; //@@ redundant?
    // quality
    gsf_nhits_ = trk->found();
    gsf_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
    gsf_chi2red_ = trk->normalizedChi2();
    // displ
    gsf_dxy_ = trk->dxy(spot);
    gsf_dxy_err_ = trk->dxyError();
    gsf_dz_ = trk->dz(spot.position());
    gsf_dz_err_ = trk->dzError();
		const auto & extra = trk->gsfExtra(); 
		gsf_ntangents_ = (extra->tangentsSize() > NHITS_MAX) ? NHITS_MAX : extra->tangentsSize();
		for(int i=0; i<gsf_ntangents_; i++) {
			gsf_hit_dpt_[i] = extra->tangents().at(i).deltaP().value();
			gsf_hit_dpt_unc_[i] = extra->tangents().at(i).deltaP().error();
		}
  } else {
    //@@ Shouldn't happen, but we take dummy values ...?
  }
}

void ElectronNtuple::fill_preid( const PreId &preid, const reco::BeamSpot &spot, const int num_gsf) {

  // Extract KF track parameters
  fill_ktf_trk( preid.trackRef(), spot );

  // ECAL/track matching parameters
  preid_e_over_p_ = preid.eopMatch();
  preid_trk_ecal_Deta_ = preid.geomMatching()[0];
  preid_trk_ecal_Dphi_ = preid.geomMatching()[1];

  // GSF tracks
  preid_gsf_success_ = false; //@@ ??
  // (p_out-p_in)/p_in from GSF track
  preid_gsf_dpt_ = preid.dpt();
  // Ratio of chi2 from GSF and KF tracks
  preid_trk_gsf_chiratio_ = preid.chi2Ratio();
  // Estimate of reduced chi2 for GSF track (assumes GSF and KF track have same d.o.f.)
  preid_gsf_chi2red_ = preid.gsfChi2();
	
  // MVA output
  preid_bdtout_ = preid.mva();

	//How many GSF it will seed
	preid_numGSF_ = num_gsf;

	//step-wise standard selection
  preid_trk_ecal_match_ = preid.ecalMatching();
	preid_trkfilter_pass_ = preid.trackFiltered();
	preid_mva_pass_ = preid.mvaSelected();
}

void ElectronNtuple::fill_ele(const reco::GsfElectronRef ele, float mvaid_v1, float mvaid_v2, float ele_conv_vtx_fit_prob, const std::vector<float>& iso_rings) {
	ele_p_			 = ele->p();
	ele_pt_			 = ele->pt();
	ele_eta_		 = ele->eta();
	ele_phi_     = ele->phi();
	ele_mvaIdV1_ = mvaid_v1;
	ele_mvaIdV2_ = mvaid_v2;
	ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob;
	ele_iso01_ = iso_rings.at(0);
	ele_iso02_ = iso_rings.at(1);
	ele_iso03_ = iso_rings.at(2);
	ele_iso04_ = iso_rings.at(3);
	fill_supercluster(ele);
}

void ElectronNtuple::fill_supercluster(const reco::GsfElectronRef ele) {

  if ( ele.isNull() ) { return; }

  // Charge //////////

  //scPixCharge();
  //isGsfCtfScPixChargeConsistent();
  //isGsfScPixChargeConsistent();
  //isGsfCtfChargeConsistent();

  // Core //////////

  core_shFracHits_ = ele->shFracInnerHits();

  // Track-Cluster matching //////////

  match_SC_EoverP_ = ele->eSuperClusterOverP();

  match_SC_dEta_ = ele->deltaEtaSuperClusterTrackAtVtx();
  match_SC_dPhi_ = ele->deltaPhiSuperClusterTrackAtVtx();

  match_seed_EoverP_     = ele->eSeedClusterOverP();
  match_seed_EoverPout_ = ele->eSeedClusterOverPout();

  match_seed_dEta_     = ele->deltaEtaSeedClusterTrackAtCalo();
  match_seed_dPhi_     = ele->deltaPhiSeedClusterTrackAtCalo();
  match_seed_dEta_vtx_ = ele->deltaEtaSeedClusterTrackAtVtx();

  match_eclu_EoverP_ = (1./ele->ecalEnergy()) - (1.0 / ele->p()); //@@ same as ele->gsfTrack()->p() ?
  match_eclu_EoverPout_ = ele->eEleClusterOverPout();

  match_eclu_dEta_ = ele->deltaEtaEleClusterTrackAtCalo();
  match_eclu_dPhi_ = ele->deltaPhiEleClusterTrackAtCalo();

  // Fiducial flags //////////

  fiducial_isEB_ = ele->isEB();
  fiducial_isEE_ = ele->isEE();

  fiducial_isGap_ = ele->isGap();
  fiducial_isEBEEGap_ = ele->isEBEEGap();
  fiducial_isEBGap_ = ele->isEBGap();
  fiducial_isEBEtaGap_ = ele->isEBEtaGap();
  fiducial_isEBPhiGap_ = ele->isEBPhiGap();
  fiducial_isEEGap_ = ele->isEEGap();
  fiducial_isEEDeeGap_ = ele->isEEDeeGap();
  fiducial_isEERingGap_ = ele->isEERingGap();

  // Shower shape //////////
  
  shape_sigmaEtaEta_   = ele->sigmaEtaEta();
  shape_sigmaIetaIeta_ = ele->sigmaIetaIeta();
  shape_sigmaIphiIphi_ = ele->sigmaIphiIphi();

  shape_e1x5_    = ele->e1x5();
  shape_e2x5Max_ = ele->e2x5Max();
  shape_e5x5_    = ele->e5x5();

  shape_r9_ = ele->r9();

  shape_HoverE_    = ele->hcalOverEcal();
  shape_HoverEBc_ = ele->hcalOverEcalBc();

  shape_hcalDepth1_    = ele->hcalDepth1OverEcal();
  shape_hcalDepth2_    = ele->hcalDepth2OverEcal();
  shape_hcalDepth1Bc_  = ele->hcalDepth1OverEcalBc();
  shape_hcalDepth2Bc_  = ele->hcalDepth2OverEcalBc();
  shape_nHcalTowersBc_ = ele->hcalTowersBehindClusters().size();

  shape_eLeft_ = ele->eLeft();
  shape_eRight_ = ele->eRight();
  shape_eTop_ = ele->eTop();
  shape_eBottom_ = ele->eBottom();

  // full 5x5

  shape_full5x5_sigmaEtaEta_   = ele->full5x5_sigmaEtaEta();
  shape_full5x5_sigmaIetaIeta_ = ele->full5x5_sigmaIetaIeta();
  shape_full5x5_sigmaIphiIphi_ = ele->full5x5_sigmaIphiIphi();
  shape_full5x5_circularity_   = 1. - ele->full5x5_e1x5() / ele->full5x5_e5x5();

  shape_full5x5_e1x5_    = ele->full5x5_e1x5();
  shape_full5x5_e2x5Max_ = ele->full5x5_e2x5Max();
  shape_full5x5_e5x5_    = ele->full5x5_e5x5();

  shape_full5x5_r9_ = ele->full5x5_r9();

  shape_full5x5_HoverE_    = ele->full5x5_hcalOverEcal();
  shape_full5x5_HoverEBc_ = ele->full5x5_hcalOverEcalBc();

  shape_full5x5_hcalDepth1_    = ele->full5x5_hcalDepth1OverEcal();
  shape_full5x5_hcalDepth2_    = ele->full5x5_hcalDepth2OverEcal();
  shape_full5x5_hcalDepth1Bc_  = ele->full5x5_hcalDepth1OverEcalBc();
  shape_full5x5_hcalDepth2Bc_  = ele->full5x5_hcalDepth2OverEcalBc();

  shape_full5x5_eLeft_ = ele->full5x5_eLeft();
  shape_full5x5_eRight_ = ele->full5x5_eRight();
  shape_full5x5_eTop_ = ele->full5x5_eTop();
  shape_full5x5_eBottom_ = ele->full5x5_eBottom();
  
  // Isolation variables //////////
  // Conversion rejection //////////
  // PFlow info //////////
  // Preselection and ambiguity //////////
  // Corrections //////////
  // ???
  
  // Brem fractions and classification //////////

  brem_frac_ = ele->fbrem();
  brem_fracTrk_ = ele->trackFbrem();
  brem_fracSC_ = ele->superClusterFbrem();
  brem_N_ = ele->numberOfBrems();

  p4kind_ = ele->candidateP4Kind();

  // SuperClusters //////////
  
  if ( ele->superCluster().isNull() ) { return; }
  const SuperClusterRef& sc = ele->superCluster();

  sc_etaWidth_ = sc->etaWidth();
  sc_phiWidth_ = sc->phiWidth();

  sc_ps_EoverEraw_ = sc->preshowerEnergy() / sc->rawEnergy();
  sc_E_ = sc->energy();
  sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );

  sc_eta_ = sc->eta();
  sc_phi_ = sc->phi();

  sc_RawE_ = sc->rawEnergy();
  sc_Nclus_ = sc->clustersSize();

  // Misc //////////

  //fixedGridRhoFastjetAll
  //1.0/ecalEnergy-1.0/trackMomentumAtVtx.R
  //ecalEnergy
  //ecalEnergyError
  //trackMomentumError
  //p4Error('P4_COMBINATION')
  //p4Error('P4_PFLOW_COMBINATION')
  //pfIsolationVariables.sumPUPt
  //classification

}

void ElectronNtuple::fill_ktf_trk( const TrackRef trk, const reco::BeamSpot &spot ) {
  if ( trk.isNonnull() ) {
    // kine
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
		trk_p_ = trk->p();
		trk_charge_ = trk->charge();
    trk_inp_ = sqrt( trk->innerMomentum().mag2() );
    trk_outp_ = sqrt( trk->outerMomentum().mag2() );
    trk_dpt_ = ( trk_inp_ > 0. ) ? fabs( trk_outp_ - trk_inp_ ) / trk_inp_ : 0.; //@@ redundant?
    // quality
    trk_nhits_ = trk->found();
    trk_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
    trk_high_purity_ = trk->quality( TrackBase::qualityByName("highPurity") );
    trk_chi2red_ = trk->normalizedChi2();
    // displ
    trk_dxy_ = trk->dxy(spot);
    trk_dxy_err_ = trk->dxyError();
    trk_dz_ = trk->dz(spot.position());
    trk_dz_err_ = trk->dzError();
  } else {
    //@@ Shouldn't happen, but we take dummy values ...?
  }
}

void ElectronNtuple::unpack_pfgsf_flags(const int flags) {
	pfgsf_gsf_has_ktf_       	= get_ith_bit(flags, 0);
	pfgsf_ktf_is_fifthStep_  	= get_ith_bit(flags, 1);
	pfgsf_gsf_ecalDriven_    	= get_ith_bit(flags, 2);
	pfgsf_gsf_trackerDriven_ 	= get_ith_bit(flags, 3);
	pfgsf_valid_gsf_brem_    	= get_ith_bit(flags, 4);
	pfgsf_passes_preselection_ = get_ith_bit(flags, 5);
	pfgsf_passes_selection_    = get_ith_bit(flags, 6);

	pfgsf_xclean_seedref_ = get_ith_bit(flags, 7);
	pfgsf_xclean_ECALDriven_too_few_hits_ = get_ith_bit(flags, 8);
	pfgsf_xclean_ECALDriven_bad_EoverP_ = get_ith_bit(flags, 9);
	pfgsf_xclean_TrkDriven_too_few_hits_ = get_ith_bit(flags, 10);
	pfgsf_xclean_TrkDriven_vs_ECALDriven_ = get_ith_bit(flags, 11);
	pfgsf_xclean_BothTrk_bad_EoverP_ = get_ith_bit(flags, 12);
	pfgsf_xclean_BothTrk_noECAL_match_ = get_ith_bit(flags, 13);
	pfgsf_xclean_AngularGsfCleaning_ = get_ith_bit(flags, 14);
	pfgsf_xclean_noECAL_match_AGAIN_ = get_ith_bit(flags, 15);
	pfgsf_xclean_FINAL_ = get_ith_bit(flags, 16);
	//@@ RB
//	pfgsf_gsf_cntr_ = get_ith_bit(flags, 31); // 31
//	pfgsf_findpfref_seednullptr_ = get_ith_bit(flags, 17);
//	pfgsf_findpfref_castnullptr_ = get_ith_bit(flags, 18);
//	pfgsf_findpfref_trknullptr_ = get_ith_bit(flags, 19);
//	pfgsf_findpfref_nosharedhits_ = get_ith_bit(flags, 20);
//	pfgsf_findpfref_nodrmatch_ = get_ith_bit(flags, 21);
//	pfgsf_findpfref_blah_ = get_ith_bit(flags, 22);
//	pfgsf_findpfref_trkmatch_ = get_ith_bit(flags, 23);
//	pfgsf_findpfref_end_ = get_ith_bit(flags, 24);

}

void ElectronNtuple::fill_GSF_ECAL_cluster_info(
	const reco::PFClusterRef cluster,
	const reco::PFTrajectoryPoint &gsf,
	noZS::EcalClusterLazyTools& tools
	) {
	gsf_ecal_cluster_e_   = cluster->energy();
	gsf_ecal_cluster_ecorr_  = cluster->correctedEnergy();
	gsf_ecal_cluster_eta_ = cluster->eta();
	gsf_ecal_cluster_deta_ = cluster->eta() - gsf.positionREP().eta();
	gsf_ecal_cluster_dphi_ = reco::deltaPhi(cluster->phi(), gsf.positionREP().phi());
	gsf_ecal_cluster_e3x3_ = tools.e3x3(*cluster);
	gsf_ecal_cluster_e5x5_ = tools.e5x5(*cluster);
	auto covs = tools.localCovariances(*cluster);
	gsf_ecal_cluster_covEtaEta_ = covs[0];
	gsf_ecal_cluster_covEtaPhi_ = covs[1];
	gsf_ecal_cluster_covPhiPhi_ = covs[2];

//	int cluster_window = (ECAL_CLUSTER_SIZE-1)/2;
//	DetId seedid = cluster->hitsAndFractions().front().first;
//	auto energies = tools.fullMatrixEnergy(
//		*cluster, seedid, -cluster_window, cluster_window,
//		-cluster_window, cluster_window);
//
//	for(size_t i=0; i<ECAL_CLUSTER_SIZE; i++) {
//		for(size_t j=0; j<ECAL_CLUSTER_SIZE; j++) {
//			gsf_ecal_cluster_ematrix_[i][j] = energies.at(i).at(j);
//		}
//	}
}

void ElectronNtuple::fill_GSF_HCAL_cluster_info(
	const reco::PFClusterRef cluster,
	const reco::PFTrajectoryPoint &gsf
	) {
	gsf_hcal_cluster_e_   = cluster->energy();
	gsf_hcal_cluster_eta_ = cluster->eta();
	gsf_hcal_cluster_deta_ = cluster->eta() - gsf.positionREP().eta();
	gsf_hcal_cluster_dphi_ = reco::deltaPhi(cluster->phi(), gsf.positionREP().phi());
}

void ElectronNtuple::fill_KTF_ECAL_cluster_info(
	const reco::PFClusterRef cluster,
	const reco::PFTrajectoryPoint &ktf,
	noZS::EcalClusterLazyTools& tools
	) {
	ktf_ecal_cluster_e_   = cluster->energy();
	ktf_ecal_cluster_ecorr_ = cluster->correctedEnergy();
	ktf_ecal_cluster_eta_ = cluster->eta();
	ktf_ecal_cluster_deta_ = cluster->eta() - ktf.positionREP().eta();
	ktf_ecal_cluster_dphi_ = reco::deltaPhi(cluster->phi(), ktf.positionREP().phi());
	ktf_ecal_cluster_e3x3_ = tools.e3x3(*cluster);
	ktf_ecal_cluster_e5x5_ = tools.e5x5(*cluster);
	auto covs = tools.localCovariances(*cluster);
	ktf_ecal_cluster_covEtaEta_ = covs[0];
	ktf_ecal_cluster_covEtaPhi_ = covs[1];
	ktf_ecal_cluster_covPhiPhi_ = covs[2];

	ktf_ecal_cluster_r9_ = ktf_ecal_cluster_e3x3_/ktf_ecal_cluster_e_;
	float e1x5 = tools.e1x5(*cluster);
	ktf_ecal_cluster_circularity_ = (ktf_ecal_cluster_e5x5_ > 0) ? 1 - e1x5/ktf_ecal_cluster_e5x5_ : -0.1;

//	int cluster_window = (ECAL_CLUSTER_SIZE-1)/2;
//	DetId seedid = cluster->hitsAndFractions().front().first;
//	auto energies = tools.fullMatrixEnergy(
//		*cluster, seedid, -cluster_window, cluster_window,
//		-cluster_window, cluster_window);
//
//	for(size_t i=0; i<ECAL_CLUSTER_SIZE; i++) {
//		for(size_t j=0; j<ECAL_CLUSTER_SIZE; j++) {
//			ktf_ecal_cluster_ematrix_[i][j] = energies.at(i).at(j);
//		}
//	}
}

void ElectronNtuple::fill_KTF_HCAL_cluster_info(
	const reco::PFClusterRef cluster,
	const reco::PFTrajectoryPoint &ktf
	) {
	ktf_hcal_cluster_e_   = cluster->energy();
	ktf_hcal_cluster_eta_ = cluster->eta();
	ktf_hcal_cluster_deta_ = cluster->eta() - ktf.positionREP().eta();
	ktf_hcal_cluster_dphi_ = reco::deltaPhi(cluster->phi(), ktf.positionREP().phi());
}


