#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDHeavyObjectCache.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TMath.h"
#include "TTree.h"
#include <iostream>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//
void IDSlimNtuple::link_tree( TTree *tree ) {

  std::cout<<"I am running IDSlimNtuple::link_tree"<<std::endl; 

  // general
  tree->Branch("evt",  &evt_ , "evt/i");  
  tree->Branch("weight", &weight_, "weight/f"); 
  tree->Branch("is_e", &is_e_, "is_e/O");
  tree->Branch("is_e_not_matched", &is_e_not_matched_, "is_e_not_matched/O"); 
  tree->Branch("is_other", &is_other_, "is_other/O");
  tree->Branch("is_egamma", &is_egamma_, "is_egamma/O");
  tree->Branch("has_trk", &has_trk_, "has_trk/O");
  tree->Branch("has_seed", &has_seed_, "has_seed/O");
  tree->Branch("has_gsf", &has_gsf_, "has_gsf/O");
  tree->Branch("has_ele", &has_ele_, "has_ele/O");
  if (largeNtuple) {
    tree->Branch("run",  &run_ , "run/i");
    tree->Branch("lumi", &lumi_, "lumi/i");
    tree->Branch("is_aod", &is_aod_, "is_aod/i");
    tree->Branch("is_mc", &is_mc_, "is_mc/i");
    tree->Branch("rho", &rho_, "rho/f");
  }

  // gen-level particles matched to reco-electron
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_tag_side", &gen_tag_side_, "gen_tag_side/I");   
  tree->Branch("gen_dR" , &gen_dR_ , "gen_dR/f" );
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_p",   &gen_p_,   "gen_p/f");
  if (largeNtuple) {
    tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");
    tree->Branch("gen_pdgid", &gen_pdgid_, "gen_pdgid/I");
    tree->Branch("gen_mom_pdgid", &gen_mom_pdgid_, "gen_mom_pdgid/I");
    tree->Branch("gen_gran_pdgid", &gen_gran_pdgid_, "gen_gran_pdgid/I");
    tree->Branch("gen_trk_dr" , &gen_trk_dr_ , "gen_trk_dr/f" );
    tree->Branch("gen_gsf_dr" , &gen_gsf_dr_ , "gen_gsf_dr/f" );
  }

  // GSF track associated to electron
  tree->Branch("gsf_dr", &gsf_dr_, "gsf_dr/f");
  tree->Branch("gsf_p",  &gsf_p_, "gsf_p/f");
  tree->Branch("gsf_pt", &gsf_pt_, "gsf_pt/f");
  tree->Branch("gsf_bdtout1", &seed_unbiased_, "gsf_bdtout1/f");
  tree->Branch("gsf_mode_p", &gsf_mode_p_, "gsf_mode_p/f");
  tree->Branch("gsf_mode_pt", &gsf_mode_pt_, "gsf_mode_pt/f");
  tree->Branch("gsf_eta", &gsf_eta_, "gsf_eta/f");
  tree->Branch("gsf_phi", &gsf_phi_, "gsf_phi/f");
  tree->Branch("gsf_mode_eta", &gsf_mode_eta_, "gsf_mode_eta/f");
  tree->Branch("gsf_mode_phi", &gsf_mode_phi_, "gsf_mode_phi/f");
  if (largeNtuple) {
    tree->Branch("gsf_bdtout2", &seed_ptbiased_, "gsf_bdtout2/f");
    tree->Branch("gsf_charge", &gsf_charge_, "gsf_charge/I");
    tree->Branch("gsf_inp", &gsf_inp_, "gsf_inp/f");
    tree->Branch("gsf_outp", &gsf_outp_, "gsf_outp/f");
    tree->Branch("gsf_missing_inner_hits", &gsf_missing_inner_hits_, "gsf_missing_inner_hits/I");
    tree->Branch("gsf_dxy",  &gsf_dxy_, "gsf_dxy/f");
    tree->Branch("gsf_dxy_err",&gsf_dxy_err_, "gsf_dxy_err/f");
    tree->Branch("gsf_dz",  &gsf_dz_, "gsf_dz/f");
    tree->Branch("gsf_dz_err",&gsf_dz_err_, "gsf_dz_err/f");
    tree->Branch("gsf_x",  &gsf_x_, "gsf_x/f");
    tree->Branch("gsf_y",  &gsf_y_, "gsf_y/f");
    tree->Branch("gsf_z",  &gsf_z_, "gsf_z/f");
  }

  // General track associated to electron
  tree->Branch("trk_dr", &trk_dr_, "trk_dr/f");
  tree->Branch("trk_p", &trk_p_, "trk_p/f");
  tree->Branch("trk_pt", &trk_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &trk_eta_, "trk_eta/f");
  tree->Branch("trk_phi", &trk_phi_, "trk_phi/f");
  if (largeNtuple) {
    tree->Branch("trk_charge", &trk_charge_, "trk_charge/I");
    tree->Branch("trk_inp", &trk_inp_, "trk_inp/f");
    tree->Branch("trk_outp", &trk_outp_, "trk_outp/f");
    tree->Branch("pdg_id", &pdg_id_, "pdg_id/I");
    tree->Branch("trk_nhits", &trk_nhits_, "trk_nhits/I");
    tree->Branch("trk_missing_inner_hits", &trk_missing_inner_hits_, "trk_missing_inner_hits/I"); 
    tree->Branch("trk_chi2red", &trk_chi2red_, "trk_chi2red/f");
    tree->Branch("trk_dxy", &trk_dxy_, "trk_dxy/f");
    tree->Branch("trk_dxy_err", &trk_dxy_err_, "trk_dxy_err/f");
    tree->Branch("trk_dz", &trk_dz_, "trk_dz/f");
    tree->Branch("trk_dz_err", &trk_dz_err_, "trk_dz_err/f");
    tree->Branch("trk_dEdx1", &trk_dEdx1_, "trk_dEdx1/f");
    tree->Branch("trk_dEdx1_Nm", &trk_dEdx1_Nm_, "trk_dEdx1_Nm/I");
    tree->Branch("trk_dEdx1_NSm", &trk_dEdx1_NSm_, "trk_dEdx1_NSm/I");
    tree->Branch("trk_high_purity", &trk_high_purity_, "trk_high_purity/I");
  }

  // Electron - kinematics
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");
  tree->Branch("core_shFracHits",&core_shFracHits_,"core_shFracHits/f");
  tree->Branch("fiducial_isEB",&fiducial_isEB_,"fiducial_isEB/I");
  tree->Branch("fiducial_isEE",&fiducial_isEE_,"fiducial_isEE/I"); 
  tree->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap_,"fiducial_isEBEEGap/I"); 
  if (largeNtuple) {
    tree->Branch("p4kind",&p4kind_,"p4kind/I"); 
    tree->Branch("ele_p_atvtx",&ele_p_atvtx_,"ele_p_atvtx/F"); 
    tree->Branch("ele_p_atcalo",&ele_p_atcalo_,"ele_p_atcalo/F"); 
  }

  // Electron - charge
  if (largeNtuple) {
    tree->Branch("chPix",&chPix_,"chPix/I");
    tree->Branch("chGCP",&chGCP_,"chGCP/I");
    tree->Branch("chGP",&chGP_,"chGP/I");
    tree->Branch("chGC",&chGC_,"chGC/I");
  }

  // Electron - id
  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_id", &ele_mva_id_, "ele_mva_id/I");
  tree->Branch("eid_ele_pt", &eid_ele_pt_, "eid_ele_pt/f");
  tree->Branch("eid_sc_eta", &eid_sc_eta_, "eid_sc_eta/f");
  tree->Branch("eid_shape_full5x5_sigmaIetaIeta", &eid_shape_full5x5_sigmaIetaIeta_, "eid_shape_full5x5_sigmaIetaIeta/f");
  tree->Branch("eid_shape_full5x5_sigmaIphiIphi", &eid_shape_full5x5_sigmaIphiIphi_, "eid_shape_full5x5_sigmaIphiIphi/f");
  tree->Branch("eid_shape_full5x5_circularity", &eid_shape_full5x5_circularity_, "eid_shape_full5x5_circularity/f");
  tree->Branch("eid_shape_full5x5_r9", &eid_shape_full5x5_r9_, "eid_shape_full5x5_r9/f");  
  tree->Branch("eid_sc_etaWidth", &eid_sc_etaWidth_, "eid_sc_etaWidth/f");
  tree->Branch("eid_sc_phiWidth", &eid_sc_phiWidth_, "eid_sc_phiWidth/f");
  tree->Branch("eid_shape_full5x5_HoverE", &eid_shape_full5x5_HoverE_, "eid_shape_full5x5_HoverE/f");
  tree->Branch("eid_trk_nhits", &eid_trk_nhits_, "eid_trk_nhits/f");
  tree->Branch("eid_trk_chi2red", &eid_trk_chi2red_, "eid_trk_chi2red/f");
  tree->Branch("eid_gsf_chi2red", &eid_gsf_chi2red_, "eid_gsf_chi2red/f");
  tree->Branch("eid_brem_frac", &eid_brem_frac_, "eid_brem_frac/f");
  tree->Branch("eid_gsf_nhits", &eid_gsf_nhits_, "eid_gsf_nhits/f");
  tree->Branch("eid_match_SC_EoverP", &eid_match_SC_EoverP_, "eid_match_SC_EoverP/f");
  tree->Branch("eid_match_eclu_EoverP", &eid_match_eclu_EoverP_, "eid_match_eclu_EoverP/f");
  tree->Branch("eid_match_SC_dEta", &eid_match_SC_dEta_, "eid_match_SC_dEta/f");
  tree->Branch("eid_match_SC_dPhi", &eid_match_SC_dPhi_, "eid_match_SC_dPhi/f");
  tree->Branch("eid_match_seed_dEta", &eid_match_seed_dEta_, "eid_match_seed_dEta/f");
  tree->Branch("eid_sc_E", &eid_sc_E_,   "eid_sc_E/f");
  tree->Branch("eid_trk_p", &eid_trk_p_, "eid_trk_p/f");
  tree->Branch("eid_rho", &eid_rho_, "eid_rho/f");
  if (largeNtuple) {
    tree->Branch("ele_conv_vtx_fit_prob", &ele_conv_vtx_fit_prob_, "ele_conv_vtx_fit_prob/f");
  }

  // Electron - isolation
  if (largeNtuple) {
    tree->Branch("ele_sumPhotonEt",        &ele_sumPhotonEt_,        "ele_sumPhotonEt/f");
    tree->Branch("ele_sumChargedHadronPt", &ele_sumChargedHadronPt_, "ele_sumChargedHadronPt/f");
    tree->Branch("ele_sumNeutralHadronEt", &ele_sumNeutralHadronEt_, "ele_sumNeutralHadronEt/f");
  }

  // Electron - further track-Cluster matching
  if (largeNtuple) {
    tree->Branch("match_seed_EoverP",&match_seed_EoverP_); 
    tree->Branch("match_seed_EoverPout",&match_seed_EoverPout_); 
    tree->Branch("match_seed_dPhi",&match_seed_dPhi_); 
    tree->Branch("match_seed_dEta_vtx",&match_seed_dEta_vtx_); 
    tree->Branch("match_eclu_EoverPout",&match_eclu_EoverPout_); 
    tree->Branch("match_eclu_dEta",&match_eclu_dEta_); 
    tree->Branch("match_eclu_dPhi",&match_eclu_dPhi_); 
  }

  // Electron energy regression                                                                                                                  
  tree->Branch("pre_ecal",&pre_ecal_);
  tree->Branch("pre_ecaltrk",&pre_ecaltrk_);
  tree->Branch("post_ecal",&post_ecal_);
  tree->Branch("post_ecaltrk",&post_ecaltrk_);

  // Electron - further full 5x5 shower shapes
  if (largeNtuple) {
    tree->Branch("shape_full5x5_e1x5",&shape_full5x5_e1x5_); 
    tree->Branch("shape_full5x5_e2x5Max",&shape_full5x5_e2x5Max_); 
    tree->Branch("shape_full5x5_e5x5",&shape_full5x5_e5x5_);
    tree->Branch("shape_full5x5_eLeft",&shape_full5x5_eLeft_); 
    tree->Branch("shape_full5x5_eRight",&shape_full5x5_eRight_); 
    tree->Branch("shape_full5x5_eTop",&shape_full5x5_eTop_); 
    tree->Branch("shape_full5x5_eBottom",&shape_full5x5_eBottom_); 
    tree->Branch("shape_full5x5_HoverEBc",&shape_full5x5_HoverEBc_); 
    tree->Branch("shape_full5x5_hcalDepth1",&shape_full5x5_hcalDepth1_); 
    tree->Branch("shape_full5x5_hcalDepth2",&shape_full5x5_hcalDepth2_); 
    tree->Branch("shape_full5x5_hcalDepth1Bc",&shape_full5x5_hcalDepth1Bc_); 
    tree->Branch("shape_full5x5_hcalDepth2Bc",&shape_full5x5_hcalDepth2Bc_); 
  }

  // Electron - brem fractions
  tree->Branch("brem_fracTrk",&brem_fracTrk_); 
  tree->Branch("brem_fracSC",&brem_fracSC_); 
  tree->Branch("brem_N",&brem_N_,"brem_N/I"); 

  // SuperCluster associated to the electron
  tree->Branch("sc_goodSeed",&sc_goodSeed_,"sc_goodSeed/O");
  tree->Branch("sc_Nclus",&sc_Nclus_,"sc_Nclus/I"); 
  tree->Branch("sc_Nclus_deta01",&sc_Nclus_deta01_,"sc_Nclus_deta01/I"); 
  tree->Branch("sc_Nclus_deta02",&sc_Nclus_deta02_,"sc_Nclus_deta02/I"); 
  tree->Branch("sc_Nclus_deta03",&sc_Nclus_deta03_,"sc_Nclus_deta03/I"); 
  if (largeNtuple) {
    tree->Branch("sc_Et",&sc_Et_); 
    tree->Branch("sc_E_ps",&sc_E_ps_,"sc_E_ps/F");
    tree->Branch("sc_E_ps1",&sc_E_ps1_,"sc_E_ps1/F");
    tree->Branch("sc_E_ps2",&sc_E_ps2_,"sc_E_ps2/F");
  }

  // Clusters making the SC 
  bool cluster_in_rootuple=false;     
  if(cluster_in_rootuple){   
    tree->Branch("sc_cluster_et",  sc_cluster_et_, "sc_cluster_et[sc_Nclus]/F");
    tree->Branch("sc_cluster_E",   sc_cluster_E_,  "sc_cluster_E[sc_Nclus]/F");
    tree->Branch("sc_cluster_eta", sc_cluster_eta_, "sc_cluster_eta[sc_Nclus]/F");
    tree->Branch("sc_cluster_phi", sc_cluster_phi_, "sc_cluster_phi[sc_Nclus]/F");
    tree->Branch("sc_cluster_nxtal", sc_cluster_nxtal_, "sc_cluster_nxtal[sc_Nclus]/I");
    tree->Branch("sc_cluster_e1x3", sc_cluster_e1x3_, "sc_cluster_e1x3[sc_Nclus]/F");
    tree->Branch("sc_cluster_e1x5", sc_cluster_e1x5_, "sc_cluster_e1x5[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2x2", sc_cluster_e2x2_, "sc_cluster_e2x2[sc_Nclus]/F");
    tree->Branch("sc_cluster_e3x3", sc_cluster_e3x3_, "sc_cluster_e3x3[sc_Nclus]/F");
    tree->Branch("sc_cluster_e5x5", sc_cluster_e5x5_, "sc_cluster_e5x5[sc_Nclus]/F");
    tree->Branch("sc_cluster_eMax", sc_cluster_eMax_, "sc_cluster_eMax[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2nd", sc_cluster_e2nd_, "sc_cluster_e2nd[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2x5Right",  sc_cluster_e2x5Right_,  "sc_cluster_e2x5Right[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2x5Left",   sc_cluster_e2x5Left_,   "sc_cluster_e2x5Left[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2x5Top",    sc_cluster_e2x5Top_,    "sc_cluster_e2x5Top[sc_Nclus]/F");
    tree->Branch("sc_cluster_e2x5Bottom", sc_cluster_e2x5Bottom_, "sc_cluster_e2x5Bottom[sc_Nclus]/F");
    tree->Branch("sc_cluster_eRight",  sc_cluster_eRight_,  "sc_cluster_eRight[sc_Nclus]/F");
    tree->Branch("sc_cluster_eLeft",   sc_cluster_eLeft_,   "sc_cluster_eLeft[sc_Nclus]/F");
    tree->Branch("sc_cluster_eTop",    sc_cluster_eTop_,    "sc_cluster_eTop[sc_Nclus]/F");
    tree->Branch("sc_cluster_eBottom", sc_cluster_eBottom_, "sc_cluster_eBottom[sc_Nclus]/F");
    tree->Branch("sc_cluster_eMaxOver2x2", sc_cluster_eMaxOver2x2_, "sc_cluster_eMaxOver2x2[sc_Nclus]/F");
    tree->Branch("sc_cluster_eMaxOver3x3", sc_cluster_eMaxOver3x3_, "sc_cluster_eMaxOver3x3[sc_Nclus]/F");
    tree->Branch("sc_cluster_eMaxOver1x3", sc_cluster_eMaxOver1x3_, "sc_cluster_eMaxOver1x3[sc_Nclus]/F");
  } else {
    tree->Branch("sc_clus1_E",      &sc_clus1_E_,      "sc_clus1_E/F");
    tree->Branch("sc_clus1_E_ov_p", &sc_clus1_E_ov_p_, "sc_clus1_E_ov_p/F");
    tree->Branch("sc_clus1_E_ov_E", &sc_clus1_E_ov_E_, "sc_clus1_E_ov_E/F");
    tree->Branch("sc_clus1_eta",    &sc_clus1_eta_,    "sc_clus1_eta/F");
    tree->Branch("sc_clus1_phi",    &sc_clus1_phi_,    "sc_clus1_phi/F");
    tree->Branch("sc_clus1_nxtal",  &sc_clus1_nxtal_,  "sc_clus1_nxtal/I");
    tree->Branch("sc_clus1_dphi",   &sc_clus1_dphi_,   "sc_clus1_dphi/F");
    tree->Branch("sc_clus1_deta",   &sc_clus1_deta_,   "sc_clus1_deta/F");
    tree->Branch("sc_clus1_ntrk_deta01",  &sc_clus1_ntrk_deta01_,  "sc_clus1_ntrk_deta01/F");
    if (largeNtuple) {
      tree->Branch("sc_clus1_et",     &sc_clus1_et_,     "sc_clus1_et/F");
    }
    //
    tree->Branch("sc_clus2_E",      &sc_clus2_E_,    "sc_clus2_E/F");
    tree->Branch("sc_clus2_E_ov_p", &sc_clus2_E_ov_p_,     "sc_clus2_E_ov_p/F");
    tree->Branch("sc_clus2_E_ov_E", &sc_clus2_E_ov_E_,     "sc_clus2_E_ov_E/F");
    tree->Branch("sc_clus2_eta",    &sc_clus2_eta_,  "sc_clus2_eta/F");
    tree->Branch("sc_clus2_phi",    &sc_clus2_phi_,  "sc_clus2_phi/F");
    tree->Branch("sc_clus2_nxtal",  &sc_clus2_nxtal_, "sc_clus2_nxtal/I");
    tree->Branch("sc_clus2_dphi",   &sc_clus2_dphi_,  "sc_clus2_dphi/F");
    tree->Branch("sc_clus2_deta",   &sc_clus2_deta_,  "sc_clus2_deta/F");
    tree->Branch("sc_clus2_ntrk_deta01",  &sc_clus2_ntrk_deta01_,  "sc_clus2_ntrk_deta01/F");
    if (largeNtuple) {
      tree->Branch("sc_clus2_et",     &sc_clus2_et_,   "sc_clus2_et/F");
    }
    //
    tree->Branch("sc_clus3_E",      &sc_clus3_E_,    "sc_clus3_E/F");
    tree->Branch("sc_clus3_E_ov_p", &sc_clus3_E_ov_p_,     "sc_clus3_E_ov_p/F");
    tree->Branch("sc_clus3_E_ov_E", &sc_clus3_E_ov_E_,     "sc_clus3_E_ov_E/F");
    tree->Branch("sc_clus3_eta",    &sc_clus3_eta_,  "sc_clus3_eta/F");
    tree->Branch("sc_clus3_phi",    &sc_clus3_phi_,  "sc_clus3_phi/F");
    tree->Branch("sc_clus3_nxtal",  &sc_clus3_nxtal_, "sc_clus3_nxtal/I");
    tree->Branch("sc_clus3_dphi",   &sc_clus3_dphi_,  "sc_clus3_dphi/F");
    tree->Branch("sc_clus3_deta",   &sc_clus3_deta_,  "sc_clus3_deta/F");
    tree->Branch("sc_clus3_ntrk_deta01",  &sc_clus3_ntrk_deta01_,  "sc_clus3_ntrk_deta01/F");
    if (largeNtuple) {
      tree->Branch("sc_clus3_et",     &sc_clus3_et_,   "sc_clus3_et/F");
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_evt( const edm::EventID& id ) {  
  if (largeNtuple) { 
    run_  = id.run();
    lumi_ = id.luminosityBlock();
  }
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_gen( const reco::GenParticlePtr genp ) {

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
  if (largeNtuple) {   
    gen_charge_ = genp->charge();
    gen_pdgid_ = genp->pdgId();
    if ( genp->mother(0) )
      gen_mom_pdgid_ = (*genp->mother(0)).pdgId();
    else 
      gen_mom_pdgid_ = 0;
    if ( (genp->mother(0))->mother(0) )
      gen_gran_pdgid_ = ( *(*genp->mother(0)).mother(0)).pdgId();
    else 
      gen_gran_pdgid_ = 0;
  }
}

void IDSlimNtuple::fill_gen( const pat::PackedGenParticleRef genp ) {  

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
  if (largeNtuple) {  
    gen_charge_ = genp->charge();
    gen_pdgid_ = genp->pdgId();
    if ( genp->mother(0) )
      gen_mom_pdgid_ = (*genp->mother(0)).pdgId();
    else 
      gen_mom_pdgid_ = 0;
    if ( (genp->mother(0))->mother(0) )
      gen_gran_pdgid_ = ( *(*genp->mother(0)).mother(0)).pdgId();
    else 
      gen_gran_pdgid_ = 0;
  }
}

void IDSlimNtuple::fill_gen_default() {

  gen_pt_  = -999.;
  gen_eta_ = -999.;
  gen_phi_ = -999.;
  gen_p_   = -999.;
  if (largeNtuple) {  
    gen_charge_ = -999; 
    gen_pdgid_  = -999;
    gen_mom_pdgid_  = -999; 
    gen_gran_pdgid_ = -999;
  }
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_trk( const reco::TrackPtr& trk,
			 const reco::BeamSpot& spot ) {       

  if ( trk.isNonnull() ) {   // should never happen
    // kine
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
    trk_p_ = trk->p();    
    if (largeNtuple) {
      trk_charge_ = trk->charge();
      if ( trk->extra().isAvailable() && trk->extra().isNonnull() ) {
	trk_inp_ = sqrt( trk->innerMomentum().mag2() );
	trk_outp_ = sqrt( trk->outerMomentum().mag2() );
      }
      // quality
      trk_nhits_ = trk->found();
      trk_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
      trk_chi2red_ = trk->normalizedChi2();
      trk_high_purity_=trk->quality( reco::TrackBase::qualityByName("highPurity") ) ;
      // displ
      trk_dxy_ = trk->dxy(spot);
      trk_dxy_err_ = trk->dxyError();
      trk_dz_ = trk->dz(spot.position());
      trk_dz_err_ = trk->dzError();        
    } 

  }
  
}

void IDSlimNtuple::fill_trk_dEdx( const reco::TrackPtr& trk,
				  std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx ) {

  if ( trk.isNonnull() ) {    // should never happen    

    if (largeNtuple) {
      const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[0]);
      const reco::DeDxData& dedx = dEdxTrack[trk];
      trk_dEdx1_=dedx.dEdx();
      trk_dEdx1_Nm_=dedx.numberOfMeasurements();
      trk_dEdx1_NSm_=dedx.numberOfSaturatedMeasurements();
    }
  }

}

void IDSlimNtuple::fill_trk_dEdx_default( ) {

  if (largeNtuple) {
    trk_dEdx1_     = -999.;
    trk_dEdx1_Nm_  = -999;
    trk_dEdx1_NSm_ = -999;
  }
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_bdt( double seed_unbiased, 
			 double seed_ptbiased ) {          
  seed_unbiased_ = seed_unbiased;
  if (largeNtuple) { 
    seed_ptbiased_ = seed_ptbiased;
  }
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_gsf( const reco::GsfTrackPtr gsf, 
			 const reco::BeamSpot& spot ) {        

  if ( gsf.isNull() ) {
    //@@ Shouldn't happen, but do we just take dummy values...? 
  } else {

    // Kinematics
    gsf_pt_ = gsf->pt();
    gsf_eta_ = gsf->eta();
    gsf_phi_ = gsf->phi();
    gsf_p_ = gsf->p();
    if (largeNtuple) {  
      gsf_charge_ = gsf->charge();
      if ( gsf->extra().isAvailable() && gsf->extra().isNonnull() ) {
	gsf_inp_ = sqrt(gsf->innerMomentum().mag2());
	gsf_outp_ = sqrt(gsf->outerMomentum().mag2());
      } else {
	gsf_inp_  = -999.;
	gsf_outp_ = -999.;
      }
    }

    // Kinematics (MODE)
    gsf_mode_pt_ = gsf->ptMode();
    gsf_mode_eta_ = gsf->etaMode();
    gsf_mode_phi_ = gsf->phiMode();
    gsf_mode_p_ = gsf->pMode();

    // Quality
    if (largeNtuple) 
      gsf_missing_inner_hits_ = gsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    
      // Displacement
    if (largeNtuple) {
      gsf_dxy_ = gsf->dxy(spot);
      gsf_dxy_err_ = gsf->dxyError();
      gsf_dz_ = gsf->dz(spot.position());
      gsf_dz_err_ = gsf->dzError();
    }
  } 
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_ele( const reco::GsfElectronPtr ele,
			 float mva_value,
			 int mva_id,
			 float ele_conv_vtx_fit_prob,
			 const double rho ) {       


  if ( ele.isNonnull() ) {   // should always be the case

    // Kinematics 
    ele_p_   = ele->p();
    ele_pt_  = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();

    // Momentum
    core_shFracHits_ = ele->shFracInnerHits();
    if (largeNtuple) {
      p4kind_  = ele->candidateP4Kind();
      ele_p_atvtx_  = sqrt(ele->trackMomentumAtVtx().mag2());
      ele_p_atcalo_ = sqrt(ele->trackMomentumAtCalo().mag2());
    }

    // Fiducial flags 
    fiducial_isEB_ = ele->isEB();
    fiducial_isEE_ = ele->isEE();
    fiducial_isEBEEGap_ = ele->isEBEEGap();

    // Charge 
    if (largeNtuple) {
      chPix_ = ele->scPixCharge();
      chGCP_ = ele->isGsfCtfScPixChargeConsistent();
      chGP_  = ele->isGsfScPixChargeConsistent();
      chGC_  = ele->isGsfCtfChargeConsistent();    
    }

    // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
    if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
    if ( mva_id > -666 )     { ele_mva_id_ = mva_id; }

    if (largeNtuple) { 
      if ( ele_conv_vtx_fit_prob > -666. ) ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; 
    }

    // ElectronID variables
    lowptgsfeleid::Features features;
    features.set(ele,rho);
    auto vfeatures = features.get();
    //@@ ORDER IS IMPORTANT!
    size_t idx = 0;
    eid_rho_ = vfeatures[idx++];
    eid_ele_pt_ = vfeatures[idx++];
    eid_sc_eta_ = vfeatures[idx++];
    eid_shape_full5x5_sigmaIetaIeta_ = vfeatures[idx++];
    eid_shape_full5x5_sigmaIphiIphi_ = vfeatures[idx++];
    eid_shape_full5x5_circularity_ = vfeatures[idx++];
    eid_shape_full5x5_r9_ = vfeatures[idx++];
    eid_sc_etaWidth_ = vfeatures[idx++];
    eid_sc_phiWidth_ = vfeatures[idx++];
    eid_shape_full5x5_HoverE_ = vfeatures[idx++];
    eid_trk_nhits_ = vfeatures[idx++];
    eid_trk_chi2red_ = vfeatures[idx++];
    eid_gsf_chi2red_ = vfeatures[idx++];
    eid_brem_frac_ = vfeatures[idx++];
    eid_gsf_nhits_ = vfeatures[idx++];
    eid_match_SC_EoverP_ = vfeatures[idx++];
    eid_match_eclu_EoverP_ = vfeatures[idx++];
    eid_match_SC_dEta_ = vfeatures[idx++];
    eid_match_SC_dPhi_ = vfeatures[idx++];
    eid_match_seed_dEta_ = vfeatures[idx++];
    eid_sc_E_ = vfeatures[idx++];
    eid_trk_p_ = vfeatures[idx++];

    // Isolation
    ele_sumPhotonEt_        = ele->pfIsolationVariables().sumPhotonEt;
    ele_sumChargedHadronPt_ = ele->pfIsolationVariables().sumChargedHadronPt;
    ele_sumNeutralHadronEt_ = ele->pfIsolationVariables().sumNeutralHadronEt;

    // Further track-Cluster matching 
    if (largeNtuple) {
      match_seed_EoverP_    = ele->eSeedClusterOverP();
      match_seed_EoverPout_ = ele->eSeedClusterOverPout();
      match_seed_dPhi_      = ele->deltaPhiSeedClusterTrackAtCalo();
      match_seed_dEta_vtx_  = ele->deltaEtaSeedClusterTrackAtVtx();
      match_eclu_EoverPout_ = ele->eEleClusterOverPout();
      match_eclu_dEta_ = ele->deltaEtaEleClusterTrackAtCalo();
      match_eclu_dPhi_ = ele->deltaPhiEleClusterTrackAtCalo();
    }

    // Further full 5x5 shower shape 
    if (largeNtuple) {
      shape_full5x5_e1x5_     = ele->full5x5_e1x5();
      shape_full5x5_e2x5Max_  = ele->full5x5_e2x5Max();
      shape_full5x5_e5x5_     = ele->full5x5_e5x5();
      shape_full5x5_eLeft_   = ele->full5x5_eLeft();
      shape_full5x5_eRight_  = ele->full5x5_eRight();
      shape_full5x5_eTop_    = ele->full5x5_eTop();
      shape_full5x5_eBottom_ = ele->full5x5_eBottom();
      shape_full5x5_HoverEBc_ = ele->full5x5_hcalOverEcalBc();
      shape_full5x5_hcalDepth1_    = ele->full5x5_hcalDepth1OverEcal();
      shape_full5x5_hcalDepth2_    = ele->full5x5_hcalDepth2OverEcal();
      shape_full5x5_hcalDepth1Bc_  = ele->full5x5_hcalDepth1OverEcalBc();
      shape_full5x5_hcalDepth2Bc_  = ele->full5x5_hcalDepth2OverEcalBc();
    }

    // Brem fractions and classification 
    brem_fracTrk_ = ele->trackFbrem();
    brem_fracSC_  = ele->superClusterFbrem();
    brem_N_ = ele->numberOfBrems();
  }    
}

void IDSlimNtuple::fill_supercluster(const reco::GsfElectronPtr ele, noZS::EcalClusterLazyTools *ecalTools_ ) {

  if ( ele.isNull() ) { return; }

  if ( ele->superCluster().isNull() ) { return; }
  const reco::SuperClusterRef& sc = ele->superCluster();

  int clusNum=0;
  float maxEne=-1;
  for(auto& cluster : sc->clusters()) {
    if (clusNum<NCLUS_MAX) {
      float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) + pow(cluster->z(),2) );
      sc_cluster_et_[clusNum]    = clusterEt;
      sc_cluster_E_[clusNum]     = cluster->energy();
      sc_cluster_eta_[clusNum]   = cluster->eta();
      sc_cluster_phi_[clusNum]   = cluster->phi();
      sc_cluster_nxtal_[clusNum] = cluster->size(); 
      sc_cluster_e1x3_[clusNum] = ecalTools_->e1x3(*cluster);
      sc_cluster_e1x5_[clusNum] = ecalTools_->e1x5(*cluster);
      sc_cluster_e2x2_[clusNum] = ecalTools_->e2x2(*cluster);
      sc_cluster_e3x3_[clusNum] = ecalTools_->e3x3(*cluster);
      sc_cluster_e5x5_[clusNum] = ecalTools_->e5x5(*cluster);
      sc_cluster_eMax_[clusNum] = ecalTools_->eMax(*cluster);
      sc_cluster_e2nd_[clusNum] = ecalTools_->e2nd(*cluster);
      sc_cluster_e2x5Right_[clusNum]  = ecalTools_->e2x5Right(*cluster);
      sc_cluster_e2x5Left_[clusNum]   = ecalTools_->e2x5Left(*cluster);
      sc_cluster_e2x5Top_[clusNum]    = ecalTools_->e2x5Top(*cluster);
      sc_cluster_e2x5Bottom_[clusNum] = ecalTools_->e2x5Bottom(*cluster);
      sc_cluster_eRight_[clusNum]  = ecalTools_->eRight(*cluster);
      sc_cluster_eLeft_[clusNum]   = ecalTools_->eLeft(*cluster);
      sc_cluster_eTop_[clusNum]    = ecalTools_->eTop(*cluster);
      sc_cluster_eBottom_[clusNum] = ecalTools_->eBottom(*cluster);
      sc_cluster_eMaxOver2x2_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e2x2(*cluster);
      sc_cluster_eMaxOver3x3_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e3x3(*cluster);
      sc_cluster_eMaxOver1x3_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e1x3(*cluster);
      if (cluster->energy() > maxEne) maxEne=cluster->energy();
      clusNum++;
    }
  }

  if (largeNtuple) {
    sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );
  }

  sc_Nclus_ = sc->clustersSize();
  float seedEne = sc->seed()->energy();
  if ( fabs(seedEne-maxEne)<0.001 ) sc_goodSeed_ = true;
}

// FC new method 
void IDSlimNtuple::fill_supercluster_miniAOD(const reco::GsfElectronPtr ele ) {

  // initialization in case of patological events
  sc_clus1_E_      = -999.;
  sc_clus1_E_ov_p_ = -999.;
  sc_clus1_E_ov_E_ = -999.;      
  sc_clus1_eta_    = -999.;
  sc_clus1_phi_    = -999.;
  sc_clus1_nxtal_  = -999;
  sc_clus1_deta_   = -999.;
  sc_clus1_dphi_   = -999.;
  sc_clus2_E_      = -999.;
  sc_clus2_E_ov_p_ = -999.;
  sc_clus2_E_ov_E_ = -999.;      
  sc_clus2_eta_    = -999.;
  sc_clus2_phi_    = -999.;
  sc_clus2_nxtal_  = -999;
  sc_clus2_deta_   = -999.;
  sc_clus2_dphi_   = -999.;
  sc_clus3_E_      = -999.;
  sc_clus3_E_ov_p_ = -999.;
  sc_clus3_E_ov_E_ = -999.;      
  sc_clus3_eta_    = -999.;
  sc_clus3_phi_    = -999.;
  sc_clus3_nxtal_  = -999;
  sc_clus3_deta_   = -999.;
  sc_clus3_dphi_   = -999.;
  sc_Nclus_deta01_ = -999;
  sc_Nclus_deta02_ = -999;   
  sc_Nclus_deta03_ = -999;  
  sc_Nclus_        = -999;
  sc_goodSeed_     = false;
  if (largeNtuple) {  
    sc_E_ps_  = -999.; 
    sc_E_ps1_ = -999.;
    sc_E_ps2_ = -999.;
    sc_Et_ = -999.;
    gsf_x_ = -999.;  
    gsf_y_ = -999.;  
    gsf_z_ = -999.;  
    sc_clus1_et_  = -999.;
    sc_clus2_et_  = -999.;
    sc_clus3_et_  = -999.;
  }

  // Analysis
  if ( ele.isNull() ) { return; }
  
  if ( ele->superCluster().isNull() ) { return; }
  const reco::SuperClusterRef& sc = ele->superCluster();

  reco::GsfTrackPtr kfTrackRef = edm::refToPtr(ele->gsfTrack());
  if (! validPtr(kfTrackRef) ) return;

  // Propagate 'electron' to ECAL surface
  double mass_=0.000511*0.000511; // ele mass 

  float p2=0;
  float px=0;
  float py=0;
  float pz=0;
  float vx=0;
  float vy=0;
  float vz=0;
  if ( kfTrackRef->extra().isAvailable() && kfTrackRef->extra().isNonnull() ) {
    p2=kfTrackRef->outerMomentum().Mag2();
    px=kfTrackRef->outerMomentum().x();
    py=kfTrackRef->outerMomentum().y();
    pz=kfTrackRef->outerMomentum().z();
    vx=kfTrackRef->outerPosition().x();
    vy=kfTrackRef->outerPosition().y();
    vz=kfTrackRef->outerPosition().z();
  } else {
    p2=pow( kfTrackRef->p() ,2 );
    px=kfTrackRef->px();
    py=kfTrackRef->py();
    pz=kfTrackRef->pz();
    vx=kfTrackRef->vx(); // must be in cm
    vy=kfTrackRef->vy();
    vz=kfTrackRef->vz();
  }

  float energy = sqrt(mass_ + p2);
  XYZTLorentzVector mom = XYZTLorentzVector(px,py,pz, energy);
  XYZTLorentzVector pos = XYZTLorentzVector(vx,vy,vz, 0.);
  float field_z=3.8;

  BaseParticlePropagator mypart(RawParticle(mom,pos), 0, 0, field_z);
  mypart.setCharge(kfTrackRef->charge());
  mypart.propagateToEcalEntrance(true); // true only first half loop , false more than one loop
  bool reach_ECAL=mypart.getSuccess(); // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps 

  // ECAL entry point for track
  // GlobalPoint ecal_pos(mypart.propagated().vertex().x(), mypart.propagated().vertex().y(), mypart.propagated().vertex().z());
  GlobalPoint ecal_pos(mypart.x(), mypart.y(), mypart.z());
  // Preshower limit
  //  bool below_ps = pow(ecal_pos.z(), 2.) > boundary_ * ecal_pos.perp2();
  // Iterate through ECAL clusters
  if (largeNtuple) {  
    gsf_x_=mypart.x();
    gsf_y_=mypart.y();
    gsf_z_=mypart.z();
  }
  
  int clusNum=0;
  float maxEne1=-1;
  float maxEne2=-1;
  float maxEne3=-1;
  int i1=-1;
  int i2=-1;
  int i3=-1;

  try{
    if(sc->clustersSize()>0 && sc->clustersBegin()!=sc->clustersEnd()){

      for(auto& cluster : sc->clusters()) {
	if (cluster->energy() > maxEne1){
	  maxEne1=cluster->energy();
	  i1=clusNum;
	}
	clusNum++;
      }
      
      if(sc->clustersSize()>1){
	clusNum=0;
	for(auto& cluster : sc->clusters()) {
	  if (clusNum!=i1) {
	    if (cluster->energy() > maxEne2){
	      maxEne2=cluster->energy();
	      i2=clusNum;
	    }
	  }
	  clusNum++;
	}
      }

      if(sc->clustersSize()>2){
	clusNum=0;
	for(auto& cluster : sc->clusters()) {
	  if (clusNum!=i1 && clusNum!=i2) {
	    if (cluster->energy() > maxEne3){
	      maxEne3=cluster->energy();
	      i3=clusNum;
	    }
	  }
	  clusNum++;
	}
      }
    }
  } catch(...){
    // std::cout<<"exception caught clusNum="<<clusNum<<" clus size"<<sc->clustersSize()<<" energy="<< sc->energy()<<std::endl;
  }

  // trovati i 3 cluster piu` energetici 
  // riempio i primi 3 cluster in E 
  // i1 i2 i3 
  clusNum=0;

  try{
    if(sc->clustersSize()>0&& sc->clustersBegin()!=sc->clustersEnd()){

      for(auto& cluster : sc->clusters()) {

	double pi_=3.1415926535;
	float deta = std::abs(ecal_pos.eta()-cluster->eta()) ;
	float dphi = std::abs(ecal_pos.phi()-cluster->phi());
	if (dphi > pi_)  dphi -= 2 * pi_;
	if(ecal_pos.phi()-cluster->phi()<0) dphi=-dphi;
	if(ecal_pos.eta()-cluster->eta()<0) deta=-deta;

	float elePmode = ele->gsfTrack()->pMode();
	float eleScEne = sc->energy();
	
	if (clusNum==i1) {
	  float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
	  sc_clus1_E_     = cluster->energy();
	  if( elePmode>0 ) sc_clus1_E_ov_p_ = cluster->energy()/elePmode;
	  if( eleScEne>0)  sc_clus1_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus1_eta_   = cluster->eta();
	  sc_clus1_phi_   = cluster->phi();
	  sc_clus1_nxtal_ =(int) cluster->size();
	  if(reach_ECAL>0){
	    sc_clus1_deta_ = deta;
	    sc_clus1_dphi_ = dphi;
	  }
	  if (largeNtuple) {
	    sc_clus1_et_  = clusterEt;
	  }

	} else if (clusNum==i2){
	  float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
	  sc_clus2_E_     = cluster->energy();
	  if( elePmode>0 ) sc_clus2_E_ov_p_ = cluster->energy()/elePmode;
	  if( eleScEne>0)  sc_clus2_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus2_eta_   = cluster->eta();
	  sc_clus2_phi_   = cluster->phi();
	  sc_clus2_nxtal_ = (int) cluster->size();
	  if(reach_ECAL>0){
	    sc_clus2_deta_ = deta;
	    sc_clus2_dphi_ = dphi;
	  }
	  if (largeNtuple) {
	    sc_clus2_et_  = clusterEt;
	  }
	  
	} else if (clusNum==i3){
	  float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
	  sc_clus3_E_     = cluster->energy();
	  if( elePmode>0 ) sc_clus3_E_ov_p_ = cluster->energy()/elePmode;
	  if( eleScEne>0)  sc_clus3_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus3_eta_   = cluster->eta();
	  sc_clus3_phi_   = cluster->phi();
	  sc_clus3_nxtal_ = (int) cluster->size();
	  if(reach_ECAL>0){
	    sc_clus3_deta_ = deta;
	    sc_clus3_dphi_ = dphi;
	  }
	  if (largeNtuple) {     
	    sc_clus3_et_  = clusterEt;
	  }
	}
	clusNum++;
      }
    }
  }catch(...){
    //
  }

  if (largeNtuple) {
    sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );
    sc_E_ps_=sc->preshowerEnergy();
    sc_E_ps1_=sc->preshowerEnergyPlane1();
    sc_E_ps2_=sc->preshowerEnergyPlane2();
  }
  sc_Nclus_ = sc->clustersSize();

  float seedEne = sc->seed()->energy();
  if ( fabs(seedEne-maxEne1)<0.001 ) sc_goodSeed_ = true;

  sc_Nclus_deta01_=0;
  sc_Nclus_deta02_=0;
  sc_Nclus_deta03_=0;
  try{
    if(sc->clustersSize()>0 && sc->clustersBegin()!=sc->clustersEnd()){
      for(auto& cluster : sc->clusters()) {
	float deta = std::abs(ecal_pos.eta()-cluster->eta()) ;
	if(deta<0.1) sc_Nclus_deta01_=sc_Nclus_deta01_+1; 
	if(deta<0.2) sc_Nclus_deta02_=sc_Nclus_deta02_+1; 
	if(deta<0.3) sc_Nclus_deta03_=sc_Nclus_deta03_+1; 
      }
    }
  }catch(...){
    //    std::cout<<"caught an exception"<<std::endl;
  }
}

// end FC

template < typename T> 
bool IDSlimNtuple::validPtr(edm::Ptr<T>& ptr){
  return (ptr.isNonnull() && ptr.isAvailable() );
}
