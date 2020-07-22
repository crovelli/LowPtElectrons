
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronFeatures.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::link_tree( TTree *tree ) {

  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");
  
  tree->Branch("weight", &weight_, "weight/f");
  tree->Branch("prescale", &prescale_, "prescale/f");
  tree->Branch("rho", &rho_, "rho/f");

  tree->Branch("is_aod", &is_aod_, "is_aod/i");
  tree->Branch("is_mc", &is_mc_, "is_mc/i");

  tree->Branch("tag_pt", &tag_pt_, "tag_pt/f");
  tree->Branch("tag_eta", &tag_eta_, "tag_eta/f");

  tree->Branch("is_e", &is_e_, "is_e/O");
  //tree->Branch("is_e_not_matched", &is_e_not_matched_, "is_e_not_matched/O");
  //tree->Branch("is_other", &is_other_, "is_other/O");
  tree->Branch("is_egamma", &is_egamma_, "is_egamma/O");

  tree->Branch("has_trk", &has_trk_, "has_trk/O");
  tree->Branch("has_seed", &has_seed_, "has_seed/O");
  tree->Branch("has_gsf", &has_gsf_, "has_gsf/O");
  tree->Branch("has_pfgsf", &has_pfgsf_, "has_pfgsf/O");
  tree->Branch("has_ele", &has_ele_, "has_ele/O");

  tree->Branch("trk_dr", &trk_dr_, "trk_dr/f");
  tree->Branch("gsf_dr", &gsf_dr_, "gsf_dr/f");
  tree->Branch("pfgsf_dr", &pfgsf_dr_, "pfgsf_dr/f");
  tree->Branch("ele_dr", &ele_dr_, "ele_dr/f");
  
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_e", &gen_e_, "gen_e/f");
  tree->Branch("gen_p", &gen_p_, "gen_p/f");
  tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");
  tree->Branch("gen_pdgid", &gen_pdgid_, "gen_pdgid/I");
  tree->Branch("gen_mom_pdgid", &gen_mom_pdgid_, "gen_mom_pdgid/I");
  tree->Branch("gen_gran_pdgid", &gen_gran_pdgid_, "gen_gran_pdgid/I");

  tree->Branch("trk_pt", &trk_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &trk_eta_, "trk_eta/f");
  tree->Branch("trk_phi", &trk_phi_, "trk_phi/f");
  tree->Branch("trk_p", &trk_p_, "trk_p/f");
  tree->Branch("trk_charge", &trk_charge_, "trk_charge/I");
  tree->Branch("trk_inp", &trk_inp_, "trk_inp/f");
  tree->Branch("trk_outp", &trk_outp_, "trk_outp/f");
  tree->Branch("trk_dpt", &trk_dpt_, "trk_dpt/f");

  tree->Branch("pdg_id", &pdg_id_, "pdg_id/I");

  tree->Branch("trk_nhits", &trk_nhits_, "trk_nhits/I");
  tree->Branch("trk_missing_inner_hits", &trk_missing_inner_hits_, "trk_missing_inner_hits/I"); 
  tree->Branch("trk_high_purity", &trk_high_purity_, "trk_high_purity/I");
  tree->Branch("trk_chi2red", &trk_chi2red_, "trk_chi2red/f");

  tree->Branch("trk_dxy", &trk_dxy_, "trk_dxy/f");
  tree->Branch("trk_dxy_err", &trk_dxy_err_, "trk_dxy_err/f");
  tree->Branch("trk_dz", &trk_dz_, "trk_dz/f");
  tree->Branch("trk_dz_err", &trk_dz_err_, "trk_dz_err/f");
  
  //tree->Branch("preid_unbiased", &preid_unbiased_, "preid_unbiased/f");
  //tree->Branch("preid_ptbiased", &preid_ptbiased_, "preid_ptbiased/f");
  tree->Branch("preid_bdtout1", &preid_unbiased_, "preid_bdtout1/f");
  tree->Branch("preid_bdtout2", &preid_ptbiased_, "preid_bdtout2/f");

  tree->Branch("seed_trk_driven", &seed_trk_driven_, "seed_trk_driven/O");
  tree->Branch("seed_ecal_driven", &seed_ecal_driven_, "seed_ecal_driven/O");

  tree->Branch("gsf_bdtout1", &seed_unbiased_, "gsf_bdtout1/f");
  tree->Branch("gsf_bdtout2", &seed_ptbiased_, "gsf_bdtout2/f");

  tree->Branch("gsf_pt", &gsf_pt_, "gsf_pt/f");
  tree->Branch("gsf_eta", &gsf_eta_, "gsf_eta/f");
  tree->Branch("gsf_phi", &gsf_phi_, "gsf_phi/f");
  tree->Branch("gsf_p", &gsf_p_, "gsf_p/f");
  tree->Branch("gsf_charge", &gsf_charge_, "gsf_charge/I");
  tree->Branch("gsf_inp", &gsf_inp_, "gsf_inp/f");
  tree->Branch("gsf_outp", &gsf_outp_, "gsf_outp/f");
  tree->Branch("gsf_dpt", &gsf_dpt_, "gsf_dpt/f");

  tree->Branch("gsf_mode_pt", &gsf_mode_pt_, "gsf_mode_pt/f");
  tree->Branch("gsf_mode_eta", &gsf_mode_eta_, "gsf_mode_eta/f");
  tree->Branch("gsf_mode_phi", &gsf_mode_phi_, "gsf_mode_phi/f");
  tree->Branch("gsf_mode_p", &gsf_mode_p_, "gsf_mode_p/f");

  tree->Branch("gsf_nhits",&gsf_nhits_, "gsf_nhits/I");
  tree->Branch("gsf_missing_inner_hits", &gsf_missing_inner_hits_, "gsf_missing_inner_hits/I");
  tree->Branch("gsf_chi2red", &gsf_chi2red_, "gsf_chi2red/f"); 
  
  tree->Branch("gsf_dxy",  &gsf_dxy_, "gsf_dxy/f");
  tree->Branch("gsf_dxy_err",&gsf_dxy_err_, "gsf_dxy_err/f");
  tree->Branch("gsf_dz",  &gsf_dz_, "gsf_dz/f");
  tree->Branch("gsf_dz_err",&gsf_dz_err_, "gsf_dz_err/f");

  //tree->Branch("gsf_ntangents", &gsf_ntangents_, "gsf_ntangents/I");
  //tree->Branch("gsf_hit_dpt", gsf_hit_dpt_, "gsf_hit_dpt[gsf_ntangents]/f");
  //tree->Branch("gsf_hit_dpt_unc", gsf_hit_dpt_unc_, "gsf_hit_dpt_unc[gsf_ntangents]/f");
  //tree->Branch("gsf_extapolated_eta", &gsf_extapolated_eta_);
  //tree->Branch("gsf_extapolated_phi", &gsf_extapolated_phi_);

  tree->Branch("pfgsf_pt", &pfgsf_pt_, "pfgsf_pt/f");
  tree->Branch("pfgsf_eta", &pfgsf_eta_, "pfgsf_eta/f");
  tree->Branch("pfgsf_phi", &pfgsf_phi_, "pfgsf_phi/f");
  tree->Branch("pfgsf_p", &pfgsf_p_, "pfgsf_p/f");
  tree->Branch("pfgsf_charge", &pfgsf_charge_, "pfgsf_charge/I");
  tree->Branch("pfgsf_inp", &pfgsf_inp_, "pfgsf_inp/f");
  tree->Branch("pfgsf_outp", &pfgsf_outp_, "pfgsf_outp/f");
  tree->Branch("pfgsf_dpt", &pfgsf_dpt_, "pfgsf_dpt/f");

  tree->Branch("pfgsf_mode_pt", &pfgsf_mode_pt_, "pfgsf_mode_pt/f");
  tree->Branch("pfgsf_mode_eta", &pfgsf_mode_eta_, "pfgsf_mode_eta/f");
  tree->Branch("pfgsf_mode_phi", &pfgsf_mode_phi_, "pfgsf_mode_phi/f");
  tree->Branch("pfgsf_mode_p", &pfgsf_mode_p_, "pfgsf_mode_p/f");

  tree->Branch("pfgsf_nhits",&pfgsf_nhits_, "pfgsf_nhits/I");
  tree->Branch("pfgsf_missing_inner_hits", &pfgsf_missing_inner_hits_, "pfgsf_missing_inner_hits/I");
  tree->Branch("pfgsf_chi2red", &pfgsf_chi2red_, "pfgsf_chi2red/f"); 
  
  tree->Branch("pfgsf_dxy",  &pfgsf_dxy_, "pfgsf_dxy/f");
  tree->Branch("pfgsf_dxy_err",&pfgsf_dxy_err_, "pfgsf_dxy_err/f");
  tree->Branch("pfgsf_dz",  &pfgsf_dz_, "pfgsf_dz/f");
  tree->Branch("pfgsf_dz_err",&pfgsf_dz_err_, "pfgsf_dz_err/f");

  //tree->Branch("pfgsf_ntangents", &pfgsf_ntangents_, "pfgsf_ntangents/I");
  //tree->Branch("pfgsf_hit_dpt", pfgsf_hit_dpt_, "pfgsf_hit_dpt[pfgsf_ntangents]/f");
  //tree->Branch("pfgsf_hit_dpt_unc", pfgsf_hit_dpt_unc_, "pfgsf_hit_dpt_unc[pfgsf_ntangents]/f");
  //tree->Branch("pfgsf_extapolated_eta", &pfgsf_extapolated_eta_);
  //tree->Branch("pfgsf_extapolated_phi", &pfgsf_extapolated_phi_);
  
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");

  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_value_retrained", &ele_mva_value_retrained_, "ele_mva_value_retrained/f");
  tree->Branch("ele_conv_vtx_fit_prob", &ele_conv_vtx_fit_prob_, "ele_conv_vtx_fit_prob/f");
  tree->Branch("ele_mva_value_depth10", &ele_mva_value_depth10_, "ele_mva_value_depth10/f");
  tree->Branch("ele_mva_value_depth15", &ele_mva_value_depth15_, "ele_mva_value_depth15/f");

  tree->Branch("eid_rho", &eid_rho_, "eid_rho/f");
  tree->Branch("eid_ele_pt", &eid_ele_pt_, "eid_ele_pt/f");

  tree->Branch("eid_trk_p", &eid_trk_p_, "eid_trk_p/f");
  tree->Branch("eid_trk_nhits", &eid_trk_nhits_, "eid_trk_nhits/f");
  tree->Branch("eid_trk_chi2red", &eid_trk_chi2red_, "eid_trk_chi2red/f");

  tree->Branch("eid_gsf_nhits", &eid_gsf_nhits_, "eid_gsf_nhits/f");
  tree->Branch("eid_gsf_chi2red", &eid_gsf_chi2red_, "eid_gsf_chi2red/f");

  tree->Branch("eid_sc_E", &eid_sc_E_, "eid_sc_E/f");
  tree->Branch("eid_sc_eta", &eid_sc_eta_, "eid_sc_eta/f");
  tree->Branch("eid_sc_etaWidth", &eid_sc_etaWidth_, "eid_sc_etaWidth/f");
  tree->Branch("eid_sc_phiWidth", &eid_sc_phiWidth_, "eid_sc_phiWidth/f");

  tree->Branch("eid_match_seed_dEta", &eid_match_seed_dEta_, "eid_match_seed_dEta/f");
  tree->Branch("eid_match_eclu_EoverP", &eid_match_eclu_EoverP_, "eid_match_eclu_EoverP/f");
  tree->Branch("eid_match_SC_EoverP", &eid_match_SC_EoverP_, "eid_match_SC_EoverP/f");
  tree->Branch("eid_match_SC_dEta", &eid_match_SC_dEta_, "eid_match_SC_dEta/f");
  tree->Branch("eid_match_SC_dPhi", &eid_match_SC_dPhi_, "eid_match_SC_dPhi/f");

  tree->Branch("eid_shape_full5x5_sigmaIetaIeta", &eid_shape_full5x5_sigmaIetaIeta_, "eid_shape_full5x5_sigmaIetaIeta/f");
  tree->Branch("eid_shape_full5x5_sigmaIphiIphi", &eid_shape_full5x5_sigmaIphiIphi_, "eid_shape_full5x5_sigmaIphiIphi/f");
  tree->Branch("eid_shape_full5x5_HoverE", &eid_shape_full5x5_HoverE_, "eid_shape_full5x5_HoverE/f");
  tree->Branch("eid_shape_full5x5_r9", &eid_shape_full5x5_r9_, "eid_shape_full5x5_r9/f");
  tree->Branch("eid_shape_full5x5_circularity", &eid_shape_full5x5_circularity_, "eid_shape_full5x5_circularity/f");

  tree->Branch("eid_brem_frac", &eid_brem_frac_, "eid_brem_frac/f");

  tree->Branch("image_gsf_ref_eta", &image_gsf_ref_eta_, "image_gsf_ref_eta/f");
  tree->Branch("image_gsf_ref_phi", &image_gsf_ref_phi_, "image_gsf_ref_phi/f");
  tree->Branch("image_gsf_ref_R", &image_gsf_ref_R_, "image_gsf_ref_R/f");
  tree->Branch("image_gsf_ref_p", &image_gsf_ref_p_, "image_gsf_ref_p/f");
  tree->Branch("image_gsf_ref_pt", &image_gsf_ref_pt_, "image_gsf_ref_pt/f");
  
  tree->Branch("image_gen_inner_eta", &image_gen_inner_eta_, "image_gen_inner_eta/f");
  tree->Branch("image_gen_inner_phi", &image_gen_inner_phi_, "image_gen_inner_phi/f");
  tree->Branch("image_gen_inner_R", &image_gen_inner_R_, "image_gen_inner_R/f");
  tree->Branch("image_gen_inner_p", &image_gen_inner_p_, "image_gen_inner_p/f");
  tree->Branch("image_gen_inner_pt", &image_gen_inner_pt_, "image_gen_inner_pt/f");

  tree->Branch("image_gsf_inner_eta", &image_gsf_inner_eta_, "image_gsf_inner_eta/f");
  tree->Branch("image_gsf_inner_phi", &image_gsf_inner_phi_, "image_gsf_inner_phi/f");
  tree->Branch("image_gsf_inner_R", &image_gsf_inner_R_, "image_gsf_inner_R/f");
  tree->Branch("image_gsf_inner_p", &image_gsf_inner_p_, "image_gsf_inner_p/f");
  tree->Branch("image_gsf_inner_pt", &image_gsf_inner_pt_, "image_gsf_inner_pt/f");
  tree->Branch("image_gsf_charge", &image_gsf_charge_, "image_gsf_charge/I");

  tree->Branch("image_gsf_proj_eta", &image_gsf_proj_eta_, "image_gsf_proj_eta/f");
  tree->Branch("image_gsf_proj_phi", &image_gsf_proj_phi_, "image_gsf_proj_phi/f");
  tree->Branch("image_gsf_proj_R", &image_gsf_proj_R_, "image_gsf_proj_R/f");
  tree->Branch("image_gsf_proj_p", &image_gsf_proj_p_, "image_gsf_proj_p/f");

  tree->Branch("image_gsf_atcalo_eta", &image_gsf_atcalo_eta_, "image_gsf_atcalo_eta/f");
  tree->Branch("image_gsf_atcalo_phi", &image_gsf_atcalo_phi_, "image_gsf_atcalo_phi/f");
  tree->Branch("image_gsf_atcalo_R", &image_gsf_atcalo_R_, "image_gsf_atcalo_R/f");
  tree->Branch("image_gsf_atcalo_p", &image_gsf_atcalo_p_, "image_gsf_atcalo_p/f");
  
  tree->Branch("image_clu_n", &image_clu_n_, "image_clu_n/I");
  tree->Branch("image_clu_eta", &image_clu_eta_, "image_clu_eta[image_clu_n]/f");
  tree->Branch("image_clu_phi", &image_clu_phi_, "image_clu_phi[image_clu_n]/f");
  tree->Branch("image_clu_e", &image_clu_e_, "image_clu_e[image_clu_n]/f");
  tree->Branch("image_clu_nhit", &image_clu_nhit_, "image_clu_nhit[image_clu_n]/I");

  tree->Branch("image_pf_n", &image_pf_n_, "image_pf_n/I");
  tree->Branch("image_pf_eta", &image_pf_eta_, "image_pf_eta[image_pf_n]/f");
  tree->Branch("image_pf_phi", &image_pf_phi_, "image_pf_phi[image_pf_n]/f");
  tree->Branch("image_pf_p", &image_pf_p_, "image_pf_p[image_pf_n]/f");
  tree->Branch("image_pf_pdgid", &image_pf_pdgid_, "image_pf_pdgid[image_pf_n]/I");
  tree->Branch("image_pf_matched", &image_pf_matched_, "image_pf_matched[image_pf_n]/I");
  tree->Branch("image_pf_lost", &image_pf_lost_, "image_pf_lost[image_pf_n]/I");

}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_evt( const edm::EventID& id ) {
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_gen( const reco::GenParticlePtr genp ) {
  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_e_ = genp->energy();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = 0;
  gen_mom_pdgid_ = 0;
  gen_gran_pdgid_ = 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_gen( const reco::CandidatePtr genp ) {
  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_e_ = genp->energy();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = 0;
  gen_mom_pdgid_ = 0;
  gen_gran_pdgid_ = 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_gen( const pat::PackedGenParticleRef genp ) {
  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_e_ = genp->energy();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = 0;
  gen_mom_pdgid_ = 0;
  gen_gran_pdgid_ = 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_trk( const reco::TrackPtr& trk,
			 const reco::BeamSpot& spot ) {
  
  if ( trk.isNonnull() ) {
    // kine
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
    trk_p_ = trk->p();
    trk_charge_ = trk->charge();
    if ( trk->extra().isAvailable() && trk->extra().isNonnull() ) {
      trk_inp_ = sqrt( trk->innerMomentum().mag2() );
      trk_outp_ = sqrt( trk->outerMomentum().mag2() );
      trk_dpt_ = ( trk_inp_ > 0. ) ? fabs( trk_outp_ - trk_inp_ ) / trk_inp_ : 0.; //@@ redundant?
    }
    // quality
    trk_nhits_ = trk->found();
    trk_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    trk_high_purity_ = trk->quality( reco::TrackBase::qualityByName("highPurity") );
    trk_chi2red_ = trk->normalizedChi2();
    // displ
    trk_dxy_ = trk->dxy(spot);
    trk_dxy_err_ = trk->dxyError();
    trk_dz_ = trk->dz(spot.position());
    trk_dz_err_ = trk->dzError();
  }
  
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_seed( bool seed_trk_driven, 
			  bool seed_ecal_driven ) {
  seed_trk_driven_ = seed_trk_driven;
  seed_ecal_driven_ = seed_ecal_driven;
}

/////////////////////////////////////////////////////////////////////////////////
//@@ to be deprecated
void IDNtuple::fill_seed( double seed_unbiased, 
			  double seed_ptbiased ) {
  std::cout << "TO BE DEPRECATED!" << std::endl;
  seed_unbiased_ = seed_unbiased;
  seed_ptbiased_ = seed_ptbiased;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_bdt( double seed_unbiased, 
			 double seed_ptbiased ) {
  seed_unbiased_ = seed_unbiased;
  seed_ptbiased_ = seed_ptbiased;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_preid( const reco::PreId& preid_ecal, 
			   const reco::PreId& preid_hcal, 
			   const reco::BeamSpot& spot, 
			   const double rho, 
                           noZS::EcalClusterLazyTools& ecal_tools ) {

  // MVA output
  preid_unbiased_ = preid_ecal.mva(0);
  preid_ptbiased_ = preid_ecal.mva(1);
  
//  // MVA decision
//  preid_unbiased_pass_ = preid_ecal.mvaSelected(0);
//  preid_ptbiased_pass_ = preid_ecal.mvaSelected(1);

  // Set seed variables
  auto vfeatures = lowptgsfeleseed::features( preid_ecal, preid_hcal, rho, spot, ecal_tools );
  
  //@@ ADD THESE PREID VARS TO THE NTUPLE???

//  //@@ ORDER IS IMPORTANT!
//  size_t idx = 0;
//  preid_trk_pt_ = vfeatures[idx++];
//  preid_trk_eta_ = vfeatures[idx++];
//  preid_trk_phi_ = vfeatures[idx++];
//  preid_trk_p_ = vfeatures[idx++];
//  preid_trk_nhits_ = vfeatures[idx++];
//  preid_trk_high_quality_ = vfeatures[idx++];
//  preid_trk_chi2red_ = vfeatures[idx++];
//  preid_rho_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_e_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_deta_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_dphi_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_e3x3_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_e5x5_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_covEtaEta_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_covEtaPhi_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_covPhiPhi_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_r9_ = vfeatures[idx++];
//  preid_ktf_ecal_cluster_circularity_ = vfeatures[idx++];
//  preid_ktf_hcal_cluster_e_ = vfeatures[idx++];
//  preid_ktf_hcal_cluster_deta_ = vfeatures[idx++];
//  preid_ktf_hcal_cluster_dphi_ = vfeatures[idx++];
//  preid_gsf_dpt_ = vfeatures[idx++];
//  preid_trk_gsf_chiratio_ = vfeatures[idx++];
//  preid_gsf_chi2red_ = vfeatures[idx++];
//  preid_trk_dxy_sig_ = vfeatures[idx++]; // must be last (not used by unbiased model)

}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_gsf( const reco::GsfTrackPtr gsf, 
			 const reco::BeamSpot& spot ) {

  if ( gsf.isNull() ) {
    //@@ Shouldn't happen, but do we just take dummy values...? 
  } else {

    // Kinematics
    gsf_pt_ = gsf->pt();
    gsf_eta_ = gsf->eta();
    gsf_phi_ = gsf->phi();
    gsf_p_ = gsf->p();
    gsf_charge_ = gsf->charge();
    if ( gsf->extra().isAvailable() && gsf->extra().isNonnull() ) {
      gsf_inp_ = sqrt(gsf->innerMomentum().mag2());
      gsf_outp_ = sqrt(gsf->outerMomentum().mag2());
      gsf_dpt_ = ( gsf_inp_ > 0. ) ? fabs( gsf_outp_ - gsf_inp_ ) / gsf_inp_ : 0.; //@@ redundant?
    }
    
    // Kinematics (MODE)
    gsf_mode_pt_ = gsf->ptMode();
    gsf_mode_eta_ = gsf->etaMode();
    gsf_mode_phi_ = gsf->phiMode();
    gsf_mode_p_ = gsf->pMode();

    // Quality
    gsf_nhits_ = gsf->found();
    gsf_missing_inner_hits_ = gsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    gsf_chi2red_ = gsf->normalizedChi2();

    // Displacement
    gsf_dxy_ = gsf->dxy(spot);
    gsf_dxy_err_ = gsf->dxyError();
    gsf_dz_ = gsf->dz(spot.position());
    gsf_dz_err_ = gsf->dzError();

    // Tangents (requires TrackExtra)
    //    const auto& extra = gsf->gsfExtra(); //@@ Collection does not exist?!
    //    if ( extra.isNonnull() ) {
    //      gsf_ntangents_ = (extra->tangentsSize() > NHITS_MAX) ? NHITS_MAX : extra->tangentsSize();
    //      for (int idx = 0; idx < gsf_ntangents_; idx++ ) {
    //	gsf_hit_dpt_[idx] = extra->tangents().at(idx).deltaP().value();
    //	gsf_hit_dpt_unc_[idx] = extra->tangents().at(idx).deltaP().error();
    //      }
    //    }
    
  } 
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_pfgsf( const reco::GsfTrackPtr pfgsf, 
			   const reco::BeamSpot& spot ) {

  if ( pfgsf.isNull() ) {
    //@@ Shouldn't happen, but do we just take dummy values...? 
  } else {

    // Kinematics
    pfgsf_pt_ = pfgsf->pt();
    pfgsf_eta_ = pfgsf->eta();
    pfgsf_phi_ = pfgsf->phi();
    pfgsf_p_ = pfgsf->p();
    pfgsf_charge_ = pfgsf->charge();
    if ( pfgsf->extra().isAvailable() && pfgsf->extra().isNonnull() ) {
      pfgsf_inp_ = sqrt(pfgsf->innerMomentum().mag2());
      pfgsf_outp_ = sqrt(pfgsf->outerMomentum().mag2());
      pfgsf_dpt_ = ( pfgsf_inp_ > 0. ) ? fabs( pfgsf_outp_ - pfgsf_inp_ ) / pfgsf_inp_ : 0.; //@@ redundant?
    }
    
    // Kinematics (MODE)
    pfgsf_mode_pt_ = pfgsf->ptMode();
    pfgsf_mode_eta_ = pfgsf->etaMode();
    pfgsf_mode_phi_ = pfgsf->phiMode();
    pfgsf_mode_p_ = pfgsf->pMode();

    // Quality
    pfgsf_nhits_ = pfgsf->found();
    pfgsf_missing_inner_hits_ = pfgsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    pfgsf_chi2red_ = pfgsf->normalizedChi2();

    // Displacement
    pfgsf_dxy_ = pfgsf->dxy(spot);
    pfgsf_dxy_err_ = pfgsf->dxyError();
    pfgsf_dz_ = pfgsf->dz(spot.position());
    pfgsf_dz_err_ = pfgsf->dzError();

    // Tangents (requires TrackExtra)
    //    const auto& extra = pfgsf->pfgsfExtra(); //@@ Collection does not exist?!
    //    if ( extra.isNonnull() ) {
    //      pfgsf_ntangents_ = (extra->tangentsSize() > NHITS_MAX) ? NHITS_MAX : extra->tangentsSize();
    //      for (int idx = 0; idx < pfgsf_ntangents_; idx++ ) {
    //	pfgsf_hit_dpt_[idx] = extra->tangents().at(idx).deltaP().value();
    //	pfgsf_hit_dpt_unc_[idx] = extra->tangents().at(idx).deltaP().error();
    //      }
    //    }
    
  } 
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_ele( const reco::GsfElectronPtr ele,
			 float mva_value,
			 float mva_value_retrained,
			 float mva_value_depth10,
			 float mva_value_depth15,
			 float ele_conv_vtx_fit_prob,
			 const double rho,
			 bool is_egamma,
			 float unbiased // Rome ...
			 ) {

  // Kinematics
  if ( is_egamma ) {
    ele_p_ = ele->p();
    ele_pt_ = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();
  } else {
    reco::GsfTrackRef gsf = ele->gsfTrack();
    if ( gsf.isNonnull() && gsf.isAvailable() ) {
      ele_p_ = gsf->p();
      ele_pt_ = gsf->pt();
      ele_eta_ = gsf->eta();
      ele_phi_ = gsf->phi();
    } else {
      std::cout << "[IDNtuple::fill_ele] ERROR: Null GsfTrackRef!" << std::endl;
    }
  }
  
  // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
  if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
  if ( mva_value_retrained > -666 ) { ele_mva_value_retrained_ = mva_value_retrained; }
  if ( ele_conv_vtx_fit_prob > -666. ) { ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; }
  if ( mva_value_depth10 > -666. ) { ele_mva_value_depth10_ = mva_value_depth10; }
  if ( mva_value_depth15 > -666. ) { ele_mva_value_depth15_ = mva_value_depth15; }

  // Set Electron variables

  //lowptgsfeleid::Features features;
  //features.set(ele,rho);
  //auto vfeatures = features.get();
  std::vector<float> vfeatures; // = getFeatures(ele,rho,unbiased); // new interface ...

  //@@ ORDER IS IMPORTANT!

  // My variables ...
//  size_t idx = 0;
//  eid_rho_ = vfeatures[idx++];
//  eid_ele_pt_ = vfeatures[idx++];
//  eid_sc_eta_ = vfeatures[idx++];
//  eid_shape_full5x5_sigmaIetaIeta_ = vfeatures[idx++];
//  eid_shape_full5x5_sigmaIphiIphi_ = vfeatures[idx++];
//  eid_shape_full5x5_circularity_ = vfeatures[idx++];
//  eid_shape_full5x5_r9_ = vfeatures[idx++];
//  eid_sc_etaWidth_ = vfeatures[idx++];
//  eid_sc_phiWidth_ = vfeatures[idx++];
//  eid_shape_full5x5_HoverE_ = vfeatures[idx++];
//  eid_trk_nhits_ = vfeatures[idx++];
//  eid_trk_chi2red_ = vfeatures[idx++];
//  eid_gsf_chi2red_ = vfeatures[idx++];
//  eid_brem_frac_ = vfeatures[idx++];
//  eid_gsf_nhits_ = vfeatures[idx++];
//  eid_match_SC_EoverP_ = vfeatures[idx++];
//  eid_match_eclu_EoverP_ = vfeatures[idx++];
//  eid_match_SC_dEta_ = vfeatures[idx++];
//  eid_match_SC_dPhi_ = vfeatures[idx++];
//  eid_match_seed_dEta_ = vfeatures[idx++];
//  eid_sc_E_ = vfeatures[idx++];
//  eid_trk_p_ = vfeatures[idx++];

  // Rome variables ...
  size_t idx = 0;
  eid_rho_ = vfeatures[idx++];
  eid_sc_eta_ = vfeatures[idx++];
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
  // following not set ...
  //eid_ele_pt_ = vfeatures[idx++];
  //eid_shape_full5x5_sigmaIetaIeta_ = vfeatures[idx++];
  //eid_shape_full5x5_sigmaIphiIphi_ = vfeatures[idx++];
  //eid_shape_full5x5_circularity_ = vfeatures[idx++];
  // following are new, not used ...
  //gsf_mode_p,core_shFracHits,gsf_bdtout1,gsf_dr,trk_dr,sc_Nclus,
  //sc_clus1_nxtal,sc_clus1_dphi,sc_clus2_dphi,sc_clus1_deta,
  //sc_clus2_deta,sc_clus1_E,sc_clus2_E,sc_clus1_E_ov_p,sc_clus2_E_ov_p  

}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_image( const float gsf_ref_eta, const float gsf_ref_phi, const float gsf_ref_R,
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
			   const std::vector<int>& pf_lost ) {
  if ( clu_eta.size() != clu_phi.size() ||
       clu_eta.size() != clu_e.size() ||
       clu_eta.size() != clu_nhit.size() ||
       pf_eta.size() != pf_phi.size() ||
       pf_eta.size() != pf_p.size() ||
       pf_eta.size() != pf_pdgid.size() ||
       pf_eta.size() != pf_matched.size() ||
       pf_eta.size() != pf_lost.size() ) { 
    throw std::runtime_error("Inconsistent vector sizes!"); 
  }
  
  image_gsf_ref_eta_  = gsf_ref_eta;
  image_gsf_ref_phi_  = gsf_ref_phi;
  image_gsf_ref_R_  = gsf_ref_R;
  image_gsf_ref_p_  = gsf_ref_p;
  image_gsf_ref_pt_ = gsf_ref_pt;

  image_gen_inner_eta_ = gen_inner_eta;
  image_gen_inner_phi_ = gen_inner_phi;
  image_gen_inner_R_ = gen_inner_R;
  image_gen_inner_p_ = gen_inner_p;
  image_gen_inner_pt_ = gen_inner_pt;

  image_gsf_inner_eta_ = gsf_inner_eta;
  image_gsf_inner_phi_ = gsf_inner_phi;
  image_gsf_inner_R_ = gsf_inner_R;
  image_gsf_inner_p_ = gsf_inner_p;
  image_gsf_inner_pt_ = gsf_inner_pt;
  image_gsf_charge_ = gsf_charge;

  image_gsf_proj_eta_ = gsf_proj_eta;
  image_gsf_proj_phi_ = gsf_proj_phi;
  image_gsf_proj_R_ = gsf_proj_R;
  image_gsf_proj_p_ = gsf_proj_p;

  image_gsf_atcalo_eta_ = gsf_atcalo_eta;
  image_gsf_atcalo_phi_ = gsf_atcalo_phi;
  image_gsf_atcalo_R_ = gsf_atcalo_R;
  image_gsf_atcalo_p_ = gsf_atcalo_p;

  image_clu_n_ = clu_eta.size() > ARRAY_SIZE ? ARRAY_SIZE : clu_eta.size();
  for ( size_t i = 0; i < image_clu_n_; ++i )  { 
    image_clu_eta_[i] = clu_eta[i];
    image_clu_phi_[i] = clu_phi[i];
    image_clu_e_[i] = clu_e[i];
    image_clu_nhit_[i] = clu_nhit[i];
  }

  image_pf_n_ = pf_eta.size() > ARRAY_SIZE ? ARRAY_SIZE : pf_eta.size();
  for ( size_t i = 0; i < image_pf_n_; ++i )  { 
    image_pf_eta_[i] = pf_eta[i];
    image_pf_phi_[i] = pf_phi[i];
    image_pf_p_[i] = pf_p[i];
    image_pf_pdgid_[i] = pf_pdgid[i];
    image_pf_matched_[i] = pf_matched[i];
    image_pf_lost_[i] = pf_lost[i];
  }

}
