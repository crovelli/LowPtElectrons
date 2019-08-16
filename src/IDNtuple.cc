#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDHeavyObjectCache.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::link_tree( TTree *tree ) {

  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");
  
  tree->Branch("weight", &weight_, "weight/f");
  tree->Branch("rho", &rho_, "rho/f");

  tree->Branch("is_aod", &is_aod_, "is_aod/i");
  tree->Branch("is_mc", &is_mc_, "is_mc/i");

  tree->Branch("is_e", &is_e_, "is_e/O");
  tree->Branch("is_e_not_matched", &is_e_not_matched_, "is_e_not_matched/O");
  tree->Branch("is_other", &is_other_, "is_other/O");
  tree->Branch("is_egamma", &is_egamma_, "is_egamma/O");

  tree->Branch("has_trk", &has_trk_, "has_trk/O");
  tree->Branch("has_seed", &has_seed_, "has_seed/O");
  tree->Branch("has_gsf", &has_gsf_, "has_gsf/O");
  tree->Branch("has_ele", &has_ele_, "has_ele/O");

  tree->Branch("trk_dr", &trk_dr_, "trk_dr/f");
  tree->Branch("gsf_dr", &gsf_dr_, "gsf_dr/f");
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
  
//  tree->Branch("preid_unbiased", &preid_unbiased_, "preid_unbiased/f");
//  tree->Branch("preid_ptbiased", &preid_ptbiased_, "preid_ptbiased/f");
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
  
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");

  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_id", &ele_mva_id_, "ele_mva_id/I");
  tree->Branch("ele_conv_vtx_fit_prob", &ele_conv_vtx_fit_prob_, "ele_conv_vtx_fit_prob/f");

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
  lowptgsfeleseed::Features features;
  features.set( preid_ecal, preid_hcal, rho, spot, ecal_tools );
  auto vfeatures = features.get();

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
void IDNtuple::fill_ele( const reco::GsfElectronPtr ele,
			 float mva_value,
			 int mva_id,
			 float ele_conv_vtx_fit_prob,
			 const double rho ) {

  // Kinematics
  ele_p_ = ele->p();
  ele_pt_ = ele->pt();
  ele_eta_ = ele->eta();
  ele_phi_ = ele->phi();
  
  // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
  if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
  if ( mva_id > -666 ) { ele_mva_id_ = mva_id; }
  if ( ele_conv_vtx_fit_prob > -666. ) { ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; }

  // Set Electron variables
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

}
