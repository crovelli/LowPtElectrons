#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDHeavyObjectCache.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::link_tree( TTree *tree ) {

  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");
  
  tree->Branch("rho", &rho_, "rho/f");
  
  tree->Branch("is_e", &is_e_, "is_e/O");
  tree->Branch("is_e_not_matched", &is_e_not_matched_, "is_e_not_matched/O");
  tree->Branch("is_other", &is_other_, "is_other/O");

  tree->Branch("has_pfgsf", &has_pfgsf_, "has_pfgsf/O");
  tree->Branch("has_pfele", &has_pfele_, "has_pfele/O");
  
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_e", &gen_e_, "gen_e/f");
  tree->Branch("gen_p", &gen_p_, "gen_p/f");
  tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");
  tree->Branch("gen_pdgid", &gen_pdgid_, "gen_pdgid/I");
  tree->Branch("gen_mom_pdgid", &gen_mom_pdgid_, "gen_mom_pdgid/I");
  tree->Branch("gen_gran_pdgid", &gen_gran_pdgid_, "gen_gran_pdgid/I");

  //tree->Branch("seed_bdt_unbiased", &seed_bdt_unbiased_, "seed_bdt_unbiased/f");
  //tree->Branch("seed_bdt_ptbiased", &seed_bdt_ptbiased_, "seed_bdt_ptbiased/f");
  tree->Branch("preid_bdtout1", &seed_bdt_unbiased_, "preid_bdtout1/f");
  tree->Branch("preid_bdtout2", &seed_bdt_ptbiased_, "preid_bdtout2/f");

  tree->Branch("trk_pt", &gsf_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &gsf_eta_, "trk_eta/f");

  tree->Branch("gsf_pt", &gsf_pt_, "gsf_pt/f");
  tree->Branch("gsf_eta", &gsf_eta_, "gsf_eta/f");
  tree->Branch("gsf_phi", &gsf_phi_, "gsf_phi/f");
  tree->Branch("gsf_p", &gsf_p_, "gsf_p/f");
  tree->Branch("gsf_charge", &gsf_charge_, "gsf_charge/I");
  //tree->Branch("gsf_inp", &gsf_inp_, "gsf_inp/f");
  //tree->Branch("gsf_outp", &gsf_outp_, "gsf_outp/f");
  //tree->Branch("gsf_dpt", &gsf_dpt_, "gsf_dpt/f");

  tree->Branch("mode_pt", &mode_pt_, "mode_pt/f");
  tree->Branch("mode_eta", &mode_eta_, "mode_eta/f");
  tree->Branch("mode_phi", &mode_phi_, "mode_phi/f");
  tree->Branch("mode_p", &mode_p_, "mode_p/f");

  tree->Branch("gsf_nhits",&gsf_nhits_, "gsf_nhits/I");
  tree->Branch("gsf_missing_inner_hits", &gsf_missing_inner_hits_, "gsf_missing_inner_hits/I");
  tree->Branch("gsf_chi2red", &gsf_chi2red_, "gsf_chi2red/f"); 
  
  tree->Branch("gsf_dxy",  &gsf_dxy_, "gsf_dxy/f");
  tree->Branch("gsf_dxy_err",&gsf_dxy_err_, "gsf_dxy_err/f");
  tree->Branch("gsf_dz",  &gsf_dz_, "gsf_dz/f");
  tree->Branch("gsf_dz_err",&gsf_dz_err_, "gsf_dz_err/f");

  //tree->Branch("gsf_ntangents", &gsf_ntangents_, "gsf_ntangents/i");
  //tree->Branch("gsf_hit_dpt", gsf_hit_dpt_, "gsf_hit_dpt[gsf_ntangents]/f");
  //tree->Branch("gsf_hit_dpt_unc", gsf_hit_dpt_unc_, "gsf_hit_dpt_unc[gsf_ntangents]/f");
  //tree->Branch("gsf_extapolated_eta", &gsf_extapolated_eta_);
  //tree->Branch("gsf_extapolated_phi", &gsf_extapolated_phi_);
  
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");

  tree->Branch("ele_mvaIdV2", &ele_mvaIdV2_, "ele_mvaIdV2/f");
  tree->Branch("ele_lowPtMva", &ele_lowPtMva_, "ele_lowPtMva/f");
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
void IDNtuple::fill_gen( const reco::GenParticleRef genp ) {
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
void IDNtuple::fill_seed( double seed_bdt_unbiased, double seed_bdt_ptbiased ) {
  seed_bdt_unbiased_ = seed_bdt_unbiased;
  seed_bdt_ptbiased_ = seed_bdt_ptbiased;
}

/////////////////////////////////////////////////////////////////////////////////
//
void IDNtuple::fill_gsf( const reco::GsfTrackRef gsf, const reco::BeamSpot& spot ) {

  if ( gsf.isNull() ) {
    //@@ Shouldn't happen, but do we just take dummy values...? 
  } else {

    // Kinematics
    gsf_pt_ = gsf->pt();
    gsf_eta_ = gsf->eta();
    gsf_phi_ = gsf->phi();
    gsf_p_ = gsf->p();
    gsf_charge_ = gsf->charge();
    //gsf_inp_ = sqrt(gsf->innerMomentum().mag2());
    //gsf_outp_ = sqrt(gsf->outerMomentum().mag2());
    //gsf_dpt_ = ( gsf_inp_ > 0. ) ? fabs( gsf_outp_ - gsf_inp_ ) / gsf_inp_ : 0.; //@@ redundant?

    // Kinematics (MODE)
    mode_pt_ = gsf->ptMode();
    mode_eta_ = gsf->etaMode();
    mode_phi_ = gsf->phiMode();
    mode_p_ = gsf->pMode();

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
void IDNtuple::fill_ele( const pat::ElectronRef ele,
			 float id_lowpt,
			 float id_v2,
			 float ele_conv_vtx_fit_prob,
			 const double rho ) {

  // Kinematics
  ele_p_ = ele->p();
  ele_pt_ = ele->pt();
  ele_eta_ = ele->eta();
  ele_phi_ = ele->phi();
  
  // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
  if ( id_lowpt > -666. ) { ele_lowPtMva_ = id_lowpt; }
  if ( id_v2 > -666. ) { ele_mvaIdV2_ = id_v2; }
  if ( ele_conv_vtx_fit_prob > -666. ) { ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; }

  // Set Electron variables
  eleid::Features features;
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace eleid {

  std::vector<float> Features::get() {
    std::vector<float> output = { 
      rho_,
      ele_pt_,
      sc_eta_,
      shape_full5x5_sigmaIetaIeta_,
      shape_full5x5_sigmaIphiIphi_,
      shape_full5x5_circularity_,
      shape_full5x5_r9_,
      sc_etaWidth_,
      sc_phiWidth_,
      shape_full5x5_HoverE_,
      trk_nhits_,
      trk_chi2red_,
      gsf_chi2red_,
      brem_frac_,
      gsf_nhits_,
      match_SC_EoverP_,
      match_eclu_EoverP_,
      match_SC_dEta_,
      match_SC_dPhi_,
      match_seed_dEta_,
      sc_E_,
      trk_p_,
    };
    return output;
  }

  void Features::set( const pat::ElectronRef& ele, double rho ) {

    // KF tracks

    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/Electron.h#L81-L84
    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/src/Electron.cc#L299-L306
    //@@ not what we want! 
    reco::TrackRef trk = ele->closestCtfTrackRef(); 
    if ( trk.isNonnull() ) {
      trk_p_ = float(trk->p());
      trk_nhits_ = float(trk->found());
      trk_chi2red_ = float(trk->normalizedChi2());
    }

    // GSF tracks
    if ( ele->core().isNonnull() ) {
      reco::GsfTrackRef gsf = ele->core()->gsfTrack();
      if ( gsf.isNonnull() ) {
	gsf_nhits_ = gsf->found();
	gsf_chi2red_ = gsf->normalizedChi2();
      }
    }
    
    // Super clusters
    if ( ele->core().isNonnull() ) {
      reco::SuperClusterRef sc = ele->core()->superCluster();
      if ( sc.isNonnull() ) {
	sc_E_ = sc->energy();
	sc_eta_ = sc->eta();
	sc_etaWidth_ = sc->etaWidth();
	sc_phiWidth_ = sc->phiWidth();
      }
    }
    
    // Track-cluster matching
    if ( ele.isNonnull() ) {
      match_seed_dEta_ = ele->deltaEtaSeedClusterTrackAtCalo();
      match_eclu_EoverP_ = (1./ele->ecalEnergy()) - (1./ele->p());
      match_SC_EoverP_ = ele->eSuperClusterOverP();
      match_SC_dEta_ = ele->deltaEtaSuperClusterTrackAtVtx();
      match_SC_dPhi_ = ele->deltaPhiSuperClusterTrackAtVtx();
    }      
    
    // Shower shape vars
    if ( ele.isNonnull() ) {
      shape_full5x5_sigmaIetaIeta_ = ele->full5x5_sigmaIetaIeta();
      shape_full5x5_sigmaIphiIphi_ = ele->full5x5_sigmaIphiIphi();
      shape_full5x5_HoverE_    = ele->full5x5_hcalOverEcal();
      shape_full5x5_r9_ = ele->full5x5_r9();
      shape_full5x5_circularity_   = 1. - ele->full5x5_e1x5() / ele->full5x5_e5x5();
    }
    
    // Misc
    rho_ = rho;
    if ( ele.isNonnull() ) {
      brem_frac_ = ele->fbrem();
      ele_pt_ = ele->pt();
    }
    
  }
  
}
