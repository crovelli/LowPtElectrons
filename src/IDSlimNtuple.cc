#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDHeavyObjectCache.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"
#include "TTree.h"
#include <iostream>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//
void IDSlimNtuple::link_tree( TTree *tree ) {

  std::cout<<"I am running IDSlimNtuple::link_tree"<<std::endl; 

  // general
  tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");  
  tree->Branch("rho", &rho_, "rho/f");
  tree->Branch("is_aod", &is_aod_, "is_aod/i");
  tree->Branch("is_mc", &is_mc_, "is_mc/i");
  tree->Branch("is_e", &is_e_, "is_e/O");
  tree->Branch("is_other", &is_other_, "is_other/O");
  tree->Branch("is_egamma", &is_egamma_, "is_egamma/O");

  // gen-level particles matched to reco-electron
  tree->Branch("gen_dR" , &gen_dR_ , "gen_dR/f" );
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_p", &gen_p_, "gen_p/f");
  tree->Branch("gen_charge", &gen_charge_, "gen_charge/I");
  tree->Branch("gen_pdgid", &gen_pdgid_, "gen_pdgid/I");
  tree->Branch("gen_mom_pdgid", &gen_mom_pdgid_, "gen_mom_pdgid/I");
  tree->Branch("gen_gran_pdgid", &gen_gran_pdgid_, "gen_gran_pdgid/I");

  // GSF track associated to electron
  tree->Branch("gsf_dr", &gsf_dr_, "gsf_dr/f");
  tree->Branch("gsf_bdtout1", &seed_unbiased_, "gsf_bdtout1/f");
  tree->Branch("gsf_bdtout2", &seed_ptbiased_, "gsf_bdtout2/f");
  tree->Branch("gsf_pt", &gsf_pt_, "gsf_pt/f");
  tree->Branch("gsf_eta", &gsf_eta_, "gsf_eta/f");
  tree->Branch("gsf_phi", &gsf_phi_, "gsf_phi/f");
  tree->Branch("gsf_p", &gsf_p_, "gsf_p/f");
  tree->Branch("gsf_charge", &gsf_charge_, "gsf_charge/I");
  tree->Branch("gsf_inp", &gsf_inp_, "gsf_inp/f");
  tree->Branch("gsf_outp", &gsf_outp_, "gsf_outp/f");
  tree->Branch("gsf_mode_pt", &gsf_mode_pt_, "gsf_mode_pt/f");
  tree->Branch("gsf_mode_eta", &gsf_mode_eta_, "gsf_mode_eta/f");
  tree->Branch("gsf_mode_phi", &gsf_mode_phi_, "gsf_mode_phi/f");
  tree->Branch("gsf_mode_p", &gsf_mode_p_, "gsf_mode_p/f");
  tree->Branch("gsf_missing_inner_hits", &gsf_missing_inner_hits_, "gsf_missing_inner_hits/I");
  tree->Branch("gsf_dxy",  &gsf_dxy_, "gsf_dxy/f");
  tree->Branch("gsf_dxy_err",&gsf_dxy_err_, "gsf_dxy_err/f");
  tree->Branch("gsf_dz",  &gsf_dz_, "gsf_dz/f");
  tree->Branch("gsf_dz_err",&gsf_dz_err_, "gsf_dz_err/f");

  //tree->Branch("gsf_ntangents", &gsf_ntangents_, "gsf_ntangents/I");
  //tree->Branch("gsf_hit_dpt", gsf_hit_dpt_, "gsf_hit_dpt[gsf_ntangents]/f");
  //tree->Branch("gsf_hit_dpt_unc", gsf_hit_dpt_unc_, "gsf_hit_dpt_unc[gsf_ntangents]/f");
  //tree->Branch("gsf_extapolated_eta", &gsf_extapolated_eta_);
  //tree->Branch("gsf_extapolated_phi", &gsf_extapolated_phi_);

  // General track associated to electron
  tree->Branch("trk_dr", &trk_dr_, "trk_dr/f");
  tree->Branch("trk_pt", &trk_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &trk_eta_, "trk_eta/f");
  tree->Branch("trk_phi", &trk_phi_, "trk_phi/f");
  tree->Branch("trk_p", &trk_p_, "trk_p/f");
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

  // Electron - kinematics
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");
  tree->Branch("p4kind",&p4kind_,"p4kind/I"); 
  tree->Branch("core_shFracHits",&core_shFracHits_,"core_shFracHits/f");
  tree->Branch("ele_p_atvtx",&ele_p_atvtx_,"ele_p_atvtx/F"); 
  tree->Branch("ele_p_atcalo",&ele_p_atcalo_,"ele_p_atcalo/F"); 

  tree->Branch("fiducial_isEB",&fiducial_isEB_,"fiducial_isEB/I");
  tree->Branch("fiducial_isEE",&fiducial_isEE_,"fiducial_isEE/I"); 
  tree->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap_,"fiducial_isEBEEGap/I"); 

  // Electron - charge
  tree->Branch("chPix",&chPix_,"chPix/I");
  tree->Branch("chGCP",&chGCP_,"chGCP/I");
  tree->Branch("chGP",&chGP_,"chGP/I");
  tree->Branch("chGC",&chGC_,"chGC/I");

  // Electron - id
  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_id", &ele_mva_id_, "ele_mva_id/I");
  tree->Branch("ele_conv_vtx_fit_prob", &ele_conv_vtx_fit_prob_, "ele_conv_vtx_fit_prob/f");
  tree->Branch("eid_rho", &eid_rho_, "eid_rho/f");
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

  // Electron - further track-Cluster matching
  tree->Branch("match_seed_EoverP",&match_seed_EoverP_); 
  tree->Branch("match_seed_EoverPout",&match_seed_EoverPout_); 
  tree->Branch("match_seed_dPhi",&match_seed_dPhi_); 
  tree->Branch("match_seed_dEta_vtx",&match_seed_dEta_vtx_); 
  tree->Branch("match_eclu_EoverPout",&match_eclu_EoverPout_); 
  tree->Branch("match_eclu_dEta",&match_eclu_dEta_); 
  tree->Branch("match_eclu_dPhi",&match_eclu_dPhi_); 

  // Electron - further full 5x5 shower shapes
  tree->Branch("shape_full5x5_e1x5",&shape_full5x5_e1x5_); 
  tree->Branch("shape_full5x5_e2x5Max",&shape_full5x5_e2x5Max_); 
  tree->Branch("shape_full5x5_e5x5",&shape_full5x5_e5x5_);
  tree->Branch("shape_full5x5_HoverEBc",&shape_full5x5_HoverEBc_); 
  tree->Branch("shape_full5x5_hcalDepth1",&shape_full5x5_hcalDepth1_); 
  tree->Branch("shape_full5x5_hcalDepth2",&shape_full5x5_hcalDepth2_); 
  tree->Branch("shape_full5x5_hcalDepth1Bc",&shape_full5x5_hcalDepth1Bc_); 
  tree->Branch("shape_full5x5_hcalDepth2Bc",&shape_full5x5_hcalDepth2Bc_); 
  tree->Branch("shape_full5x5_eLeft",&shape_full5x5_eLeft_); 
  tree->Branch("shape_full5x5_eRight",&shape_full5x5_eRight_); 
  tree->Branch("shape_full5x5_eTop",&shape_full5x5_eTop_); 
  tree->Branch("shape_full5x5_eBottom",&shape_full5x5_eBottom_); 

  // Electron - brem fractions
  tree->Branch("brem_fracTrk",&brem_fracTrk_); 
  tree->Branch("brem_fracSC",&brem_fracSC_); 
  tree->Branch("brem_N",&brem_N_,"brem_N/I"); 

  // SuperCluster associated to the electron
  tree->Branch("sc_Et",&sc_Et_); 
  tree->Branch("sc_Nclus",&sc_Nclus_,"sc_Nclus/I"); 
  tree->Branch("sc_goodSeed",&sc_goodSeed_,"sc_goodSeed/O");

  // Clusters making the SC
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
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_evt( const edm::EventID& id ) {  
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_gen( const reco::GenParticlePtr genp ) {

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
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

void IDSlimNtuple::fill_gen( const pat::PackedGenParticleRef genp ) {  

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
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

void IDSlimNtuple::fill_gen_default() {

  gen_pt_  = -999.;
  gen_eta_ = -999.;
  gen_phi_ = -999.;
  gen_p_   = -999.;
  gen_charge_ = -999; 
  gen_pdgid_  = -999;
  gen_mom_pdgid_  = -999; 
  gen_gran_pdgid_ = -999;
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_trk( const reco::TrackPtr& trk,
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
    }
    // quality
    trk_nhits_ = trk->found();
    trk_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    trk_chi2red_ = trk->normalizedChi2();
    // displ
    trk_dxy_ = trk->dxy(spot);
    trk_dxy_err_ = trk->dxyError();
    trk_dz_ = trk->dz(spot.position());
    trk_dz_err_ = trk->dzError();    
  }
  
}

void IDSlimNtuple::fill_trk_dEdx( const reco::TrackPtr& trk,
				  std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx ) {

  if ( trk.isNonnull() ) {
    const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[0]);
    const reco::DeDxData& dedx = dEdxTrack[trk];
    trk_dEdx1_=dedx.dEdx();
    trk_dEdx1_Nm_=dedx.numberOfMeasurements();
    trk_dEdx1_NSm_=dedx.numberOfSaturatedMeasurements();
  }

}

void IDSlimNtuple::fill_trk_dEdx_default( ) {

  trk_dEdx1_     = -999.;
  trk_dEdx1_Nm_  = -999;
  trk_dEdx1_NSm_ = -999;
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_bdt( double seed_unbiased, 
			 double seed_ptbiased ) {          
  seed_unbiased_ = seed_unbiased;
  seed_ptbiased_ = seed_ptbiased;
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
    gsf_charge_ = gsf->charge();
    if ( gsf->extra().isAvailable() && gsf->extra().isNonnull() ) {
      gsf_inp_ = sqrt(gsf->innerMomentum().mag2());
      gsf_outp_ = sqrt(gsf->outerMomentum().mag2());
    }
    
    // Kinematics (MODE)
    gsf_mode_pt_ = gsf->ptMode();
    gsf_mode_eta_ = gsf->etaMode();
    gsf_mode_phi_ = gsf->phiMode();
    gsf_mode_p_ = gsf->pMode();

    // Quality
    gsf_missing_inner_hits_ = gsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

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
void IDSlimNtuple::fill_ele( const reco::GsfElectronPtr ele,
			 float mva_value,
			 int mva_id,
			 float ele_conv_vtx_fit_prob,
			 const double rho ) {       


  if ( ele.isNonnull() ) {  

    // Kinematics 
    ele_p_   = ele->p();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();

    // Momentum
    p4kind_  = ele->candidateP4Kind();
    core_shFracHits_ = ele->shFracInnerHits();
    ele_p_atvtx_  = sqrt(ele->trackMomentumAtVtx().mag2());
    ele_p_atcalo_ = sqrt(ele->trackMomentumAtCalo().mag2());
    // Fiducial flags 
    fiducial_isEB_ = ele->isEB();
    fiducial_isEE_ = ele->isEE();
    fiducial_isEBEEGap_ = ele->isEBEEGap();

    // Charge 
    chPix_ = ele->scPixCharge();
    chGCP_ = ele->isGsfCtfScPixChargeConsistent();
    chGP_  = ele->isGsfScPixChargeConsistent();
    chGC_  = ele->isGsfCtfChargeConsistent();    

    // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
    if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
    if ( mva_id > -666 ) { ele_mva_id_ = mva_id; }
    if ( ele_conv_vtx_fit_prob > -666. ) { ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; }
    
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

    // Further track-Cluster matching 
    match_seed_EoverP_    = ele->eSeedClusterOverP();
    match_seed_EoverPout_ = ele->eSeedClusterOverPout();
    match_seed_dPhi_      = ele->deltaPhiSeedClusterTrackAtCalo();
    match_seed_dEta_vtx_  = ele->deltaEtaSeedClusterTrackAtVtx();
    match_eclu_EoverPout_ = ele->eEleClusterOverPout();
    match_eclu_dEta_ = ele->deltaEtaEleClusterTrackAtCalo();
    match_eclu_dPhi_ = ele->deltaPhiEleClusterTrackAtCalo();

    // Further full 5x5 shower shape 
    shape_full5x5_e1x5_     = ele->full5x5_e1x5();
    shape_full5x5_e2x5Max_  = ele->full5x5_e2x5Max();
    shape_full5x5_e5x5_     = ele->full5x5_e5x5();
    shape_full5x5_HoverEBc_ = ele->full5x5_hcalOverEcalBc();
    shape_full5x5_hcalDepth1_    = ele->full5x5_hcalDepth1OverEcal();
    shape_full5x5_hcalDepth2_    = ele->full5x5_hcalDepth2OverEcal();
    shape_full5x5_hcalDepth1Bc_  = ele->full5x5_hcalDepth1OverEcalBc();
    shape_full5x5_hcalDepth2Bc_  = ele->full5x5_hcalDepth2OverEcalBc();
    shape_full5x5_eLeft_   = ele->full5x5_eLeft();
    shape_full5x5_eRight_  = ele->full5x5_eRight();
    shape_full5x5_eTop_    = ele->full5x5_eTop();
    shape_full5x5_eBottom_ = ele->full5x5_eBottom();

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
  sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );
  sc_Nclus_ = sc->clustersSize();

  float seedEne = sc->seed()->energy();
  if ( fabs(seedEne-maxEne)<0.001 ) sc_goodSeed_ = true;
}

