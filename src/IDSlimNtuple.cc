#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronFeatures.h"
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"
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

  // gen-level particles matched to reco-electron
  tree->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
  tree->Branch("gen_eta", &gen_eta_, "gen_eta/f");
  tree->Branch("gen_tag_side", &gen_tag_side_, "gen_tag_side/I");   
  tree->Branch("gen_dR" , &gen_dR_ , "gen_dR/f" );
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_p",   &gen_p_,   "gen_p/f");

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

  // General track associated to electron
  tree->Branch("trk_dr", &trk_dr_, "trk_dr/f");
  tree->Branch("trk_p", &trk_p_, "trk_p/f");
  tree->Branch("trk_pt", &trk_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &trk_eta_, "trk_eta/f");
  tree->Branch("trk_phi", &trk_phi_, "trk_phi/f");

  // Electron - kinematics
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");
  tree->Branch("core_shFracHits",&core_shFracHits_,"core_shFracHits/f");
  tree->Branch("fiducial_isEB",&fiducial_isEB_,"fiducial_isEB/I");
  tree->Branch("fiducial_isEE",&fiducial_isEE_,"fiducial_isEE/I"); 
  tree->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap_,"fiducial_isEBEEGap/I"); 

  // Electron - id
  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
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

  /*
  // Electron energy regression                                                                                                                  
  tree->Branch("pre_ecal",&pre_ecal_);
  tree->Branch("pre_ecaltrk",&pre_ecaltrk_);
  tree->Branch("post_ecal",&post_ecal_);
  tree->Branch("post_ecaltrk",&post_ecaltrk_);
  tree->Branch("sc_raw_energy",&sc_raw_energy_);
  tree->Branch("sc_energy",&sc_energy_);
  */

  // Electron - brem fractions
  tree->Branch("brem_fracTrk",&brem_fracTrk_); 
  tree->Branch("brem_fracSC",&brem_fracSC_); 
  tree->Branch("brem_N",&brem_N_,"brem_N/I"); 

  // SuperCluster associated to the electron
  tree->Branch("sc_Nclus",&sc_Nclus_,"sc_Nclus/I"); 

  // Clusters making the SC 
  tree->Branch("sc_clus1_E",      &sc_clus1_E_,      "sc_clus1_E/F");
  tree->Branch("sc_clus1_E_ov_p", &sc_clus1_E_ov_p_, "sc_clus1_E_ov_p/F");
  tree->Branch("sc_clus1_E_ov_E", &sc_clus1_E_ov_E_, "sc_clus1_E_ov_E/F");
  tree->Branch("sc_clus1_eta",    &sc_clus1_eta_,    "sc_clus1_eta/F");
  tree->Branch("sc_clus1_phi",    &sc_clus1_phi_,    "sc_clus1_phi/F");
  tree->Branch("sc_clus1_nxtal",  &sc_clus1_nxtal_,  "sc_clus1_nxtal/I");
  tree->Branch("sc_clus1_dphi",   &sc_clus1_dphi_,   "sc_clus1_dphi/F");
  tree->Branch("sc_clus1_deta",   &sc_clus1_deta_,   "sc_clus1_deta/F");
  //
  tree->Branch("sc_clus2_E",      &sc_clus2_E_,    "sc_clus2_E/F");
  tree->Branch("sc_clus2_E_ov_p", &sc_clus2_E_ov_p_,     "sc_clus2_E_ov_p/F");
  tree->Branch("sc_clus2_E_ov_E", &sc_clus2_E_ov_E_,     "sc_clus2_E_ov_E/F");
  tree->Branch("sc_clus2_eta",    &sc_clus2_eta_,  "sc_clus2_eta/F");
  tree->Branch("sc_clus2_phi",    &sc_clus2_phi_,  "sc_clus2_phi/F");
  tree->Branch("sc_clus2_nxtal",  &sc_clus2_nxtal_, "sc_clus2_nxtal/I");
  tree->Branch("sc_clus2_dphi",   &sc_clus2_dphi_,  "sc_clus2_dphi/F");
  tree->Branch("sc_clus2_deta",   &sc_clus2_deta_,  "sc_clus2_deta/F");
  //
  tree->Branch("sc_clus3_E",      &sc_clus3_E_,    "sc_clus3_E/F");
  tree->Branch("sc_clus3_E_ov_p", &sc_clus3_E_ov_p_,     "sc_clus3_E_ov_p/F");
  tree->Branch("sc_clus3_E_ov_E", &sc_clus3_E_ov_E_,     "sc_clus3_E_ov_E/F");
  tree->Branch("sc_clus3_eta",    &sc_clus3_eta_,  "sc_clus3_eta/F");
  tree->Branch("sc_clus3_phi",    &sc_clus3_phi_,  "sc_clus3_phi/F");
  tree->Branch("sc_clus3_nxtal",  &sc_clus3_nxtal_, "sc_clus3_nxtal/I");
  tree->Branch("sc_clus3_dphi",   &sc_clus3_dphi_,  "sc_clus3_dphi/F");
  tree->Branch("sc_clus3_deta",   &sc_clus3_deta_,  "sc_clus3_deta/F");
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_evt( const edm::EventID& id ) {  
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_gen( const reco::GenParticlePtr genp ) {

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
}

void IDSlimNtuple::fill_gen( const pat::PackedGenParticleRef genp ) {  

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
}

void IDSlimNtuple::fill_gen_default() {

  gen_pt_  = -999.;
  gen_eta_ = -999.;
  gen_phi_ = -999.;
  gen_p_   = -999.;
}

/////////////////////////////////////////////////////////////////////////////////
/*
void IDSlimNtuple::fill_trk( const reco::TrackPtr& trk,
			 const reco::BeamSpot& spot ) {       

  if ( trk.isNonnull() ) {   // should never happen
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
    trk_p_ = trk->p();    
  }
}
*/

void IDSlimNtuple::fill_trk_default() {

  trk_pt_  = -999.;
  trk_eta_ = -999.;
  trk_phi_ = -999.;
  trk_p_   = -999.;
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

    // Kinematics (MODE)
    gsf_mode_pt_ = gsf->ptMode();
    gsf_mode_eta_ = gsf->etaMode();
    gsf_mode_phi_ = gsf->phiMode();
    // gsf_mode_p_ = gsf->pMode();
    // std::cout << "HAND: gsf_mode_p_ = " << gsf->pMode() << std::endl;

  } 
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_ele( const reco::GsfElectronPtr ele,
			     float mva_value,
			     float ele_conv_vtx_fit_prob,
			     const double rho, float unbSeed ) {       


  if ( ele.isNonnull() ) {   // should always be the case

    // Kinematics 
    ele_p_   = ele->p();
    ele_pt_  = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();

    // Momentum
    // core_shFracHits_ = ele->shFracInnerHits();
    // std::cout << "HAND: core_shFracHits_ = " << ele->shFracInnerHits() << std::endl;

    // Fiducial flags 
    fiducial_isEB_ = ele->isEB();
    fiducial_isEE_ = ele->isEE();
    fiducial_isEBEEGap_ = ele->isEBEEGap();

    // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
    if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }

    // ElectronID variables
    std::vector<float> vfeatures = lowptgsfeleid::features(ele, rho, unbSeed);
    //@@ ORDER IS IMPORTANT!
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
    gsf_mode_p_ = vfeatures[idx++]; 
    core_shFracHits_ = vfeatures[idx++];   
    seed_unbiased_ = vfeatures[idx++]; 
    gsf_dr_ = vfeatures[idx++]; 
    trk_dr_ = vfeatures[idx++]; 
    sc_Nclus_ = vfeatures[idx++];   
    sc_clus1_nxtal_ = vfeatures[idx++];  
    sc_clus1_dphi_ = vfeatures[idx++]; 
    sc_clus2_dphi_ = vfeatures[idx++]; 
    sc_clus1_deta_ = vfeatures[idx++];  
    sc_clus2_deta_ = vfeatures[idx++];  
    sc_clus1_E_ = vfeatures[idx++];  
    sc_clus2_E_ = vfeatures[idx++];  
    sc_clus1_E_ov_p_ = vfeatures[idx++];
    sc_clus2_E_ov_p_ = vfeatures[idx++];

    // Brem fractions and classification 
    brem_fracTrk_ = ele->trackFbrem();
    brem_fracSC_  = ele->superClusterFbrem();
    brem_N_ = ele->numberOfBrems();

    // Other ID variables (removed from features)
    eid_ele_pt_ = ele->pt();
    eid_shape_full5x5_sigmaIetaIeta_ = ele->full5x5_sigmaIetaIeta();
    eid_shape_full5x5_sigmaIphiIphi_ = ele->full5x5_sigmaIphiIphi();
    eid_shape_full5x5_circularity_ = 1. - ele->full5x5_e1x5() / ele->full5x5_e5x5();
  }    
}

// FC new method 
void IDSlimNtuple::fill_supercluster_miniAOD(const reco::GsfElectronPtr ele ) {

  // initialization in case of patological events
  // sc_clus1_E_      = -999.;
  // sc_clus1_E_ov_p_ = -999.;
  sc_clus1_E_ov_E_ = -999.;      
  sc_clus1_eta_    = -999.;
  sc_clus1_phi_    = -999.;
  // sc_clus1_nxtal_  = -999;
  // sc_clus1_deta_   = -999.;
  // sc_clus1_dphi_   = -999.;
  // sc_clus2_E_      = -999.;
  // sc_clus2_E_ov_p_ = -999.;
  sc_clus2_E_ov_E_ = -999.;      
  sc_clus2_eta_    = -999.;
  sc_clus2_phi_    = -999.;
  sc_clus2_nxtal_  = -999;
  //sc_clus2_deta_   = -999.;
  //sc_clus2_dphi_   = -999.;
  sc_clus3_E_      = -999.;
  sc_clus3_E_ov_p_ = -999.;
  sc_clus3_E_ov_E_ = -999.;      
  sc_clus3_eta_    = -999.;
  sc_clus3_phi_    = -999.;
  sc_clus3_nxtal_  = -999;
  sc_clus3_deta_   = -999.;
  sc_clus3_dphi_   = -999.;
  // sc_Nclus_        = -999;

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

  BaseParticlePropagator mypart(RawParticle(mom, pos, kfTrackRef->charge()), 0, 0, field_z);
  mypart.propagateToEcalEntrance(true); // true only first half loop , false more than one loop
  bool reach_ECAL=mypart.getSuccess(); // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps 

  // ECAL entry point for track
  GlobalPoint ecal_pos(mypart.particle().vertex().x(), mypart.particle().vertex().y(), mypart.particle().vertex().z());
  // Preshower limit
  //  bool below_ps = pow(ecal_pos.z(), 2.) > boundary_ * ecal_pos.perp2();
  // Iterate through ECAL clusters
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
	  // sc_clus1_E_     = cluster->energy();
	  // std::cout << "HAND: sc_clus1_E_ = " << cluster->energy() << std::endl;
	  // if( elePmode>0 ) sc_clus1_E_ov_p_ = cluster->energy()/elePmode;
	  // std::cout << "HAND: sc_clus1_E_ov_p_ = " << cluster->energy()/elePmode << std::endl;
	  if( eleScEne>0)  sc_clus1_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus1_eta_   = cluster->eta();
	  sc_clus1_phi_   = cluster->phi();
	  // sc_clus1_nxtal_ =(int) cluster->size();
	  // std::cout << "HAND: sc_clus1_nxtal_ = " << sc_clus1_nxtal_ << std::endl;
	  // if(reach_ECAL>0){
	  //sc_clus1_deta_ = deta;
	  //sc_clus1_dphi_ = dphi;
	  //}
	  // std::cout << "HAND: sc_clus1_dphi_ = " << dphi 
	  //	       << ", sc_clus1_deta_ = " << deta << std::endl;
	  
	} else if (clusNum==i2){
	  // sc_clus2_E_     = cluster->energy();
	  // std::cout << "HAND: sc_clus2_E_ = " << cluster->energy() << std::endl;
	  // if( elePmode>0 ) sc_clus2_E_ov_p_ = cluster->energy()/elePmode;
	  // std::cout << "HAND: sc_clus2_E_ov_p_ = " << cluster->energy()/elePmode << std::endl;
	  if( eleScEne>0)  sc_clus2_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus2_eta_   = cluster->eta();
	  sc_clus2_phi_   = cluster->phi();
	  sc_clus2_nxtal_ = (int) cluster->size();
	  //if(reach_ECAL>0){
	  //sc_clus2_deta_ = deta;
	  //sc_clus2_dphi_ = dphi;
	  //}
	  // std::cout << "HAND: sc_clus2_dphi_ = " << dphi 
	  // << ", sc_clus2_deta_ = " << deta << std::endl;
	  
	} else if (clusNum==i3){
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
	}
	clusNum++;
      }
    }
  }catch(...){
    //
  }

  // sc_Nclus_ = sc->clustersSize();
  // std::cout << "HAND: sc_Nclus_ = " << sc->clustersSize() << std::endl;

}

// end FC

template < typename T> 
bool IDSlimNtuple::validPtr(edm::Ptr<T>& ptr){
  return (ptr.isNonnull() && ptr.isAvailable() );
}
