#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronFeatures.h"
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
  tree->Branch("genOther_dR" , &genOther_dR_ , "genOther_dR/f" );
  tree->Branch("gen_phi", &gen_phi_, "gen_phi/f");
  tree->Branch("gen_p",   &gen_p_,   "gen_p/f");

  // GSF track associated to electron
  tree->Branch("gsf_p",  &gsf_p_, "gsf_p/f");
  tree->Branch("gsf_pt", &gsf_pt_, "gsf_pt/f");
  tree->Branch("gsf_bdtout2", &seed_ptbiased_, "gsf_bdtout2/f");
  tree->Branch("gsf_mode_pt", &gsf_mode_pt_, "gsf_mode_pt/f");
  tree->Branch("gsf_eta", &gsf_eta_, "gsf_eta/f");
  tree->Branch("gsf_phi", &gsf_phi_, "gsf_phi/f");
  tree->Branch("gsf_mode_eta", &gsf_mode_eta_, "gsf_mode_eta/f");
  tree->Branch("gsf_mode_phi", &gsf_mode_phi_, "gsf_mode_phi/f");

  // General track associated to electron
  tree->Branch("trk_p", &trk_p_, "trk_p/f");
  tree->Branch("trk_pt", &trk_pt_, "trk_pt/f");
  tree->Branch("trk_eta", &trk_eta_, "trk_eta/f");
  tree->Branch("trk_phi", &trk_phi_, "trk_phi/f");
  tree->Branch("trk_dEdx1", &trk_dEdx1_, "trk_dEdx1/f");   
  tree->Branch("trk_dEdx1_Nm", &trk_dEdx1_Nm_, "trk_dEdx1_Nm/I");  
  tree->Branch("trk_dEdx1_NSm", &trk_dEdx1_NSm_, "trk_dEdx1_NSm/I");     

  // Electron - kinematics
  tree->Branch("ele_p", &ele_p_, "ele_p/f");
  tree->Branch("ele_pt", &ele_pt_, "ele_pt/f");
  tree->Branch("ele_eta", &ele_eta_, "ele_eta/f");
  tree->Branch("ele_phi", &ele_phi_, "ele_phi/f");
  tree->Branch("fiducial_isEB",&fiducial_isEB_,"fiducial_isEB/I");
  tree->Branch("fiducial_isEE",&fiducial_isEE_,"fiducial_isEE/I"); 
  tree->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap_,"fiducial_isEBEEGap/I"); 

  // Electron - id
  tree->Branch("ele_mva_value", &ele_mva_value_, "ele_mva_value/f");
  tree->Branch("ele_mva_id", &ele_mva_id_, "ele_mva_id/I");
  tree->Branch("eid_rho", &eid_rho_, "eid_rho/f");
  tree->Branch("eid_sc_eta", &eid_sc_eta_, "eid_sc_eta/f");
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
  tree->Branch("eid_gsf_mode_p", &eid_gsf_mode_p_, "eid_gsf_mode_p/f");
  tree->Branch("eid_core_shFracHits", &eid_core_shFracHits_, "eid_core_shFracHits/f");
  tree->Branch("eid_gsf_bdtout1", &eid_gsf_bdtout1_, "eid_gsf_bdtout1/f");
  tree->Branch("eid_gsf_dr", &eid_gsf_dr_, "eid_gsf_dr/f");
  tree->Branch("eid_trk_dr", &eid_trk_dr_, "eid_trk_dr/f");
  tree->Branch("eid_sc_Nclus", &eid_sc_Nclus_, "eid_sc_Nclus/f");
  tree->Branch("eid_sc_clus1_nxtal", &eid_sc_clus1_nxtal_, "eid_sc_clus1_nxtal/f");
  tree->Branch("eid_sc_clus1_dphi", &eid_sc_clus1_dphi_, "eid_sc_clus1_dphi/f");
  tree->Branch("eid_sc_clus2_dphi", &eid_sc_clus2_dphi_, "eid_sc_clus2_dphi/f");
  tree->Branch("eid_sc_clus1_deta", &eid_sc_clus1_deta_, "eid_sc_clus1_deta/f");
  tree->Branch("eid_sc_clus2_deta", &eid_sc_clus2_deta_, "eid_sc_clus2_deta/f");
  tree->Branch("eid_sc_clus1_E", &eid_sc_clus1_E_, "eid_sc_clus1_E/f");
  tree->Branch("eid_sc_clus2_E", &eid_sc_clus2_E_, "eid_sc_clus2_E/f");
  tree->Branch("eid_sc_clus1_E_ov_p", &eid_sc_clus1_E_ov_p_, "eid_sc_clus1_E_ov_p/f");
  tree->Branch("eid_sc_clus2_E_ov_p", &eid_sc_clus2_E_ov_p_, "eid_sc_clus2_E_ov_p/f");

  // Electron energy regression 
  tree->Branch("pre_ecal",&pre_ecal_);
  tree->Branch("pre_ecaltrk",&pre_ecaltrk_);
  tree->Branch("post_ecal",&post_ecal_);
  tree->Branch("post_ecaltrk",&post_ecaltrk_);
  tree->Branch("sc_raw_energy",&sc_raw_energy_);
  tree->Branch("sc_energy",&sc_energy_);

  // Electron - brem fractions
  tree->Branch("brem_fracTrk",&brem_fracTrk_); 
  tree->Branch("brem_fracSC",&brem_fracSC_); 
  tree->Branch("brem_N",&brem_N_,"brem_N/I"); 

  // SuperCluster associated to the electron
  tree->Branch("sc_goodSeed",&sc_goodSeed_,"sc_goodSeed/O");
  tree->Branch("sc_Nclus_deta01",&sc_Nclus_deta01_,"sc_Nclus_deta01/I"); 
  tree->Branch("sc_Nclus_deta02",&sc_Nclus_deta02_,"sc_Nclus_deta02/I"); 
  tree->Branch("sc_Nclus_deta03",&sc_Nclus_deta03_,"sc_Nclus_deta03/I"); 

  // Clusters making the SC 
  tree->Branch("sc_clus1_E_ov_E", &sc_clus1_E_ov_E_, "sc_clus1_E_ov_E/F");
  tree->Branch("sc_clus1_eta",    &sc_clus1_eta_,    "sc_clus1_eta/F");
  tree->Branch("sc_clus1_phi",    &sc_clus1_phi_,    "sc_clus1_phi/F");
  //
  tree->Branch("sc_clus2_E_ov_E", &sc_clus2_E_ov_E_, "sc_clus2_E_ov_E/F");
  tree->Branch("sc_clus2_eta",    &sc_clus2_eta_,    "sc_clus2_eta/F");
  tree->Branch("sc_clus2_phi",    &sc_clus2_phi_,    "sc_clus2_phi/F");
  tree->Branch("sc_clus2_nxtal",  &sc_clus2_nxtal_,  "sc_clus2_nxtal/I");
  //
  tree->Branch("sc_clus3_E",      &sc_clus3_E_,    "sc_clus3_E/F");
  tree->Branch("sc_clus3_E_ov_p", &sc_clus3_E_ov_p_,     "sc_clus3_E_ov_p/F");
  tree->Branch("sc_clus3_E_ov_E", &sc_clus3_E_ov_E_,     "sc_clus3_E_ov_E/F");
  tree->Branch("sc_clus3_eta",    &sc_clus3_eta_,  "sc_clus3_eta/F");
  tree->Branch("sc_clus3_phi",    &sc_clus3_phi_,  "sc_clus3_phi/F");
  tree->Branch("sc_clus3_nxtal",  &sc_clus3_nxtal_, "sc_clus3_nxtal/I");
  tree->Branch("sc_clus3_dphi",   &sc_clus3_dphi_,  "sc_clus3_dphi/F");
  tree->Branch("sc_clus3_deta",   &sc_clus3_deta_,  "sc_clus3_deta/F");

  // Distance wrt muons
  tree->Branch("minDrWithMu", &minDrWithMu_);
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
void IDSlimNtuple::fill_trk( const reco::TrackPtr& trk ) {

  if ( trk.isNonnull() ) {   // should never happen
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
    trk_p_ = trk->p();    
  }
}

void IDSlimNtuple::fill_trk_dEdx( const reco::TrackRef& trk,
				  std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx ) {

  if ( trk.isNonnull() ) {    // should never happen    
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
void IDSlimNtuple::fill_gsf( const reco::GsfTrackPtr gsf ){ 

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
  } 
}

/////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuple::fill_ele( const reco::GsfElectronPtr ele,
			     float mva_value,
			     int mva_id,
			     const double rho, 
			     float seed_unbiased,
			     float field_z) {       
  

  if ( ele.isNonnull() ) {   // should always be the case

    // Kinematics 
    ele_p_   = ele->p();
    ele_pt_  = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();

    // Fiducial flags 
    fiducial_isEB_ = ele->isEB();
    fiducial_isEE_ = ele->isEE();
    fiducial_isEBEEGap_ = ele->isEBEEGap();

    if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
    if ( mva_id > -666 )     { ele_mva_id_ = mva_id; }

    // ElectronID variables: order is important!
    std::vector<float> vfeatures;
    vfeatures = lowptgsfeleid::features_V1(*ele, rho, seed_unbiased, field_z);
    
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
    eid_gsf_mode_p_ = vfeatures[idx++];  
    eid_core_shFracHits_ = vfeatures[idx++]; 
    eid_gsf_bdtout1_ = vfeatures[idx++];  
    eid_gsf_dr_ = vfeatures[idx++];    
    eid_trk_dr_ = vfeatures[idx++];    
    eid_sc_Nclus_ = vfeatures[idx++]; 
    eid_sc_clus1_nxtal_ = vfeatures[idx++]; 
    eid_sc_clus1_dphi_ = vfeatures[idx++]; 
    eid_sc_clus2_dphi_ = vfeatures[idx++]; 
    eid_sc_clus1_deta_ = vfeatures[idx++]; 
    eid_sc_clus2_deta_ = vfeatures[idx++]; 
    eid_sc_clus1_E_ = vfeatures[idx++];  
    eid_sc_clus2_E_ = vfeatures[idx++];  
    eid_sc_clus1_E_ov_p_ = vfeatures[idx++]; 
    eid_sc_clus2_E_ov_p_ = vfeatures[idx++];   

    // Brem fractions and classification 
    brem_fracTrk_ = ele->trackFbrem();
    brem_fracSC_  = ele->superClusterFbrem();
    brem_N_ = ele->numberOfBrems();
  }    
}

// FC new method 
void IDSlimNtuple::fill_supercluster_miniAOD(const reco::GsfElectronPtr ele ) {

  // initialization in case of patological events
  sc_clus1_E_ov_E_ = -999.;      
  sc_clus1_eta_    = -999.;
  sc_clus1_phi_    = -999.;
  sc_clus2_E_ov_E_ = -999.;      
  sc_clus2_eta_    = -999.;
  sc_clus2_phi_    = -999.;
  sc_clus2_nxtal_  = -999;
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
  sc_goodSeed_     = false;

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
  float field_z=3.8112;

  BaseParticlePropagator mypart(RawParticle(mom,pos,kfTrackRef->charge()), 0, 0, field_z);
  mypart.propagateToEcalEntrance(true); // true only first half loop , false more than one loop
  bool reach_ECAL=mypart.getSuccess(); // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps 

  // ECAL entry point for track
  GlobalPoint ecal_pos(mypart.particle().x(), mypart.particle().y(), mypart.particle().z());
  
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
	  if( eleScEne>0)  sc_clus1_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus1_eta_   = cluster->eta();
	  sc_clus1_phi_   = cluster->phi();

	} else if (clusNum==i2){
	  if( eleScEne>0)  sc_clus2_E_ov_E_ = cluster->energy()/eleScEne;
	  sc_clus2_eta_   = cluster->eta();
	  sc_clus2_phi_   = cluster->phi();
	  sc_clus2_nxtal_ = (int) cluster->size();
	  
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
