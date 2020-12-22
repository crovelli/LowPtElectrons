#ifndef reduceSize_h
#define reduceSize_h
#include "TROOT.h"
#include "TRint.h"
#include <TChain.h>
#include <TFile.h>
#include <iostream>

class reduceSize{
  
 public: 
  reduceSize();
  void train();
  virtual ~reduceSize();
  
  UInt_t          evt;
  Float_t         weight;
  Bool_t          is_e;
  Bool_t          is_e_not_matched;
  Bool_t          is_other;
  Bool_t          is_egamma;
  Bool_t          has_trk;
  Bool_t          has_seed;
  Bool_t          has_gsf;
  Bool_t          has_ele;
  Float_t         gen_pt;
  Float_t         gen_eta;
  Int_t           gen_tag_side;
  Float_t         gen_dR;
  Float_t         gen_phi;
  Float_t         gen_p;
  Float_t         gsf_p;
  Float_t         gsf_pt;
  Float_t         gsf_bdtout2;
  Float_t         gsf_mode_pt;
  Float_t         gsf_eta;
  Float_t         gsf_phi;
  Float_t         gsf_mode_eta;
  Float_t         gsf_mode_phi;
  Float_t         trk_p;
  Float_t         trk_pt;
  Float_t         trk_eta;
  Float_t         trk_phi;
  Float_t         ele_p;
  Float_t         ele_pt;
  Float_t         ele_eta;
  Float_t         ele_phi;
  Int_t           fiducial_isEB;
  Int_t           fiducial_isEE;
  Int_t           fiducial_isEBEEGap;
  Float_t         ele_mva_value;
  Int_t           ele_mva_id;
  Float_t         eid_rho;
  Float_t         eid_sc_eta;
  Float_t         eid_shape_full5x5_r9;
  Float_t         eid_sc_etaWidth;
  Float_t         eid_sc_phiWidth;
  Float_t         eid_shape_full5x5_HoverE;
  Float_t         eid_trk_nhits;
  Float_t         eid_trk_chi2red;
  Float_t         eid_gsf_chi2red;
  Float_t         eid_brem_frac;
  Float_t         eid_gsf_nhits;
  Float_t         eid_match_SC_EoverP;
  Float_t         eid_match_eclu_EoverP;
  Float_t         eid_match_SC_dEta;
  Float_t         eid_match_SC_dPhi;
  Float_t         eid_match_seed_dEta;
  Float_t         eid_sc_E;
  Float_t         eid_trk_p;
  Float_t         eid_gsf_mode_p;
  Float_t         eid_core_shFracHits;
  Float_t         eid_gsf_bdtout1;
  Float_t         eid_gsf_dr;
  Float_t         eid_trk_dr;
  Float_t         eid_sc_Nclus;
  Float_t         eid_sc_clus1_nxtal;
  Float_t         eid_sc_clus1_dphi;
  Float_t         eid_sc_clus2_dphi;
  Float_t         eid_sc_clus1_deta;
  Float_t         eid_sc_clus2_deta;
  Float_t         eid_sc_clus1_E;
  Float_t         eid_sc_clus2_E;
  Float_t         eid_sc_clus1_E_ov_p;
  Float_t         eid_sc_clus2_E_ov_p;
  Float_t         pre_ecal;
  Float_t         pre_ecaltrk;
  Float_t         post_ecal;
  Float_t         post_ecaltrk;
  Float_t         sc_raw_energy;
  Float_t         sc_energy;
  Float_t         brem_fracTrk;
  Float_t         brem_fracSC;
  Int_t           brem_N;
  Bool_t          sc_goodSeed;
  Int_t           sc_Nclus_deta01;
  Int_t           sc_Nclus_deta02;
  Int_t           sc_Nclus_deta03;
  Float_t         sc_clus1_E_ov_E;
  Float_t         sc_clus1_eta;
  Float_t         sc_clus1_phi;
  Float_t         sc_clus2_E_ov_E;
  Float_t         sc_clus2_eta;
  Float_t         sc_clus2_phi;
  Int_t           sc_clus2_nxtal;
  Float_t         sc_clus3_E;
  Float_t         sc_clus3_E_ov_p;
  Float_t         sc_clus3_E_ov_E;
  Float_t         sc_clus3_eta;
  Float_t         sc_clus3_phi;
  Int_t           sc_clus3_nxtal;
  Float_t         sc_clus3_dphi;
  Float_t         sc_clus3_deta;
  Int_t type;

};

#endif 
#ifdef reduceSize_cxx

reduceSize::reduceSize(){

}


reduceSize::~reduceSize()
{

}

#endif 

