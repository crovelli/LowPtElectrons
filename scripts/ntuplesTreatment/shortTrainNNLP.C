#define shortTrainNNLP_cxx
#include "shortTrainNNLP.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMultiLayerPerceptron.h>
#include <TMLPAnalyzer.h>
#include <TNeuron.h>
#include <TSynapse.h>
#include <TLegend.h>

#include <stdio.h>

using namespace std; 

void shortTrainNNLP::train(Int_t ntrain=100) {

  std::cout<<"hello ... welcome to shortTrain"<<std::endl; 
   const char *fname = "/eos/cms/store/user/crovelli/LowPtEle/Batch3/BuToKJpsiToee_all.root";

   TFile *input = 0;
      input = TFile::Open(fname);

   if (!input) return;
   TDirectory * dir = (TDirectory*)input->Get("/eos/cms/store/user/crovelli/LowPtEle/Batch3/BuToKJpsiToee_all.root:/ntuplizer");
   TTree *fChain ;
   dir->GetObject("tree",fChain);


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
   Float_t         gsf_dr;
   Float_t         gsf_pt;
   Float_t         gsf_bdtout1;
   Float_t         gsf_mode_p;
   Float_t         gsf_eta;
   Float_t         trk_dr;
   Float_t         trk_pt;
   Float_t         trk_eta;
   Float_t         ele_p;
   Float_t         ele_pt;
   Float_t         ele_eta;
   Float_t         core_shFracHits;
   Int_t           fiducial_isEB;
   Int_t           fiducial_isEE;
   Int_t           fiducial_isEBEEGap;
   Float_t         ele_mva_value;
   Int_t           ele_mva_id;
   Float_t         eid_ele_pt;
   Float_t         eid_sc_eta;
   Float_t         eid_shape_full5x5_sigmaIetaIeta;
   Float_t         eid_shape_full5x5_sigmaIphiIphi;
   Float_t         eid_shape_full5x5_circularity;
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
   Float_t         eid_rho;
   Float_t         brem_fracTrk;
   Float_t         brem_fracSC;
   Int_t           brem_N;
   Bool_t          sc_goodSeed;
   Int_t           sc_Nclus;
   Int_t           sc_Nclus_deta01;
   Int_t           sc_Nclus_deta02;
   Int_t           sc_Nclus_deta03;
   Float_t         sc_clus1_E;
   Float_t         sc_clus1_E_ov_p;
   Float_t         sc_clus1_E_ov_E;
   Float_t         sc_clus1_eta;
   Float_t         sc_clus1_phi;
   Int_t           sc_clus1_nxtal;
   Float_t         sc_clus1_dphi;
   Float_t         sc_clus1_deta;
   Float_t         sc_clus1_ntrk_deta01;
   Float_t         sc_clus2_E;
   Float_t         sc_clus2_E_ov_p;
   Float_t         sc_clus2_E_ov_E;
   Float_t         sc_clus2_eta;
   Float_t         sc_clus2_phi;
   Int_t           sc_clus2_nxtal;
   Float_t         sc_clus2_dphi;
   Float_t         sc_clus2_deta;
   Float_t         sc_clus2_ntrk_deta01;
   Float_t         sc_clus3_E;
   Float_t         sc_clus3_E_ov_p;
   Float_t         sc_clus3_E_ov_E;
   Float_t         sc_clus3_eta;
   Float_t         sc_clus3_phi;
   Int_t           sc_clus3_nxtal;
   Float_t         sc_clus3_dphi;
   Float_t         sc_clus3_deta;
   Float_t         sc_clus3_ntrk_deta01;
   Int_t type;
   //
 fChain->SetBranchAddress("evt", &evt);
 fChain->SetBranchAddress("weight", &weight);
 fChain->SetBranchAddress("is_e", &is_e);
 fChain->SetBranchAddress("is_e_not_matched", &is_e_not_matched);
 fChain->SetBranchAddress("is_other", &is_other);
 fChain->SetBranchAddress("is_egamma", &is_egamma);
 fChain->SetBranchAddress("has_trk", &has_trk);
 fChain->SetBranchAddress("has_seed", &has_seed);
 fChain->SetBranchAddress("has_gsf", &has_gsf);
 fChain->SetBranchAddress("has_ele", &has_ele);
 fChain->SetBranchAddress("gen_pt", &gen_pt);
 fChain->SetBranchAddress("gen_eta", &gen_eta);
 fChain->SetBranchAddress("gen_tag_side", &gen_tag_side);
 fChain->SetBranchAddress("gen_dR", &gen_dR);
 fChain->SetBranchAddress("gsf_dr", &gsf_dr);
 fChain->SetBranchAddress("gsf_pt", &gsf_pt);
 fChain->SetBranchAddress("gsf_bdtout1", &gsf_bdtout1);
 fChain->SetBranchAddress("gsf_mode_p", &gsf_mode_p);
 fChain->SetBranchAddress("gsf_eta", &gsf_eta);
 fChain->SetBranchAddress("trk_dr", &trk_dr);
 fChain->SetBranchAddress("trk_pt", &trk_pt);
 fChain->SetBranchAddress("trk_eta", &trk_eta);
 fChain->SetBranchAddress("ele_p", &ele_p);
 fChain->SetBranchAddress("ele_pt", &ele_pt);
 fChain->SetBranchAddress("ele_eta", &ele_eta);
 fChain->SetBranchAddress("core_shFracHits", &core_shFracHits);
 fChain->SetBranchAddress("fiducial_isEB", &fiducial_isEB);
 fChain->SetBranchAddress("fiducial_isEE", &fiducial_isEE);
 fChain->SetBranchAddress("fiducial_isEBEEGap", &fiducial_isEBEEGap);
 fChain->SetBranchAddress("ele_mva_value", &ele_mva_value);
 fChain->SetBranchAddress("ele_mva_id", &ele_mva_id);
 fChain->SetBranchAddress("eid_ele_pt", &eid_ele_pt);
 fChain->SetBranchAddress("eid_sc_eta", &eid_sc_eta);
 fChain->SetBranchAddress("eid_shape_full5x5_sigmaIetaIeta", &eid_shape_full5x5_sigmaIetaIeta);
 fChain->SetBranchAddress("eid_shape_full5x5_sigmaIphiIphi", &eid_shape_full5x5_sigmaIphiIphi);
 fChain->SetBranchAddress("eid_shape_full5x5_circularity", &eid_shape_full5x5_circularity);
 fChain->SetBranchAddress("eid_shape_full5x5_r9", &eid_shape_full5x5_r9);
 fChain->SetBranchAddress("eid_sc_etaWidth", &eid_sc_etaWidth);
 fChain->SetBranchAddress("eid_sc_phiWidth", &eid_sc_phiWidth);
 fChain->SetBranchAddress("eid_shape_full5x5_HoverE", &eid_shape_full5x5_HoverE);
 fChain->SetBranchAddress("eid_trk_nhits", &eid_trk_nhits);
 fChain->SetBranchAddress("eid_trk_chi2red", &eid_trk_chi2red);
 fChain->SetBranchAddress("eid_gsf_chi2red", &eid_gsf_chi2red);
 fChain->SetBranchAddress("eid_brem_frac", &eid_brem_frac);
 fChain->SetBranchAddress("eid_gsf_nhits", &eid_gsf_nhits);
 fChain->SetBranchAddress("eid_match_SC_EoverP", &eid_match_SC_EoverP);
 fChain->SetBranchAddress("eid_match_eclu_EoverP", &eid_match_eclu_EoverP);
 fChain->SetBranchAddress("eid_match_SC_dEta", &eid_match_SC_dEta);
 fChain->SetBranchAddress("eid_match_SC_dPhi", &eid_match_SC_dPhi);
 fChain->SetBranchAddress("eid_match_seed_dEta", &eid_match_seed_dEta);
 fChain->SetBranchAddress("eid_sc_E", &eid_sc_E);
 fChain->SetBranchAddress("eid_trk_p", &eid_trk_p);
 fChain->SetBranchAddress("eid_rho", &eid_rho);
 fChain->SetBranchAddress("brem_fracTrk", &brem_fracTrk);
 fChain->SetBranchAddress("brem_fracSC", &brem_fracSC);
 fChain->SetBranchAddress("brem_N", &brem_N);
 fChain->SetBranchAddress("sc_goodSeed", &sc_goodSeed);
 fChain->SetBranchAddress("sc_Nclus", &sc_Nclus);
 fChain->SetBranchAddress("sc_Nclus_deta01", &sc_Nclus_deta01);
 fChain->SetBranchAddress("sc_Nclus_deta02", &sc_Nclus_deta02);
 fChain->SetBranchAddress("sc_Nclus_deta03", &sc_Nclus_deta03);
 fChain->SetBranchAddress("sc_clus1_E", &sc_clus1_E);
 fChain->SetBranchAddress("sc_clus1_E_ov_p", &sc_clus1_E_ov_p);
 fChain->SetBranchAddress("sc_clus1_E_ov_E", &sc_clus1_E_ov_E);
 fChain->SetBranchAddress("sc_clus1_eta", &sc_clus1_eta);
 fChain->SetBranchAddress("sc_clus1_phi", &sc_clus1_phi);
 fChain->SetBranchAddress("sc_clus1_nxtal", &sc_clus1_nxtal);
 fChain->SetBranchAddress("sc_clus1_dphi", &sc_clus1_dphi);
 fChain->SetBranchAddress("sc_clus1_deta", &sc_clus1_deta);
 fChain->SetBranchAddress("sc_clus1_ntrk_deta01", &sc_clus1_ntrk_deta01);
 fChain->SetBranchAddress("sc_clus2_E", &sc_clus2_E);
 fChain->SetBranchAddress("sc_clus2_E_ov_p", &sc_clus2_E_ov_p);
 fChain->SetBranchAddress("sc_clus2_E_ov_E", &sc_clus2_E_ov_E);
 fChain->SetBranchAddress("sc_clus2_eta", &sc_clus2_eta);
 fChain->SetBranchAddress("sc_clus2_phi", &sc_clus2_phi);
 fChain->SetBranchAddress("sc_clus2_nxtal", &sc_clus2_nxtal);
 fChain->SetBranchAddress("sc_clus2_dphi", &sc_clus2_dphi);
 fChain->SetBranchAddress("sc_clus2_deta", &sc_clus2_deta);
 fChain->SetBranchAddress("sc_clus2_ntrk_deta01", &sc_clus2_ntrk_deta01);
 fChain->SetBranchAddress("sc_clus3_E", &sc_clus3_E);
 fChain->SetBranchAddress("sc_clus3_E_ov_p", &sc_clus3_E_ov_p);
 fChain->SetBranchAddress("sc_clus3_E_ov_E", &sc_clus3_E_ov_E);
 fChain->SetBranchAddress("sc_clus3_eta", &sc_clus3_eta);
 fChain->SetBranchAddress("sc_clus3_phi", &sc_clus3_phi);
 fChain->SetBranchAddress("sc_clus3_nxtal", &sc_clus3_nxtal);
 fChain->SetBranchAddress("sc_clus3_dphi", &sc_clus3_dphi);
 fChain->SetBranchAddress("sc_clus3_deta", &sc_clus3_deta);
 fChain->SetBranchAddress("sc_clus3_ntrk_deta01", &sc_clus3_ntrk_deta01);

 
 TFile *fileNew = TFile::Open("/tmp/simu_shortLP.root","RECREATE");     
 fileNew->ls();
 fileNew->cd();
 TDirectory *myDir = (TDirectory*)fileNew->mkdir("ntuplizer");
 myDir->cd();
 TTree *simu = new TTree("tree", "Filtered Monte Carlo Events");  


// general                                                                                                                                
simu->Branch("evt",  &evt, "evt/i");
simu->Branch("weight", &weight, "weight/f");
simu->Branch("is_e", &is_e, "is_e/O");
simu->Branch("is_e_not_matched", &is_e_not_matched, "is_e_not_matched/O");
simu->Branch("is_other", &is_other, "is_other/O");
simu->Branch("is_egamma", &is_egamma, "is_egamma/O");
simu->Branch("has_trk", &has_trk, "has_trk/O");
simu->Branch("has_seed", &has_seed, "has_seed/O");
simu->Branch("has_gsf", &has_gsf, "has_gsf/O");
simu->Branch("has_ele", &has_ele, "has_ele/O");
simu->Branch("gen_pt" , &gen_pt , "gen_pt/f" );
simu->Branch("gen_eta", &gen_eta, "gen_eta/f");
simu->Branch("gen_tag_side", &gen_tag_side, "gen_tag_side/I");
simu->Branch("gen_dR" , &gen_dR , "gen_dR/f" );
simu->Branch("gsf_dr", &gsf_dr, "gsf_dr/f");
simu->Branch("gsf_pt", &gsf_pt, "gsf_pt/f");
simu->Branch("gsf_bdtout1", &gsf_bdtout1, "gsf_bdtout1/f");
simu->Branch("gsf_mode_p", &gsf_mode_p, "gsf_mode_p/f");
simu->Branch("gsf_eta", &gsf_eta, "gsf_eta/f");
simu->Branch("trk_dr", &trk_dr, "trk_dr/f");
simu->Branch("trk_pt", &trk_pt, "trk_pt/f");
simu->Branch("trk_eta", &trk_eta, "trk_eta/f");
simu->Branch("ele_p", &ele_p, "ele_p/f");
simu->Branch("ele_pt", &ele_pt, "ele_pt/f");
simu->Branch("ele_eta", &ele_eta, "ele_eta/f");
simu->Branch("core_shFracHits",&core_shFracHits,"core_shFracHits/f");
simu->Branch("fiducial_isEB",&fiducial_isEB,"fiducial_isEB/I");
simu->Branch("fiducial_isEE",&fiducial_isEE,"fiducial_isEE/I");
simu->Branch("fiducial_isEBEEGap",&fiducial_isEBEEGap,"fiducial_isEBEEGap/I");
simu->Branch("ele_mva_value", &ele_mva_value, "ele_mva_value/f");
simu->Branch("ele_mva_id", &ele_mva_id, "ele_mva_id/I");
simu->Branch("eid_rho", &eid_rho, "eid_rho/f");
simu->Branch("eid_ele_pt", &eid_ele_pt, "eid_ele_pt/f");
simu->Branch("eid_sc_eta", &eid_sc_eta, "eid_sc_eta/f");
simu->Branch("eid_shape_full5x5_sigmaIetaIeta", &eid_shape_full5x5_sigmaIetaIeta, "eid_shape_full5x5_sigmaIetaIeta/f");
simu->Branch("eid_shape_full5x5_sigmaIphiIphi", &eid_shape_full5x5_sigmaIphiIphi, "eid_shape_full5x5_sigmaIphiIphi/f");
simu->Branch("eid_shape_full5x5_circularity", &eid_shape_full5x5_circularity, "eid_shape_full5x5_circularity/f");
simu->Branch("eid_shape_full5x5_r9", &eid_shape_full5x5_r9, "eid_shape_full5x5_r9/f");
simu->Branch("eid_sc_etaWidth", &eid_sc_etaWidth, "eid_sc_etaWidth/f");
simu->Branch("eid_sc_phiWidth", &eid_sc_phiWidth, "eid_sc_phiWidth/f");
simu->Branch("eid_shape_full5x5_HoverE", &eid_shape_full5x5_HoverE, "eid_shape_full5x5_HoverE/f");
simu->Branch("eid_trk_nhits", &eid_trk_nhits, "eid_trk_nhits/f");
simu->Branch("eid_trk_chi2red", &eid_trk_chi2red, "eid_trk_chi2red/f");
simu->Branch("eid_gsf_chi2red", &eid_gsf_chi2red, "eid_gsf_chi2red/f");
simu->Branch("eid_brem_frac", &eid_brem_frac, "eid_brem_frac/f");
simu->Branch("eid_gsf_nhits", &eid_gsf_nhits, "eid_gsf_nhits/f");
simu->Branch("eid_match_SC_EoverP", &eid_match_SC_EoverP, "eid_match_SC_EoverP/f");
simu->Branch("eid_match_eclu_EoverP", &eid_match_eclu_EoverP, "eid_match_eclu_EoverP/f");
simu->Branch("eid_match_SC_dEta", &eid_match_SC_dEta, "eid_match_SC_dEta/f");
simu->Branch("eid_match_SC_dPhi", &eid_match_SC_dPhi, "eid_match_SC_dPhi/f");
simu->Branch("eid_match_seed_dEta", &eid_match_seed_dEta, "eid_match_seed_dEta/f");
simu->Branch("eid_sc_E", &eid_sc_E,   "eid_sc_E/f");
simu->Branch("eid_trk_p", &eid_trk_p, "eid_trk_p/f");
 simu->Branch("brem_fracTrk",&brem_fracTrk,"brem_fracTrk/f");
 simu->Branch("brem_fracSC",&brem_fracSC,"brem_fracSC/f");
simu->Branch("brem_N",&brem_N,"brem_N/I");
simu->Branch("sc_goodSeed",&sc_goodSeed,"sc_goodSeed/O");
simu->Branch("sc_Nclus",&sc_Nclus,"sc_Nclus/I");
simu->Branch("sc_Nclus_deta01",&sc_Nclus_deta01,"sc_Nclus_deta01/I");
simu->Branch("sc_Nclus_deta02",&sc_Nclus_deta02,"sc_Nclus_deta02/I");
simu->Branch("sc_Nclus_deta03",&sc_Nclus_deta03,"sc_Nclus_deta03/I");
simu->Branch("sc_clus1_E",     &sc_clus1_E,     "sc_clus1_E/F");
simu->Branch("sc_clus1_E_ov_p",     &sc_clus1_E_ov_p,     "sc_clus1_E_ov_p/F");
simu->Branch("sc_clus1_E_ov_E",     &sc_clus1_E_ov_E,     "sc_clus1_E_ov_E/F");
simu->Branch("sc_clus1_eta",   &sc_clus1_eta,   "sc_clus1_eta/F");
simu->Branch("sc_clus1_phi",   &sc_clus1_phi,   "sc_clus1_phi/F");
simu->Branch("sc_clus1_nxtal", &sc_clus1_nxtal, "sc_clus1_nxtal/I");
simu->Branch("sc_clus1_dphi",  &sc_clus1_dphi,  "sc_clus1_dphi/F");
simu->Branch("sc_clus1_deta",  &sc_clus1_deta,  "sc_clus1_deta/F");
simu->Branch("sc_clus1_ntrk_deta01",  &sc_clus1_ntrk_deta01,  "sc_clus1_ntrk_deta01/F");
simu->Branch("sc_clus2_E",     &sc_clus2_E,    "sc_clus2_E/F");
simu->Branch("sc_clus2_E_ov_p",     &sc_clus2_E_ov_p,     "sc_clus2_E_ov_p/F");
simu->Branch("sc_clus2_E_ov_E",     &sc_clus2_E_ov_E,     "sc_clus2_E_ov_E/F");
simu->Branch("sc_clus2_eta",   &sc_clus2_eta,  "sc_clus2_eta/F");
simu->Branch("sc_clus2_phi",   &sc_clus2_phi,  "sc_clus2_phi/F");
simu->Branch("sc_clus2_nxtal", &sc_clus2_nxtal, "sc_clus2_nxtal/I");
simu->Branch("sc_clus2_dphi",  &sc_clus2_dphi,  "sc_clus2_dphi/F");
simu->Branch("sc_clus2_deta",  &sc_clus2_deta,  "sc_clus2_deta/F");
simu->Branch("sc_clus2_ntrk_deta01",  &sc_clus2_ntrk_deta01,  "sc_clus2_ntrk_deta01/F");
simu->Branch("sc_clus3_E",     &sc_clus3_E,    "sc_clus3_E/F");
simu->Branch("sc_clus3_E_ov_p",     &sc_clus3_E_ov_p,     "sc_clus3_E_ov_p/F");
simu->Branch("sc_clus3_E_ov_E",     &sc_clus3_E_ov_E,     "sc_clus3_E_ov_E/F");
simu->Branch("sc_clus3_eta",   &sc_clus3_eta,  "sc_clus3_eta/F");
simu->Branch("sc_clus3_phi",   &sc_clus3_phi,  "sc_clus3_phi/F");
simu->Branch("sc_clus3_nxtal", &sc_clus3_nxtal, "sc_clus3_nxtal/I");
simu->Branch("sc_clus3_dphi",  &sc_clus3_dphi,  "sc_clus3_dphi/F");
simu->Branch("sc_clus3_deta",  &sc_clus3_deta,  "sc_clus3_deta/F");
simu->Branch("sc_clus3_ntrk_deta01",  &sc_clus3_ntrk_deta01,  "sc_clus3_ntrk_deta01/F");
simu->Branch("type",   &type,   "type/I");

//
 const int MaxEle=20; 
 int is_egammaV[MaxEle]={0};
 int is_eV[MaxEle]={0};
 float gsf_bdtout1V[MaxEle]={0};
 float gen_ptV[MaxEle]={0}; 
 float ele_mva_valueV[MaxEle]={0};
 float ele_mva_idV[MaxEle]={0};
 float gsf_mode_pV[MaxEle]={0};
 float eid_rhoV[MaxEle]={0};
 float eid_ele_ptV[MaxEle]={0};
 float eid_sc_etaV[MaxEle]={0};
 float eid_shape_full5x5_sigmaIetaIetaV[MaxEle]={0};
 float eid_shape_full5x5_sigmaIphiIphiV[MaxEle]={0};
 float eid_shape_full5x5_circularityV[MaxEle]={0};
 float eid_shape_full5x5_r9V[MaxEle]={0};
 float eid_sc_etaWidthV[MaxEle]={0};
 float eid_sc_phiWidthV[MaxEle]={0};
 float eid_shape_full5x5_HoverEV[MaxEle]={0};
 float eid_trk_nhitsV[MaxEle]={0};
 float eid_trk_chi2redV[MaxEle]={0};
 float eid_gsf_chi2redV[MaxEle]={0};
 float eid_brem_fracV[MaxEle]={0};
 float eid_gsf_nhitsV[MaxEle]={0};
 float eid_match_SC_EoverPV[MaxEle]={0};
 float eid_match_eclu_EoverPV[MaxEle]={0};
 float eid_match_SC_dEtaV[MaxEle]={0};
 float eid_match_SC_dPhiV[MaxEle]={0};
 float eid_match_seed_dEtaV[MaxEle]={0};
 float eid_sc_EV[MaxEle]={0};
 float eid_trk_pV[MaxEle]={0};
 float sc_goodSeedV[MaxEle]={0};
 float core_shFracHitsV[MaxEle]={0};
 float gsf_drV[MaxEle]={0};
 float trk_drV[MaxEle]={0};
 float sc_NclusV[MaxEle]={0};
 int sc_clus1_nxtalV[MaxEle]={0};
 int sc_clus2_nxtalV[MaxEle]={0};
 int sc_clus3_nxtalV[MaxEle]={0};
 float sc_clus1_dphiV[MaxEle]={0};
 float sc_clus2_dphiV[MaxEle]={0};
 float sc_clus3_dphiV[MaxEle]={0};
 float sc_clus1_detaV[MaxEle]={0};
 float sc_clus2_detaV[MaxEle]={0};
 float sc_clus3_detaV[MaxEle]={0};
 float sc_clus1_EV[MaxEle]={0};
 float sc_clus2_EV[MaxEle]={0};
 float sc_clus3_EV[MaxEle]={0};
 int sc_clus1_ntrk_deta01V[MaxEle]={0};
 int sc_clus2_ntrk_deta01V[MaxEle]={0};
 int sc_clus3_ntrk_deta01V[MaxEle]={0};
 //new
 int evtV[MaxEle]={0};
 float weightV[MaxEle]={0};
 int is_e_not_matchedV[MaxEle]={0}; 
 int is_otherV[MaxEle]={0}; 
 int has_trkV[MaxEle]={0}; 
 int has_seedV[MaxEle]={0}; 
 int has_gsfV[MaxEle]={0}; 
 int has_eleV[MaxEle]={0}; 
 float gen_etaV[MaxEle]={0}; 
 float gen_tag_sideV[MaxEle]={0};
 float gen_dRV[MaxEle]={0};
 float gsf_ptV[MaxEle]={0};
 float gsf_etaV[MaxEle]={0};
 float trk_ptV[MaxEle]={0};
 float trk_etaV[MaxEle]={0};
 float ele_pV[MaxEle]={0};
 float ele_ptV[MaxEle]={0};
 float ele_etaV[MaxEle]={0};
 int   fiducial_isEBV[MaxEle]={0};
 int   fiducial_isEEV[MaxEle]={0};
 int   fiducial_isEBEEGapV[MaxEle]={0};
 float brem_fracTrkV[MaxEle]={0}; 
 float brem_fracSCV[MaxEle]={0}; 
 float brem_NV[MaxEle]={0}; 
 float sc_clus1_E_ov_pV[MaxEle]={0};
 float sc_clus2_E_ov_pV[MaxEle]={0};
 float sc_clus3_E_ov_pV[MaxEle]={0};
 float sc_clus1_E_ov_EV[MaxEle]={0};
 float sc_clus2_E_ov_EV[MaxEle]={0};
 float sc_clus3_E_ov_EV[MaxEle]={0};
 float sc_clus1_etaV[MaxEle]={0};
 float sc_clus2_etaV[MaxEle]={0};
 float sc_clus3_etaV[MaxEle]={0};
 float sc_clus1_phiV[MaxEle]={0};
 float sc_clus2_phiV[MaxEle]={0};
 float sc_clus3_phiV[MaxEle]={0};
 int sc_Nclus_deta01V[MaxEle]={0};
 int sc_Nclus_deta02V[MaxEle]={0};
 int sc_Nclus_deta03V[MaxEle]={0};


 //
 int typeV[MaxEle]={0};



   type = 1;

   int nloop=fChain->GetEntries();
   int nmax=200000000;
   if(nloop > nmax) nloop=nmax; 

   uint evt_old=0; 
   
   int ele_counter=0; 

   TH1F *bg_pt = new TH1F("bg_pt", "bg_pt", 100, 0., 20.);
   TH1F *sig_pt = new TH1F("sig_pt", "sig_pt", 100, 0., 20.);
   bg_pt->SetDirectory(0);
   sig_pt->SetDirectory(0);


   for (int i = 0; i < nloop; i++) {
     fChain->GetEntry(i);

     if(evt_old==0) evt_old=evt; // first entry 

     if(evt_old!=evt) {
       // fill var for old event                                                                                                                           
       int nele_eg0=0;
       int nele_eg1=0;
       int nele_eg2=0;
       int nele_eg0f=0;
       int nele_eg1f=0;
       int nele_eg2f=0;


       if(ele_counter>0){
	 for (int j=0; j<ele_counter; j++){

	   if(is_egammaV[j]==0 ){  // x CHIARA: only LP are kept,  comment this line if you want to keep all electrons 

	     // if(typeV[j]==1 || (typeV[j]==0 && j%20==0 )) { // x CHIARA: change here if you want to save 1/N LP background entries  
	     if(typeV[j]==1 || (typeV[j]==0 )) {               // x CHIARA: use this line to save all LP background entries
		is_egamma=is_egammaV[j];
		is_e=is_eV[j];
		gen_pt=gen_ptV[j];
		ele_mva_value=ele_mva_valueV[j];
		ele_mva_id=ele_mva_idV[j];
		gsf_mode_p=gsf_mode_pV[j];
		eid_rho=eid_rhoV[j];
		eid_ele_pt=eid_ele_ptV[j];
		eid_sc_eta=eid_sc_etaV[j];
		eid_shape_full5x5_sigmaIetaIeta=eid_shape_full5x5_sigmaIetaIetaV[j];
		eid_shape_full5x5_sigmaIphiIphi=eid_shape_full5x5_sigmaIphiIphiV[j];
		eid_shape_full5x5_circularity=eid_shape_full5x5_circularityV[j];
		eid_shape_full5x5_r9=eid_shape_full5x5_r9V[j];
		eid_sc_etaWidth=eid_sc_etaWidthV[j];
		eid_sc_phiWidth=eid_sc_phiWidthV[j];
		eid_shape_full5x5_HoverE=eid_shape_full5x5_HoverEV[j];
		eid_trk_nhits=eid_trk_nhitsV[j];
		eid_trk_chi2red=eid_trk_chi2redV[j];
		eid_gsf_chi2red=eid_gsf_chi2redV[j];
		eid_brem_frac=eid_brem_fracV[j];
		eid_gsf_nhits=eid_gsf_nhitsV[j];
		eid_match_SC_EoverP=eid_match_SC_EoverPV[j];
		eid_match_eclu_EoverP=eid_match_eclu_EoverPV[j];
		eid_match_SC_dEta=eid_match_SC_dEtaV[j];
		eid_match_SC_dPhi=eid_match_SC_dPhiV[j];
		eid_match_seed_dEta=eid_match_seed_dEtaV[j];
		eid_sc_E=eid_sc_EV[j];
		eid_trk_p=eid_trk_pV[j];
		gsf_bdtout1=gsf_bdtout1V[j];
		sc_goodSeed=sc_goodSeedV[j];
		core_shFracHits=core_shFracHitsV[j];
		gsf_dr=gsf_drV[j];
		trk_dr=trk_drV[j];
		sc_Nclus=sc_NclusV[j];
		sc_clus1_nxtal=sc_clus1_nxtalV[j];
		sc_clus2_nxtal=sc_clus2_nxtalV[j];
		sc_clus3_nxtal=sc_clus3_nxtalV[j];
		sc_clus1_dphi=sc_clus1_dphiV[j];
		sc_clus2_dphi=sc_clus2_dphiV[j];
		sc_clus3_dphi=sc_clus3_dphiV[j];
		sc_clus1_deta=sc_clus1_detaV[j];
		sc_clus2_deta=sc_clus2_detaV[j];
		sc_clus3_deta=sc_clus3_detaV[j];
		sc_clus1_E=sc_clus1_EV[j];
		sc_clus2_E=sc_clus2_EV[j];
		sc_clus3_E=sc_clus3_EV[j];
		sc_clus1_ntrk_deta01=sc_clus1_ntrk_deta01V[j];
		sc_clus2_ntrk_deta01=sc_clus2_ntrk_deta01V[j];
		sc_clus3_ntrk_deta01=sc_clus3_ntrk_deta01V[j];
		type=typeV[j];

 //new
		evt=evtV[j];			 
		weight=weightV[j];			 
		is_e_not_matched=is_e_not_matchedV[j]; 		 
		is_other=is_otherV[j]; 			 
		has_trk=has_trkV[j]; 			 
		has_seed=has_seedV[j]; 			 
		has_gsf=has_gsfV[j];
		has_ele=has_eleV[j];
		gen_eta=gen_etaV[j];
		gen_tag_side=gen_tag_sideV[j];
		gen_dR=gen_dRV[j];
		gsf_pt=gsf_ptV[j];
		gsf_eta=gsf_etaV[j];
		trk_pt=trk_ptV[j];
		trk_eta=trk_etaV[j];
		ele_p=ele_pV[j];
		ele_pt=ele_ptV[j];
		ele_eta=ele_etaV[j];
		fiducial_isEB=fiducial_isEBV[j];
		fiducial_isEE=fiducial_isEEV[j];
		fiducial_isEBEEGap=fiducial_isEBEEGapV[j];
		brem_fracTrk=brem_fracTrkV[j];
		brem_fracSC=brem_fracSCV[j];
		brem_N=brem_NV[j];
		sc_clus1_E_ov_p=sc_clus1_E_ov_pV[j];
		sc_clus2_E_ov_p=sc_clus2_E_ov_pV[j];
		sc_clus3_E_ov_p=sc_clus3_E_ov_pV[j];
		sc_clus1_E_ov_E=sc_clus1_E_ov_EV[j];
		sc_clus2_E_ov_E=sc_clus2_E_ov_EV[j];
		sc_clus3_E_ov_E=sc_clus3_E_ov_EV[j];
		sc_clus1_eta=sc_clus1_etaV[j];
		sc_clus2_eta=sc_clus2_etaV[j];
		sc_clus3_eta=sc_clus3_etaV[j];
		sc_clus1_phi=sc_clus1_phiV[j];
		sc_clus2_phi=sc_clus2_phiV[j];
		sc_clus3_phi=sc_clus3_phiV[j];
		sc_Nclus_deta01=sc_Nclus_deta01V[j];
		sc_Nclus_deta02=sc_Nclus_deta02V[j];
		sc_Nclus_deta03=sc_Nclus_deta03V[j];    

		evt=evt_old;
		simu->Fill();
		if(typeV[j]==1){
		  nele_eg0=nele_eg0+1;
		  sig_pt->Fill(gsf_mode_pV[j]);
		}
		if(typeV[j]==0){
		  nele_eg0f=nele_eg0f+1;
		  bg_pt->Fill(gsf_mode_pV[j]);
		}

	      }
	   }
	   if(is_egammaV[j]==1){
	     if(typeV[j]==1)nele_eg1=nele_eg1+1;
	     if(typeV[j]==0)nele_eg1f=nele_eg1f+1;
	   } else if(is_egammaV[j]==2){
	     if(typeV[j]==1)nele_eg2=nele_eg2+1;
	     if(typeV[j]==0)nele_eg2f=nele_eg2f+1;
	   }


	 }
	 
       }
       
       // reset everything 
       ele_counter=0; 
       for (int j=0; j<MaxEle; j++){
	 is_egammaV[j]=0;
	 is_eV[j]=0;
	 gen_ptV[j]=0;
	 ele_mva_valueV[j]=0;
	 ele_mva_idV[j]=0;
	 gsf_mode_pV[j]=0;
	 eid_rhoV[j]=0;
	 eid_ele_ptV[j]=0;
	 eid_sc_etaV[j]=0;
	 eid_shape_full5x5_sigmaIetaIetaV[j]=0;
	 eid_shape_full5x5_sigmaIphiIphiV[j]=0;
	 eid_shape_full5x5_circularityV[j]=0;
	 eid_shape_full5x5_r9V[j]=0;
	 eid_sc_etaWidthV[j]=0;
	 eid_sc_phiWidthV[j]=0;
	 eid_shape_full5x5_HoverEV[j]=0;
	 eid_trk_nhitsV[j]=0;
	 eid_trk_chi2redV[j]=0;
	 eid_gsf_chi2redV[j]=0;
	 eid_brem_fracV[j]=0;
	 eid_gsf_nhitsV[j]=0;
	 eid_match_SC_EoverPV[j]=0;
	 eid_match_eclu_EoverPV[j]=0;
	 eid_match_SC_dEtaV[j]=0;
	 eid_match_SC_dPhiV[j]=0;
	 eid_match_seed_dEtaV[j]=0;
	 eid_sc_EV[j]=0;
	 eid_trk_pV[j]=0;
	 gsf_bdtout1V[j]=0;
	 sc_goodSeedV[j]=0;
	 core_shFracHitsV[j]=0;
	 gsf_drV[j]=0;
	 trk_drV[j]=0;
	 sc_NclusV[j]=0;
	 sc_clus1_nxtalV[j]=0;
	 sc_clus2_nxtalV[j]=0;
	 sc_clus3_nxtalV[j]=0;
	 sc_clus1_dphiV[j]=0;
	 sc_clus2_dphiV[j]=0;
	 sc_clus3_dphiV[j]=0;
	 sc_clus1_detaV[j]=0;
	 sc_clus2_detaV[j]=0;
	 sc_clus3_detaV[j]=0;
	 sc_clus1_EV[j]=0;
	 sc_clus2_EV[j]=0;
	 sc_clus3_EV[j]=0;
	 sc_clus1_ntrk_deta01V[j]=0;
	 sc_clus2_ntrk_deta01V[j]=0;
	 sc_clus3_ntrk_deta01V[j]=0;
	 typeV[j]=0;
	 // new
  evtV[j]=0;
  weightV[j]=0;
  is_e_not_matchedV[j]=0; 
  is_otherV[j]=0; 
  has_trkV[j]=0; 
  has_seedV[j]=0; 
  has_gsfV[j]=0; 
  has_eleV[j]=0; 
  gen_etaV[j]=0; 
  gen_tag_sideV[j]=0;
  gen_dRV[j]=0;
  gsf_ptV[j]=0;
  gsf_etaV[j]=0;
  trk_ptV[j]=0;
  trk_etaV[j]=0;
  ele_pV[j]=0;
  ele_ptV[j]=0;
  ele_etaV[j]=0;
    fiducial_isEBV[j]=0;
    fiducial_isEEV[j]=0;
    fiducial_isEBEEGapV[j]=0;
  brem_fracTrkV[j]=0; 
  brem_fracSCV[j]=0; 
  brem_NV[j]=0; 
  sc_clus1_E_ov_pV[j]=0;
  sc_clus2_E_ov_pV[j]=0;
  sc_clus3_E_ov_pV[j]=0;
  sc_clus1_E_ov_EV[j]=0;
  sc_clus2_E_ov_EV[j]=0;
  sc_clus3_E_ov_EV[j]=0;
  sc_clus1_etaV[j]=0;
  sc_clus2_etaV[j]=0;
  sc_clus3_etaV[j]=0;
  sc_clus1_phiV[j]=0;
  sc_clus2_phiV[j]=0;
  sc_clus3_phiV[j]=0;
  sc_Nclus_deta01V[j]=0;
  sc_Nclus_deta02V[j]=0;
  sc_Nclus_deta03V[j]=0;

	 evt=evt_old;

       }

       //       std::cout<<"Event "<<evt_old<<" "<<nele_eg0<<" "<<nele_eg0f<<" "<<nele_eg1<<" "<<nele_eg1f<<" "<<nele_eg2<<" "<<nele_eg2f<<std::endl;



       // reload the event i 
       fChain->GetEntry(i);
       evt_old=evt; 
     } 

     // process a normal entry
   
     if(evt_old==evt ){


       //if (1) { // chiaraaaa
       if(is_e==1&& gen_pt>0 && TMath::Abs(gen_eta)<2.5 &&gen_dR<0.03 ){
       // is matched to a good MC electron 

       bool pcheck=false; // is already in the list ? 
       int    jele=-1;
       for(int j=0; j<ele_counter; j++){
	 if(gen_ptV[j]>0&&TMath::Abs(gen_pt-gen_ptV[j])<0.001){
	   pcheck=true;
	   jele=j;
	 }
       }
       if(!pcheck&&ele_counter<MaxEle){
	 // ele new  
	 is_egammaV[ele_counter]=is_egamma;
	 is_eV[ele_counter]=is_e;
	 gen_ptV[ele_counter]=gen_pt;
	 ele_mva_valueV[ele_counter]=ele_mva_value;
	 ele_mva_idV[ele_counter]=ele_mva_id;
	 gsf_mode_pV[ele_counter]=gsf_mode_p;
	 eid_rhoV[ele_counter]=eid_rho;
	 eid_ele_ptV[ele_counter]=eid_ele_pt;
	 eid_sc_etaV[ele_counter]=eid_sc_eta;
	 eid_shape_full5x5_sigmaIetaIetaV[ele_counter]=eid_shape_full5x5_sigmaIetaIeta;
	 eid_shape_full5x5_sigmaIphiIphiV[ele_counter]=eid_shape_full5x5_sigmaIphiIphi;
	 eid_shape_full5x5_circularityV[ele_counter]=eid_shape_full5x5_circularity;
	 eid_shape_full5x5_r9V[ele_counter]=eid_shape_full5x5_r9;
	 eid_sc_etaWidthV[ele_counter]=eid_sc_etaWidth;
	 eid_sc_phiWidthV[ele_counter]=eid_sc_phiWidth;
	 eid_shape_full5x5_HoverEV[ele_counter]=eid_shape_full5x5_HoverE;
	 eid_trk_nhitsV[ele_counter]=eid_trk_nhits;
	 eid_trk_chi2redV[ele_counter]=eid_trk_chi2red;
	 eid_gsf_chi2redV[ele_counter]=eid_gsf_chi2red;
	 eid_brem_fracV[ele_counter]=eid_brem_frac;
	 eid_gsf_nhitsV[ele_counter]=eid_gsf_nhits;
	 eid_match_SC_EoverPV[ele_counter]=eid_match_SC_EoverP;
	 eid_match_eclu_EoverPV[ele_counter]=eid_match_eclu_EoverP;
	 eid_match_SC_dEtaV[ele_counter]=eid_match_SC_dEta;
	 eid_match_SC_dPhiV[ele_counter]=eid_match_SC_dPhi;
	 eid_match_seed_dEtaV[ele_counter]=eid_match_seed_dEta;
	 eid_sc_EV[ele_counter]=eid_sc_E;
	 eid_trk_pV[ele_counter]=eid_trk_p;
	 gsf_bdtout1V[ele_counter]=gsf_bdtout1;
	 sc_goodSeedV[ele_counter]=sc_goodSeed;
	 core_shFracHitsV[ele_counter]=core_shFracHits;
	 gsf_drV[ele_counter]=gsf_dr;
	 trk_drV[ele_counter]=trk_dr;
	 sc_NclusV[ele_counter]=sc_Nclus;
	 sc_clus1_nxtalV[ele_counter]=sc_clus1_nxtal;
	 sc_clus2_nxtalV[ele_counter]=sc_clus2_nxtal;
	 sc_clus3_nxtalV[ele_counter]=sc_clus3_nxtal;
	 sc_clus1_dphiV[ele_counter]=sc_clus1_dphi;
	 sc_clus2_dphiV[ele_counter]=sc_clus2_dphi;
	 sc_clus3_dphiV[ele_counter]=sc_clus3_dphi;
	 sc_clus1_detaV[ele_counter]=sc_clus1_deta;
	 sc_clus2_detaV[ele_counter]=sc_clus2_deta;
	 sc_clus3_detaV[ele_counter]=sc_clus3_deta;
	 sc_clus1_EV[ele_counter]=sc_clus1_E;
	 sc_clus2_EV[ele_counter]=sc_clus2_E;
	 sc_clus3_EV[ele_counter]=sc_clus3_E;
	 sc_clus1_ntrk_deta01V[ele_counter]=sc_clus1_ntrk_deta01;
	 sc_clus2_ntrk_deta01V[ele_counter]=sc_clus2_ntrk_deta01;
	 sc_clus3_ntrk_deta01V[ele_counter]=sc_clus3_ntrk_deta01;
	 typeV[ele_counter]=1; 
 //new
	 evtV[ele_counter]=evt;			 
	 weightV[ele_counter]=weight;			 
	 is_e_not_matchedV[ele_counter]=is_e_not_matched; 		 
	 is_otherV[ele_counter]=is_other; 			 
	 has_trkV[ele_counter]=has_trk; 			 
	 has_seedV[ele_counter]=has_seed; 			 
	 has_gsfV[ele_counter]=has_gsf;
	 has_eleV[ele_counter]=has_ele;
	 gen_etaV[ele_counter]=gen_eta;
	 gen_tag_sideV[ele_counter]=gen_tag_side;
	 gen_dRV[ele_counter]=gen_dR;
	 gsf_ptV[ele_counter]=gsf_pt;
	 gsf_etaV[ele_counter]=gsf_eta;
	 trk_ptV[ele_counter]=trk_pt;
	 trk_etaV[ele_counter]=trk_eta;
	 ele_pV[ele_counter]=ele_p;
	 ele_ptV[ele_counter]=ele_pt;
	 ele_etaV[ele_counter]=ele_eta;
	 fiducial_isEBV[ele_counter]=fiducial_isEB;
	 fiducial_isEEV[ele_counter]=fiducial_isEE;
	 fiducial_isEBEEGapV[ele_counter]=fiducial_isEBEEGap;
	 brem_fracTrkV[ele_counter]=brem_fracTrk;
	 brem_fracSCV[ele_counter]=brem_fracSC;
	 brem_NV[ele_counter]=brem_N;
	 sc_clus1_E_ov_pV[ele_counter]=sc_clus1_E_ov_p;
	 sc_clus2_E_ov_pV[ele_counter]=sc_clus2_E_ov_p;
	 sc_clus3_E_ov_pV[ele_counter]=sc_clus3_E_ov_p;
	 sc_clus1_E_ov_EV[ele_counter]=sc_clus1_E_ov_E;
	 sc_clus2_E_ov_EV[ele_counter]=sc_clus2_E_ov_E;
	 sc_clus3_E_ov_EV[ele_counter]=sc_clus3_E_ov_E;
	 sc_clus1_etaV[ele_counter]=sc_clus1_eta;
	 sc_clus2_etaV[ele_counter]=sc_clus2_eta;
	 sc_clus3_etaV[ele_counter]=sc_clus3_eta;
	 sc_clus1_phiV[ele_counter]=sc_clus1_phi;
	 sc_clus2_phiV[ele_counter]=sc_clus2_phi;
	 sc_clus3_phiV[ele_counter]=sc_clus3_phi;
	 sc_Nclus_deta01V[ele_counter]=sc_Nclus_deta01;
	 sc_Nclus_deta02V[ele_counter]=sc_Nclus_deta02;
	 sc_Nclus_deta03V[ele_counter]=sc_Nclus_deta03;

	 ele_counter=ele_counter+1;

       } else if(pcheck){ 
	 is_egammaV[jele]=2; // matched to both LP ele and GSF, data are for the LP 
       }
       
       } else if(is_e==0) {
       // background 

	 bool pcheck=false; // is already in the list ? 
	 int jele=-1;
	 if(is_egamma==1){ // is a fake GSF
	   for(int j=0; j<ele_counter; j++){
	     if(is_egammaV[j]==0){  // loop on LP only 
	       if(TMath::Abs(eid_trk_pV[j]-eid_trk_p)<0.0001){ // use the same trk: GSF matched to LP
		 pcheck=true;
		 jele=j;
	       }
	     }
	   }
	 }
	 if(!pcheck&&ele_counter<MaxEle){
	   // fake ele new  
	   is_egammaV[ele_counter]=is_egamma;
	   is_eV[ele_counter]=is_e;
	   gen_ptV[ele_counter]=gen_pt;
	   ele_mva_valueV[ele_counter]=ele_mva_value;
	   ele_mva_idV[ele_counter]=ele_mva_id;
	   gsf_mode_pV[ele_counter]=gsf_mode_p;
	   eid_rhoV[ele_counter]=eid_rho;
	   eid_ele_ptV[ele_counter]=eid_ele_pt;
	   eid_sc_etaV[ele_counter]=eid_sc_eta;
	   eid_shape_full5x5_sigmaIetaIetaV[ele_counter]=eid_shape_full5x5_sigmaIetaIeta;
	   eid_shape_full5x5_sigmaIphiIphiV[ele_counter]=eid_shape_full5x5_sigmaIphiIphi;
	   eid_shape_full5x5_circularityV[ele_counter]=eid_shape_full5x5_circularity;
	   eid_shape_full5x5_r9V[ele_counter]=eid_shape_full5x5_r9;
	   eid_sc_etaWidthV[ele_counter]=eid_sc_etaWidth;
	   eid_sc_phiWidthV[ele_counter]=eid_sc_phiWidth;
	   eid_shape_full5x5_HoverEV[ele_counter]=eid_shape_full5x5_HoverE;
	   eid_trk_nhitsV[ele_counter]=eid_trk_nhits;
	   eid_trk_chi2redV[ele_counter]=eid_trk_chi2red;
	   eid_gsf_chi2redV[ele_counter]=eid_gsf_chi2red;
	   eid_brem_fracV[ele_counter]=eid_brem_frac;
	   eid_gsf_nhitsV[ele_counter]=eid_gsf_nhits;
	   eid_match_SC_EoverPV[ele_counter]=eid_match_SC_EoverP;
	   eid_match_eclu_EoverPV[ele_counter]=eid_match_eclu_EoverP;
	   eid_match_SC_dEtaV[ele_counter]=eid_match_SC_dEta;
	   eid_match_SC_dPhiV[ele_counter]=eid_match_SC_dPhi;
	   eid_match_seed_dEtaV[ele_counter]=eid_match_seed_dEta;
	   eid_sc_EV[ele_counter]=eid_sc_E;
	   eid_trk_pV[ele_counter]=eid_trk_p;
	   gsf_bdtout1V[ele_counter]=gsf_bdtout1;
	   sc_goodSeedV[ele_counter]=sc_goodSeed;
	   core_shFracHitsV[ele_counter]=core_shFracHits;
	   gsf_drV[ele_counter]=gsf_dr;
	   trk_drV[ele_counter]=trk_dr;
	   sc_NclusV[ele_counter]=sc_Nclus;
	   sc_clus1_nxtalV[ele_counter]=sc_clus1_nxtal;
	   sc_clus2_nxtalV[ele_counter]=sc_clus2_nxtal;
	   sc_clus3_nxtalV[ele_counter]=sc_clus3_nxtal;
	   sc_clus1_dphiV[ele_counter]=sc_clus1_dphi;
	   sc_clus2_dphiV[ele_counter]=sc_clus2_dphi;
	   sc_clus3_dphiV[ele_counter]=sc_clus3_dphi;
	   sc_clus1_detaV[ele_counter]=sc_clus1_deta;
	   sc_clus2_detaV[ele_counter]=sc_clus2_deta;
	   sc_clus3_detaV[ele_counter]=sc_clus3_deta;
	   sc_clus1_EV[ele_counter]=sc_clus1_E;
	   sc_clus2_EV[ele_counter]=sc_clus2_E;
	   sc_clus3_EV[ele_counter]=sc_clus3_E;
	   sc_clus1_ntrk_deta01V[ele_counter]=sc_clus1_ntrk_deta01;
	   sc_clus2_ntrk_deta01V[ele_counter]=sc_clus2_ntrk_deta01;
	   sc_clus3_ntrk_deta01V[ele_counter]=sc_clus3_ntrk_deta01;
	   typeV[ele_counter]=0; // is background
//new
	 evtV[ele_counter]=evt;			 
	 weightV[ele_counter]=weight;			 
	 is_e_not_matchedV[ele_counter]=is_e_not_matched; 		 
	 is_otherV[ele_counter]=is_other; 			 
	 has_trkV[ele_counter]=has_trk; 			 
	 has_seedV[ele_counter]=has_seed; 			 
	 has_gsfV[ele_counter]=has_gsf;
	 has_eleV[ele_counter]=has_ele;
	 gen_etaV[ele_counter]=gen_eta;
	 gen_tag_sideV[ele_counter]=gen_tag_side;
	 gen_dRV[ele_counter]=gen_dR;
	 gsf_ptV[ele_counter]=gsf_pt;
	 gsf_etaV[ele_counter]=gsf_eta;
	 trk_ptV[ele_counter]=trk_pt;
	 trk_etaV[ele_counter]=trk_eta;
	 ele_pV[ele_counter]=ele_p;
	 ele_ptV[ele_counter]=ele_pt;
	 ele_etaV[ele_counter]=ele_eta;
	 fiducial_isEBV[ele_counter]=fiducial_isEB;
	 fiducial_isEEV[ele_counter]=fiducial_isEE;
	 fiducial_isEBEEGapV[ele_counter]=fiducial_isEBEEGap;
	 brem_fracTrkV[ele_counter]=brem_fracTrk;
	 brem_fracSCV[ele_counter]=brem_fracSC;
	 brem_NV[ele_counter]=brem_N;
	 sc_clus1_E_ov_pV[ele_counter]=sc_clus1_E_ov_p;
	 sc_clus2_E_ov_pV[ele_counter]=sc_clus2_E_ov_p;
	 sc_clus3_E_ov_pV[ele_counter]=sc_clus3_E_ov_p;
	 sc_clus1_E_ov_EV[ele_counter]=sc_clus1_E_ov_E;
	 sc_clus2_E_ov_EV[ele_counter]=sc_clus2_E_ov_E;
	 sc_clus3_E_ov_EV[ele_counter]=sc_clus3_E_ov_E;
	 sc_clus1_etaV[ele_counter]=sc_clus1_eta;
	 sc_clus2_etaV[ele_counter]=sc_clus2_eta;
	 sc_clus3_etaV[ele_counter]=sc_clus3_eta;
	 sc_clus1_phiV[ele_counter]=sc_clus1_phi;
	 sc_clus2_phiV[ele_counter]=sc_clus2_phi;
	 sc_clus3_phiV[ele_counter]=sc_clus3_phi;
	 sc_Nclus_deta01V[ele_counter]=sc_Nclus_deta01;
	 sc_Nclus_deta02V[ele_counter]=sc_Nclus_deta02;
	 sc_Nclus_deta03V[ele_counter]=sc_Nclus_deta03;
	   ele_counter=ele_counter+1;
	 } else if(pcheck){ 
	   if(is_egammaV[jele]==0)is_egammaV[jele]=2; // matched to both LP ele and GSF 
	 }
       }

     } 

     

   }

   delete input;
   simu->Write();
   delete fileNew;
  std::cout<<"hello ... ending shortTrain"<<std::endl; 

}
int main(int argc, char **argv){

  shortTrainNNLP t;
  t.train(100);

  return 0;

}
