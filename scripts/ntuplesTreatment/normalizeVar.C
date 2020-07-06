#define normalizeVar_cxx
#include "normalizeVar.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <stdio.h>
using namespace std; 

void normalizeVar::train() {

  const char *origfname = "/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/DoubleElectronGun_withRegression.root";
  
  TFile *origInput = 0;
  origInput = TFile::Open(origfname);
  if (!origInput) return;
  std::cout << "just opened your file" << std::endl;

  TDirectory *dir = (TDirectory*)origInput->Get("/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/DoubleElectronGun_withRegression.root:/ntuplizer");
  TTree *fChain ; 
  dir->GetObject("tree",fChain);
  std::cout << "just taken your tree" << std::endl;

  // branches
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
  fChain->SetBranchAddress("pre_ecal",&pre_ecal);
  fChain->SetBranchAddress("pre_ecaltrk",&pre_ecaltrk);
  fChain->SetBranchAddress("post_ecal",&post_ecal);
  fChain->SetBranchAddress("post_ecaltrk",&post_ecaltrk);
  
  // New tree
  TFile *fileNew = TFile::Open("/tmp/simu_shortLP.root","RECREATE");     
  fileNew->ls();
  fileNew->cd();
  TDirectory *myDir = (TDirectory*)fileNew->mkdir("ntuplizer");
  myDir->cd();
  TTree *simu = new TTree("tree", "Normalized Monte Carlo Events");  

  // new tree branches   
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
  simu->Branch("pre_ecal", &pre_ecal, "pre_ecal/F");
  simu->Branch("pre_ecaltrk", &pre_ecaltrk, "pre_ecaltrk/F");
  simu->Branch("post_ecal", &post_ecal, "post_ecal/F");
  simu->Branch("post_ecaltrk", &post_ecaltrk, "post_ecaltrk/F");

  simu->Branch("type",   &type,   "type/I");

  // Loop over entries
  int nloop=fChain->GetEntries();
  for (int i = 0; i < nloop; i++) {

    if(i%100000==0)   std::cout<<"processing entry ... "<<i<<std::endl;
    fChain->GetEntry(i);
    
    // check variables and normalize outlayers
    normalize();
	      
    simu->Fill();
  }

  simu->Write();
  delete fileNew;
  std::cout<<"hello ... ending shortTrain"<<std::endl; 
}

void normalizeVar::normalize(){

  if(gsf_mode_p<0) gsf_mode_p=0;
  if(gsf_mode_p>1000) gsf_mode_p=1000;
  if(eid_rho<0) eid_rho=0;
  if(eid_rho>100) eid_rho=100;
  if(eid_ele_pt<0) eid_ele_pt=0;
  if(eid_ele_pt>1000) eid_ele_pt=1000;
  if(eid_sc_eta<-5) eid_sc_eta=-5;
  if(eid_sc_eta>5) eid_sc_eta=5;
  if(eid_shape_full5x5_sigmaIetaIeta<0)eid_shape_full5x5_sigmaIetaIeta=0;
  if(eid_shape_full5x5_sigmaIetaIeta>1)eid_shape_full5x5_sigmaIetaIeta=1;
  if(eid_shape_full5x5_sigmaIphiIphi<0) eid_shape_full5x5_sigmaIphiIphi=0;
  if(eid_shape_full5x5_sigmaIphiIphi>1) eid_shape_full5x5_sigmaIphiIphi=1;
  if(eid_shape_full5x5_circularity<0) eid_shape_full5x5_circularity=0;
  if(eid_shape_full5x5_circularity>1) eid_shape_full5x5_circularity=1;
  if(eid_shape_full5x5_r9<0) eid_shape_full5x5_r9=0;
  if(eid_shape_full5x5_r9>2) eid_shape_full5x5_r9=2;
  if(eid_sc_etaWidth<0) eid_sc_etaWidth=0;
  if(eid_sc_etaWidth>3.14) eid_sc_etaWidth=3.14;
  if(eid_sc_phiWidth<0) eid_sc_phiWidth=0;
  if(eid_sc_phiWidth>3.14) eid_sc_phiWidth=3.14;
  if(eid_shape_full5x5_HoverE<0) eid_shape_full5x5_HoverE=0;
  if(eid_shape_full5x5_HoverE>50) eid_shape_full5x5_HoverE=50;
  if(eid_trk_nhits<-1) eid_trk_nhits=-1;
  if(eid_trk_nhits>50) eid_trk_nhits=50;
  if(eid_trk_chi2red<-1) eid_trk_chi2red=-1;
  if(eid_trk_chi2red>50) eid_trk_chi2red=50;
  if(eid_gsf_chi2red<-1) eid_gsf_chi2red=-1;
  if(eid_gsf_chi2red>100) eid_gsf_chi2red=100;
  if(eid_brem_frac<0) eid_brem_frac=-1;
  if(eid_brem_frac>1) eid_brem_frac=1;
  if(eid_gsf_nhits<-1) eid_gsf_nhits=-1;
  if(eid_gsf_nhits>50) eid_gsf_nhits=50;
  if(eid_match_SC_EoverP<0) eid_match_SC_EoverP=0;
  if(eid_match_SC_EoverP>100) eid_match_SC_EoverP=100;
  if(eid_match_eclu_EoverP<-0.001) eid_match_eclu_EoverP=-0.001;
  if(eid_match_eclu_EoverP>0.001) eid_match_eclu_EoverP=0.001;
  eid_match_eclu_EoverP=eid_match_eclu_EoverP*1.E7;
  if(eid_match_SC_dEta<-10)eid_match_SC_dEta=-10;
  if(eid_match_SC_dEta>10)eid_match_SC_dEta=10;
  if(eid_match_SC_dPhi<-3.14)eid_match_SC_dPhi=-3.14;
  if(eid_match_SC_dPhi>3.14)eid_match_SC_dPhi=3.14;
  if(eid_match_seed_dEta<-10)eid_match_seed_dEta=-10;
  if(eid_match_seed_dEta>10)eid_match_seed_dEta=10;
  if(eid_sc_E<0) eid_sc_E=0;
  if(eid_sc_E>1000) eid_sc_E=1000;
  if(eid_trk_p<-1) eid_trk_p=-1;
  if(eid_trk_p>1000) eid_trk_p=1000;
  if(gsf_bdtout1<-20) gsf_bdtout1=-20;
  if(gsf_bdtout1>20) gsf_bdtout1=20;

  if(core_shFracHits<0) core_shFracHits=0;
  if(core_shFracHits>1) core_shFracHits=1;
  if(gsf_dr<0) gsf_dr=5;
  if(gsf_dr>5) gsf_dr=5;
  if(trk_dr<0) trk_dr=5;
  if(trk_dr>5) trk_dr=5;
  if(sc_Nclus<0) sc_Nclus=0;
  if(sc_Nclus>20) sc_Nclus=20;
  if(sc_clus1_nxtal<0) sc_clus1_nxtal=0;
  if(sc_clus2_nxtal<0) sc_clus2_nxtal=0;
  if(sc_clus3_nxtal<0) sc_clus3_nxtal=0;
  if(sc_clus1_nxtal>100) sc_clus1_nxtal=100;
  if(sc_clus2_nxtal>100) sc_clus2_nxtal=100;
  if(sc_clus3_nxtal>100) sc_clus3_nxtal=100;

  if(sc_clus1_dphi<-3.14) sc_clus1_dphi=-5;
  if(sc_clus1_dphi>3.14) sc_clus1_dphi=5;
  if(sc_clus2_dphi<-3.14) sc_clus2_dphi=-5;
  if(sc_clus2_dphi>3.14)  sc_clus2_dphi=5;
  if(sc_clus3_dphi<-3.14) sc_clus3_dphi=-5;
  if(sc_clus3_dphi>3.14)  sc_clus3_dphi=5;
  if(sc_clus1_deta<-5) sc_clus1_deta=-5;
  if(sc_clus1_deta>5)  sc_clus1_deta=5;
  if(sc_clus2_deta<-5) sc_clus2_deta=-5;
  if(sc_clus2_deta>5)  sc_clus2_deta=5;
  if(sc_clus3_deta<-5) sc_clus3_deta=-5;
  if(sc_clus3_deta>5)  sc_clus3_deta=5;
  if(sc_clus1_E<0)    sc_clus1_E=0;
  if(sc_clus1_E>1000) sc_clus1_E=1000;
  if(sc_clus2_E<0)    sc_clus2_E=0;
  if(sc_clus2_E>1000) sc_clus2_E=1000;
  if(sc_clus3_E<0)    sc_clus3_E=0;
  if(sc_clus3_E>1000) sc_clus3_E=1000;

  if(sc_clus1_ntrk_deta01<0|| sc_clus1_E==0) sc_clus1_ntrk_deta01=-1;
  if(sc_clus2_ntrk_deta01<0|| sc_clus2_E==0) sc_clus2_ntrk_deta01=-1;
  if(sc_clus3_ntrk_deta01<0|| sc_clus3_E==0) sc_clus3_ntrk_deta01=-1;

  if(sc_clus1_E_ov_p<0) sc_clus1_E_ov_p=-1;
  if(sc_clus2_E_ov_p<0) sc_clus2_E_ov_p=-1;
  if(sc_clus3_E_ov_p<0) sc_clus3_E_ov_p=-1;

  if(pre_ecal<0) pre_ecal=0;
  if(pre_ecal>1000) pre_ecal=1000;
  if(pre_ecaltrk<0) pre_ecaltrk=0;
  if(pre_ecaltrk>1000) pre_ecaltrk=1000;
  if(post_ecal<0) post_ecal=0;
  if(post_ecal>1000) post_ecal=1000;
  if(post_ecaltrk<0) post_ecaltrk=0;
  if(post_ecaltrk>1000) post_ecaltrk=1000;

}

int main(int argc, char **argv){

  normalizeVar t;
  t.train();

  return 0;

}
