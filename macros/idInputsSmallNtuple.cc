// C++ includes
#include <iostream>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>

using namespace std;

void idInputsSmallNtuple() {
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile *theFile = TFile::Open("/eos/cms/store/user/crovelli/LowPtEle/Batch2/BuToKJpsiToee_slim.root");
  TTree *theTree = (TTree*)theFile->Get("ntuplizer/tree");  

  // Variables
  const int nVar_trk  = 3;
  const int nVar_gsf  = 6;
  const int nVar_ele  = 4;
  const int nVar_eid  = 13;
  const int nVar_shapeFull = 7;
  const int nVar_brem = 3;
  const int nVar_sc = 3;
  const int nVar_cluster1 = 6;     // Miniaod
  const int nVar_cluster2 = 6;
  const int nVar_cluster3 = 6;
  const int nVar_iso = 3;

  TString variables_trk[nVar_trk] = { "trk_dr", "trk_pt", "trk_eta" };
  TString variables_gsf[nVar_gsf] = { "gsf_dr", "gsf_pt", "gsf_eta", "eid_gsf_nhits", "eid_gsf_chi2red", "gsf_bdtout1" };
  TString variables_ele[nVar_ele] = { "eid_ele_pt", "ele_eta", "ele_mva_value", "core_shFracHits" };
  TString variables_eid[nVar_eid] = { "eid_sc_etaWidth", "eid_sc_phiWidth", "eid_match_seed_dEta", "eid_match_eclu_EoverP", "eid_match_SC_EoverP", "eid_match_SC_dEta", "eid_match_SC_dPhi", "eid_shape_full5x5_sigmaIetaIeta", "eid_shape_full5x5_sigmaIphiIphi", "eid_shape_full5x5_HoverE", "eid_shape_full5x5_r9", "eid_brem_frac", "eid_shape_full5x5_circularity" };
  TString variables_shapeFull[nVar_shapeFull] = { "shape_full5x5_e1x5", "shape_full5x5_e2x5Max", "shape_full5x5_e5x5", "shape_full5x5_eLeft", "shape_full5x5_eRight", "shape_full5x5_eTop", "shape_full5x5_eBottom" };
  TString variables_brem[nVar_brem] = { "brem_fracTrk", "brem_fracSC", "brem_N" };
  TString variables_sc[nVar_sc] = { "eid_sc_eta", "sc_Nclus", "sc_goodSeed" };
  TString variables_cluster1[nVar_cluster1] = { "sc_clus1_nxtal", "sc_clus1_dphi", "sc_clus1_deta", "sc_clus1_E_ov_p", "sc_clus1_E_ov_E", "sc_clus1_ntrk_deta01" };
  TString variables_cluster2[nVar_cluster2] = { "sc_clus2_nxtal", "sc_clus2_dphi", "sc_clus2_deta", "sc_clus2_E_ov_p", "sc_clus2_E_ov_E", "sc_clus2_ntrk_deta01" };
  TString variables_cluster3[nVar_cluster3] = { "sc_clus3_nxtal", "sc_clus3_dphi", "sc_clus3_deta", "sc_clus3_E_ov_p", "sc_clus3_E_ov_E", "sc_clus3_ntrk_deta01" };
  TString variables_iso[nVar_iso] = { "ele_sumPhotonEt", "ele_sumChargedHadronPt", "ele_sumNeutralHadronEt" };

  // Binning
  int bins_trk[nVar_trk]   = { 20, 50, 25 };
  int bins_gsf[nVar_gsf]   = { 20, 50, 25, 30, 50, 100 };
  int bins_ele[nVar_ele]   = { 50, 30, 100, 10 };
  int bins_eid_EB[nVar_eid] = { 50, 50, 50, 50, 50, 50, 50, 100, 100, 50, 50, 100, 100 };
  int bins_eid_EE[nVar_eid] = { 50, 50, 50, 50, 50, 50, 50, 100, 100, 50, 50, 100, 100 };
  int bins_shapeFull_EB[nVar_shapeFull] = { 50, 50, 50, 50, 50, 50, 50 };
  int bins_shapeFull_EE[nVar_shapeFull] = { 50, 50, 50, 50, 50, 50, 50 };
  int bins_brem_EB[nVar_brem] = { 50, 50, 15 };
  int bins_brem_EE[nVar_brem] = { 50, 50, 15 };
  int bins_sc[nVar_sc] = { 30, 50, 2 };
  int bins_cluster1[nVar_cluster1] = { 10, 50, 40, 50, 50, 10 };
  int bins_cluster2[nVar_cluster2] = { 10, 50, 40, 50, 50, 10 };
  int bins_cluster3[nVar_cluster3] = { 10, 50, 40, 50, 50, 10 };
  int bins_iso_EB[nVar_iso]  = { 50, 50, 50 };
  int bins_iso_EE[nVar_iso]  = { 50, 50, 50 };
    
  // Low edge
  float lowedge_trk[nVar_trk]    = { 0., 0., -3. };
  float lowedge_gsf[nVar_gsf]    = { 0., 0., -3., 0, 0, -6. };
  float lowedge_ele[nVar_ele]    = { 0., -3., -10., 0. };
  float lowedge_eid_EB[nVar_eid] = { 0., 0., -1.5, -0.0000001, 0., -1., -2., 0., 0., 0., 0., -0.1, 0. };     
  float lowedge_eid_EE[nVar_eid] = { 0., 0., -1.5, -0.0000001, 0., -1., -2., 0., 0., 0., 0., -0.1, 0. };  		   
  float lowedge_shapeFull_EB[nVar_shapeFull] = { 0., 0., 0., 0., 0., 0., 0. };
  float lowedge_shapeFull_EE[nVar_shapeFull] = { 0., 0., 0., 0., 0., 0., 0. };
  float lowedge_brem_EB[nVar_brem] = { 0., 0., 0. };
  float lowedge_brem_EE[nVar_brem] = { 0., 0., 0. };
  float lowedge_sc[nVar_sc] = { -3., 0., -0.5 };
  float lowedge_cluster1[nVar_cluster1] = { 0, -1.5, -1.0, 0, 0, 0 };
  float lowedge_cluster2[nVar_cluster2] = { 0, -1.5, -1.0, 0, 0, 0 };
  float lowedge_cluster3[nVar_cluster3] = { 0, -1.5, -1.0, 0, 0, 0 };
  float lowedge_iso_EB[nVar_iso] = { 0., 0., 0. };
  float lowedge_iso_EE[nVar_iso] = { 0., 0., 0. };

  // High edge
  float highedge_trk[nVar_trk]    = { 1., 50., 3. };
  float highedge_gsf[nVar_gsf]    = { 1., 50., 3., 30, 150, 10. };
  float highedge_ele[nVar_ele]    = { 100., 3., 10., 1. };
  float highedge_eid_EB[nVar_eid] = { 2., 2., 1.5, 0.0000001, 4., 1., 2., 0.012, 0.012, 3., 1.2, 1., 1. }; 	   
  float highedge_eid_EE[nVar_eid] = { 2., 2., 1.5, 0.0000001, 4., 1., 2., 0.020, 0.020, 3., 1.2, 1., 1. }; 	
  float highedge_shapeFull_EB[nVar_shapeFull] = { 10., 10., 10., 3., 3., 3., 3. };
  float highedge_shapeFull_EE[nVar_shapeFull] = { 10., 10., 10., 3., 3., 3., 3. };
  float highedge_brem_EB[nVar_brem] = { 1., 1., 15. };
  float highedge_brem_EE[nVar_brem] = { 1., 1., 15. };
  float highedge_sc[nVar_sc] = { 3., 20., 1.5 };
  float highedge_cluster1[nVar_cluster1] = { 10, 1.5, 1.0, 1., 1., 10 };
  float highedge_cluster2[nVar_cluster2] = { 10, 1.5, 1.0, 1., 1., 10 };
  float highedge_cluster3[nVar_cluster3] = { 10, 1.5, 1.0, 1., 1., 10 };    
  float highedge_iso_EB[nVar_iso] = { 1., 1., 1. };
  float highedge_iso_EE[nVar_iso] = { 1., 1., 1. };

  // Log scale
  bool logY_trk[nVar_trk]  = { 1, 1, 0 };
  bool logY_gsf[nVar_gsf]  = { 1, 1, 0, 0, 1, 0 };
  bool logY_ele[nVar_ele]  = { 1, 0, 0 };
  bool logY_eid_EB[nVar_eid] = { 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1 };
  bool logY_eid_EE[nVar_eid] = { 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1 };
  bool logY_shapeFull_EB[nVar_shapeFull] = { 0, 0, 0, 1, 1, 1, 1 };
  bool logY_shapeFull_EE[nVar_shapeFull] = { 0, 0, 0, 1, 1, 1, 1 };
  bool logY_brem_EB[nVar_brem] = { 0, 0, 0 };
  bool logY_brem_EE[nVar_brem] = { 0, 0, 0 };
  bool logY_sc[nVar_sc] = { 0, 0, 0 };
  bool logY_cluster1[nVar_cluster1]  = { 0, 0, 0, 0, 0, 1 };
  bool logY_cluster2[nVar_cluster2]  = { 0, 0, 0, 0, 0, 1 };
  bool logY_cluster3[nVar_cluster3]  = { 0, 0, 0, 0, 0, 1 };
  bool logY_iso_EB[nVar_iso] = { 1, 1, 1 };
  bool logY_iso_EE[nVar_iso] = { 1, 1, 1 };

  // Fakes first
  bool fakesfirst_trk[nVar_trk]  = { 1, 0, 0 };
  bool fakesfirst_gsf[nVar_gsf]  = { 0, 1, 0, 0, 0, 1 };
  bool fakesfirst_ele[nVar_ele]  = { 1, 0, 1 };
  bool fakesfirst_eid_EB[nVar_eid] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }; 
  bool fakesfirst_eid_EE[nVar_eid] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }; 
  bool fakesfirst_shapeFull_EB[nVar_shapeFull] = { 1, 1, 1, 0, 0, 0, 0 };
  bool fakesfirst_shapeFull_EE[nVar_shapeFull] = { 0, 0, 0, 0, 0, 0, 0 };
  bool fakesfirst_brem_EB[nVar_brem] = { 0, 0, 0 };
  bool fakesfirst_brem_EE[nVar_brem] = { 0, 0, 0 };
  bool fakesfirst_sc[nVar_sc] = { 0, 1, 0 };
  bool fakesfirst_cluster1[nVar_cluster1]  = { 1, 0, 0, 0, 0, 0 };
  bool fakesfirst_cluster2[nVar_cluster2]  = { 1, 0, 0, 0, 0, 0 };
  bool fakesfirst_cluster3[nVar_cluster3]  = { 1, 0, 0, 0, 0, 0 };
  bool fakesfirst_iso_EB[nVar_iso] = { 0, 0, 0 };
  bool fakesfirst_iso_EE[nVar_iso] = { 0, 0, 0 };

  // Selections
  TString signal_comm    = "is_e && gen_dR<0.1 && trk_eta>-2.4 && trk_eta<2.4";
  TString fakes_comm     = "!is_e && gen_dR>0.1 && trk_eta>-2.4 && trk_eta<2.4";
  TString signal_commEB  = signal_comm+" && fiducial_isEB";
  TString fakes_commEB   = fakes_comm +" && fiducial_isEB";
  TString signal_commEE  = signal_comm+" && fiducial_isEE";
  TString fakes_commEE   = fakes_comm +" && fiducial_isEE";
  TString signal_goodSeed_EB = signal_commEB + " && sc_goodSeed";
  TString fakes_goodSeed_EB  = fakes_commEB  + " && sc_goodSeed";
  TString signal_goodSeed_EE = signal_commEE + " && sc_goodSeed";
  TString fakes_goodSeed_EE  = fakes_commEE  + " && sc_goodSeed";
  TString signal_badSeed_EB  = signal_commEB + " && !sc_goodSeed";
  TString fakes_badSeed_EB   = fakes_commEB  + " && !sc_goodSeed";
  TString signal_badSeed_EE  = signal_commEE + " && !sc_goodSeed";
  TString fakes_badSeed_EE   = fakes_commEE  + " && !sc_goodSeed";
  TString signal_cluster1 = signal_comm + " && sc_clus1_E>0 && sc_clus1_dphi>-50";
  TString fakes_cluster1  = fakes_comm  + " && sc_clus1_E>0 && sc_clus1_dphi>-50";
  TString signal_cluster2 = signal_comm + " && sc_clus2_E>0 && sc_clus2_dphi>-50";
  TString fakes_cluster2  = fakes_comm  + " && sc_clus2_E>0 && sc_clus2_dphi>-50";
  TString signal_cluster3 = signal_comm + " && sc_clus3_E>0 && sc_clus3_dphi>-50";
  TString fakes_cluster3  = fakes_comm  + " && sc_clus3_E>0 && sc_clus3_dphi>-50";

  // Histos, tracks
  TH1F *H_Fakes_trk[nVar_trk];
  TH1F *H_Signal_trk[nVar_trk]; 
  for (int ii=0; ii<nVar_trk; ii++) {
    H_Fakes_trk[ii]  = new TH1F("H_Fakes_"+variables_trk[ii],  variables_trk[ii], bins_trk[ii], lowedge_trk[ii], highedge_trk[ii]); 
    H_Signal_trk[ii] = new TH1F("H_Signal_"+variables_trk[ii], variables_trk[ii], bins_trk[ii], lowedge_trk[ii], highedge_trk[ii]); 
    H_Fakes_trk[ii]->Sumw2();
    H_Signal_trk[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_trk[ii],  variables_trk[ii], fakes_comm);
    theTree->Project("H_Signal_"+variables_trk[ii], variables_trk[ii], signal_comm);
    cout << H_Fakes_trk[ii]->Integral() << " " << H_Signal_trk[ii]->Integral() << endl;
    H_Fakes_trk[ii] ->Scale(1./H_Fakes_trk[ii]->Integral());
    H_Signal_trk[ii]->Scale(1./H_Signal_trk[ii]->Integral());
  }

  // Histos, gsf tracks
  TH1F *H_Fakes_gsf[nVar_gsf];
  TH1F *H_Signal_gsf[nVar_gsf]; 
  for (int ii=0; ii<nVar_gsf; ii++) {
    H_Fakes_gsf[ii]  = new TH1F("H_Fakes_"+variables_gsf[ii],  variables_gsf[ii], bins_gsf[ii], lowedge_gsf[ii], highedge_gsf[ii]); 
    H_Signal_gsf[ii] = new TH1F("H_Signal_"+variables_gsf[ii], variables_gsf[ii], bins_gsf[ii], lowedge_gsf[ii], highedge_gsf[ii]); 
    H_Fakes_gsf[ii]->Sumw2();
    H_Signal_gsf[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_gsf[ii],  variables_gsf[ii], fakes_comm);
    theTree->Project("H_Signal_"+variables_gsf[ii], variables_gsf[ii], signal_comm);
    cout << H_Fakes_gsf[ii]->Integral() << " " << H_Signal_gsf[ii]->Integral() << endl;
    H_Fakes_gsf[ii] ->Scale(1./H_Fakes_gsf[ii]->Integral());
    H_Signal_gsf[ii]->Scale(1./H_Signal_gsf[ii]->Integral());
  }

  // Histos, low pT electrons
  TH1F *H_Fakes_ele[nVar_ele];
  TH1F *H_Signal_ele[nVar_ele]; 
  for (int ii=0; ii<nVar_ele; ii++) {
    H_Fakes_ele[ii]  = new TH1F("H_Fakes_"+variables_ele[ii],  variables_ele[ii], bins_ele[ii], lowedge_ele[ii], highedge_ele[ii]); 
    H_Signal_ele[ii] = new TH1F("H_Signal_"+variables_ele[ii], variables_ele[ii], bins_ele[ii], lowedge_ele[ii], highedge_ele[ii]); 
    H_Fakes_ele[ii]->Sumw2();
    H_Signal_ele[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_ele[ii],  variables_ele[ii], fakes_comm);
    theTree->Project("H_Signal_"+variables_ele[ii], variables_ele[ii], signal_comm);
    cout << H_Fakes_ele[ii]->Integral() << " " << H_Signal_ele[ii]->Integral() << endl;
    H_Fakes_ele[ii] ->Scale(1./H_Fakes_ele[ii]->Integral());
    H_Signal_ele[ii]->Scale(1./H_Signal_ele[ii]->Integral());
  }

  // Histos, id for low pT electrons (barrel)
  TH1F *H_Fakes_eid_EB[nVar_eid];
  TH1F *H_Signal_eid_EB[nVar_eid]; 
  for (int ii=0; ii<nVar_eid; ii++) {
    H_Fakes_eid_EB[ii]  = new TH1F("H_Fakes_"+variables_eid[ii]+"_EB",  variables_eid[ii]+"_EB", bins_eid_EB[ii], lowedge_eid_EB[ii], highedge_eid_EB[ii]); 
    H_Signal_eid_EB[ii] = new TH1F("H_Signal_"+variables_eid[ii]+"_EB", variables_eid[ii]+"_EB", bins_eid_EB[ii], lowedge_eid_EB[ii], highedge_eid_EB[ii]); 
    H_Fakes_eid_EB[ii]->Sumw2();
    H_Signal_eid_EB[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_eid[ii]+"_EB",  variables_eid[ii], fakes_commEB);
    theTree->Project("H_Signal_"+variables_eid[ii]+"_EB", variables_eid[ii], signal_commEB);
    cout << H_Fakes_eid_EB[ii]->Integral() << " " << H_Signal_eid_EB[ii]->Integral() << endl;
    H_Fakes_eid_EB[ii] ->Scale(1./H_Fakes_eid_EB[ii]->Integral());
    H_Signal_eid_EB[ii]->Scale(1./H_Signal_eid_EB[ii]->Integral());
  }

  // Histos, id for low pT electrons (endcap)
  TH1F *H_Fakes_eid_EE[nVar_eid];
  TH1F *H_Signal_eid_EE[nVar_eid]; 
  for (int ii=0; ii<nVar_eid; ii++) {
    H_Fakes_eid_EE[ii]  = new TH1F("H_Fakes_"+variables_eid[ii]+"_EE",  variables_eid[ii]+"_EE", bins_eid_EE[ii], lowedge_eid_EE[ii], highedge_eid_EE[ii]); 
    H_Signal_eid_EE[ii] = new TH1F("H_Signal_"+variables_eid[ii]+"_EE", variables_eid[ii]+"_EE", bins_eid_EE[ii], lowedge_eid_EE[ii], highedge_eid_EE[ii]); 
    H_Fakes_eid_EE[ii]->Sumw2();
    H_Signal_eid_EE[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_eid[ii]+"_EE",  variables_eid[ii], fakes_commEE);
    theTree->Project("H_Signal_"+variables_eid[ii]+"_EE", variables_eid[ii], signal_commEE);
    cout << H_Fakes_eid_EE[ii]->Integral() << " " << H_Signal_eid_EE[ii]->Integral() << endl;
    H_Fakes_eid_EE[ii] ->Scale(1./H_Fakes_eid_EE[ii]->Integral());
    H_Signal_eid_EE[ii]->Scale(1./H_Signal_eid_EE[ii]->Integral());
  }

  // Histos, further id for low pT electrons : shower shapes Full5x5 (barrel)
  TH1F *H_Fakes_shapeFull_EB[nVar_shapeFull];
  TH1F *H_Signal_shapeFull_EB[nVar_shapeFull]; 
  for (int ii=0; ii<nVar_shapeFull; ii++) {
    H_Fakes_shapeFull_EB[ii]  = new TH1F("H_Fakes_"+variables_shapeFull[ii]+"_EB",  variables_shapeFull[ii]+"_EB", bins_shapeFull_EB[ii], lowedge_shapeFull_EB[ii], highedge_shapeFull_EB[ii]); 
    H_Signal_shapeFull_EB[ii] = new TH1F("H_Signal_"+variables_shapeFull[ii]+"_EB", variables_shapeFull[ii]+"_EB", bins_shapeFull_EB[ii], lowedge_shapeFull_EB[ii], highedge_shapeFull_EB[ii]); 
    H_Fakes_shapeFull_EB[ii]->Sumw2();
    H_Signal_shapeFull_EB[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_shapeFull[ii]+"_EB",  variables_shapeFull[ii], fakes_commEB);
    theTree->Project("H_Signal_"+variables_shapeFull[ii]+"_EB", variables_shapeFull[ii], signal_commEB);
    cout << H_Fakes_shapeFull_EB[ii]->Integral() << " " << H_Signal_shapeFull_EB[ii]->Integral() << endl;
    H_Fakes_shapeFull_EB[ii] ->Scale(1./H_Fakes_shapeFull_EB[ii]->Integral());
    H_Signal_shapeFull_EB[ii]->Scale(1./H_Signal_shapeFull_EB[ii]->Integral());
  }

  // Histos, further id for low pT electrons : shower shapes full 5x5 (endcap)
  TH1F *H_Fakes_shapeFull_EE[nVar_shapeFull];
  TH1F *H_Signal_shapeFull_EE[nVar_shapeFull]; 
  for (int ii=0; ii<nVar_shapeFull; ii++) {
    H_Fakes_shapeFull_EE[ii]  = new TH1F("H_Fakes_"+variables_shapeFull[ii]+"_EE",  variables_shapeFull[ii]+"_EE", bins_shapeFull_EE[ii], lowedge_shapeFull_EE[ii], highedge_shapeFull_EE[ii]); 
    H_Signal_shapeFull_EE[ii] = new TH1F("H_Signal_"+variables_shapeFull[ii]+"_EE", variables_shapeFull[ii]+"_EE", bins_shapeFull_EE[ii], lowedge_shapeFull_EE[ii], highedge_shapeFull_EE[ii]); 
    H_Fakes_shapeFull_EE[ii]->Sumw2();
    H_Signal_shapeFull_EE[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_shapeFull[ii]+"_EE",  variables_shapeFull[ii], fakes_commEE);
    theTree->Project("H_Signal_"+variables_shapeFull[ii]+"_EE", variables_shapeFull[ii], signal_commEE);
    cout << H_Fakes_shapeFull_EE[ii]->Integral() << " " << H_Signal_shapeFull_EE[ii]->Integral() << endl;
    H_Fakes_shapeFull_EE[ii] ->Scale(1./H_Fakes_shapeFull_EE[ii]->Integral());
    H_Signal_shapeFull_EE[ii]->Scale(1./H_Signal_shapeFull_EE[ii]->Integral());
  }

  // Histos, further id for low pT electrons : brem-related (barrel)
  TH1F *H_Fakes_brem_EB[nVar_brem];
  TH1F *H_Signal_brem_EB[nVar_brem]; 
  for (int ii=0; ii<nVar_brem; ii++) {
    H_Fakes_brem_EB[ii]  = new TH1F("H_Fakes_"+variables_brem[ii]+"_EB",  variables_brem[ii]+"_EB", bins_brem_EB[ii], lowedge_brem_EB[ii], highedge_brem_EB[ii]); 
    H_Signal_brem_EB[ii] = new TH1F("H_Signal_"+variables_brem[ii]+"_EB", variables_brem[ii]+"_EB", bins_brem_EB[ii], lowedge_brem_EB[ii], highedge_brem_EB[ii]); 
    H_Fakes_brem_EB[ii]->Sumw2();
    H_Signal_brem_EB[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_brem[ii]+"_EB",  variables_brem[ii], fakes_commEB);
    theTree->Project("H_Signal_"+variables_brem[ii]+"_EB", variables_brem[ii], signal_commEB);
    cout << H_Fakes_brem_EB[ii]->Integral() << " " << H_Signal_brem_EB[ii]->Integral() << endl;
    H_Fakes_brem_EB[ii] ->Scale(1./H_Fakes_brem_EB[ii]->Integral());
    H_Signal_brem_EB[ii]->Scale(1./H_Signal_brem_EB[ii]->Integral());
  }

  // Histos, further id for low pT electrons : brem-related (endcap)
  TH1F *H_Fakes_brem_EE[nVar_brem];
  TH1F *H_Signal_brem_EE[nVar_brem]; 
  for (int ii=0; ii<nVar_brem; ii++) {
    H_Fakes_brem_EE[ii]  = new TH1F("H_Fakes_"+variables_brem[ii]+"_EE",  variables_brem[ii]+"_EE", bins_brem_EE[ii], lowedge_brem_EE[ii], highedge_brem_EE[ii]); 
    H_Signal_brem_EE[ii] = new TH1F("H_Signal_"+variables_brem[ii]+"_EE", variables_brem[ii]+"_EE", bins_brem_EE[ii], lowedge_brem_EE[ii], highedge_brem_EE[ii]); 
    H_Fakes_brem_EE[ii]->Sumw2();
    H_Signal_brem_EE[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_brem[ii]+"_EE",  variables_brem[ii], fakes_commEE);
    theTree->Project("H_Signal_"+variables_brem[ii]+"_EE", variables_brem[ii], signal_commEE);
    cout << H_Fakes_brem_EE[ii]->Integral() << " " << H_Signal_brem_EE[ii]->Integral() << endl;
    H_Fakes_brem_EE[ii] ->Scale(1./H_Fakes_brem_EE[ii]->Integral());
    H_Signal_brem_EE[ii]->Scale(1./H_Signal_brem_EE[ii]->Integral());
  }

  // Histos, superclusters
  TH1F *H_Fakes_sc[nVar_sc];
  TH1F *H_Signal_sc[nVar_sc]; 
  for (int ii=0; ii<nVar_sc; ii++) {
    H_Fakes_sc[ii]  = new TH1F("H_Fakes_"+variables_sc[ii],  variables_sc[ii], bins_sc[ii], lowedge_sc[ii], highedge_sc[ii]); 
    H_Signal_sc[ii] = new TH1F("H_Signal_"+variables_sc[ii], variables_sc[ii], bins_sc[ii], lowedge_sc[ii], highedge_sc[ii]); 
    H_Fakes_sc[ii]->Sumw2();
    H_Signal_sc[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_sc[ii],  variables_sc[ii], fakes_comm);
    theTree->Project("H_Signal_"+variables_sc[ii], variables_sc[ii], signal_comm);
    cout << variables_sc[ii] << ", " << H_Fakes_sc[ii]->Integral() << ", " << H_Signal_sc[ii]->Integral() << endl;
    H_Fakes_sc[ii] ->Scale(1./H_Fakes_sc[ii]->Integral());
    H_Signal_sc[ii]->Scale(1./H_Signal_sc[ii]->Integral());
    cout << variables_sc[ii] << ", " << H_Fakes_sc[ii]->Integral() << ", " << H_Signal_sc[ii]->Integral() << endl;
  }

  // Histos, clusters
  TH1F *H_Fakes_cluster1[nVar_cluster1];
  TH1F *H_Signal_cluster1[nVar_cluster1]; 
  for (int ii=0; ii<nVar_cluster1; ii++) {
    H_Fakes_cluster1[ii]  = new TH1F("H_Fakes_"+variables_cluster1[ii],  variables_cluster1[ii], bins_cluster1[ii], lowedge_cluster1[ii], highedge_cluster1[ii]); 
    H_Signal_cluster1[ii] = new TH1F("H_Signal_"+variables_cluster1[ii], variables_cluster1[ii], bins_cluster1[ii], lowedge_cluster1[ii], highedge_cluster1[ii]); 
    H_Fakes_cluster1[ii]->Sumw2();
    H_Signal_cluster1[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_cluster1[ii],  variables_cluster1[ii], fakes_cluster1);
    theTree->Project("H_Signal_"+variables_cluster1[ii], variables_cluster1[ii], signal_cluster1);
    cout << variables_cluster1[ii] << ", " << H_Fakes_cluster1[ii]->Integral() << ", " << H_Signal_cluster1[ii]->Integral() << endl;
    H_Fakes_cluster1[ii] ->Scale(1./H_Fakes_cluster1[ii]->Integral());
    H_Signal_cluster1[ii]->Scale(1./H_Signal_cluster1[ii]->Integral());
    cout << variables_cluster1[ii] << ", " << H_Fakes_cluster1[ii]->Integral() << ", " << H_Signal_cluster1[ii]->Integral() << endl;
  }
  //
  TH1F *H_Fakes_cluster2[nVar_cluster2];
  TH1F *H_Signal_cluster2[nVar_cluster2]; 
  for (int ii=0; ii<nVar_cluster2; ii++) {
    H_Fakes_cluster2[ii]  = new TH1F("H_Fakes_"+variables_cluster2[ii],  variables_cluster2[ii], bins_cluster2[ii], lowedge_cluster2[ii], highedge_cluster2[ii]); 
    H_Signal_cluster2[ii] = new TH1F("H_Signal_"+variables_cluster2[ii], variables_cluster2[ii], bins_cluster2[ii], lowedge_cluster2[ii], highedge_cluster2[ii]); 
    H_Fakes_cluster2[ii]->Sumw2();
    H_Signal_cluster2[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_cluster2[ii],  variables_cluster2[ii], fakes_cluster2);
    theTree->Project("H_Signal_"+variables_cluster2[ii], variables_cluster2[ii], signal_cluster2);
    cout << variables_cluster2[ii] << ", " << H_Fakes_cluster2[ii]->Integral() << ", " << H_Signal_cluster2[ii]->Integral() << endl;
    H_Fakes_cluster2[ii] ->Scale(1./H_Fakes_cluster2[ii]->Integral());
    H_Signal_cluster2[ii]->Scale(1./H_Signal_cluster2[ii]->Integral());
    cout << variables_cluster2[ii] << ", " << H_Fakes_cluster2[ii]->Integral() << ", " << H_Signal_cluster2[ii]->Integral() << endl;
  }
  //
  TH1F *H_Fakes_cluster3[nVar_cluster3];
  TH1F *H_Signal_cluster3[nVar_cluster3]; 
  for (int ii=0; ii<nVar_cluster3; ii++) {
    H_Fakes_cluster3[ii]  = new TH1F("H_Fakes_"+variables_cluster3[ii],  variables_cluster3[ii], bins_cluster3[ii], lowedge_cluster3[ii], highedge_cluster3[ii]); 
    H_Signal_cluster3[ii] = new TH1F("H_Signal_"+variables_cluster3[ii], variables_cluster3[ii], bins_cluster3[ii], lowedge_cluster3[ii], highedge_cluster3[ii]); 
    H_Fakes_cluster3[ii]->Sumw2();
    H_Signal_cluster3[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_cluster3[ii],  variables_cluster3[ii], fakes_cluster3);
    theTree->Project("H_Signal_"+variables_cluster3[ii], variables_cluster3[ii], signal_cluster3);
    cout << variables_cluster3[ii] << ", " << H_Fakes_cluster3[ii]->Integral() << ", " << H_Signal_cluster3[ii]->Integral() << endl;
    H_Fakes_cluster3[ii] ->Scale(1./H_Fakes_cluster3[ii]->Integral());
    H_Signal_cluster3[ii]->Scale(1./H_Signal_cluster3[ii]->Integral());
    cout << variables_cluster3[ii] << ", " << H_Fakes_cluster3[ii]->Integral() << ", " << H_Signal_cluster3[ii]->Integral() << endl;
  }
  //
  // Histos, isolation for low pT electrons (barrel)
  TH1F *H_Fakes_iso_EB[nVar_iso];
  TH1F *H_Signal_iso_EB[nVar_iso]; 
  for (int ii=0; ii<nVar_iso; ii++) {
    H_Fakes_iso_EB[ii]  = new TH1F("H_Fakes_"+variables_iso[ii]+"_EB",  variables_iso[ii]+"_EB", bins_iso_EB[ii], lowedge_iso_EB[ii], highedge_iso_EB[ii]); 
    H_Signal_iso_EB[ii] = new TH1F("H_Signal_"+variables_iso[ii]+"_EB", variables_iso[ii]+"_EB", bins_iso_EB[ii], lowedge_iso_EB[ii], highedge_iso_EB[ii]); 
    H_Fakes_iso_EB[ii]->Sumw2();
    H_Signal_iso_EB[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_iso[ii]+"_EB",  variables_iso[ii], fakes_commEB);
    theTree->Project("H_Signal_"+variables_iso[ii]+"_EB", variables_iso[ii], signal_commEB);
    cout << H_Fakes_iso_EB[ii]->Integral() << " " << H_Signal_iso_EB[ii]->Integral() << endl;
    H_Fakes_iso_EB[ii] ->Scale(1./H_Fakes_iso_EB[ii]->Integral());
    H_Signal_iso_EB[ii]->Scale(1./H_Signal_iso_EB[ii]->Integral());
  }

  // Histos, isolation for low pT electrons (endcap)
  TH1F *H_Fakes_iso_EE[nVar_iso];
  TH1F *H_Signal_iso_EE[nVar_iso]; 
  for (int ii=0; ii<nVar_iso; ii++) {
    H_Fakes_iso_EE[ii]  = new TH1F("H_Fakes_"+variables_iso[ii]+"_EE",  variables_iso[ii]+"_EE", bins_iso_EE[ii], lowedge_iso_EE[ii], highedge_iso_EE[ii]); 
    H_Signal_iso_EE[ii] = new TH1F("H_Signal_"+variables_iso[ii]+"_EE", variables_iso[ii]+"_EE", bins_iso_EE[ii], lowedge_iso_EE[ii], highedge_iso_EE[ii]); 
    H_Fakes_iso_EE[ii]->Sumw2();
    H_Signal_iso_EE[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_iso[ii]+"_EE",  variables_iso[ii], fakes_commEE);
    theTree->Project("H_Signal_"+variables_iso[ii]+"_EE", variables_iso[ii], signal_commEE);
    cout << H_Fakes_iso_EE[ii]->Integral() << " " << H_Signal_iso_EE[ii]->Integral() << endl;
    H_Fakes_iso_EE[ii] ->Scale(1./H_Fakes_iso_EE[ii]->Integral());
    H_Signal_iso_EE[ii]->Scale(1./H_Signal_iso_EE[ii]->Integral());
  }


  // Cosmetics
  for (int ii=0; ii<nVar_trk; ii++) {   
    H_Fakes_trk[ii]  -> SetLineColor(2); 
    H_Signal_trk[ii] -> SetLineColor(4); 
    H_Fakes_trk[ii]  -> SetLineWidth(2);     
    H_Signal_trk[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_gsf; ii++) {   
    H_Fakes_gsf[ii]  -> SetLineColor(2); 
    H_Signal_gsf[ii] -> SetLineColor(4); 
    H_Fakes_gsf[ii]  -> SetLineWidth(2);     
    H_Signal_gsf[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_ele; ii++) {   
    H_Fakes_ele[ii]  -> SetLineColor(2); 
    H_Signal_ele[ii] -> SetLineColor(4); 
    H_Fakes_ele[ii]  -> SetLineWidth(2);     
    H_Signal_ele[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_eid; ii++) {   
    H_Fakes_eid_EB[ii]  -> SetLineColor(2); 
    H_Signal_eid_EB[ii] -> SetLineColor(4); 
    H_Fakes_eid_EB[ii]  -> SetLineWidth(2);     
    H_Signal_eid_EB[ii] -> SetLineWidth(2);   
    H_Fakes_eid_EE[ii]  -> SetLineColor(2); 
    H_Signal_eid_EE[ii] -> SetLineColor(4); 
    H_Fakes_eid_EE[ii]  -> SetLineWidth(2);     
    H_Signal_eid_EE[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_shapeFull; ii++) {   
    H_Fakes_shapeFull_EB[ii]  -> SetLineColor(2); 
    H_Signal_shapeFull_EB[ii] -> SetLineColor(4); 
    H_Fakes_shapeFull_EB[ii]  -> SetLineWidth(2);     
    H_Signal_shapeFull_EB[ii] -> SetLineWidth(2);   
    H_Fakes_shapeFull_EE[ii]  -> SetLineColor(2); 
    H_Signal_shapeFull_EE[ii] -> SetLineColor(4); 
    H_Fakes_shapeFull_EE[ii]  -> SetLineWidth(2);     
    H_Signal_shapeFull_EE[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_brem; ii++) {   
    H_Fakes_brem_EB[ii]  -> SetLineColor(2); 
    H_Signal_brem_EB[ii] -> SetLineColor(4); 
    H_Fakes_brem_EB[ii]  -> SetLineWidth(2);     
    H_Signal_brem_EB[ii] -> SetLineWidth(2);   
    H_Fakes_brem_EE[ii]  -> SetLineColor(2); 
    H_Signal_brem_EE[ii] -> SetLineColor(4); 
    H_Fakes_brem_EE[ii]  -> SetLineWidth(2);     
    H_Signal_brem_EE[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_sc; ii++) {   
    H_Fakes_sc[ii]  -> SetLineColor(2); 
    H_Signal_sc[ii] -> SetLineColor(4); 
    H_Fakes_sc[ii]  -> SetLineWidth(2);     
    H_Signal_sc[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_cluster1; ii++) {   
    H_Fakes_cluster1[ii]  -> SetLineColor(2); 
    H_Signal_cluster1[ii] -> SetLineColor(4); 
    H_Fakes_cluster1[ii]  -> SetLineWidth(2);     
    H_Signal_cluster1[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_cluster2; ii++) {   
    H_Fakes_cluster2[ii]  -> SetLineColor(2); 
    H_Signal_cluster2[ii] -> SetLineColor(4); 
    H_Fakes_cluster2[ii]  -> SetLineWidth(2);     
    H_Signal_cluster2[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_cluster3; ii++) {   
    H_Fakes_cluster3[ii]  -> SetLineColor(2); 
    H_Signal_cluster3[ii] -> SetLineColor(4); 
    H_Fakes_cluster3[ii]  -> SetLineWidth(2);     
    H_Signal_cluster3[ii] -> SetLineWidth(2);   
  }
  for (int ii=0; ii<nVar_iso; ii++) {   
    H_Fakes_iso_EB[ii]  -> SetLineColor(2); 
    H_Signal_iso_EB[ii] -> SetLineColor(4); 
    H_Fakes_iso_EB[ii]  -> SetLineWidth(2);     
    H_Signal_iso_EB[ii] -> SetLineWidth(2);   
    H_Fakes_iso_EE[ii]  -> SetLineColor(2); 
    H_Signal_iso_EE[ii] -> SetLineColor(4); 
    H_Fakes_iso_EE[ii]  -> SetLineWidth(2);     
    H_Signal_iso_EE[ii] -> SetLineWidth(2);   
  }


  // Plots
  TLegend *leg;
  leg = new TLegend(0.55,0.60,0.80,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(H_Signal_trk[0], "Signal", "lp");
  leg->AddEntry(H_Fakes_trk[0],  "Fakes",  "lp");

  for (int ii=0; ii<nVar_trk; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_trk[ii]) H_Fakes_trk[ii] -> Draw();
    else H_Signal_trk[ii]->Draw();
    H_Signal_trk[ii]->Draw("same");
    H_Fakes_trk[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_trk[ii]);
    c.SaveAs(variables_trk[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_gsf; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_gsf[ii]) H_Fakes_gsf[ii] -> Draw();
    else H_Signal_gsf[ii]->Draw();
    H_Signal_gsf[ii]->Draw("same");
    H_Fakes_gsf[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_gsf[ii]);
    c.SaveAs(variables_gsf[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_ele; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_ele[ii]) H_Fakes_ele[ii] -> Draw();
    else H_Signal_ele[ii]->Draw();
    H_Signal_ele[ii]->Draw("same");
    H_Fakes_ele[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_ele[ii]);
    c.SaveAs(variables_ele[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_eid; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_eid_EB[ii]) H_Fakes_eid_EB[ii] -> Draw();
    else H_Signal_eid_EB[ii]->Draw();
    H_Signal_eid_EB[ii]->Draw("same");
    H_Fakes_eid_EB[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_eid_EB[ii]);
    c.SaveAs(variables_eid[ii]+"_EB.png");
  }
  //
  for (int ii=0; ii<nVar_eid; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_eid_EE[ii]) H_Fakes_eid_EE[ii] -> Draw();
    else H_Signal_eid_EE[ii]->Draw();
    H_Signal_eid_EE[ii]->Draw("same");
    H_Fakes_eid_EE[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_eid_EE[ii]);
    c.SaveAs(variables_eid[ii]+"_EE.png");
  }
  //
  for (int ii=0; ii<nVar_shapeFull; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_shapeFull_EB[ii]) H_Fakes_shapeFull_EB[ii] -> Draw();
    else H_Signal_shapeFull_EB[ii]->Draw();
    H_Signal_shapeFull_EB[ii]->Draw("same");
    H_Fakes_shapeFull_EB[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_shapeFull_EB[ii]);
    c.SaveAs(variables_shapeFull[ii]+"_EB.png");
  }
  //
  for (int ii=0; ii<nVar_shapeFull; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_shapeFull_EE[ii]) H_Fakes_shapeFull_EE[ii] -> Draw();
    else H_Signal_shapeFull_EE[ii]->Draw();
    H_Signal_shapeFull_EE[ii]->Draw("same");
    H_Fakes_shapeFull_EE[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_shapeFull_EE[ii]);
    c.SaveAs(variables_shapeFull[ii]+"_EE.png");
  }
  //
  for (int ii=0; ii<nVar_brem; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_brem_EB[ii]) H_Fakes_brem_EB[ii] -> Draw();
    else H_Signal_brem_EB[ii]->Draw();
    H_Signal_brem_EB[ii]->Draw("same");
    H_Fakes_brem_EB[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_brem_EB[ii]);
    c.SaveAs(variables_brem[ii]+"_EB.png");
  }
  //
  for (int ii=0; ii<nVar_brem; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_brem_EE[ii]) H_Fakes_brem_EE[ii] -> Draw();
    else H_Signal_brem_EE[ii]->Draw();
    H_Signal_brem_EE[ii]->Draw("same");
    H_Fakes_brem_EE[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_brem_EE[ii]);
    c.SaveAs(variables_brem[ii]+"_EE.png");
  }
  //
  for (int ii=0; ii<nVar_sc; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_sc[ii]) H_Fakes_sc[ii] -> Draw();
    else H_Signal_sc[ii]->Draw();
    H_Signal_sc[ii]->Draw("same");
    H_Fakes_sc[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_sc[ii]);
    c.SaveAs(variables_sc[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_cluster1; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_cluster1[ii]) H_Fakes_cluster1[ii] -> Draw();
    else H_Signal_cluster1[ii]->Draw();
    H_Signal_cluster1[ii]->Draw("same");
    H_Fakes_cluster1[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_cluster1[ii]);
    c.SaveAs(variables_cluster1[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_cluster2; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_cluster2[ii]) H_Fakes_cluster2[ii] -> Draw();
    else H_Signal_cluster2[ii]->Draw();
    H_Signal_cluster2[ii]->Draw("same");
    H_Fakes_cluster2[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_cluster2[ii]);
    c.SaveAs(variables_cluster2[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_cluster3; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_cluster3[ii]) H_Fakes_cluster3[ii] -> Draw();
    else H_Signal_cluster3[ii]->Draw();
    H_Signal_cluster3[ii]->Draw("same");
    H_Fakes_cluster3[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_cluster3[ii]);
    c.SaveAs(variables_cluster3[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_iso; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_iso_EB[ii]) H_Fakes_iso_EB[ii] -> Draw();
    else H_Signal_iso_EB[ii]->Draw();
    H_Signal_iso_EB[ii]->Draw("same");
    H_Fakes_iso_EB[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_iso_EB[ii]);
    c.SaveAs(variables_iso[ii]+"_EB.png");
  }
  //
  for (int ii=0; ii<nVar_iso; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_iso_EE[ii]) H_Fakes_iso_EE[ii] -> Draw();
    else H_Signal_iso_EE[ii]->Draw();
    H_Signal_iso_EE[ii]->Draw("same");
    H_Fakes_iso_EE[ii] ->Draw("same");
    leg->Draw();
    c.SetLogy(logY_iso_EE[ii]);
    c.SaveAs(variables_iso[ii]+"_EE.png");
  }

}


