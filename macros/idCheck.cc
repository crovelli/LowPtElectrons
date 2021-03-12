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

void idCheck() {
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile *theFile = TFile::Open("/eos/cms/store/user/crovelli/LowPtEle/Batch1_UL/BuToKJpsiToee_UL18_AOD__prescale0d30.root");
  //TFile *theFile = TFile::Open("/eos/cms/store/user/crovelli/LowPtEle/Batch1_UL_newBkg/BuToKJpsiToee_UL18_AOD__prescale0d30.root");
  TTree *theTree = (TTree*)theFile->Get("ntuplizer/tree");  
  
  // Variables
  const int nVar_trk  = 4;
  const int nVar_gsf  = 5;
  const int nVar_ele  = 3;
  const int nVar_eid  = 12;
  const int nVar_cluster1 = 5;  
  const int nVar_cluster2 = 4;

  TString variables_trk[nVar_trk] = { "eid_trk_dr", "eid_trk_p", "eid_trk_nhits", "eid_trk_chi2red" };
  TString variables_gsf[nVar_gsf] = { "eid_gsf_dr", "eid_gsf_mode_p", "eid_gsf_nhits", "eid_gsf_chi2red", "eid_gsf_bdtout1" };
  TString variables_ele[nVar_ele] = { "eid_rho", "eid_sc_eta", "eid_core_shFracHits" };
  TString variables_eid[nVar_eid] = { "eid_shape_full5x5_r9","eid_sc_etaWidth", "eid_sc_phiWidth", "eid_shape_full5x5_HoverE","eid_brem_frac","eid_match_SC_EoverP","eid_match_eclu_EoverP","eid_match_SC_dEta","eid_match_SC_dPhi","eid_match_seed_dEta","eid_sc_E","eid_sc_Nclus"};
  TString variables_cluster1[nVar_cluster1] = { "eid_sc_clus1_nxtal", "eid_sc_clus1_dphi", "eid_sc_clus1_deta",   "eid_sc_clus1_E_ov_p", "eid_sc_clus1_E" };
  TString variables_cluster2[nVar_cluster2] = { "eid_sc_clus2_dphi",  "eid_sc_clus2_deta", "eid_sc_clus2_E_ov_p", "eid_sc_clus2_E" };

  // Binning and range
  int bins_trk[nVar_trk]       = { 50,   60,  30,  50 };
  float lowedge_trk[nVar_trk]  = { 0.,   0.,  0., 0. };
  float highedge_trk[nVar_trk] = { 0.01, 30., 30., 5. };

  int bins_gsf[nVar_gsf]       = { 50,   60,  30,  50,  90  };
  float lowedge_gsf[nVar_gsf]  = { 0.,   0.,   0., 0., -10. };
  float highedge_gsf[nVar_gsf] = { 0.003, 30., 30., 5.,  15. };

  int bins_ele[nVar_ele]       = {100,   48, 10  };
  float lowedge_ele[nVar_ele]  = {  0.,-2.4,  0. };
  float highedge_ele[nVar_ele] = { 50., 2.4,  1. }; 

  int bins_eid[nVar_eid]       = { 50, 50, 50, 50, 50, 50,  50,   50,   50,   50,  100,  30 };
  float lowedge_eid[nVar_eid]  = { 0., 0., 0., 0., 0.,  0., 0., -0.5, -0.5, -0.15,   0., 0. };
  float highedge_eid[nVar_eid] = { 1., 1., 1., 0.5, 1., 3., 0.5,  0.5,  0.5, 0.15,  50., 30.};

  int bins_cluster1[nVar_cluster1]       = { 10,  50,   50,  50,  80 };  
  float lowedge_cluster1[nVar_cluster1]  = {  0, -0.5, -0.5,  0.,  0.};
  float highedge_cluster1[nVar_cluster1] = { 10,  0.5,  0.5,  1., 20.};

  int bins_cluster2[nVar_cluster2]       = {   50,  50,   50, 40 };  
  float lowedge_cluster2[nVar_cluster2]  = { -0.5, -0.5,  0.,  0.};
  float highedge_cluster2[nVar_cluster2] = {  0.5,  0.5,  1., 10.};

  // Log scale
  bool logY_trk[nVar_trk] = { 0, 0, 0, 0 };
  bool logY_gsf[nVar_gsf] = { 0, 0, 0, 0, 0 };
  bool logY_ele[nVar_ele] = { 0, 0, 0 };
  bool logY_eid[nVar_eid] = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
  bool logY_cluster1[nVar_cluster1] = { 0, 0, 0, 0, 0 };
  bool logY_cluster2[nVar_cluster2] = { 0, 0, 0, 0 };

  // Fakes first
  bool fakesfirst_trk[nVar_trk] = { 0, 1, 1, 1 };
  bool fakesfirst_gsf[nVar_gsf] = { 0, 1, 0, 0, 1 };
  bool fakesfirst_ele[nVar_ele] = { 0, 0, 0 };
  bool fakesfirst_eid[nVar_eid] = { 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0 }; 
  bool fakesfirst_cluster1[nVar_cluster1] = { 1, 0, 0, 1, 1 };
  bool fakesfirst_cluster2[nVar_cluster2] = { 0, 0, 0, 1 };

  // Selections
  TString signal_comm    = "is_e";
  TString fakes_comm     = "!is_e";
  TString signal_cluster1 = signal_comm + " && eid_sc_clus1_E>0 && eid_sc_clus1_dphi>-50";
  TString fakes_cluster1  = fakes_comm  + " && eid_sc_clus1_E>0 && eid_sc_clus1_dphi>-50";
  TString signal_cluster2 = signal_comm + " && eid_sc_clus2_E>0 && eid_sc_clus2_dphi>-50";
  TString fakes_cluster2  = fakes_comm  + " && eid_sc_clus2_E>0 && eid_sc_clus2_dphi>-50";

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

  // Histos, id for low pT electrons
  TH1F *H_Fakes_eid[nVar_eid];
  TH1F *H_Signal_eid[nVar_eid]; 
  for (int ii=0; ii<nVar_eid; ii++) {
    H_Fakes_eid[ii]  = new TH1F("H_Fakes_"+variables_eid[ii],  variables_eid[ii]+"", bins_eid[ii], lowedge_eid[ii], highedge_eid[ii]); 
    H_Signal_eid[ii] = new TH1F("H_Signal_"+variables_eid[ii], variables_eid[ii]+"", bins_eid[ii], lowedge_eid[ii], highedge_eid[ii]); 
    H_Fakes_eid[ii]->Sumw2();
    H_Signal_eid[ii]->Sumw2();
    theTree->Project("H_Fakes_"+variables_eid[ii],  variables_eid[ii], fakes_comm);
    theTree->Project("H_Signal_"+variables_eid[ii], variables_eid[ii], signal_comm);
    cout << H_Fakes_eid[ii]->Integral() << " " << H_Signal_eid[ii]->Integral() << endl;
    H_Fakes_eid[ii] ->Scale(1./H_Fakes_eid[ii]->Integral());
    H_Signal_eid[ii]->Scale(1./H_Signal_eid[ii]->Integral());
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
    H_Fakes_eid[ii]  -> SetLineColor(2); 
    H_Signal_eid[ii] -> SetLineColor(4); 
    H_Fakes_eid[ii]  -> SetLineWidth(2);     
    H_Signal_eid[ii] -> SetLineWidth(2);   
    H_Fakes_eid[ii]  -> SetLineColor(2); 
    H_Signal_eid[ii] -> SetLineColor(4); 
    H_Fakes_eid[ii]  -> SetLineWidth(2);     
    H_Signal_eid[ii] -> SetLineWidth(2);   
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
    if (fakesfirst_trk[ii]) H_Fakes_trk[ii] -> Draw("hist");
    else H_Signal_trk[ii]->Draw("hist");
    H_Signal_trk[ii]->Draw("samehist");
    H_Fakes_trk[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_trk[ii]);
    c.SaveAs(variables_trk[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_gsf; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_gsf[ii]) H_Fakes_gsf[ii] -> Draw("hist");
    else H_Signal_gsf[ii]->Draw("hist");
    H_Signal_gsf[ii]->Draw("samehist");
    H_Fakes_gsf[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_gsf[ii]);
    c.SaveAs(variables_gsf[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_ele; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_ele[ii]) H_Fakes_ele[ii] -> Draw("hist");
    else H_Signal_ele[ii]->Draw("hist");
    H_Signal_ele[ii]->Draw("samehist");
    H_Fakes_ele[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_ele[ii]);
    c.SaveAs(variables_ele[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_eid; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_eid[ii]) H_Fakes_eid[ii] -> Draw("hist");
    else H_Signal_eid[ii]->Draw("hist");
    H_Signal_eid[ii]->Draw("samehist");
    H_Fakes_eid[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_eid[ii]);
    c.SaveAs(variables_eid[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_cluster1; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_cluster1[ii]) H_Fakes_cluster1[ii] -> Draw("hist");
    else H_Signal_cluster1[ii]->Draw("hist");
    H_Signal_cluster1[ii]->Draw("samehist");
    H_Fakes_cluster1[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_cluster1[ii]);
    c.SaveAs(variables_cluster1[ii]+".png");
  }
  //
  for (int ii=0; ii<nVar_cluster2; ii++) {  
    TCanvas c("c","",1);
    if (fakesfirst_cluster2[ii]) H_Fakes_cluster2[ii] -> Draw("hist");
    else H_Signal_cluster2[ii]->Draw("hist");
    H_Signal_cluster2[ii]->Draw("samehist");
    H_Fakes_cluster2[ii] ->Draw("samehist");
    leg->Draw();
    c.SetLogy(logY_cluster2[ii]);
    c.SaveAs(variables_cluster2[ii]+".png");
  }
}


