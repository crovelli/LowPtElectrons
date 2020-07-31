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

void checkRegression() {
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // Files
  TFile *theFileNo   = TFile::Open("/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL__normalized.root");
  TFile *theFileWith = TFile::Open("/eos/cms/store/user/crovelli/LowPtEle/Batch1_Aug22/miniaod/BuToKJpsiToeeALL_withRegression__normalized.root");
  TTree *theTreeNo   = (TTree*)theFileNo->Get("ntuplizer/tree");  
  TTree *theTreeWith = (TTree*)theFileWith->Get("ntuplizer/tree");  

  // Variables
  const int nVar = 3;
  TString variables[nVar] = { "eid_match_SC_EoverP", "eid_match_eclu_EoverP", "eid_sc_E" };

  // Binning
  int bins[nVar]       = { 50, 50, 25 };
  float lowedge[nVar]  = { 0., -0.5, 0. };
  float highedge[nVar] = { 2., 0.5, 15. };

  // Log scale
  bool logY[nVar] = { 0, 0, 0 };

  // Fakes first
  bool fakesfirst[nVar] = { 0, 0, 0 };

  // Selections
  TString signal_commEB = "is_e && fiducial_isEB";
  TString fakes_commEB  = "!is_e && fiducial_isEB";
  TString signal_commEE = "is_e && fiducial_isEE";
  TString fakes_commEE  = "!is_e && fiducial_isEE";


  // Histos, no regression, barrel
  TH1F *H_Fakes_NoReg_EB[nVar];
  TH1F *H_Signal_NoReg_EB[nVar]; 
  for (int ii=0; ii<nVar; ii++) {
    H_Fakes_NoReg_EB[ii]  = new TH1F("H_FakesNoReg_"+variables[ii]+"_EB",  "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Signal_NoReg_EB[ii] = new TH1F("H_SignalNoReg_"+variables[ii]+"_EB", "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Fakes_NoReg_EB[ii]->Sumw2();
    H_Signal_NoReg_EB[ii]->Sumw2();
    theTreeNo->Project("H_FakesNoReg_"+variables[ii]+"_EB",  variables[ii], fakes_commEB);
    theTreeNo->Project("H_SignalNoReg_"+variables[ii]+"_EB", variables[ii], signal_commEB);
    cout << "No regression, EB, var " << ii << " => " << H_Fakes_NoReg_EB[ii]->Integral() << " " << H_Signal_NoReg_EB[ii]->Integral() << endl;
    H_Fakes_NoReg_EB[ii] ->Scale(1./H_Fakes_NoReg_EB[ii]->Integral());
    H_Signal_NoReg_EB[ii]->Scale(1./H_Signal_NoReg_EB[ii]->Integral());
  }

  // Histos, no regression, endcap
  TH1F *H_Fakes_NoReg_EE[nVar];
  TH1F *H_Signal_NoReg_EE[nVar]; 
  for (int ii=0; ii<nVar; ii++) {
    H_Fakes_NoReg_EE[ii]  = new TH1F("H_FakesNoReg_"+variables[ii]+"_EE",  "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Signal_NoReg_EE[ii] = new TH1F("H_SignalNoReg_"+variables[ii]+"_EE", "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Fakes_NoReg_EE[ii]->Sumw2();
    H_Signal_NoReg_EE[ii]->Sumw2();
    theTreeNo->Project("H_FakesNoReg_"+variables[ii]+"_EE",  variables[ii], fakes_commEE);
    theTreeNo->Project("H_SignalNoReg_"+variables[ii]+"_EE", variables[ii], signal_commEE);
    cout << "No regression, EE, var " << ii << " => " << H_Fakes_NoReg_EE[ii]->Integral() << " " << H_Signal_NoReg_EE[ii]->Integral() << endl;
    H_Fakes_NoReg_EE[ii] ->Scale(1./H_Fakes_NoReg_EE[ii]->Integral());
    H_Signal_NoReg_EE[ii]->Scale(1./H_Signal_NoReg_EE[ii]->Integral());
  }

  // Histos, with regression, barrel
  TH1F *H_Fakes_WithReg_EB[nVar];
  TH1F *H_Signal_WithReg_EB[nVar]; 
  for (int ii=0; ii<nVar; ii++) {
    H_Fakes_WithReg_EB[ii]  = new TH1F("H_FakesWithReg_"+variables[ii]+"_EB",  "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Signal_WithReg_EB[ii] = new TH1F("H_SignalWithReg_"+variables[ii]+"_EB", "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Fakes_WithReg_EB[ii]->Sumw2();
    H_Signal_WithReg_EB[ii]->Sumw2();
    theTreeWith->Project("H_FakesWithReg_"+variables[ii]+"_EB",  variables[ii], fakes_commEB);
    theTreeWith->Project("H_SignalWithReg_"+variables[ii]+"_EB", variables[ii], signal_commEB);
    cout << "With regression, EB, var " << ii << " => " << H_Fakes_WithReg_EB[ii]->Integral() << " " << H_Signal_WithReg_EB[ii]->Integral() << endl;
    H_Fakes_WithReg_EB[ii] ->Scale(1./H_Fakes_WithReg_EB[ii]->Integral());
    H_Signal_WithReg_EB[ii]->Scale(1./H_Signal_WithReg_EB[ii]->Integral());
  }

  // Histos, with regression, endcap
  TH1F *H_Fakes_WithReg_EE[nVar];
  TH1F *H_Signal_WithReg_EE[nVar]; 
  for (int ii=0; ii<nVar; ii++) {
    H_Fakes_WithReg_EE[ii]  = new TH1F("H_FakesWithReg_"+variables[ii]+"_EE",  "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Signal_WithReg_EE[ii] = new TH1F("H_SignalWithReg_"+variables[ii]+"_EE", "", bins[ii], lowedge[ii], highedge[ii]); 
    H_Fakes_WithReg_EE[ii]->Sumw2();
    H_Signal_WithReg_EE[ii]->Sumw2();
    theTreeWith->Project("H_FakesWithReg_"+variables[ii]+"_EE",  variables[ii], fakes_commEE);
    theTreeWith->Project("H_SignalWithReg_"+variables[ii]+"_EE", variables[ii], signal_commEE);
    cout << "With regression, EE, var " << ii << " => " << H_Fakes_WithReg_EE[ii]->Integral() << " " << H_Signal_WithReg_EE[ii]->Integral() << endl;
    H_Fakes_WithReg_EE[ii] ->Scale(1./H_Fakes_WithReg_EE[ii]->Integral());
    H_Signal_WithReg_EE[ii]->Scale(1./H_Signal_WithReg_EE[ii]->Integral());
  }


  // Cosmetics
  for (int ii=0; ii<nVar; ii++) {   
    H_Fakes_NoReg_EB[ii]  -> SetLineColor(2); 
    H_Signal_NoReg_EB[ii] -> SetLineColor(2); 
    H_Fakes_NoReg_EB[ii]  -> SetLineWidth(2);     
    H_Signal_NoReg_EB[ii] -> SetLineWidth(2);   

    H_Fakes_NoReg_EE[ii]  -> SetLineColor(2); 
    H_Signal_NoReg_EE[ii] -> SetLineColor(2); 
    H_Fakes_NoReg_EE[ii]  -> SetLineWidth(2);     
    H_Signal_NoReg_EE[ii] -> SetLineWidth(2);   

    H_Fakes_WithReg_EB[ii]  -> SetLineColor(4); 
    H_Signal_WithReg_EB[ii] -> SetLineColor(4); 
    H_Fakes_WithReg_EB[ii]  -> SetLineWidth(2);     
    H_Signal_WithReg_EB[ii] -> SetLineWidth(2);   

    H_Fakes_WithReg_EE[ii]  -> SetLineColor(4); 
    H_Signal_WithReg_EE[ii] -> SetLineColor(4); 
    H_Fakes_WithReg_EE[ii]  -> SetLineWidth(2);     
    H_Signal_WithReg_EE[ii] -> SetLineWidth(2);   
  }


  // Plots
  TLegend *leg2;
  leg2 = new TLegend(0.55,0.60,0.80,0.85);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(H_Signal_NoReg_EB[0],  "No regression",   "lp");
  leg2->AddEntry(H_Fakes_WithReg_EB[0], "With regression", "lp");

  for (int ii=0; ii<nVar; ii++) {  

    TCanvas c5("c5","",1);
    H_Fakes_WithReg_EB[ii] -> Draw();
    H_Fakes_NoReg_EB[ii]   -> Draw("same");
    leg2->Draw();
    c5.SetLogy(logY[ii]);
    c5.SaveAs(variables[ii]+"_Fakes_EB.png");

    TCanvas c6("c6","",1);
    H_Fakes_WithReg_EE[ii] -> Draw();
    H_Fakes_NoReg_EE[ii]   -> Draw("same");
    leg2->Draw();
    c6.SetLogy(logY[ii]);
    c6.SaveAs(variables[ii]+"_Fakes_EE.png");

    TCanvas c7("c7","",1);
    H_Signal_WithReg_EB[ii] -> Draw();
    H_Signal_NoReg_EB[ii]   -> Draw("same");
    leg2->Draw();
    c7.SetLogy(logY[ii]);
    c7.SaveAs(variables[ii]+"_Signal_EB.png");

    TCanvas c8("c8","",1);
    H_Signal_WithReg_EE[ii] -> Draw();
    H_Signal_NoReg_EE[ii]   -> Draw("same");
    leg2->Draw();
    c8.SetLogy(logY[ii]);
    c8.SaveAs(variables[ii]+"_Signal_EE.png");
  }


}


