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

// Utilities
#include "../src/FiguresOfMeritEvaluator.cc"

using namespace std;

void idRocs() {
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile *theFile = TFile::Open("../run/output_fromJune28.root");
  TTree *theTree = (TTree*)theFile->Get("ntuplizer/tree");  

  // BDT output
  TH1F *H_FakesPFEle  = new TH1F("H_FakesPFEle", "",100,-1.,1.);                      
  TH1F *H_SignalPFEle = new TH1F("H_SignalPFEle","",100,-1.,1.);                      
  H_FakesPFEle->Sumw2();
  H_SignalPFEle->Sumw2();    
  theTree->Project("H_FakesPFEle",  "ele_mva_value","ele_pt>0.5 && fabs(ele_eta)<2.4 && has_ele && is_egamma && !is_e");
  theTree->Project("H_SignalPFEle", "ele_mva_value","ele_pt>0.5 && fabs(ele_eta)<2.4 && has_ele && is_egamma && is_e");
  H_FakesPFEle ->Scale(1./H_FakesPFEle->Integral());
  H_SignalPFEle->Scale(1./H_SignalPFEle->Integral()); 
  //
  TH1F *H_FakesLowPtEle  = new TH1F("H_FakesLowPtEle", "",100,-10.,10.);                      
  TH1F *H_SignalLowPtEle = new TH1F("H_SignalLowPtEle","",100,-10.,10.);                      
  H_FakesLowPtEle->Sumw2();
  H_SignalLowPtEle->Sumw2();    
  theTree->Project("H_FakesLowPtEle",  "ele_mva_value","ele_pt>0.5 && fabs(ele_eta)<2.4 && has_ele && !is_egamma && !is_e");
  theTree->Project("H_SignalLowPtEle", "ele_mva_value","ele_pt>0.5 && fabs(ele_eta)<2.4 && has_ele && !is_egamma && is_e");
  H_FakesLowPtEle ->Scale(1./H_FakesLowPtEle->Integral()); 
  H_SignalLowPtEle->Scale(1./H_SignalLowPtEle->Integral()); 
  
  // Cosmetics
  H_FakesPFEle  -> SetLineColor(2);   
  H_FakesPFEle  -> SetLineWidth(2);
  H_SignalPFEle -> SetLineColor(4);   
  H_SignalPFEle -> SetLineWidth(2);
  //
  H_FakesLowPtEle  -> SetLineColor(2);   
  H_FakesLowPtEle  -> SetLineWidth(2);
  H_SignalLowPtEle -> SetLineColor(4);   
  H_SignalLowPtEle -> SetLineWidth(2);


  // Plots
  TLegend *leg;
  leg = new TLegend(0.15,0.55,0.50,0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(H_SignalPFEle, "Signal", "lp");
  leg->AddEntry(H_FakesPFEle, "Fakes", "lp");

  TCanvas c1("c1","",1);
  H_FakesPFEle  -> GetXaxis()->SetTitle("BDT output");
  H_SignalPFEle -> GetXaxis()->SetTitle("BDT output");
  H_FakesPFEle->Draw();
  H_SignalPFEle->Draw("same");
  leg->Draw();
  c1.SetLogy();
  c1.SaveAs("bdtPFEle.png");

  TCanvas c2("c2","",1);
  H_FakesLowPtEle  -> GetXaxis()->SetTitle("BDT output");
  H_SignalLowPtEle -> GetXaxis()->SetTitle("BDT output");
  H_FakesLowPtEle->Draw();
  H_SignalLowPtEle->Draw("same");
  leg->Draw();
  c2.SaveAs("bdtLowPtEle.png");


  // Producing the FOMs
  FiguresOfMeritEvaluator myEval_lowPtEle;
  myEval_lowPtEle.addSignal("lowPtEle", H_SignalLowPtEle);
  myEval_lowPtEle.addBackgrounds(H_FakesLowPtEle);
  myEval_lowPtEle.setCutDirection(">");
  TGraph *myGraph_lowPtEle= myEval_lowPtEle.getFOM("lowPtEle",2);
  myGraph_lowPtEle->SetTitle("LowPtElectrons");
  myGraph_lowPtEle->SetName("lowPtEle_graph");
  myGraph_lowPtEle->GetXaxis()->SetTitle("Mistag Rate");
  myGraph_lowPtEle->GetYaxis()->SetTitle("Efficiency");

  TCanvas c3("c3","",1);
  c3.SetGrid();
  myGraph_lowPtEle->SetMarkerColor(4);
  myGraph_lowPtEle->SetMarkerStyle(20);
  myGraph_lowPtEle->SetMarkerSize(1);
  myGraph_lowPtEle->Draw("AP");
  c3.SetLogx();
  c3.SaveAs("rocLowPtEle.png");
}


