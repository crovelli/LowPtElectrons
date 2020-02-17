#include "../interface/FiguresOfMeritEvaluator.h"

#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

FiguresOfMeritEvaluator::FiguresOfMeritEvaluator() {

  m_xmin = 0.0;
  m_xmax = 1.0;
  m_ymin = 0.0;
  m_ymax = 1.0;

  m_signalTitle     = "signal";
  m_backgroundTitle = "background";
}

void FiguresOfMeritEvaluator::addSignal(const char *nameVar, TH1F* sig) {
  
  m_signalHisto.push_back(sig);
  m_names.push_back(TString(nameVar));
}

void FiguresOfMeritEvaluator::addBackgrounds(TH1F* bkg0, TH1F* bkg1) {
  
  if(bkg1) bkg0->Add(bkg1);
  m_bkgHisto.push_back(bkg0);
}

TGraph* FiguresOfMeritEvaluator::getFOM(const char *nameVar, int option) {
  
  TGraph *outGraph = new TGraph();

  int indexVar = -1;
  
  for(unsigned int ivar=0; ivar<m_signalHisto.size(); ivar++) {
    
    if (m_names[ivar].Contains(nameVar) && 
	TString(nameVar).Contains(m_names[ivar])) indexVar=ivar;
  }
  
  if ( indexVar==-1 ) {
    
    std::cout << "ERROR! The requested variable ( "
	      << nameVar << " ) is not in the list of known variables!" << std::endl;
    return 0;
  }

  TH1F *signal       = m_signalHisto[indexVar];
  TH1F *background   = m_bkgHisto[indexVar];
  const char *cutDir = m_direction[indexVar];
  
  if( signal && background ) {
    TAxis *axisS = signal->GetXaxis();
    TAxis *axisB = background->GetXaxis();
    int nBinsSig = axisS->GetNbins();
    int nBinsBkg = axisB->GetNbins();
    
    if( nBinsSig!=nBinsBkg ) {
      std::cout << "ERROR! signal and background histograms have different binning." 
		<< std::endl;
      return 0;
    }

    // needed also overflow + underflow
    outGraph->Set(nBinsSig+2);

    double signalIntegral     = signal->Integral(0,nBinsSig+1);
    double backgroundIntegral = background->Integral(0,nBinsSig+1);
    
    double tmpSignalIntegral     =0.0;
    double tmpBackgroundIntegral =0.0;

    // cout << "signalIntegral = "     << signalIntegral     << endl;
    // cout << "backgroundIntegral = " << backgroundIntegral << endl;

    for ( int ibin=0; ibin<=nBinsSig+1; ibin++) {

      if( strcmp(cutDir,"<")==0 ) {
	tmpSignalIntegral     = signal->Integral(0,ibin);
	tmpBackgroundIntegral = background->Integral(0,ibin);
      }
      else if( strcmp(cutDir,">")==0 ) {
	tmpSignalIntegral     = signal->Integral(ibin,nBinsSig+1);
	tmpBackgroundIntegral = background->Integral(ibin,nBinsSig+1);
      }
      else {
	std::cout << "CONFIGURATION ERROR! direction of the cut not set." << std::endl
		  << "Please use: \">\" for var>x0 or  \"<\" for var<x0" << std::endl;
	return 0;
      }

      // cout << "ibin = " << ibin 
      //	   << ": signal = " << tmpSignalIntegral << ", background = " << tmpBackgroundIntegral << endl;

      double signalEff     = tmpSignalIntegral / signalIntegral;
      double backgroundEff = tmpBackgroundIntegral / backgroundIntegral;

      if( option == 0 ) {
	outGraph->SetPoint(ibin,signalEff,1-backgroundEff);
      }
      else if( option == 1 ) {
	outGraph->SetPoint(ibin,signalEff,backgroundEff);
      }
      else if( option == 2 ) {
	outGraph->SetPoint(ibin,backgroundEff,signalEff);
      }
      else {
	std::cout << "unrecognized option" << std::endl;
	return 0;
      }
    }
  }

  else {
    std::cout << "ERROR! Cannot find signal or background histogram for variable "
	      << nameVar << std::endl;
    return 0;
  }
  
  return outGraph;
}

void FiguresOfMeritEvaluator:: drawResults(const char *fileName, int option) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  if( m_signalHisto.size()!=m_bkgHisto.size() ) {
    std::cout << "ERROR! for some variable signal or background missing. Exiting!" << std::endl;
    return;
  }
  
  TCanvas c1("c1","--c1--",356,0,656,800);
  TLegend* leg = new TLegend(0.50,0.75,0.70,0.90);
  leg->SetFillStyle(0); 
  leg->SetBorderSize(0); 
  leg->SetTextSize(0.03); 
  leg->SetFillColor(0);

  for( unsigned int ivar=0; ivar<m_signalHisto.size(); ivar++) {
    
    const char *name = m_names[ivar].Data();
    
    TGraph *graph = getFOM(name);

    if( graph ) {
      
      std::cout << "---> processing " << name << "..."; 
      
      leg->AddEntry(graph,name,"p");
      
      graph->SetTitle("");
      graph->SetMarkerStyle(8);
      graph->SetMarkerColor(ivar+1);
      graph->SetLineColor(ivar+1);
      graph->GetXaxis()->SetRangeUser(m_xmin,m_xmax);
      graph->GetYaxis()->SetRangeUser(m_ymin,m_ymax);
      
      std::string sigSuffix = " efficiency";
      std::string xAxisName = std::string(m_signalTitle) + sigSuffix;
      
      std::string bkgSuffix = (option==0) ? " rejection" : " efficiency";
      std::string yAxisName = std::string(m_backgroundTitle) + bkgSuffix;
      
      graph->GetXaxis()->SetTitle(xAxisName.c_str());
      graph->GetYaxis()->SetTitle(yAxisName.c_str());
      
      if(ivar==0) graph->Draw("AP");
      else  graph->Draw("P");
      
      std::cout << " done." << std::endl;
    }
    else {
      std::cout << "WARNING! The requested variable " << name 
		<< " could not have been drawn" << std::endl;
    }
  }
  
  leg->Draw();
  c1.SaveAs(fileName);
}

void FiguresOfMeritEvaluator::setRange(double xmin, double xmax, double ymin, double ymax) {
  
  m_xmin = xmin;
  m_xmax = xmax;
  m_ymin = ymin;
  m_ymax = ymax;
}
