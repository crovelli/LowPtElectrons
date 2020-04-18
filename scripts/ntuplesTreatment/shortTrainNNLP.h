#ifndef shortTrainNNLP_h
#define shortTrainNNLP_h
#include "TROOT.h"
#include "TRint.h"


#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "TH2.h"
#include "TF1.h"
#include "TH1.h"

class shortTrainNNLP{

 public: 
  shortTrainNNLP();
  void train(Int_t ntrains);
  virtual ~shortTrainNNLP();
};

#endif 
#ifdef shortTrainNNLP_cxx

shortTrainNNLP::shortTrainNNLP(){

}


shortTrainNNLP::~shortTrainNNLP()
{

}
#endif 

