#include "LowPtElectrons/LowPtElectrons/interface/ElectronNtuple.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/Point3D.h"

void ElectronNtuple::link_tree(TTree *tree) {
	

	tree->Branch("run",  &run_ , "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt",  &evt_ , "evt/i");

	tree_->Branch("gen_pt" , &gen_pt_ , "gen_pt/f" );
	tree_->Branch("gen_eta", &gen_eta_, "gen_eta/f");
	tree_->Branch("gen_phi", &gen_phi_, "gen_phi/f");

	tree_->Branch("trk_pt",				 	&trk_pt_				   , "trk_pt/f");
	tree_->Branch("trk_eta",		 	  &trk_eta_		       , "trk_eta/f");
	tree_->Branch("trk_phi",		 	  &trk_phi_     		 , "trk_phi/f");
	tree_->Branch("trk_nhits",			&trk_nhits_        , "trk_nhits/f");
	tree_->Branch("trk_high_purity",&trk_high_purity_	 , "trk_high_purity/i");
	tree_->Branch("trk_dxy",			  &trk_dxy_		  		 , "trk_dxy/f");
	tree_->Branch("trk_dxy_err",		&trk_dxy_err_			 , "trk_dxy_err/f");
	tree_->Branch("trk_inp",  			&trk_inp_  	    	 , "trk_inp/f");
	tree_->Branch("trk_outp",	  		&trk_outp_	  	   , "trk_outp/f");
	tree_->Branch("trk_chi2red",    &trk_chi2red_      , "trk_chi2red/f"); 

	tree_->Branch("preid_ibin", 					&preid_ibin_          , "preid_ibin/I");
	tree_->Branch("preid_trk_ecal_match", &preid_trk_ecal_match_, "preid_trk_ecal_match/O");
	tree_->Branch("preid_bdtout",				 	&preid_bdtout_		  	, "preid_bdtout/f");
	tree_->Branch("preid_trk_ecal_Deta",	&preid_trk_ecal_Deta_ , "preid_trk_ecal_Deta/f");
	tree_->Branch("preid_trk_ecal_Dphi",	&preid_trk_ecal_Dphi_ , "preid_trk_ecal_Dphi/f");
	tree_->Branch("preid_e_over_p",			 	&preid_e_over_p_			, "preid_e_over_p/f");
	//stage 2, with GSF
	tree_->Branch("preid_gsf_success"			, &preid_gsf_success_     , "preid_gsf_success/O");
	tree_->Branch("preid_gsf_dpt"					, &preid_gsf_dpt_		  		, "preid_gsf_dpt/f");					
	tree_->Branch("preid_trk_gsf_chiratio", &preid_trk_gsf_chiratio_, "preid_trk_gsf_chiratio/f");
	tree_->Branch("preid_gsf_chi2red"     , &preid_gsf_chi2red_     , "preid_gsf_chi2red/f");     
	
	tree_->Branch("gsf_pt",				 	&gsf_pt_				   , "gsf_pt/f");
	tree_->Branch("gsf_eta",		 	  &gsf_eta_		       , "gsf_eta/f");
	tree_->Branch("gsf_phi",		 	  &gsf_phi_     		 , "gsf_phi/f");
	tree_->Branch("gsf_nhits",			&gsf_nhits_        , "gsf_nhits/f");
	tree_->Branch("gsf_dxy",			  &gsf_dxy_		  		 , "gsf_dxy/f");
	tree_->Branch("gsf_dxy_err",		&gsf_dxy_err_			 , "gsf_dxy_err/f");
	tree_->Branch("gsf_inp",  			&gsf_inp_  	    	 , "gsf_inp/f");
	tree_->Branch("gsf_outp",	  		&gsf_outp_	  	   , "gsf_outp/f");
	tree_->Branch("gsf_chi2red",    &gsf_chi2red_      , "gsf_chi2red/f"); 

	tree_->Branch("ele_pt",				 	&ele_pt_				   , "ele_pt/f");
	tree_->Branch("ele_eta",		 	  &ele_eta_		       , "ele_eta/f");
	tree_->Branch("ele_phi",		 	  &ele_phi_     		 , "ele_phi/f");
}

//fillers
void fill_evt(/* TODO */) {
	//TODO
}

void fill_gen(const GenParticleRef genp) {
	gen_pt_  = genp->pt();
	gen_eta_ = genp->eta();
	gen_phi_ = genp->phi();
}

void fill_gsf_trk(const GsfTrackRef trk, const reco::BeamSpot &spot) {
	gsf_pt_			 = trk->pt();
	gsf_eta_		 = trk->eta();
	gsf_phi_     = trk->phi();
	gsf_nhits_   = trk->found();
	gsf_dxy_		 = trk->dxy(spot);
	gsf_dxy_err_ = trk->dxyError();
	gsf_inp_  	 = sqrt(trk->innerMomentum().mag2());
	gsf_outp_	   = sqrt(trk->outerMomentum().mag2());
	gsf_chi2red_ = trk->normalizedChi2();
}

void fill_preid(const PreId &preid) {
}

void fill_ele(const GsfElectronRef ele) {
	ele_pt_			 = ele->pt();
	ele_eta_		 = ele->eta();
	ele_phi_     = ele->phi();
}

void fill_ktf_trk(const TrackRef trk) {
	trk_pt_			 = trk->pt();
	trk_eta_		 = trk->eta();
	trk_phi_     = trk->phi();
	trk_nhits_   = trk->found();
	trk_dxy_		 = trk->dxy(spot);
	trk_dxy_err_ = trk->dxyError();
	trk_inp_  	 = sqrt(trk->innerMomentum().mag2());
	trk_outp_	   = sqrt(trk->outerMomentum().mag2());
	trk_chi2red_ = trk->normalizedChi2();
	trk_high_purity_ = trackRef->quality(
		TrackBase::qualityByName("highPurity")
		);
}

