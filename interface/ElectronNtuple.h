#ifndef LowPtElectrons_LowPtElectrons_ElectronNtuple
#define LowPtElectrons_LowPtElectrons_ElectronNtuple

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "FWCore/Framework/interface/Event.h"
class TTree;
/*namespace reco {
	class GsfTrackRef;
	class GsfElectronRef;
	class GenParticleRef;
	class PreId;
	class TrackRef;
	class BeamSpot;
}

namespace edm {
	class EventID;
	}*/

class ElectronNtuple {
	/*Small class to provide fillers and hide tree I/O*/
public:
	ElectronNtuple() {}
	
	void reset() {
		ElectronNtuple dummy; //create a new one
		*this = dummy; //use assignment to reset
	}

	void link_tree(TTree *tree);

	//fillers
	void fill_evt(const edm::EventID &id);
	void fill_gen(const reco::GenParticleRef genp);
	void fill_gsf_trk(const reco::GsfTrackRef trk, const reco::BeamSpot &spot);
	void fill_preid(const reco::PreId &preid, const reco::BeamSpot &spot);
	void fill_ele(const reco::GsfElectronRef ele);
	void fill_ktf_trk(const reco::TrackRef trk, const reco::BeamSpot &spot);

	void is_e(bool t=true) {is_e_=t;}
	void is_e_not_matched(bool t=true) {is_e_not_matched_=t;}
	void is_other(bool t=true) {is_other_=t;}

private:
	//only simple types (no arrays allowed, otherwise the reset() method fails;
	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	unsigned long long evt_ = 0;

	bool is_e_ = false;
	bool is_e_not_matched_ = false;
	bool is_other_ = false;
	
	// GEN electrons
	float gen_pt_ = -1;
	float gen_eta_ = -1;
	float gen_phi_ = -1;

	// KF tracks: kine
	float trk_pt_ = -1.;
	float trk_eta_ = -1.;
	float trk_phi_ = -1.;
	float trk_inp_ = -1.;
	float trk_outp_ = -1.;
	float trk_dpt_ = -1.;
	// KF tracks: quality
	float trk_nhits_ = -1;
	int   trk_high_purity_ = 0;
	float trk_chi2red_ = -1.;
	// KF tracks: displ
	float trk_dxy_ = -1;
	float trk_dxy_err_ = -1;
	float trk_dz_ = -1;
	float trk_dz_err_ = -1;

	// PreId: ECAL/track matching
	bool  preid_trk_ecal_match_ = false;
	float preid_e_over_p_ = -1.;
	float preid_trk_ecal_Deta_ = -1.;
	float preid_trk_ecal_Dphi_ = -1.;
	// PreId: GSF track parameters
	bool  preid_gsf_success_ = false;
	float preid_gsf_dpt_ = -1.;
	float preid_trk_gsf_chiratio_ = -1.;
	float preid_gsf_chi2red_ = -1.;
	// PreId: MVA output
	float preid_bdtout_ = -1.;
	int   preid_ibin_ = -1;
	
	// GSF tracks: kine
	float gsf_pt_ = -1.;
	float gsf_eta_ = -1.;
	float gsf_phi_ = -1.;
	float gsf_inp_ = -1.;
	float gsf_outp_ = -1.;
	float gsf_dpt_ = -1.;
	// GSF tracks: quality
	float gsf_nhits_ = -1;
	float gsf_chi2red_ = -1.;
	// GSF tracks: displ
	float gsf_dxy_ = -1;
	float gsf_dxy_err_ = -1;
	float gsf_dz_ = -1;
	float gsf_dz_err_ = -1;

	// GSF electrons
	float ele_pt_ = -1.;
	float ele_eta_ = -1.;
	float ele_phi_ = -1.;

};

#endif
