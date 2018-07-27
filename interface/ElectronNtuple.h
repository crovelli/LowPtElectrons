#ifndef LowPtElectrons_LowPtElectrons_ElectronNtuple
#define LowPtElectrons_LowPtElectrons_ElectronNtuple

class TTree;
namespace reco {
	class GsfTrackRef;
	class GsfElectronRef;
	class GenParticleRef;
	class PreId;
	class TrackRef;
	class BeamSpot;
}

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
	void fill_evt(/* TODO */);
	void fill_gen(const GenParticleRef genp);
	void fill_gsf_trk(const GsfTrackRef trk, const reco::BeamSpot &spot);
	void fill_preid(const PreId &preid);
	void fill_ele(const GsfElectronRef ele);
	void fill_ktf_trk(const TrackRef trk, const reco::BeamSpot &spot);

private:
	//only simple types (no arrays allowed, otherwise the reset() method fails;
	unsigned int lumi_ = 0;
	unsigned int run_ = 0;
	unsigned long long evt_ = 0;
	
	float gen_pt_ = -1;
	float gen_eta_ = -1;
	float gen_phi_ = -1;
	
	float trk_pt_ = -1.;
	float trk_eta_ = -1.;
	float trk_phi_ = -1.;
	float trk_nhits_ = -1;
	int   trk_high_purity_ = 0;
	float trk_dxy_ = -1;
	float trk_dxy_err_ = -1;
	float trk_inp_ = -1.;
	float trk_outp_ = -1.;
	float trk_chi2red_ = -1.;

	int   preid_ibin_ = -1;
	bool  preid_trk_ecal_match_ = false;
	float preid_bdtout_ = -1.;
	float preid_trk_ecal_Deta_ = -1.;
	float preid_trk_ecal_Dphi_ = -1.;
	float preid_e_over_p_ = -1.;
	//stage 2, with GSF
	bool  preid_gsf_success_ = false;
	float preid_gsf_dpt_ = -1.;
	float preid_trk_gsf_chiratio_ = -1.;
	float preid_gsf_chi2red_ = -1.;
	
	float gsf_pt_ = -1.;
	float gsf_eta_ = -1.;
	float gsf_phi_ = -1.;
	float gsf_nhits_ = -1;
	float gsf_dxy_ = -1;
	float gsf_dxy_err_ = -1;
	float gsf_inp_ = -1.;
	float gsf_outp_ = -1.;
	float gsf_chi2red_ = -1.;

	float ele_pt_ = -1.;
	float ele_eta_ = -1.;
	float ele_phi_ = -1.;
};

#endif
