#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDNtuple.h"
#include "TRandom3.h"
#include "TTree.h"
#include <set>
#include <vector>

namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }

////////////////////////////////////////////////////////////////////////////////
//
class IDFeatures : public edm::EDAnalyzer {
  
public:
  
  explicit IDFeatures( const edm::ParameterSet& );
  ~IDFeatures() {}
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
private:
  
  void electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& prunedGenParticles,
		       const edm::Handle< edm::View<reco::Candidate> >& packedGenParticles,
		       std::set<reco::CandidatePtr>& electrons_from_B );

  void electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& genParticles,
		       std::set<reco::CandidatePtr>& electrons_from_B );
  
  bool isAncestor( const reco::Candidate* ancestor, 
		   const reco::Candidate* particle );

  void matchGenToTrk( const std::set<reco::CandidatePtr>& electrons_from_B,
		      const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
		      std::vector<reco::TrackRef>& other_trk );

  void matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		       const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
		       std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed );

  void matchSeedToGsf( const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
		       const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		       std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf );

  void matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< edm::View<reco::GsfElectron> >& electrons,
		      std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele );

  void matchGenToGsf( const std::set<reco::CandidatePtr>& electrons_from_B,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
		      std::vector<reco::GsfTrackRef>& other_gsf );

  void matchTrkToGsf( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf );
  
  reco::GsfTrackRef matchGsfToEgammaGsf( const reco::GsfTrackRef gsfTrack,
					 const edm::Handle< std::vector<reco::GsfTrack> >& egammaGsfTracks );

  inline bool filterTrack( reco::TrackRef trk ) {
    if ( !(trk->quality(reco::TrackBase::qualityByName("highPurity"))) ) { return false; }
    if ( trk->pt() < minTrackPt_ ) { return false; }
    return true;
  }
    
private:
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtuple ntuple_;
  bool check_from_B_;
  double dr_max_; // DeltaR matching
  double fakes_multiplier_;
  int isAOD_;
  bool hasGEN_;
  double minTrackPt_;

  // Available in AOD and MINIAOD
  const edm::EDGetTokenT<double> rho_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaSeedUnbiased_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaSeedPtbiased_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaIDLowPt_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaIDv2_;
  // Available in AOD only
  const edm::EDGetTokenT< std::vector<reco::Track> > ctfTracks_;
  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeeds_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > egammaGsfTracks_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > electrons_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > egammaElectrons_;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_;
  // Available in MINIAOD only
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > egammaGsfTracksMAOD_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > electronsMAOD_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > egammaElectronsMAOD_;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_;
  const edm::EDGetTokenT< edm::View<reco::Candidate> > packedGenParticles_;

  //const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;
  
};

////////////////////////////////////////////////////////////////////////////////
//
IDFeatures::IDFeatures( const edm::ParameterSet& cfg ) :
  tree_(0),
  ntuple_{},
  check_from_B_(cfg.getParameter<bool>("checkFromB")),
  dr_max_(cfg.getParameter<double>("drMax")),
  fakes_multiplier_(cfg.getParameter<double>("fakesMultiplier")),
  isAOD_(-1),
  hasGEN_(true),
  minTrackPt_(cfg.getParameter<double>("minTrackPt")),
  // Available in AOD and MINIAOD
  rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
  beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
  gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  mvaSeedUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVASeedUnbiased"))),
  mvaSeedPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVASeedPtbiased"))),
  mvaIDLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDLowPt"))),
  mvaIDv2_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDV2"))),
  // Available in AOD only
  ctfTracks_(consumes< std::vector<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
  eleSeeds_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))),
  egammaGsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("egammaGsfTracks"))),
  electrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons"))),
  egammaElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("egammaElectrons"))),
  genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
  // Available in MINIAOD only
  egammaGsfTracksMAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("egammaGsfTracks_mAOD"))),
  electronsMAOD_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons_mAOD"))),
  egammaElectronsMAOD_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("egammaElectrons_mAOD"))),
  prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
  packedGenParticles_(consumes< edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("packedGenParticles")))
  //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
  // Available in MINIAOD only
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
  }

////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDFeatures::beginRun( const edm::Run& run,
			   const edm::EventSetup& es ) {
  
}

////////////////////////////////////////////////////////////////////////////////
//
void IDFeatures::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // Get electrons, identify if data or MC and RECO/AOD or MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > electrons;
  if ( isAOD_ == -1 ) {
    event.getByToken(electrons_, electrons); // std::vector<reco::GsfElectron>
    if ( electrons.isValid() ) {
      isAOD_ = 1;
      std::cout << "File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(electronsMAOD_,electrons); // std::vector<pat::Electron>
      if ( electrons.isValid() ) { 
	isAOD_ = 0;
	std::cout << "File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << " failed to find a standard AOD or miniAOD particle collection " 
	  << std::endl;
      }
    }
  } else if ( isAOD_ == 1 ) {
    event.getByToken(electrons_, electrons); // std::vector<reco::GsfElectron>
  } else if ( isAOD_ == 0 ) {
    event.getByToken(electronsMAOD_,electrons); // std::vector<pat::Electron>
  } else {
    throw cms::Exception(" Invalid value for isAOD: ") 
      << isAOD_ 
      << std::endl;
  }

  // Get reco::Tracks
  edm::Handle< std::vector<reco::Track> > ctfTracks;
  if ( isAOD_ == 1 ) { event.getByToken(ctfTracks_, ctfTracks); }

  // Get ElectronSeeds
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeeds;
  if ( isAOD_ == 1 ) { event.getByToken(eleSeeds_, eleSeeds); }

  // Get GsfTracks for low pT
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracks;
  event.getByToken(gsfTracks_, gsfTracks);

  // Get GsfTracks for Egamma
  edm::Handle< std::vector<reco::GsfTrack> > egammaGsfTracks;
  if      ( isAOD_ == 1 ) { event.getByToken(egammaGsfTracks_, egammaGsfTracks); }
  else if ( isAOD_ == 0 ) { event.getByToken(egammaGsfTracksMAOD_, egammaGsfTracks); }

  // Get PF electrons
  edm::Handle< edm::View<reco::GsfElectron> > egammaElectrons;
  if      ( isAOD_ == 1 ) { event.getByToken(egammaElectrons_, egammaElectrons); }
  else if ( isAOD_ == 0 ) { event.getByToken(egammaElectronsMAOD_, egammaElectrons); }

  // Get GenParticles
  edm::Handle< edm::View<reco::GenParticle> > genParticles;
  if ( hasGEN_ ) {
    if ( isAOD_ == 1 ) { 
      event.getByToken(genParticles_, genParticles);
      if ( !(genParticles.isValid()) ) { 
	hasGEN_ = false;
	std::cout << "No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
      event.getByToken(prunedGenParticles_, genParticles);
      if ( !(genParticles.isValid()) ) { 
	hasGEN_ = false;
	std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
      }
    }
  }
    
  edm::Handle< edm::View<reco::Candidate> > packedGenParticles;
  if ( hasGEN_ && isAOD_ == 0 ) { event.getByToken(packedGenParticles_, packedGenParticles); }

  edm::Handle<double> rho;
  event.getByToken(rho_, rho);
  
  edm::Handle<reco::BeamSpot> beamspot;
  event.getByToken(beamspot_, beamspot);

  edm::Handle< edm::ValueMap<float> > mvaSeedUnbiased;
  event.getByToken(mvaSeedUnbiased_, mvaSeedUnbiased);

  edm::Handle< edm::ValueMap<float> > mvaSeedPtbiased;
  event.getByToken(mvaSeedPtbiased_, mvaSeedPtbiased);

  edm::Handle< edm::ValueMap<float> > mvaIDLowPt;
  event.getByToken(mvaIDLowPt_, mvaIDLowPt);
  
  edm::Handle< edm::ValueMap<float> > mvaIDv2;
  event.getByToken(mvaIDv2_, mvaIDv2);

  //edm::Handle< edm::ValueMap<float> > convVtxFitProb;
  //event.getByToken(convVtxFitProb_, convVtxFitProb);

  // Find GEN electrons from B decays
  std::set<reco::CandidatePtr> electrons_from_B;
  if ( hasGEN_ ) {
    if      ( isAOD_ == 1 ) { electronsFromB( genParticles, electrons_from_B ); }
    else if ( isAOD_ == 0 ) { electronsFromB( genParticles, packedGenParticles, electrons_from_B ); }
  }

  // Match GEN electrons to reco::Tracks
  std::map<reco::CandidatePtr, reco::TrackRef> gen2trk;
  std::vector<reco::TrackRef> other_trk;
  matchGenToTrk( electrons_from_B, ctfTracks, gen2trk, other_trk );

  // Match reco::Tracks to low pT ElectronSeeds
  std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed;
  matchTrkToSeed( ctfTracks, eleSeeds, trk2seed );

  // Match reco::Tracks to low pT ElectronSeeds
  std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf;
  matchSeedToGsf( eleSeeds, gsfTracks, seed2gsf );

  // Match reco::Tracks to low pT GsfTracks
  std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf;
  matchTrkToGsf( ctfTracks, gsfTracks, trk2gsf );

  // Match GEN electrons to low pT GsfTracks
  std::map<reco::CandidatePtr, reco::GsfTrackRef> gen2gsf;
  std::vector<reco::GsfTrackRef> other_gsf;
  matchGenToGsf( electrons_from_B, gsfTracks, gen2gsf, other_gsf );

  // Match GEN electrons to EGamma GsfTracks
  std::map<reco::CandidatePtr, reco::GsfTrackRef> gen2egammagsf;
  std::vector<reco::GsfTrackRef> other_egammagsf;
  matchGenToGsf( electrons_from_B, egammaGsfTracks, gen2egammagsf, other_egammagsf );
  
  // Match GsfTracks to low pT electrons
  std::map<reco::GsfTrackRef, reco::GsfElectronPtr> gsf2ele;
  matchGsfToEle( gsfTracks, electrons, gsf2ele );

  // Match GsfTracks to EGamma PF electrons
  std::map<reco::GsfTrackRef, reco::GsfElectronPtr> gsf2egammaele;
  matchGsfToEle( egammaGsfTracks, egammaElectrons, gsf2egammaele );
  
  //////////
  // Fill ntuple with signal electrons
  //////////

  for ( auto gen_electron : electrons_from_B ) {

    // Init and set truth label
    ntuple_.reset();
    ntuple_.is_e(true);
    ntuple_.is_other(false);

    // Fill Rho, Event, and GEN branches
    ntuple_.set_rho(*rho);
    ntuple_.fill_evt(event.id());
    ntuple_.fill_gen(gen_electron);

    // Check if GEN electron is matched to a reco::Track
    const auto& matched_trk = gen2trk.find(gen_electron);
    if ( matched_trk != gen2trk.end() && matched_trk->second.isNonnull() ) { 
      
      // Check if reco::Track is matched to an ElectronSeed
      const auto& matched_seed = trk2seed.find(matched_trk->second);
      if ( matched_seed != trk2seed.end() && matched_seed->second.isNonnull() ) { 
	
	// Check if ElectronSeed is matched to a reco::GsfTrack
	const auto& matched_gsf = seed2gsf.find(matched_seed->second);
	if ( matched_gsf != seed2gsf.end() && matched_gsf->second.isNonnull() ) { 
	  
	  ntuple_.fill_gsf(matched_gsf->second, *beamspot);

	  // Store Seed BDT discrimator values
	  float unbiased = (*mvaSeedUnbiased)[matched_gsf->second];
	  float ptbiased = (*mvaSeedPtbiased)[matched_gsf->second];
	  ntuple_.fill_seed( unbiased, ptbiased );
	  
	  // Check if GsfTrack is matched to an electron 
	  const auto& matched_ele = gsf2ele.find(matched_gsf->second);
	  if ( matched_ele != gsf2ele.end() ) {
	    
	    //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
	    float id_lowpt = -999.;
	    if ( mvaIDLowPt.isValid() && mvaIDLowPt->size() == electrons->size() ) {
	      id_lowpt = mvaIDLowPt->get( matched_ele->second.key() );
	    }

	    //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
	    float id_v2 = -999.;
	    if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
	      id_v2 = mvaIDv2->get( matched_ele->second.key() );
	    }
	    
	    //@@ dirty hack as is not in Event nor embedded in pat::Electron
	    float conv_vtx_fit_prob = -999.;
	    //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == electrons->size() ) {
	    //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
	    //}

	    ntuple_.fill_ele( matched_ele->second, id_lowpt, id_v2, conv_vtx_fit_prob, *rho );

	    //@@ Add SuperCluster vars?
	    //ntuple_.fill_supercluster(matched_ele->second);

	  } // Find reco::Track
	} // Find ElectronSeed
      } // Find GsfTrack
    } // Find GenParticle

    // Check if GEN electron is matched to a EGamma GsfTrack
    const auto& matched_egammagsf = gen2egammagsf.find(gen_electron);
    if ( matched_egammagsf != gen2egammagsf.end() ) { 

      // Record presence of EGAMMA EGamma GsfTrack
      ntuple_.has_egamma_gsf(true);

      // Check if GsfTrack is matched to an electron 
      const auto& matched_egammaele = gsf2egammaele.find(matched_egammagsf->second);
      if ( matched_egammaele != gsf2egammaele.end() ) {

	ntuple_.has_egamma_ele(true);

	//@@ dirty hack as ID is not in Event nor embedded in pat::Electron
	//float id_v2 = -999.;
	//if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
	//id_v2 = mvaIDv2->get( matched_ele->second.key() );
	//}
	
      } // Find GsfTrack

    } // Find GenParticle

    // Fill tree
    tree_->Fill();

  } // Fill ntuple with signal electrons

//  for ( auto gen : electrons_from_B ) {
//
//    // Init and set truth label
//    ntuple_.reset();
//    ntuple_.is_e(true);
//    ntuple_.is_other(false);
//
//    // Fill Rho, Event, and GEN branches
//    ntuple_.set_rho(*rho);
//    ntuple_.fill_evt(event.id());
//    ntuple_.fill_gen(gen);
//
//    // Check if GEN electron is matched to a GsfTrack
//    const auto& matched_gsf = gen2gsf.find(gen);
//    if ( matched_gsf != gen2gsf.end() ) { 
//
//      ntuple_.fill_gsf(matched_gsf->second, *beamspot);
//
//      // Store Seed BDT discrimator values
//      float unbiased = (*mvaSeedUnbiased)[matched_gsf->second];
//      float ptbiased = (*mvaSeedPtbiased)[matched_gsf->second];
//      ntuple_.fill_seed( unbiased, ptbiased );
//
//      // Check if GsfTrack is matched to an electron 
//      const auto& matched_ele = gsf2ele.find(matched_gsf->second);
//      if ( matched_ele != gsf2ele.end() ) {
//
//	//@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	float id_lowpt = -999.;
//	if ( mvaIDLowPt.isValid() && mvaIDLowPt->size() == electrons->size() ) {
//	  id_lowpt = mvaIDLowPt->get( matched_ele->second.key() );
//	}
//
//	//@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	float id_v2 = -999.;
//	if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
//	  id_v2 = mvaIDv2->get( matched_ele->second.key() );
//	}
//	
//	//@@ dirty hack as is not in Event nor embedded in pat::Electron
//	float conv_vtx_fit_prob = -999.;
//	//if ( convVtxFitProb.isValid() && convVtxFitProb->size() == electrons->size() ) {
//	//  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
//	//}
//
//	ntuple_.fill_ele( matched_ele->second, id_lowpt, id_v2, conv_vtx_fit_prob, *rho );
//
//	//@@ Add SuperCluster vars?
//	//ntuple_.fill_supercluster(matched_ele->second);
//
//      } // Find GsfTrack
//
//    } // Find GenParticle
//
//    // Check if GEN electron is matched to a EGamma GsfTrack
//    const auto& matched_egammagsf = gen2egammagsf.find(gen);
//    if ( matched_egammagsf != gen2egammagsf.end() ) { 
//
//      // Record presence of EGAMMA EGamma GsfTrack
//      ntuple_.has_egamma_gsf(true);
//
//      // Check if GsfTrack is matched to an electron 
//      const auto& matched_egammaele = gsf2egammaele.find(matched_egammagsf->second);
//      if ( matched_egammaele != gsf2egammaele.end() ) {
//
//	ntuple_.has_egamma_ele(true);
//
//	//@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	//float id_v2 = -999.;
//	//if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
//	//id_v2 = mvaIDv2->get( matched_ele->second.key() );
//	//}
//	
//      } // Find GsfTrack
//
//    } // Find GenParticle
//
//    // Fill tree
//    tree_->Fill();
//
//  } // Fill ntuple with signal electrons
//
  //////////
  // Fill ntuple with background pT electrons (prescaled by 'fakesMultiplier' configurable)
  //////////

  std::vector<int> indices;
  unsigned int nfakes = 0;
  while ( indices.size() < other_gsf.size() && // stop when all tracks considered
	  ( nfakes < fakes_multiplier_ || fakes_multiplier_ < 0 ) ) { // stop when fakesMultiplier is satisfied
    int index = int( gRandom->Rndm() * other_gsf.size() ); // pick a track at random
    if ( std::find( indices.begin(),
		    indices.end(),
		    index ) != indices.end() ) { continue; } // consider each track only once
    indices.push_back(index); // record tracks used
    nfakes++;
    reco::GsfTrackRef other = other_gsf.at(index);

    // Init and set truth label
    ntuple_.reset();
    ntuple_.is_e(false);
    ntuple_.is_other(true);

    // Fill Rho, Event, and GEN branches
    ntuple_.set_rho(*rho);
    ntuple_.fill_evt(event.id());
    
    ntuple_.fill_gsf(other, *beamspot);

    // Store Seed BDT discrimator values
    float unbiased = (*mvaSeedUnbiased)[other];
    float ptbiased = (*mvaSeedPtbiased)[other];
    ntuple_.fill_seed( unbiased, ptbiased );
    
    // Check if GsfTrack is matched to an electron 
    const auto& matched_ele = gsf2ele.find(other);
    if ( matched_ele != gsf2ele.end() ) {
      
      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float id_lowpt = -999.;
      if ( mvaIDLowPt.isValid() && mvaIDLowPt->size() == electrons->size() ) {
	id_lowpt = mvaIDLowPt->get( matched_ele->second.key() );
      }
      
      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float id_v2 = -999.;
      if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
        id_v2 = mvaIDv2->get( matched_ele->second.key() );
      }
      
      //@@ dirty hack as is not in Event nor embedded in pat::Electron
      float conv_vtx_fit_prob = -999.;
      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == electrons->size() ) {
      //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
      //}
      
      ntuple_.fill_ele( matched_ele->second, id_lowpt, id_v2, conv_vtx_fit_prob, *rho );
      
      //@@ Add SuperCluster vars?
      //ntuple_.fill_supercluster(matched_ele->second);

    } // Find GsfTrack

    // Find closest EGamma GsfTrack 
    reco::GsfTrackRef matched_egammagsf = matchGsfToEgammaGsf( other, egammaGsfTracks );
    if ( matched_egammagsf.isNonnull() ) {

      ntuple_.has_egamma_gsf(true);

      // Check if GsfTrack is matched to a EGAMMA electron 
      const auto& matched_egammaele = gsf2egammaele.find(matched_egammagsf);
      if ( matched_egammaele != gsf2egammaele.end() ) {
	ntuple_.has_egamma_ele(true);
      }

    } // Find GsfTrack

    // Fill tree
    tree_->Fill();

  } // Fill ntuple with background low pT electrons
  
}

////////////////////////////////////////////////////////////////////////////////
// Assumes pruned and packed GenParticles from MINIAOD
void IDFeatures::electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& prunedGenParticles,
				 const edm::Handle< edm::View<reco::Candidate> >& packedGenParticles,
				 std::set<reco::CandidatePtr>& electrons_from_B ) {

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Examples
  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/packedGenParticles_cfi.py
  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py

  electrons_from_B.clear();
  
  // Iterate through ("status 1") packed candidates
  for ( size_t ipa = 0; ipa < packedGenParticles->size(); ipa++ ) {
    reco::CandidatePtr packed(packedGenParticles, ipa);
    
    // Find GEN electrons
    bool is_ele = std::abs(packed->pdgId()) == 11; // && packed->isLastCopy() 

    // Search for ancestors only for GEN electrons ... 
    if ( !is_ele ) { continue; }
    
    // Does GEN ele come from B decay?
    bool comes_from_B = false;
    for ( size_t ipr = 0; ipr < prunedGenParticles->size(); ipr++ ) {
      reco::GenParticlePtr pruned(prunedGenParticles, ipr);
      comes_from_B = 
	std::abs(pruned->pdgId()) > 510 &&
	std::abs(pruned->pdgId()) < 546 &&
	packed->mother(0) != nullptr && 
	isAncestor( &*pruned, &*packed );

      // Is coming from a B
      if ( is_ele && ( comes_from_B || !check_from_B_ ) ) {
	//electrons_from_B.insert(packed);
      }
      
    }
    
    //std::cout << "is_ele: " << is_ele << " comes_from_B: " << comes_from_B << std::endl;
    
  } // packedGenParticles loop
  
}

////////////////////////////////////////////////////////////////////////////////
// Assumes GenParticles from RECO/AOD
void IDFeatures::electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& genParticles,
				 std::set<reco::CandidatePtr>& electrons_from_B ) {
  
  electrons_from_B.clear();
  for ( size_t idx = 0; idx < genParticles->size(); idx++ ) {
    reco::CandidatePtr genp(genParticles, idx);
    
    // Last copy of GEN electron 
    bool is_ele = std::abs(genp->pdgId()) == 11; // && genp->isLastCopy();
    
    // Does GEN ele comes from B decay?
    bool comes_from_B = 
      genp->numberOfMothers() >= 1 &&
      std::abs(genp->mother()->pdgId()) > 510 &&
      std::abs(genp->mother()->pdgId()) < 546;
    
    if (!comes_from_B &&                                                      // from B decay
	genp->numberOfMothers() >= 1 && genp->mother() &&                     // has mother
	std::abs(genp->mother()->pdgId()) == 443 &&                           // mother is J/psi
	genp->mother()->numberOfMothers() >= 1 && genp->mother()->mother() && // has grandmother
	std::abs(genp->mother()->mother()->pdgId()) > 510 &&                  // grandmother is B
	std::abs(genp->mother()->mother()->pdgId()) < 546 ) {
      comes_from_B = true;
    }
    
    // is coming from a B
    if ( is_ele && ( comes_from_B || !check_from_B_ ) ) {
      electrons_from_B.insert(genp);
    }
    
  } // genParticles loop

}

////////////////////////////////////////////////////////////////////////////////
// 
bool IDFeatures::isAncestor( const reco::Candidate* ancestor, 
			     const reco::Candidate* particle ) {
  
  // particle is already the ancestor
  if ( ancestor == particle ) return true;
  
  // otherwise loop on mothers, if any and return true if the ancestor is found
  for ( size_t i = 0; i < particle->numberOfMothers(); i++ ) {
    if ( isAncestor(ancestor,particle->mother(i)) ) return true;
  }
  
  // if we did not return yet, then particle and ancestor are not relatives
  return false;
  
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchGenToTrk( const std::set<reco::CandidatePtr>& electrons_from_B,
				const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
				std::vector<reco::TrackRef>& other_trk ) {
  
  gen2trk.clear();
  other_trk.clear();
  std::set<reco::TrackRef> matched;

  // Match ctfTracks to GEN ele from B decays
  for ( const auto& gen : electrons_from_B ) {
    size_t best_idx = 666;
    double min_dr2 = dr_max_*dr_max_;
    for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
      reco::TrackRef trk(ctfTracks,idx);
      if ( !filterTrack(trk) ) { continue; }
      double dr2 = reco::deltaR2(trk->eta(),
				 trk->phi(),
				 gen->eta(), 
				 gen->phi() );
      if ( dr2 < min_dr2 ) {
	min_dr2 = dr2;
	best_idx = idx;
      }
    }
    if ( min_dr2 < dr_max_*dr_max_ ) {
      reco::TrackRef trk(ctfTracks,best_idx);
      gen2trk.insert( std::pair<reco::CandidatePtr, reco::TrackRef>( gen, trk ) );
      matched.insert(trk);
    }
  }

  // Identify ctfTracks from background, i.e. not matched to a GEN electron
  other_trk.reserve( ctfTracks->size() - matched.size() );
  for ( size_t idx = 0; idx < ctfTracks->size(); idx++ ) {
    reco::TrackRef trk(ctfTracks,idx);
    if ( !filterTrack(trk) ) { continue; }
    if ( matched.find(trk) == matched.end() ) {
      other_trk.push_back(trk);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				 const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
				 std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed ) {
  
  trk2seed.clear();
  std::set<reco::TrackRef> matched;

  // Store ctfTrack -> eleSeed link
  for ( size_t idx = 0; idx < eleSeeds->size(); ++idx ) {
    reco::ElectronSeedRef ele_seed(eleSeeds, idx);
    if ( ele_seed.isNull() ) { continue; }
    reco::TrackRef trk = ele_seed->ctfTrack();
    if ( trk.isNull() || trk.id() != ctfTracks.id() ) { 
      std::cout << "Collections do not match, or null track Ref!!" << std::endl;
      continue;
    }
    trk2seed[trk] = ele_seed;
    matched.insert(trk);
  }

  // Add null eleSeed Ref for unmatched ctfTracks
  for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
    reco::TrackRef trk(ctfTracks, idx);
    if ( std::find( matched.begin(), matched.end(), trk ) != matched.end() ) { continue; }
    trk2seed[trk] = reco::ElectronSeedRef();
  }

//  int cntr = 0;
//  typedef std::map<reco::TrackRef, reco::ElectronSeedRef>::const_iterator Iter;
//  for ( Iter iter = trk2seed.begin(); iter != trk2seed.end(); ++iter ) {
//    if ( iter->second.isNonnull() ) { cntr++; }
//  }
//  std::cout << " # ctfTracks: " << ctfTracks->size()
//	    << " # eleSeeds: " << eleSeeds->size()
//	    << " # map entries: " << trk2seed.size()
//	    << " # nonnull ptrs: " << cntr
//	    << " # matched: " << matched.size()
//	    << std::endl;
  
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchSeedToGsf( const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
				 const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				 std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf ) {
  
  seed2gsf.clear();
  std::set<reco::ElectronSeedRef> matched;

  // Store eleSeed -> gsfTrack link
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackRef gsf(gsfTracks, idx);
    if ( gsf.isNull() ) { continue; }
    edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
    if ( seed.isNull() ) { continue; }
    reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
    if ( ele_seed.isNull() || ele_seed.id() != eleSeeds.id() ) { 
      std::cout << "Collections do not match, or null seed Ref!!" << std::endl;
      continue;
    }
    seed2gsf[ele_seed] = gsf;
    matched.insert(ele_seed);
  }

  // Add null gsfTrack Ref for unmatched eleSeeds
  for ( size_t idx = 0; idx < eleSeeds->size(); ++idx ) {
    reco::ElectronSeedRef seed(eleSeeds, idx);
    if ( std::find( matched.begin(), matched.end(), seed ) != matched.end() ) { continue; }
    seed2gsf[seed] = reco::GsfTrackRef();
  }

//  int cntr = 0;
//  typedef std::map<reco::ElectronSeedRef, reco::GsfTrackRef>::const_iterator Iter;
//  for ( Iter iter = seed2gsf.begin(); iter != seed2gsf.end(); ++iter ) {
//    if ( iter->second.isNonnull() ) { cntr++; }
//  }
//  std::cout << " # eleSeeds: " << eleSeeds->size()
//	    << " # gsfTracks: " << gsfTracks->size()
//	    << " # map entries: " << seed2gsf.size()
//	    << " # nonnull ptrs: " << cntr
//	    << " # matched: " << matched.size()
//	    << std::endl;
  
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				const edm::Handle< edm::View<reco::GsfElectron> >& electrons,
				std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele ) {
  gsf2ele.clear();
  for ( size_t idx = 0; idx < electrons->size(); ++idx ) {
    reco::GsfElectronPtr ele(electrons, idx);
    reco::GsfTrackRef gsf = ele->gsfTrack();
    if ( gsf2ele.find(gsf) != gsf2ele.end() ) {
      std::cout << "THIS SHOULD NEVER HAPPEN! Multiple low pT electrons matched to the same GSFTrack?!"
		<< std::endl;
    } else {
      gsf2ele.insert( std::pair<reco::GsfTrackRef, reco::GsfElectronPtr>(gsf, ele) );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchTrkToGsf( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf ) {
  
  trk2gsf.clear();
  std::set<reco::TrackRef> matched;

  // Store ctfTrack -> gsfTrack link
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackRef gsf(gsfTracks, idx);
    if ( gsf.isNull() ) { continue; }
    edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
    if ( seed.isNull() ) { continue; }
    reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
    if ( ele_seed.isNull() ) { continue; }
    reco::TrackRef trk = ele_seed->ctfTrack();
    if ( trk.isNull() || trk.id() != ctfTracks.id() ) { 
      std::cout << "Collections do not match, or null track Ref!!" << std::endl;
      continue;
    }
    trk2gsf[trk] = gsf;
    matched.insert(trk);
  }

  // Add null gsfTrack Ref for unmatched ctfTracks
  for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
    reco::TrackRef trk(ctfTracks, idx);
    if ( std::find( matched.begin(), matched.end(), trk ) != matched.end() ) { continue; }
    trk2gsf[trk] = reco::GsfTrackRef();
  }

//  int cntr = 0;
//  typedef std::map<reco::TrackRef, reco::GsfTrackRef>::const_iterator Iter;
//  for ( Iter iter = trk2gsf.begin(); iter != trk2gsf.end(); ++iter ) {
//    if ( iter->second.isNonnull() ) { cntr++; }
//  }
//  std::cout << " # ctfTracks: " << ctfTracks->size()
//	    << " # gsfTracks: " << gsfTracks->size()
//	    << " # map entries: " << trk2gsf.size()
//	    << " # nonnull ptrs: " << cntr
//	    << " # matched: " << matched.size()
//	    << std::endl;
  
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchGenToGsf( const std::set<reco::CandidatePtr>& electrons_from_B,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
				std::vector<reco::GsfTrackRef>& other_gsf ) {
  
  gen2gsf.clear();
  other_gsf.clear();
  std::set<reco::GsfTrackRef> matched;

  // Match GsfTracks to GEN ele from B decays
  for ( const auto& gen : electrons_from_B ) {
    size_t best_idx = 666;
    double min_dr2 = dr_max_*dr_max_;
    for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
      reco::GsfTrackRef gsf(gsfTracks,idx);
      double dr2 = reco::deltaR2(gsf->etaMode(), //@@ use mode value for eta ?!
				 gsf->phiMode(), //@@ use mode value for phi ?!
				 gen->eta(), 
				 gen->phi() );
      if ( dr2 < min_dr2 ) {
	min_dr2 = dr2;
	best_idx = idx;
      }
    }
    if ( min_dr2 < dr_max_*dr_max_ ) {
      reco::GsfTrackRef gsf(gsfTracks,best_idx);
      gen2gsf.insert( std::pair<reco::CandidatePtr, reco::GsfTrackRef>( gen, gsf ) );
      matched.insert(gsf);
    }
  }

  // Identify GsfTracks from background, i.e. not matched to a GEN electron
  other_gsf.reserve( gsfTracks->size() - matched.size() );
  for ( size_t idx = 0; idx < gsfTracks->size(); idx++ ) {
    reco::GsfTrackRef gsf(gsfTracks,idx);
    if ( matched.find(gsf) == matched.end() ) {
      other_gsf.push_back(gsf);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// 
reco::GsfTrackRef IDFeatures::matchGsfToEgammaGsf( const reco::GsfTrackRef gsfTrack,
						   const edm::Handle< std::vector<reco::GsfTrack> >& egammaGsfTracks ) {
  size_t best_idx = -1;
  double min_dr2 = dr_max_*dr_max_;
  for ( size_t idx = 0; idx < egammaGsfTracks->size(); ++idx ) {
    reco::GsfTrackRef egammagsf(egammaGsfTracks,idx);
    double dr2 = reco::deltaR2(egammagsf->etaMode(), //@@ use mode value for eta ?!
			       egammagsf->phiMode(), //@@ use mode value for phi ?!
			       gsfTrack->etaMode(),  //@@ use mode value for eta ?!
			       gsfTrack->phiMode() ); //@@ use mode value for phi ?!
    if ( dr2 < min_dr2 ) {
      min_dr2 = dr2;
      best_idx = idx;
    }
  }
  if ( min_dr2 < dr_max_*dr_max_ ) {
    return reco::GsfTrackRef(egammaGsfTracks,best_idx);
  } else {
    return reco::GsfTrackRef(egammaGsfTracks.id()); // null Ref
  }  
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDFeatures);
