#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
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
#include <math.h>

namespace reco { typedef edm::RefToBase<CaloCluster> CaloClusterRef; }
namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

////////////////////////////////////////////////////////////////////////////////
//
class IDFeatures : public edm::EDAnalyzer {
  
public:
  
  explicit IDFeatures( const edm::ParameterSet& );
  ~IDFeatures() {}
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
private:

  // Add GEN info for signal electrons
  void fill( const edm::Event& event,
	     const edm::EventSetup& setup,
	     std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
	     std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
	     std::map<reco::CandidatePtr, reco::GsfElectronPtr>& gen2ele,
	     std::map<reco::CandidatePtr, float>& gen2trk_dr2,
	     std::map<reco::CandidatePtr, float>& gen2gsf_dr2_mode,
	     std::map<reco::CandidatePtr, float>& gen2gsf_dr2,
	     std::map<reco::CandidatePtr, float>& gen2ele_dr2,
	     const std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed,
	     const std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf,
	     const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
	     bool is_egamma = false );

  // Add signal electrons or fake, low pT or PF
  void fill( bool is_signal_ele,
	     const edm::Event& event,
	     const edm::EventSetup& setup,
	     const std::vector<reco::TrackRef>& tracks,
	     const std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed,
	     const std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf,
	     const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
	     bool is_egamma = false,
	     bool reset_and_fill = true );

  void fill( bool is_signal_ele,
	     const edm::Event& event,
	     const edm::EventSetup& setup,
	     const std::vector<reco::GsfTrackRef>& tracks,
	     const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
	     bool is_egamma,
	     bool reset_and_fill );
  
  //////////

//  void electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& genParticles,
//		       std::set<reco::CandidatePtr>& electrons_from_B );

  void electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& genParticles,
		       std::set<reco::CandidatePtr>& electrons_from_B );

//  void electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& prunedGenParticles,
//		       const edm::Handle< edm::View<reco::Candidate> >& packedGenParticles,
//		       std::set<reco::CandidatePtr>& electrons_from_B );
  
  bool isAncestor( const reco::Candidate* ancestor, 
		   const reco::Candidate* particle );

  //////////

  void matchGenToTrk( const std::set<reco::CandidatePtr>& electrons_from_B,
		      const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
		      std::map<reco::CandidatePtr, float>& gen2trk_dr2,
		      std::vector<reco::TrackRef>& other_trk );

//  void matchGenToTrk( const std::set<reco::CandidatePtr>& electrons_from_B,
//		      const edm::Handle<pat::PackedCandidateCollection>& packedCands,
//		      const edm::Handle<pat::PackedCandidateCollection>& lostTracks,
//		      std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
//		      std::vector<reco::TrackRef>& other_trk );

  void matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		       const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
		       std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed );

  void matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		       const std::vector<reco::ElectronSeedRef>& eleSeeds,
		       std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed );

  void matchSeedToGsf( const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
		       const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		       std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf );

  void matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
		      std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele );

  //////////

  void matchGenToGsf( const std::set<reco::CandidatePtr>& electrons_from_B,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
		      std::map<reco::CandidatePtr, float>& gen2gsf_dr2_mode,
		      std::map<reco::CandidatePtr, float>& gen2gsf_dr2,
		      std::vector<reco::GsfTrackRef>& other_gsf );

  void matchGenToEle( const std::set<reco::CandidatePtr>& electrons_from_B,
		      const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
		      std::map<reco::CandidatePtr, reco::GsfElectronPtr>& gen2ele,
		      std::map<reco::CandidatePtr, float>& gen2ele_dr2,
		      std::vector<reco::GsfElectronPtr>& other_ele );

  void matchTrkToGsf( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::map<reco::TrackRef, reco::GsfTrackRef>& trk2gsf );

  void checkGsfToTrk( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< edm::Association< std::vector<reco::Track> > >& gsf2trk,
		      const std::map<reco::TrackRef, reco::GsfTrackRef>& trk2gsf );
  
  //////////

  reco::GsfTrackRef matchGsfToEgammaGsf( const reco::GsfTrackRef gsfTrack,
					 const edm::Handle< std::vector<reco::GsfTrack> >& egammaGsfTracks );

  inline bool filterTrack( reco::TrackRef trk ) {
    if ( !(trk->quality(reco::TrackBase::qualityByName("highPurity"))) ) { return false; }
    if ( trk->pt() < minTrackPt_ ) { return false; }
    return true;
  }

  class Match {
  public:
    Match( reco::CandidatePtr gen, int itrk, double dr2 ) {
      gen_ = gen; idx_ = itrk; dr2_ = dr2;
    };
    reco::CandidatePtr gen_;
    int idx_;
    double dr2_;
    static bool compare_by_dr2( const Match& a, const Match& b ) {
      return a.dr2_ < b.dr2_;
    };
  };
  
private:

  // Misc

  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double dr_max_; // Max DeltaR value considered
  double dr_threshold_; // Threshold for DeltaR matching
  double prescale_;
  int isAOD_;
  bool hasGEN_;
  double minTrackPt_;

  // Generic collections

  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticles_; // AOD
  edm::Handle< std::vector<reco::GenParticle> > genParticlesH_;

//  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // mAOD
//  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;
//
//  const edm::EDGetTokenT< edm::View<reco::Candidate> > packedGenParticles_; // mAOD
//  edm::Handle< edm::View<reco::Candidate> > packedGenParticlesH_;

  const edm::EDGetTokenT< std::vector<reco::Track> > ctfTracks_; // AOD
  edm::Handle< std::vector<reco::Track> > ctfTracksH_;

  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;

  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;

  //noZS::EcalClusterLazyTools ecalTools_;
  
  const edm::EDGetTokenT<reco::SuperClusterCollection> barrelSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> barrelSCsH_;

  const edm::EDGetTokenT<reco::SuperClusterCollection> endcapSCs_; // AOD
  edm::Handle<reco::SuperClusterCollection> endcapSCsH_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> packedCands_; // mAOD
  edm::Handle<pat::PackedCandidateCollection> packedCandsH_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracks_; // mAOD
  edm::Handle<pat::PackedCandidateCollection> lostTracksH_;

  // Low pT collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeeds_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsEcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsEcalH_;

  const edm::EDGetTokenT< std::vector<reco::PreId> > preIdsHcal_; // AOD
  edm::Handle< std::vector<reco::PreId> > preIdsHcalH_;

  const edm::EDGetTokenT< edm::ValueMap<reco::PreIdRef> > preIdRefs_; // AOD
  edm::Handle< edm::ValueMap<reco::PreIdRef> > preIdRefsH_;

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; // AOD and mAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_MAOD_; // mAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;

  const edm::EDGetTokenT< edm::Association<reco::TrackCollection> > gsfTrackLinks_; // AOD
  edm::Handle<edm::Association<reco::TrackCollection> > gsfTrackLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; // mAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; // mAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; //  on the fly?
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;

//  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValue_; // on the fly?
//  edm::Handle< edm::ValueMap<float> > mvaValueH_;

//  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaId_; // on the fly?
//  edm::Handle< edm::ValueMap<bool> > mvaIdH_;

  // EGamma collections

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeedsEGamma_; // AOD
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsEGammaH_; // AOD

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_; // AOD
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_MAOD_; // mAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksEGammaH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectronsEGamma_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectronsEGamma_MAOD_; // mAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;

  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaIdEGamma_; // on the fly?
  edm::Handle< edm::ValueMap<bool> > mvaIdEGammaH_;

  // Conversions

  //@@ const edm::EDGetTokenT< edm::ValueMap<float> > convVtxFitProb_;
  
};

////////////////////////////////////////////////////////////////////////////////
//
IDFeatures::IDFeatures( const edm::ParameterSet& cfg ) :
  tree_(0),
  ntuple_{},
  verbose_(cfg.getParameter<int>("verbose")),
  check_from_B_(cfg.getParameter<bool>("checkFromB")),
  dr_max_(cfg.getParameter<double>("drMax")),
  dr_threshold_(cfg.getParameter<double>("drThreshold")),
  prescale_(cfg.getParameter<double>("prescale")),
  isAOD_(-1),
  hasGEN_(true),
  minTrackPt_(cfg.getParameter<double>("minTrackPt")),
  // Generic collections
  rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
  rhoH_(),
  beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
  beamspotH_(),
  genParticles_(consumes< std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
  genParticlesH_(),
  //genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
  //prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
  //genParticlesH_(),
  //packedGenParticles_(consumes< edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("packedGenParticles"))),
  //packedGenParticlesH_(),
  ctfTracks_(consumes< std::vector<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
  ctfTracksH_(),
  ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
  ebRecHitsH_(),
  eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
  eeRecHitsH_(),
  //ecalTools_(),
  barrelSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("barrelSuperClusters"))),
  barrelSCsH_(),
  endcapSCs_(consumes<reco::SuperClusterCollection>(cfg.getParameter<edm::InputTag>("endcapSuperClusters"))),
  endcapSCsH_(),
  packedCands_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("packedCandidates"))),
  packedCandsH_(),
  lostTracks_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
  lostTracksH_(),
  // Low pT collections
  eleSeeds_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeeds"))),
  eleSeedsH_(),
  preIdsEcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsEcal"))),
  preIdsEcalH_(),
  preIdsHcal_(consumes< std::vector<reco::PreId> >(cfg.getParameter<edm::InputTag>("preIdsHcal"))),
  preIdsHcalH_(),
  preIdRefs_(consumes< edm::ValueMap<reco::PreIdRef> >(cfg.getParameter<edm::InputTag>("preIdRefs"))),
  preIdRefsH_(),
  gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  gsfTracksH_(),
  gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
  gsfElectrons_MAOD_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons_MAOD"))),
  gsfElectronsH_(),
  gsfTrackLinks_(consumes< edm::Association<reco::TrackCollection> >(cfg.getParameter<edm::InputTag>("gsfTrackLinks"))),
  gsfTrackLinksH_(),
  packedCandLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("packedCandLinks"))),
  packedCandLinksH_(),
  lostTrackLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("lostTrackLinks"))),
  lostTrackLinksH_(),
  mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
  mvaUnbiasedH_(),
  mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
  mvaPtbiasedH_(),
  mvaValueLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPt"))),
  mvaValueLowPtH_(),
  //mvaValue_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValue"))),
  //mvaValueH_(),
  //mvaId_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("mvaId"))),
  //mvaIdH_(),
  // EGamma collections
  eleSeedsEGamma_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeedsEGamma"))),
  eleSeedsEGammaH_(),
  gsfTracksEGamma_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma"))),
  gsfTracksEGamma_MAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma_MAOD"))),
  gsfTracksEGammaH_(),
  gsfElectronsEGamma_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma"))),
  gsfElectronsEGamma_MAOD_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma_MAOD"))),
  gsfElectronsEGammaH_(),
  mvaValueEGamma_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueEGamma"))),
  mvaValueEGammaH_(),
  mvaIdEGamma_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("mvaIdEGamma"))),
  mvaIdEGammaH_()
  // Conversions
  //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
  }

////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDFeatures::beginRun( const edm::Run& run,
			   const edm::EventSetup& es ) {
  //@@ ?
}

////////////////////////////////////////////////////////////////////////////////
//
void IDFeatures::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_); // std::vector<reco::GsfElectron>
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(gsfElectrons_MAOD_,gsfElectronsH_); // std::vector<pat::Electron>
      if ( gsfElectronsH_.isValid() ) { 
	isAOD_ = 0;
	std::cout << "File contains MINIAOD data tier!" << std::endl;
      } else {
	throw cms::Exception(" Collection not found: ") 
	  << " failed to find a standard AOD or miniAOD particle collection " 
	  << std::endl;
      }
    }
  } else if ( isAOD_ == 1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_); // std::vector<reco::GsfElectron>
  } else if ( isAOD_ == 0 ) {
    event.getByToken(gsfElectrons_MAOD_,gsfElectronsH_); // std::vector<pat::Electron>
  } else {
    throw cms::Exception(" Invalid value for isAOD: ") 
      << isAOD_ 
      << std::endl;
  }

  // Generic collections 
  event.getByToken(rho_, rhoH_);
  event.getByToken(beamspot_, beamspotH_);
  
  // GEN particles
  if ( hasGEN_ ) {
    if ( isAOD_ == 1 ) { 
      event.getByToken(genParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	hasGEN_ = false;
	std::cout << "No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
//      event.getByToken(prunedGenParticles_, genParticlesH_);
//      if ( !(genParticlesH_.isValid()) ) { 
//	hasGEN_ = false;
//	std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
//      }
    }
  }
//  if ( hasGEN_ && isAOD_ == 0 ) { event.getByToken(packedGenParticles_, packedGenParticlesH_); }

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
  }

  // RecHits and SuperClusters
  if ( isAOD_ == 1 ) { 
    event.getByToken(ebRecHits_, ebRecHitsH_);
    event.getByToken(eeRecHits_, eeRecHitsH_);
    event.getByToken(barrelSCs_, barrelSCsH_);
    event.getByToken(endcapSCs_, endcapSCsH_);
    //ecalTools_ = noZS::EcalClusterLazyTools(event, setup, ebRecHitsH_, eeRecHitsH_);
  }

  // ElectronSeeds and PreIds
  if ( isAOD_ == 1 ) { 
    event.getByToken(eleSeeds_, eleSeedsH_); 
    event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); 
    event.getByToken(preIdsEcal_, preIdsEcalH_); 
    event.getByToken(preIdsHcal_, preIdsHcalH_); 
    event.getByToken(preIdRefs_, preIdRefsH_); 
  }

  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  if      ( isAOD_ == 1 ) { 
    event.getByToken(gsfTrackLinks_, gsfTrackLinksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCandLinks_, packedCandLinksH_); 
    event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  }

  // EGamma collections 
  if      ( isAOD_ == 1 ) { event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfTracksEGamma_, gsfTracksEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(gsfTracksEGamma_MAOD_, gsfTracksEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfElectronsEGamma_, gsfElectronsEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(gsfElectronsEGamma_MAOD_, gsfElectronsEGammaH_); }

  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
  //event.getByToken(mvaValue_, mvaValueH_);
  //event.getByToken(mvaId_, mvaIdH_);
  event.getByToken(mvaValueEGamma_, mvaValueEGammaH_);
  event.getByToken(mvaIdEGamma_, mvaIdEGammaH_);

  // Conversions
  //event.getByToken(convVtxFitProb_, convVtxFitProbH_);

  //////////

  // Find GEN electrons from B decays
  std::set<reco::CandidatePtr> electrons_from_B;
  if ( hasGEN_ ) {
    if      ( isAOD_ == 1 ) { electronsFromB( genParticlesH_, electrons_from_B ); }
    //else if ( isAOD_ == 0 ) { electronsFromB( genParticlesH_, packedGenParticlesH_, electrons_from_B ); }
  }

  // Match GEN electrons to reco::Tracks
  std::map<reco::CandidatePtr, reco::TrackRef> gen2trk;
  std::map<reco::CandidatePtr, float> gen2trk_dr2;
  std::vector<reco::TrackRef> other_trk;
  if ( isAOD_ == 1 ) matchGenToTrk( electrons_from_B, ctfTracksH_, 
				    gen2trk, gen2trk_dr2, other_trk );

//  // Match GEN electrons to reco::SuperClusters
//  std::map<reco::CandidatePtr, reco::TrackRef> gen2sc;
//  std::vector<reco::TrackRef> other_sc;
//  if ( isAOD_ == 1 ) matchGenToSC( electrons_from_B, barrelSCsH_, endcapSCsH_, gen2sc, other_sc );

  //////////
  // Low pT links

  // Match GEN electrons to low pT GsfTracks
  std::map<reco::CandidatePtr, reco::GsfTrackRef> gen2gsf;
  std::map<reco::CandidatePtr, float> gen2gsf_dr2_mode;
  std::map<reco::CandidatePtr, float> gen2gsf_dr2;
  std::vector<reco::GsfTrackRef> other_gsf;
  if ( isAOD_ == 1 ) matchGenToGsf( electrons_from_B, gsfTracksH_, 
				    gen2gsf, 
				    gen2gsf_dr2_mode, gen2gsf_dr2, 
				    other_gsf );

  // Match GEN electrons to EGamma (PF) GsfElectrons
  std::map<reco::CandidatePtr, reco::GsfElectronPtr> gen2ele;
  std::map<reco::CandidatePtr, float> gen2ele_dr2;
  std::vector<reco::GsfElectronPtr> other_ele;
  if ( isAOD_ == 1 ) matchGenToEle( electrons_from_B, gsfElectronsH_, 
				    gen2ele, 
				    gen2ele_dr2, 
				    other_ele );
  
  //@@
  // Match reco::Tracks to low pT ElectronSeeds
  std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed;
  if ( isAOD_ == 1 ) matchTrkToSeed( ctfTracksH_, eleSeedsH_, trk2seed );

  // Match low pT ElectronSeeds to low pT GsfTracks
  std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf;
  if ( isAOD_ == 1 ) matchSeedToGsf( eleSeedsH_, gsfTracksH_, seed2gsf );
  
  // Match low pT GsfTracks to low pT electrons
  std::map<reco::GsfTrackRef, reco::GsfElectronPtr> gsf2ele;
  if ( isAOD_ == 1 ) matchGsfToEle( gsfTracksH_, gsfElectronsH_, gsf2ele );

  if ( verbose_ > 0 ) {
    std::cout << "Collection sizes: " << std::endl
	      << "  ctfTracksH_->size(): " << int(ctfTracksH_.isValid() ? ctfTracksH_->size() : -1) << std::endl
	      << "  eleSeedsH_->size(): " << int(eleSeedsH_.isValid() ? eleSeedsH_->size() : -1) << std::endl
	      << "  preIdsEcalH_->size(): " << int(preIdsEcalH_.isValid() ? preIdsEcalH_->size() : -1) << std::endl
	      << "  preIdsHcalH_->size(): " << int(preIdsHcalH_.isValid() ? preIdsHcalH_->size() : -1) << std::endl
	      << "  preIdRefsH_->size(): " << int(preIdRefsH_.isValid() ? preIdRefsH_->size() : -1) << std::endl
	      << "  gsfTracksH_->size(): " << int(gsfTracksH_.isValid() ? gsfTracksH_->size() : -1) << std::endl
	      << "  gsfElectronsH_->size(): " << int(gsfElectronsH_.isValid() ? gsfElectronsH_->size() : -1) << std::endl
	      << "  electrons_from_B.size(): " << electrons_from_B.size() << std::endl
	      << "  gen2trk.size(): " << gen2trk.size() << std::endl
	      << "  trk2seed.size(): " << trk2seed.size() << std::endl
	      << "  seed2gsf.size(): " << seed2gsf.size() << std::endl
	      << "  gsf2ele.size(): " << gsf2ele.size() << std::endl;
  }
  
  //////////
  // EGamma links
  // (More convoluted) 

  // Match GEN electrons to EGamma GsfTracks
  std::map<reco::CandidatePtr, reco::GsfTrackRef> gen2gsf_egamma;
  std::map<reco::CandidatePtr, float> gen2gsf_dr2_mode_egamma;
  std::map<reco::CandidatePtr, float> gen2gsf_dr2_egamma;
  std::vector<reco::GsfTrackRef> other_gsf_egamma;
  if ( isAOD_ == 1 ) matchGenToGsf( electrons_from_B, gsfTracksEGammaH_, 
				    gen2gsf_egamma, 
				    gen2gsf_dr2_mode_egamma, gen2gsf_dr2_egamma, 
				    other_gsf_egamma );

  // Match GEN electrons to EGamma (PF) GsfElectrons
  std::map<reco::CandidatePtr, reco::GsfElectronPtr> gen2ele_egamma;
  std::map<reco::CandidatePtr, float> gen2ele_dr2_egamma;
  std::vector<reco::GsfElectronPtr> other_ele_egamma;
  if ( isAOD_ == 1 ) matchGenToEle( electrons_from_B, gsfElectronsEGammaH_, 
				    gen2ele_egamma, 
				    gen2ele_dr2_egamma, 
				    other_ele_egamma );
  
  // Match EGamma ElectronSeeds to EGamma GsfTracks
  std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf_egamma;
  if ( isAOD_ == 1 ) matchSeedToGsf( eleSeedsEGammaH_, gsfTracksEGammaH_, seed2gsf_egamma );

  // Remove unmatched seeds (i.e. map entries with null GsfTrackRef)
  std::vector<reco::ElectronSeedRef> seeds;
  for ( const auto& iter : seed2gsf_egamma ) {
    if ( iter.second.isNonnull() ) { seeds.push_back(iter.first); }
  }

  // Match reco::Tracks to EGamma ElectronSeeds
  std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed_egamma;
  if ( isAOD_ == 1 ) matchTrkToSeed( ctfTracksH_, seeds, trk2seed_egamma );
  
  // Match EGamma GsfTracks to EGamma (PF) electrons
  std::map<reco::GsfTrackRef, reco::GsfElectronPtr> gsf2ele_egamma;
  if ( isAOD_ == 1 ) matchGsfToEle( gsfTracksEGammaH_, gsfElectronsEGammaH_, gsf2ele_egamma );

  if ( verbose_ > 0 ) {
    std::cout << "Collection sizes: " << std::endl
	      << "  ctfTracksH_->size(): " << ctfTracksH_->size() << std::endl
	      << "  eleSeedsEGammaH_->size(): " << eleSeedsEGammaH_->size() << std::endl
	      << "  seeds->size(): " << seeds.size() << std::endl
	      << "  gsfTracksEGammaH_->size(): " << gsfTracksEGammaH_->size() << std::endl
	      << "  gsfElectronsEGammaH_->size(): " << gsfElectronsEGammaH_->size() << std::endl
	      << "  electrons_from_B.size(): " << electrons_from_B.size() << std::endl
	      << "  gen2trk_egamma.size(): " << gen2trk.size() << std::endl
	      << "  trk2seed_egamma.size(): " << trk2seed_egamma.size() << std::endl
	      << "  seed2gsf_egamma.size(): " << seed2gsf_egamma.size() << std::endl
	      << "  gsf2ele_egamma.size(): " << gsf2ele_egamma.size() << std::endl;
  }

  //////////
  // Misc

//  // Match reco::Tracks to EGamma GsfTracks
//  std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf_egamma;
//  if ( isAOD_ == 1 ) matchTrkToGsf( ctfTracksH_, gsfTracksEGammaH_, trk2gsf_egamma );

//  // Match reco::Tracks to EGamma GsfTracks via EGamma ElectronSeeds
//  std::map<reco::TrackRef, reco::ElectronSeedRef> trk2seed_egamma;
//  std::map<reco::ElectronSeedRef, reco::GsfTrackRef> seed2gsf_egamma;
//  if ( isAOD_ == 1 ) matchTrkToGsf( ctfTracksH_, eleSeedsEGammaH_, gsfTracksEGammaH_, 
//				    trk2seed_egamma, seed2gsf_egamma );
  
//  // Match reco::Tracks to low pT GsfTracks
//  std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf;
//  if ( isAOD_ == 1 ) matchTrkToGsf( ctfTracksH_, gsfTracksH_, trk2gsf );

//  // Sanity checks
//  if ( isAOD_ == 1 ) { checkGsfToTrk( ctfTracksH_, gsfTracksH_, gsfTrackLinksH_, trk2gsf ); } 
//  //else if ( isAOD_ == 0 ) { checkGsfToTrk( packedCandLinksH_, lostTrackLinksH_ );

//  // Match GEN electrons to EGamma GsfTracks
//  std::map<reco::CandidatePtr, reco::GsfTrackRef> gen2gsf_egamma;
//  std::vector<reco::GsfTrackRef> other_gsf_egamma;
//  if ( isAOD_ == 1 ) 
//    matchGenToGsf( electrons_from_B, gsfTracksEGammaH_, gen2gsf_egamma, other_gsf_egamma );

  //////////

  //@@
  // Fill ntuple with signal and fake electrons from low pT reconstruction 
  //std::cout << "!!!LOW PT!!!" << std::endl;
  fill(event, setup,
       gen2trk, gen2gsf, gen2ele,
       gen2trk_dr2,
       gen2gsf_dr2_mode, gen2gsf_dr2,
       gen2ele_dr2,
       trk2seed, seed2gsf, gsf2ele );
  fill( false, event, setup,
	other_trk,
	trk2seed, seed2gsf, gsf2ele );
  
  // Fill ntuple with signal and fake electrons from EGamma reconstruction  
  //std::cout << "!!!EGAMMA!!!" << std::endl;
  fill( event, setup,
	gen2trk, gen2gsf_egamma, gen2ele_egamma,
	gen2trk_dr2,
	gen2gsf_dr2_mode_egamma, gen2gsf_dr2_egamma,
	gen2ele_dr2_egamma,
	trk2seed_egamma, seed2gsf_egamma, gsf2ele_egamma, true );
  fill( false, event, setup,
	other_trk,
	trk2seed_egamma, seed2gsf_egamma, gsf2ele_egamma, true );
  
}

////////////////////////////////////////////////////////////////////////////////
/*
  
 */ 
void IDFeatures::fill( const edm::Event& event,
		       const edm::EventSetup& setup,
		       std::map<reco::CandidatePtr, reco::TrackRef>& gen2trk,
		       std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
		       std::map<reco::CandidatePtr, reco::GsfElectronPtr>& gen2ele,
		       std::map<reco::CandidatePtr, float>& gen2trk_dr2,
		       std::map<reco::CandidatePtr, float>& gen2gsf_dr2_mode,
		       std::map<reco::CandidatePtr, float>& gen2gsf_dr2,
		       std::map<reco::CandidatePtr, float>& gen2ele_dr2,
		       const std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed,
		       const std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf,
		       const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
		       bool is_egamma ) {

  if ( verbose_ > 1 ) {
    std::cout << "#gen2trk: " << gen2trk.size() << " " << std::endl;
  }

//  // Add outer loop of GEN particles around nested method
//  for ( const auto& iter : gen2trk ) {
//    reco::CandidatePtr gen = iter.first;
//    if ( gen.isNull() ) { continue; }
//    ntuple_.reset(); // Init tree here
//    ntuple_.is_e(true);
//    ntuple_.fill_gen(gen); // Add gen info
//    ntuple_.fill_evt(event.id()); // Set Event
//    ntuple_.set_rho(*rhoH_); // Fill Rho
//    ntuple_.is_egamma(is_egamma); // Define if EGamma (or low pT)
//    reco::TrackRef trk = iter.second;
//    if ( trk.isNonnull() ) { 
//      ntuple_.has_trk(true); // Has reco::Track
////      std::cout << "Matched GenParticle #" << matched_trk->first.key()
////		<< " to reco::Track #" << matched_trk->second.key()
////		<< std::endl;
//      fill( true, // is_signal_ele
//	    event,
//	    setup,
//	    std::vector<reco::TrackRef>{trk}, // single TrackRef for this GEN particle
//	    trk2seed,
//	    seed2gsf,
//	    gsf2ele,
//	    is_egamma,
//	    false ); // Do not reset_and_fill in nested method (do it here!)
//    }
//    tree_->Fill(); // Fill tree here
//  }

  // Add outer loop of GEN particles around nested method
  for ( auto iter : gen2gsf ) {

    reco::CandidatePtr gen = iter.first;
    if ( gen.isNull() ) { continue; }

    ntuple_.reset(); // Init tree here
    ntuple_.is_e(true);
    ntuple_.fill_gen(gen); // Add gen info
    ntuple_.fill_evt(event.id()); // Set Event
    ntuple_.set_rho(*rhoH_); // Fill Rho
    ntuple_.is_egamma(is_egamma); // Define if EGamma (or low pT)

    // Record dr2, then set "has_gsf" if dr < dr_threshold_
    reco::GsfTrackRef gsf = iter.second;
    if ( gsf.isNonnull() ) {
      //@@ntuple_.gsf_dr_mode(gen2gsf_dr2_mode[gen]); 
      ntuple_.gsf_dr(gen2gsf_dr2_mode[gen]); //@@gen2gsf_dr2[gen]); 

      // Find matched electron
      if ( gen2ele.find(gen) != gen2ele.end() ) {
	reco::GsfElectronPtr ele = gen2ele[gen];
	if ( ele.isNonnull() ) { 
	  ntuple_.ele_dr(gen2ele_dr2[gen]); 
	}
      }

      if ( gen2gsf_dr2_mode[gen] < dr_threshold_*dr_threshold_ ) { 
	if ( verbose_ > 3 ) {
	  std::cout << "Matched GenParticle #" << gen.key()
		    << " to reco::GsfTrack #" << gsf.key()
		    << std::endl;
	}

	// Locate reco::Track that seeds this GsfTrack (if tracker-driven)
	std::vector<reco::TrackRef> one_trk;
	if ( gsf.isNonnull() ) { 
	  edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
	  if ( seed.isNonnull() ) { 
	    reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
	    if ( ele_seed.isNonnull() ) { 
	      reco::TrackRef trk = ele_seed->ctfTrack();
	      if ( trk.isNonnull() ) { 
		one_trk.push_back(trk);
		if ( gen2trk[gen] == trk ) { 
		  ntuple_.trk_dr(gen2trk_dr2[gen]); 
		  if ( verbose_ > 3 ) {
		    std::cout << "Matched GenParticle #" << gen.key()
			      << " to reco::GsfTrack #" << gsf.key()
			      << " (via reco::GsfTrack)"
			      << std::endl;
		  }
		} else {
		  std::cout << "WARNING! Track mismatch id/key: " 
			    << gen2trk[gen].id() << "/" << gen2trk[gen].key()
			    << " vs " << trk.id() << "/" << trk.key()
			    << std::endl;
		}
	      }
	    }
	  }
	}

	if ( !one_trk.empty() ) { 
	  fill( true, // is_signal_ele
		event,
		setup,
		one_trk, // single TrackRef for this GEN particle
		trk2seed,
		seed2gsf,
		gsf2ele,
		is_egamma,
		false ); // Do not reset_and_fill in nested method (do it here!)
	} else {
	  // ECAL-driven? Then set null Ref and handle in nested fill() method
	  fill( true, // is_signal_ele
		event,
		setup,
		std::vector<reco::GsfTrackRef>{gsf},
		gsf2ele,
		is_egamma,
		false ); // Do not reset_and_fill in nested method (do it here!)
	}

      } else { // dR(gen,gsf) < dr_threshold_ 

	// Attempt to match to electron
	if ( gen2ele.find(gen) != gen2ele.end() ) { 
	  reco::GsfElectronPtr ele = gen2ele[gen];

	  // Locate reco::Track that seeds this GsfTrack (if tracker-driven)
	  std::vector<reco::TrackRef> one_trk;
	  if ( ele.isNonnull() ) { 
	    reco::GsfTrackRef gsf = ele->gsfTrack();
	    if ( gsf.isNonnull() ) { 
	      edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
	      if ( seed.isNonnull() ) { 
		reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
		if ( ele_seed.isNonnull() ) { 
		  reco::TrackRef trk = ele_seed->ctfTrack();
		  if ( trk.isNonnull() ) { 
		    one_trk.push_back(trk);
		    if ( gen2trk[gen] == trk ) { 
		      ntuple_.trk_dr(gen2trk_dr2[gen]); 
		      if ( verbose_ > 3 ) {
			std::cout << "Matched GenParticle #" << gen.key()
				  << " to reco::TrackRef #" << trk.key()
				  << " via reco::GsfElectron #" << ele.key()
				  << std::endl;
		      }
		    } else {
		      std::cout << "WARNING! Track mismatch id/key: " 
				<< gen2trk[gen].id() << "/" << gen2trk[gen].key()
				<< " vs " << trk.id() << "/" << trk.key()
				<< std::endl;
		    }
		  }
		}
	      }
	    }
	  }

	  if ( !one_trk.empty() ) { 
	    fill( true, // is_signal_ele
		  event,
		  setup,
		  one_trk, // single TrackRef for this GEN particle
		  trk2seed,
		  seed2gsf,
		  gsf2ele,
		  is_egamma,
		  false ); // Do not reset_and_fill in nested method (do it here!)
	  }

	} else if ( gen2trk.find(gen) != gen2trk.end() ) { // fill at least reco::Track
	  reco::TrackRef trk = gen2trk[gen];
	  if ( trk.isNonnull() ) { 
	    ntuple_.trk_dr(gen2trk_dr2[gen]); 
	    if ( verbose_ > 3 ) {
	      std::cout << "Matched GenParticle #" << gen.key()
			<< " to reco::Track #" << gsf.key()
			<< " (in absence of reco::GsfTrack)"
			<< std::endl;
	    }
	    fill( true, // is_signal_ele
		  event,
		  setup,
		  std::vector<reco::TrackRef>{trk}, // single TrackRef for this GEN particle
		  trk2seed,
		  seed2gsf,
		  gsf2ele,
		  is_egamma,
		  false ); // Do not reset_and_fill in nested method (do it here!)
	  }
	}

      }
 
    } else {
      
      //std::cout << "ERROR! GsfTrackRef is null!!" << std::endl;
      
      // Find matched electron
      if ( gen2ele.find(gen) != gen2ele.end() ) {
	reco::GsfElectronPtr ele = gen2ele[gen];
	ntuple_.ele_dr(gen2ele_dr2[gen]); 
	std::vector<reco::TrackRef> one_trk;
	if ( ele.isNonnull() ) { 
	  ntuple_.ele_dr(gen2ele_dr2[gen]); 
	  reco::GsfTrackRef gsf = ele->gsfTrack();
	  if ( gsf.isNonnull() ) { 
	    edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
	    if ( seed.isNonnull() ) { 
	      reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
	      if ( ele_seed.isNonnull() ) { 
		reco::TrackRef trk = ele_seed->ctfTrack();
		if ( trk.isNonnull() ) { 
		  one_trk.push_back(trk);
		  if ( gen2trk[gen] == trk ) { 
		    ntuple_.trk_dr(gen2trk_dr2[gen]); 
		    if ( verbose_ > 3 ) {
		      std::cout << "Matched GenParticle #" << gen.key()
				<< " to reco::TrackRef #" << trk.key()
				<< " via reco::GsfElectron #" << ele.key()
				<< std::endl;
		    }
		  } else {
		    std::cout << "WARNING! Track mismatch id/key: " 
			      << gen2trk[gen].id() << "/" << gen2trk[gen].key()
			      << " vs " << trk.id() << "/" << trk.key()
			      << std::endl;
		  }
		}
	      }
	    }
	  }
	}
	if ( !one_trk.empty() ) { 
	  fill( true, // is_signal_ele
		event,
		setup,
		one_trk, // single TrackRef for this GEN particle
		trk2seed,
		seed2gsf,
		gsf2ele,
		is_egamma,
		false ); // Do not reset_and_fill in nested method (do it here!)
	}
      }

    }
    
    tree_->Fill(); // Fill tree here
    
  } // GEN ele loop
  
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::fill( bool is_signal_ele,
		       const edm::Event& event,
		       const edm::EventSetup& setup,
		       const std::vector<reco::TrackRef>& tracks,
		       const std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed,
		       const std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf,
		       const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
		       bool is_egamma,
		       bool reset_and_fill ) {

  if ( verbose_ > 1 ) {
    if ( reset_and_fill ) { std::cout << "#gen2trk: 0 "; }
    std::cout << "#tracks: " << tracks.size()
	      << " #trk2seed: " << trk2seed.size()
	      << " #seed2gsf: " << seed2gsf.size()
	      << " #gsf2ele: " << gsf2ele.size()
	      << std::endl;
  }

  // Outermost loop is over reco::Tracks
  for ( const auto& trk : tracks ) {
    if ( trk.isNull() || ( !is_signal_ele && gRandom->Rndm() > prescale_ ) ) { continue; }
    
    // Init tree
    if ( reset_and_fill ) { ntuple_.reset(); }
    
    ntuple_.is_e(is_signal_ele); // Set truth label
    ntuple_.is_other(!is_signal_ele); // Set truth label
    ntuple_.is_egamma(is_egamma); // Define if EGamma (or low pT)
    if ( !is_signal_ele ) { ntuple_.set_weight(prescale_); } // Set background weight

    ntuple_.fill_evt(event.id()); // Set Event ID
    ntuple_.set_rho(*rhoH_); // Fill Rho

    // Store reco::Track info
    ntuple_.has_trk(true); // Has reco::Track
    ntuple_.fill_trk( edm::refToPtr(trk), *beamspotH_ );
    const auto& matched_seed = trk2seed.find(trk);
    if ( matched_seed != trk2seed.end() && matched_seed->second.isNonnull() ) { 
      ntuple_.has_seed(true);
      if ( verbose_ > 3 ) {
	std::cout << "Matched reco::Track #" << trk.key()
		  << " to ElectronSeed #" << matched_seed->second.key()
		  << std::endl;
      }

      // Store PreId info and seed BDT discrimator values
      if ( !is_egamma && isAOD_ == 1 
	   && preIdRefsH_.isValid() && preIdsEcalH_.isValid() && preIdsHcalH_.isValid() ) { 
	const reco::PreIdRef preid_ecal = (*preIdRefsH_)[matched_seed->second->ctfTrack()];
	const reco::PreIdRef preid_hcal(preIdsHcalH_,preid_ecal.index());
	if ( preid_ecal.isNonnull() && preid_hcal.isNonnull() ) {
	  noZS::EcalClusterLazyTools ecal_tools(event, setup, ebRecHits_, eeRecHits_);
	  ntuple_.fill_preid( *preid_ecal,
			      *preid_hcal,
			      *beamspotH_,
			      *rhoH_,
	                      ecal_tools );
	}
      }
      
      // Check if ElectronSeed is matched to a reco::GsfTrack
      const auto& matched_gsf = seed2gsf.find(matched_seed->second);
      if ( matched_gsf != seed2gsf.end() && matched_gsf->second.isNonnull() ) { 

	if ( verbose_ > 3 ) {
	  std::cout << "Matched ElectronSeed #" << matched_gsf->first.key()
		    << " to GsfTrack #" << matched_gsf->second.key()
		    << std::endl;
	}

	fill( is_signal_ele,
	      event,
	      setup,
	      std::vector<reco::GsfTrackRef>{matched_gsf->second},
	      gsf2ele,
	      is_egamma,
	      reset_and_fill );

//	ntuple_.has_gsf(true);
////	std::cout << "Matched ElectronSeed #" << matched_gsf->first.key()
////		  << " to GsfTrack #" << matched_gsf->second.key()
////		  << std::endl;
//	
//	ntuple_.fill_gsf(matched_gsf->second, *beamspotH_);
//	
//	// Store Seed BDT discrimator values (at GsfTrack-level)
//	if ( !is_egamma ) { 
//	  float unbiased = (*mvaUnbiasedH_)[matched_gsf->second];
//	  float ptbiased = (*mvaPtbiasedH_)[matched_gsf->second];
//	  ntuple_.fill_seed( unbiased, ptbiased );
//	}
//	
//	// Check if GsfTrack is matched to an electron 
//	const auto& matched_ele = gsf2ele.find(matched_gsf->second);
//	if ( matched_ele != gsf2ele.end() && matched_ele->second.isNonnull() ) {
//	  ntuple_.has_ele(true);
////	  std::cout << "Matched GsfTrack #" << matched_ele->first.key()
////		    << " to GsfElectron #" << matched_ele->second.key()
////		    << std::endl;
//	  
//	  //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	  float mva_value = -999.;
//	  int mva_id = -999;
//	  if ( !is_egamma ) {
//	    if ( mvaValueLowPtH_.isValid() && mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
//	      mva_value = mvaValueLowPtH_->get( matched_ele->second.key() );
//	    } else {
//	      std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
//	    }
//	  } else {
//	    if ( mvaValueEGammaH_.isValid() && mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
//	      mva_value = mvaValueEGammaH_->get( matched_ele->second.key() );
//	    } else {
//	      std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
//	    }
//	    if ( mvaIdEGammaH_.isValid() && mvaIdEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
//	      mva_id = mvaIdEGammaH_->get( matched_ele->second.key() ) ? 1 : 0;
//	    } else {
//	      std::cout << "ERROR! Issue matching MVA ID to GsfElectrons!" << std::endl;
//	    }
//	  }
//	  
//	  //@@ dirty hack as is not in Event nor embedded in pat::Electron
//	  float conv_vtx_fit_prob = -999.;
//	  //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
//	  //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
//	  //}
//	  
//	  ntuple_.fill_ele( matched_ele->second, mva_value, mva_id, conv_vtx_fit_prob, *rhoH_ );
//	  
//	  //@@ Add SuperCluster vars?
//	  //ntuple_.fill_supercluster(matched_ele->second);
//	  
//	} // Find GsfElectron

      } // Find GsfTrack
    } // Find ElectronSeed
  
    // Fill tree
    if ( reset_and_fill ) { tree_->Fill(); }

  } // Loop over tracks

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::fill( bool is_signal_ele,
		       const edm::Event& event,
		       const edm::EventSetup& setup,
		       const std::vector<reco::GsfTrackRef>& tracks,
		       const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
		       bool is_egamma,
		       bool reset_and_fill ) {
  
  if ( verbose_ > 1 ) {
    if ( reset_and_fill ) { 
      std::cout << "#tracks: 0"
		<< " #gsf2trk: 0"
		<< " #trk2seed: 0"
		<< " #seed2gsf: 0";
    }
    std::cout << "#gsftracks: " << tracks.size()
	      << " #gsf2ele: " << gsf2ele.size()
	      << std::endl;
  }

  // Outermost loop is over reco::GsfTracks
  for ( const auto& gsf : tracks ) {

    if ( gsf.isNull() /*@@|| ( !is_signal_ele && gRandom->Rndm() > prescale_ )*/ ) { continue; }
    
    // Init tree
    if ( reset_and_fill ) { ntuple_.reset(); }
    
    ntuple_.is_e(is_signal_ele); // Set truth label
    ntuple_.is_other(!is_signal_ele); // Set truth label
    ntuple_.is_egamma(is_egamma); // Define if EGamma (or low pT)
    if ( !is_signal_ele ) { ntuple_.set_weight(prescale_); } // Set background weight

    ntuple_.fill_evt(event.id()); // Set Event ID
    ntuple_.set_rho(*rhoH_); // Fill Rho

    // Store reco::GsfTrack info
    ntuple_.has_gsf(true);
    ntuple_.fill_gsf(edm::refToPtr(gsf), *beamspotH_);
	
    // Store Seed BDT discrimator values (at GsfTrack-level)
    if ( !is_egamma ) { 
      float unbiased = (*mvaUnbiasedH_)[gsf];
      float ptbiased = (*mvaPtbiasedH_)[gsf];
      ntuple_.fill_seed( unbiased, ptbiased );
    }
    
    // Check if GsfTrack is matched to an electron 
    const auto& matched_ele = gsf2ele.find(gsf);
    if ( matched_ele != gsf2ele.end() && matched_ele->second.isNonnull() ) {
      ntuple_.has_ele(true);
      if ( verbose_ > 3 ) {
	std::cout << "Matched GsfTrack #" << matched_ele->first.key()
		  << " to GsfElectron #" << matched_ele->second.key()
		  << std::endl;
      }

      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
      float mva_value = -999.;
      int mva_id = -999;
      if ( !is_egamma ) {
	if ( mvaValueLowPtH_.isValid() && mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
	  mva_value = mvaValueLowPtH_->get( matched_ele->second.key() );
	} else {
	  std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
      } else {
	if ( mvaValueEGammaH_.isValid() && mvaValueEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_value = mvaValueEGammaH_->get( matched_ele->second.key() );
	} else {
	  std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
	}
	if ( mvaIdEGammaH_.isValid() && mvaIdEGammaH_->size() == gsfElectronsEGammaH_->size() ) {
	  mva_id = mvaIdEGammaH_->get( matched_ele->second.key() ) ? 1 : 0;
	} else {
	  std::cout << "ERROR! Issue matching MVA ID to GsfElectrons!" << std::endl;
	}
      }
      
      //@@ dirty hack as is not in Event nor embedded in pat::Electron
      float conv_vtx_fit_prob = -999.;
      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
      //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
      //}
      
      ntuple_.fill_ele( matched_ele->second, mva_value, mva_id, conv_vtx_fit_prob, *rhoH_ );
      
      //@@ Add SuperCluster vars?
      //ntuple_.fill_supercluster(matched_ele->second);
      
    } // Find GsfElectron
    
    // Fill tree
    if ( reset_and_fill ) { tree_->Fill(); }
    
  } // Loop over tracks
  
}

////////////////////////////////////////////////////////////////////////////////
// Assumes GenParticles from RECO/AOD
//void IDFeatures::electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& genParticles,
//				 std::set<reco::CandidatePtr>& electrons_from_B ) {
void IDFeatures::electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& genParticles,
				 std::set<reco::CandidatePtr>& electrons_from_B ) {
  
  bool tag_side_muon = false;

  electrons_from_B.clear();
  for ( size_t idx = 0; idx < genParticles->size(); idx++ ) {
    //reco::CandidatePtr genp(genParticles, idx);
    reco::GenParticlePtr genp(genParticles, idx);
    
    // Last copy of GEN electron 
    bool is_ele = std::abs(genp->pdgId()) == 11 && genp->isLastCopy(); //@@ not a method of Candidate

    // Does GEN ele comes from B decay?
    bool non_resonant = 
      genp->numberOfMothers() >= 1 && genp->mother() &&                     // has mother
      std::abs(genp->mother()->pdgId()) > 510 &&                            // mother is B
      std::abs(genp->mother()->pdgId()) < 546;                              // mother is B
    bool resonant = 
      genp->numberOfMothers() >= 1 && genp->mother() &&                     // has mother
      std::abs(genp->mother()->pdgId()) == 443 &&                           // mother is J/psi
      genp->mother()->numberOfMothers() >= 1 && genp->mother()->mother() && // has grandmother
      std::abs(genp->mother()->mother()->pdgId()) > 510 &&                  // grandmother is B
      std::abs(genp->mother()->mother()->pdgId()) < 546;                    // grandmother is B
    
    //  Check for tag side muon
    float muon_pt = 7.;
    float muon_eta = 1.5;
    tag_side_muon |= ( std::abs(genp->pdgId()) == 13 && genp->isLastCopy() 
		       && 
		       ( ( genp->numberOfMothers() >= 1 &&
			   genp->mother() &&
			   std::abs(genp->mother()->pdgId()) > 510 &&
			   std::abs(genp->mother()->pdgId()) < 546 && 
			   genp->mother()->pt() > muon_pt && 
			   std::abs(genp->mother()->eta()) < muon_eta ) 
			 ||
			 ( genp->numberOfMothers() >= 1 && 
			   genp->mother() &&
			   genp->mother()->numberOfMothers() >= 1 && 
			   genp->mother()->mother() &&
			   std::abs(genp->mother()->mother()->pdgId()) > 510 &&
			   std::abs(genp->mother()->mother()->pdgId()) < 546 && 
			   genp->mother()->mother()->pt() > muon_pt &&
			   std::abs(genp->mother()->mother()->eta()) < muon_eta ) ) );
    
    // is coming from a B
    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
      electrons_from_B.insert(genp);
//      std::cout << "electronsFromB: "
//		<< " #gen_electrons: " << electrons_from_B.size()
//		<< " resonant? " << resonant
//		<< " non resonant? " << non_resonant
//		<< " tag-side muon? " << tag_side_muon
//		<< std::endl;
    }
    
  } // genParticles loop

  if ( !tag_side_muon ) { electrons_from_B.clear(); }

}

//////////////////////////////////////////////////////////////////////////////////
//// Assumes pruned and packed GenParticles from MINIAOD
//void IDFeatures::electronsFromB( const edm::Handle< edm::View<reco::GenParticle> >& prunedGenParticles,
//				 const edm::Handle< edm::View<reco::Candidate> >& packedGenParticles,
//				 std::set<reco::CandidatePtr>& electrons_from_B ) {
//
//  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Examples
//  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/packedGenParticles_cfi.py
//  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py
//
//  electrons_from_B.clear();
//  
//  // Iterate through ("status 1") packed candidates
//  for ( size_t ipa = 0; ipa < packedGenParticles->size(); ipa++ ) {
//    reco::CandidatePtr packed(packedGenParticles, ipa);
//    
//    // Find GEN electrons
//    bool is_ele = std::abs(packed->pdgId()) == 11; // && packed->isLastCopy() 
//
//    // Search for ancestors only for GEN electrons ... 
//    if ( !is_ele ) { continue; }
//    
//    // Does GEN ele come from B decay?
//    bool comes_from_B = false;
//    for ( size_t ipr = 0; ipr < prunedGenParticles->size(); ipr++ ) {
//      reco::GenParticlePtr pruned(prunedGenParticles, ipr);
//      comes_from_B = 
//	std::abs(pruned->pdgId()) > 510 &&
//	std::abs(pruned->pdgId()) < 546 &&
//	packed->mother(0) != nullptr && 
//	isAncestor( &*pruned, &*packed );
//
//      // Is coming from a B
//      if ( is_ele && ( comes_from_B || !check_from_B_ ) ) {
//	//electrons_from_B.insert(packed);
//      }
//      
//    }
//    
//    //std::cout << "is_ele: " << is_ele << " comes_from_B: " << comes_from_B << std::endl;
//    
//  } // packedGenParticles loop
//  
//}

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
				std::map<reco::CandidatePtr, float>& gen2trk_dr2,
				std::vector<reco::TrackRef>& other_trk ) {

  // Identify all possible matches b/w GEN ele and tracks, store indices and dr2
  std::vector<Match> matches;
  matches.reserve( electrons_from_B.size()*ctfTracks->size() );
  for ( auto gen : electrons_from_B ) {
    for ( size_t itrk = 0; itrk < ctfTracks->size(); ++itrk ) {
      reco::TrackRef trk(ctfTracks,itrk);
      if ( !filterTrack(trk) ) { continue; }
      matches.emplace_back(gen,itrk,reco::deltaR2(trk->eta(),
						  trk->phi(),
						  gen->eta(), 
						  gen->phi()) );
    }
  }

  // Sort all matches by dr2
  std::sort(matches.begin(), matches.end(), Match::compare_by_dr2);

  // Store best matches by dr2
  std::set<reco::CandidatePtr> matched_gen;
  std::set<reco::TrackRef> matched_trk;
  for ( auto iter : matches ) {
    reco::CandidatePtr gen = iter.gen_;
    reco::TrackRef trk(ctfTracks,iter.idx_);
    float dr2 = iter.dr2_;
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) == matched_gen.end() &&
	 std::find( matched_trk.begin(), matched_trk.end(), trk ) == matched_trk.end() && 
	 dr2 < 0.01 ) { //@@ max dR is 0.1
      gen2trk.insert( std::pair<reco::CandidatePtr, reco::TrackRef>( gen, trk ) );
      gen2trk_dr2.insert( std::pair<reco::CandidatePtr, float>( gen, dr2 ) );
      matched_gen.insert(gen);
      matched_trk.insert(trk);
      if ( matched_gen.size() >= electrons_from_B.size() &&
	   matched_trk.size() >= electrons_from_B.size() ) {
	break;
      }
    }
  }

  if ( verbose_ > 0 ) {
    if ( gen2trk.size() <  electrons_from_B.size() ) {
      std::cout << "WARNING: Fewer matches than GEN electrons!" << std::endl;
    }
  }

//  std::cout << "gen2trk size: " << gen2trk.size() << std::endl;
//  size_t idx = 0; 
//  for ( auto iter : gen2trk ) {
//    reco::CandidatePtr gen = iter.first;
//    reco::TrackRef trk = iter.second;
//    std::cout << "  idx: " << idx
//      	      << " gen key: " << gen.key()
//	      << " trk key: " << trk.key()
//	      << " dr2: " << gen2trk_dr2[iter.first]
//	      << std::endl
//	      << "  gen pt: " << gen->pt()
//	      << " gen eta: " << gen->eta()
//	      << " gen phi: " << gen->phi()
//	      << std::endl
//	      << "  trk pt: " << trk->pt()
//	      << " trk eta: " << trk->eta()
//	      << " trk phi: " << trk->phi()
//	      << std::endl;
//    idx++;
//  }

//  gen2trk.clear();
//  gen2trk_dr2.clear();
//  other_trk.clear();
//  //std::set<reco::CandidatePtr> matched_gen;
//  //std::set<reco::TrackRef> matched_trk;
//  matched_gen.clear();
//  matched_trk.clear();
//
//  // Match ctfTracks to GEN ele from B decays
//  for ( const auto& gen : electrons_from_B ) {
//    size_t best_idx = 666;
//    //double min_dr2 = dr_max_*dr_max_;
//    double min_dr2 = 0.01;
//    for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
//      reco::TrackRef trk(ctfTracks,idx);
//      if ( !filterTrack(trk) ) { continue; }
//      double dr2 = reco::deltaR2(trk->eta(),
//				 trk->phi(),
//				 gen->eta(), 
//				 gen->phi() );
//      if ( dr2 < min_dr2 ) {
//	min_dr2 = dr2;
//	best_idx = idx;
//      }
//    }
//    //if ( min_dr2 < dr_max_*dr_max_ ) {
//    if ( best_idx < 666 && min_dr2 < 0.01 ) {
//      reco::TrackRef trk(ctfTracks,best_idx);
//      gen2trk.insert( std::pair<reco::CandidatePtr, reco::TrackRef>( gen, trk ) );
//      gen2trk_dr2.insert( std::pair<reco::CandidatePtr, float>( gen, min_dr2 ) );
//      matched_gen.insert(gen);
//      matched_trk.insert(trk);
//    }
//  }
//
//  std::cout << "gen2trk size: " << gen2trk.size() << std::endl;
//  idx = 0; 
//  for ( auto iter : gen2trk ) {
//    reco::CandidatePtr gen = iter.first;
//    reco::TrackRef trk = iter.second;
//    std::cout << "  idx: " << idx
//      	      << " gen key: " << gen.key()
//	      << " trk key: " << trk.key()
//	      << " dr2: " << gen2trk_dr2[iter.first]
//	      << std::endl
//	      << "  gen pt: " << gen->pt()
//	      << " gen eta: " << gen->eta()
//	      << " gen phi: " << gen->phi()
//	      << std::endl
//	      << "  trk pt: " << trk->pt()
//	      << " trk eta: " << trk->eta()
//	      << " trk phi: " << trk->phi()
//	      << std::endl;
//    idx++;
//  }

  // Add null reco::Track Ref for unmatched GEN particles
  for ( const auto& gen : electrons_from_B ) {
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) != matched_gen.end() ) { continue; }
    gen2trk[gen] = reco::TrackRef();
    gen2trk_dr2[gen] = -1.;
  }

  // Identify ctfTracks from background, i.e. not matched to a GEN electron
  other_trk.reserve( ctfTracks->size() - matched_trk.size() );
  for ( size_t idx = 0; idx < ctfTracks->size(); idx++ ) {
    reco::TrackRef trk(ctfTracks,idx);
    if ( !filterTrack(trk) ) { continue; }
    if ( matched_trk.find(trk) == matched_trk.end() ) {
      other_trk.push_back(trk);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				 const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
				 std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed ) {
  std::vector<reco::ElectronSeedRef> seeds;
  for ( size_t idx = 0; idx < eleSeeds->size(); ++idx ) {
    reco::ElectronSeedRef ele_seed(eleSeeds,idx);
    seeds.push_back(ele_seed); 
  }
  matchTrkToSeed( ctfTracks, seeds, trk2seed );
}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchTrkToSeed( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				 const std::vector<reco::ElectronSeedRef>& eleSeeds,
				 std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed ) {
  
  //std::cout << "IDFeatures::matchTrkToSeed()" << std::endl;

  trk2seed.clear();
  std::set<reco::TrackRef> matched;
  typedef std::map<reco::TrackRef, reco::ElectronSeedRef>::iterator Iter;

  // Store ctfTrack -> eleSeed link
  for ( size_t idx = 0; idx < eleSeeds.size(); ++idx ) {
    reco::ElectronSeedRef ele_seed(eleSeeds[idx]);
    if ( ele_seed.isNull() ) { 
      std::cout << "ERROR! Null ElectronSeedRef!" << std::endl;
      continue;
    }
    if ( ele_seed->isEcalDriven() ) { 
      //@@ if not used, then see "ERROR! Single TrackRef mapped to many ElectronSeedRefs!"
      if ( verbose_ > 2 ) { std::cout << "INFO: ECAL-driven seed!" << std::endl; }
      //continue;
    }
    if ( verbose_ > 1 ) { 
      std::cout << "INFO: ECAL-driven seed: " << ele_seed->isEcalDriven() << std::endl
		<< "       TRK-driven seed: " << ele_seed->isTrackerDriven() << std::endl;
    }
    reco::CaloClusterRef calo = ele_seed->caloCluster();
    if ( calo.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null CaloClusterRef!" << std::endl; }
    }
    reco::TrackRef trk = ele_seed->ctfTrack();
    if ( trk.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null TrackRef!" << std::endl; }
      continue;
    } 
    if ( trk.id() != ctfTracks.id() ) { 
      std::cout << "ERROR! Collections do not match!" << std::endl;
      continue;
    } 
    Iter iter = trk2seed.find(trk);
    if ( iter == trk2seed.end() ) { 
      trk2seed[trk] = ele_seed;
      matched.insert(trk);
      if ( verbose_ > 2 ) {
	std::cout << " Match! " 
		  << "  Trk key: " << trk.key()
		  << "  Seed key: " << ele_seed.key()
		  << "  Calo key: " << ( calo.isNull() ? -1 : calo.key() )
		  << std::endl 
		  << "  trk pT:  " << trk->pt()
		  << "  eta " << trk->eta()
		  << "  phi " << trk->phi() 
		  << std::endl 
		  << "  calo Et: " << ( calo.isNull() ? -1. : calo->energy() / std::cosh(calo->eta()) )
		  << "  eta " << ( calo.isNull() ? -1. : calo->eta() )
		  << "  phi " << ( calo.isNull() ? -1. : calo->phi() )
		  << std::endl 
		  << "  dR(trk,calo):  "
		  << ( calo.isNull() ? -1. : sqrt( reco::deltaR2( trk->eta(),
								  trk->phi(),
								  calo->eta(),
								  calo->phi() ) ) )
		  << std::endl;
      }
      
    } else { // if ( iter == trk2seed.end() )
      
      std::cout << "ERROR! Single TrackRef mapped to many ElectronSeedRefs!" << std::endl;
      if ( iter->first.isNonnull() && iter->second.isNonnull() ) {
	std::cout << " old: seed key: " << std::setw(3) << iter->second.key() 
		  << " trk key: " << std::setw(3) << iter->first.key() 
		  << " pT: " << iter->first->pt() 
		  << " eta: " << iter->first->eta() 
		  << " phi: " << iter->first->phi() 
		  << std::endl;
	if ( iter->second->caloCluster().isNonnull() ) {
	  std::cout << " calo:                           Et: " 
		    << iter->second->caloCluster()->energy() 
	    / std::cosh(iter->second->caloCluster()->eta())
		    << " eta: " << iter->second->caloCluster()->eta() 
		    << " phi: " << iter->second->caloCluster()->phi() 
		    << ", dR: " << sqrt( reco::deltaR2( iter->first->eta(),
							iter->first->phi(),
							iter->second->caloCluster()->eta(), 
							iter->second->caloCluster()->phi() ) )
		    << std::endl;
	}
      }
      
      if ( trk.isNonnull() && ele_seed.isNonnull() ) {
	std::cout << " new: seed key: " << std::setw(3) << ele_seed.key() 
		  << " trk key: " << std::setw(3) << trk.key() 
		  << " pT: " << trk->pt() 
		  << " eta: " << trk->eta() 
		  << " phi: " << trk->phi() 
		  << std::endl;
	if ( calo.isNonnull() ) {
	  std::cout << " calo:                           Et: " 
			<< calo->energy() 
	    / std::cosh(calo->eta())
		    << " eta: " << calo->eta() 
		    << " phi: " << calo->phi() 
		    << ", dR: " << sqrt( reco::deltaR2( trk->eta(),
							trk->phi(),
							calo->eta(), 
							calo->phi() ) )
		    << std::endl;
	}
      }
      
      if ( trk.isNonnull() && iter->first.isNonnull() ) {
	if ( reco::deltaR2( trk->eta(), 
			    trk->phi(), 
			    iter->second->ctfTrack()->eta(), 
			    iter->second->ctfTrack()->phi() ) < 
	     reco::deltaR2( iter->first->eta(),
			    iter->first->phi(),
			    iter->second->ctfTrack()->eta(), 
			    iter->second->ctfTrack()->phi() ) ) {
	  std::cout << "INFO: Overwriting ctfTrackRef->electronSeedRef entry!" << std::endl;
	  trk2seed[trk] = ele_seed; 
	  matched.insert(trk);
	}
      }
      
    } // if ( iter == trk2seed.end() ) 
  } // eleSeeds iterator

  // Add null eleSeed Ref for unmatched ctfTracks
  for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
    reco::TrackRef trk(ctfTracks, idx);
    if ( std::find( matched.begin(), matched.end(), trk ) != matched.end() ) { continue; }
    trk2seed[trk] = reco::ElectronSeedRef();
  }

//  int cntr = 0;
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
				 std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf ) {
  
  seed2gsf.clear();
  std::set<reco::ElectronSeedRef> matched;
  typedef std::map<reco::ElectronSeedRef, reco::GsfTrackRef>::const_iterator Iter;

  // Store eleSeed -> gsfTrack link
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackRef gsf(gsfTracks, idx);
    if ( gsf.isNull() ) { 
      std::cout << "ERROR! Null GsfTrackRef!" << std::endl;
      continue;
    }
    edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
    if ( seed.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null TrajectorySeedRef!" << std::endl; }
      continue;
    }
    reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
    if ( ele_seed.isNull() ) { 
      std::cout << "ERROR! Null ElectronSeedRef!" << std::endl;
      continue;
    }
    if ( ele_seed.id() != eleSeeds.id() ) { 
      std::cout << "ERROR! Collections do not match!" << std::endl;
      continue;
    }
    reco::TrackRef trk = ele_seed->ctfTrack();
    if ( trk.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null TrackRef!" << std::endl; }
    }
    reco::CaloClusterRef calo = ele_seed->caloCluster();
    if ( calo.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null CaloClusterRef!" << std::endl; }
    }
    Iter iter = seed2gsf.find(ele_seed);
    if ( iter == seed2gsf.end() ) { 
      seed2gsf[ele_seed] = gsf;
      matched.insert(ele_seed);
      if ( verbose_ > 2 ) {
	std::cout << " Match! " 
		  << "  Seed key: " << ele_seed.key()
		  << "  Trk key: " << ( trk.isNull() ? -1 : trk.key() )
		  << "  Calo key: " << ( calo.isNull() ? -1 : calo.key() )
		  << "  Gsf key: " << gsf.key() 
		  << std::endl 
		  << "  trk pT:  " << ( trk.isNull() ? -1. : trk->pt() )
		  << "  eta " << ( trk.isNull() ? -1. : trk->eta() )
		  << "  phi " << ( trk.isNull() ? -1. : trk->phi() )
		  << std::endl 
		  << "  calo Et: " << ( calo.isNull() ? -1. : calo->energy() 
					/ std::cosh(calo->eta()) )
		  << "  eta " << ( calo.isNull() ? -1. : calo->eta() )
		  << "  phi " << ( calo.isNull() ? -1. : calo->phi() )
		  << std::endl 
		  << "  gsf pT:  " << gsf->ptMode()
		  << "  eta " << gsf->etaMode()
		  << "  phi " << gsf->phiMode() 
		  << std::endl 
		  << "  dR(trk,gsf):  "
		  << ( trk.isNull() ? -1. : sqrt( reco::deltaR2( trk->eta(),
								 trk->phi(),
								 gsf->etaMode(),
								 gsf->phiMode() ) ) )
		  << std::endl 
		  << "  dR(calo,gsf): "
		  << ( calo.isNull() ? -1. : sqrt( reco::deltaR2( calo->eta(),
								  calo->phi(),
								  gsf->etaMode(),
								  gsf->phiMode() ) ) )
		  << std::endl;
      }
    } else {
      std::cout << "ERROR! Single ElectronSeedRef mapped to many GsfTrackRefs!" << std::endl;
      if ( trk.isNonnull() && iter->first->ctfTrack().isNonnull() ) {
	std::cout << " old seed key, pT, eta, phi: " 
		  << std::setw(3) 
		  << iter->first.key() << " "
		  << iter->first->ctfTrack()->pt() << " "
		  << iter->first->ctfTrack()->eta() << " "
		  << iter->first->ctfTrack()->phi() 
		  << std::endl 
		  << " old gsf  key, pT, eta, phi: " 
		  << std::setw(3) 
		  << iter->second.key() << " "
		  << iter->second->pt() << " "       //@@ Mode?
		  << iter->second->eta() << " "      //@@ Mode?
		  << iter->second->phi() << ", dR: " //@@ Mode?
		  << sqrt( reco::deltaR2( iter->first->ctfTrack()->eta(),
					  iter->first->ctfTrack()->phi(),
					  iter->second->eta(),    //@@ Mode?
					  iter->second->phi() ) ) //@@ Mode?
		  << std::endl 
		  << " new seed key, pT, eta, phi: " 
		  << std::setw(3) 
		  << ele_seed.key() << " "
		  << trk->pt() << " "
		  << trk->eta() << " "
		  << trk->phi()
		  << std::endl 
		  << " new gsf  key, pT, eta, phi: " 
		  << std::setw(3) 
		  << gsf.key() << " "
		  << gsf->pt() << " "          //@@ Mode?
		  << gsf->eta() << " "         //@@ Mode?
		  << gsf->phi() << ", dR: "    //@@ Mode?
		  << sqrt( reco::deltaR2( trk->eta(), 
					  trk->phi(), 
					  gsf->eta(),      //@@ Mode?
					  gsf->phi() ) )   //@@ Mode?
		  << std::endl; 
	if ( reco::deltaR2( trk->eta(), 
			    trk->phi(), 
			    gsf->etaMode(),
			    gsf->phiMode() ) <
	     reco::deltaR2( iter->first->ctfTrack()->eta(),
			    iter->first->ctfTrack()->phi(),
			    iter->second->etaMode(),
			    iter->second->phiMode() ) ) {
	  std::cout << "INFO: Overwriting electronSeedRef->gsfTrackRef entry!" << std::endl;
	  seed2gsf[ele_seed] = gsf; 
	  matched.insert(ele_seed);
	}
      }
    }
  }

  // Add null gsfTrack Ref for unmatched eleSeeds
  for ( size_t idx = 0; idx < eleSeeds->size(); ++idx ) {
    reco::ElectronSeedRef seed(eleSeeds, idx);
    if ( std::find( matched.begin(), matched.end(), seed ) != matched.end() ) { continue; }
    seed2gsf[seed] = reco::GsfTrackRef();
  }

//  int cntr = 0;
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
				const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
				std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele ) {
  
  gsf2ele.clear();
  std::set<reco::GsfTrackRef> matched;
  typedef std::map<reco::GsfTrackRef, reco::GsfElectronPtr>::const_iterator Iter;

  // Store gsfTrack -> gsfElectron link
  for ( size_t idx = 0; idx < gsfElectrons->size(); ++idx ) {
    reco::GsfElectronPtr ele(gsfElectrons, idx);
    if ( ele.isNull() ) { 
      std::cout << "ERROR! Null GsfElectronPtr!" << std::endl;
      continue;
    }
    reco::GsfTrackRef gsf = ele->gsfTrack();
    if ( gsf.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO! Null GsfTrackRef!" << std::endl; }
      continue;
    }
    if ( gsf.id() != gsfTracks.id() ) { 
      std::cout << "ERROR! Collections do not match!!" << std::endl;
      continue;
    }
    Iter iter = gsf2ele.find(gsf);
    if ( iter == gsf2ele.end() ) { 
      gsf2ele[gsf] = ele;
      matched.insert(gsf);
    } else {
      std::cout << "ERROR! Single GsfTrackRef mapped to many GsfElectronPtrs!" << std::endl;
    }
  }

  // Add null gsfTrack Ref for unmatched eleSeeds
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackRef gsf(gsfTracks, idx);
    if ( std::find( matched.begin(), matched.end(), gsf ) != matched.end() ) { continue; }
    gsf2ele[gsf] = reco::GsfElectronPtr();
  }

//  int cntr = 0;
//  for ( Iter iter = gsf2ele.begin(); iter != gsf2ele.end(); ++iter ) {
//    if ( iter->second.isNonnull() ) { cntr++; }
//  }
//  std::cout << " # gsfTracks: " << gsfTracks->size()
//	    << " # gsfElectrons: " << gsfElectrons->size()
//	    << " # map entries: " << gsf2ele.size()
//	    << " # nonnull ptrs: " << cntr
//	    << " # matched: " << matched.size()
//	    << std::endl;
  
}

//////////////////////////////////////////////////////////////////////////////////
//// 
//void IDFeatures::matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
//				const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
//				std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele ) {
//  gsf2ele.clear();
//  for ( size_t idx = 0; idx < gsfElectrons->size(); ++idx ) {
//    reco::GsfElectronPtr ele(gsfElectrons, idx);
//    reco::GsfTrackRef gsf = ele->gsfTrack();
//    if ( gsf2ele.find(gsf) != gsf2ele.end() ) {
//      std::cout << "THIS SHOULD NEVER HAPPEN! Multiple low pT electrons matched to the same GSFTrack?!"
//		<< std::endl;
//    } else {
//      gsf2ele.insert( std::pair<reco::GsfTrackRef, reco::GsfElectronPtr>(gsf, ele) );
//    }
//  }
//}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchGenToGsf( const std::set<reco::CandidatePtr>& electrons_from_B,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				std::map<reco::CandidatePtr, reco::GsfTrackRef>& gen2gsf,
				std::map<reco::CandidatePtr, float>& gen2gsf_dr2_mode,
				std::map<reco::CandidatePtr, float>& gen2gsf_dr2,
				std::vector<reco::GsfTrackRef>& other_gsf ) {

  gen2gsf.clear();
  gen2gsf_dr2_mode.clear();
  gen2gsf_dr2.clear();
  
  // Identify all possible matches b/w GEN ele and tracks, store indices and dr2
  std::vector<Match> matches_mode;
  std::vector<Match> matches;
  matches_mode.reserve( electrons_from_B.size()*gsfTracks->size() );
  matches.reserve( electrons_from_B.size()*gsfTracks->size() );
  for ( auto gen : electrons_from_B ) {
    for ( size_t itrk = 0; itrk < gsfTracks->size(); ++itrk ) {
      reco::GsfTrackRef trk(gsfTracks,itrk);
      matches_mode.emplace_back(gen,itrk,reco::deltaR2(trk->etaMode(),
						       trk->phiMode(),
						       gen->eta(), 
						       gen->phi()) );
      matches.emplace_back(gen,itrk,reco::deltaR2(trk->eta(),
						  trk->phi(),
						  gen->eta(), 
						  gen->phi()) );
    }
  }
  
  // Sort all matches by dr2
  std::sort(matches_mode.begin(), matches_mode.end(), Match::compare_by_dr2);
  std::sort(matches.begin(), matches.end(), Match::compare_by_dr2);

  // Store best matches by dr2
  std::set<reco::CandidatePtr> matched_gen;
  std::set<reco::GsfTrackRef> matched_trk;
  for ( auto iter : matches_mode ) {
    reco::CandidatePtr gen = iter.gen_;
    reco::GsfTrackRef trk(gsfTracks,iter.idx_);
    float dr2 = iter.dr2_;
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) == matched_gen.end() &&
	 std::find( matched_trk.begin(), matched_trk.end(), trk ) == matched_trk.end() && 
	 dr2 < dr_max_*dr_max_ ) { //@@ max deltaR considered
      gen2gsf.insert( std::pair<reco::CandidatePtr, reco::GsfTrackRef>( gen, trk ) );
      gen2gsf_dr2_mode.insert( std::pair<reco::CandidatePtr, float>( gen, dr2 ) );
      matched_gen.insert(gen);
      matched_trk.insert(trk);
      // store dr2 for "non-mode" eta/phi values 
      auto match = std::find_if(matches.begin(), 
				matches.end(), 
				[gen](const Match& m) {return m.gen_ == gen;}
				);
      if ( match != matches.end() ) {
	gen2gsf_dr2.insert( std::pair<reco::CandidatePtr, float>( match->gen_, match->dr2_ ) );
      } else { 
	gen2gsf_dr2.insert( std::pair<reco::CandidatePtr, float>( gen, -1. ) );
      }
      // Check if match found for all gen particles
      if ( matched_gen.size() >= electrons_from_B.size() &&
	   matched_trk.size() >= electrons_from_B.size() ) {
	break;
      }
    }
  }

  if ( verbose_ > 0 ) {
    if ( gen2gsf.size() <  electrons_from_B.size() ) {
      std::cout << "WARNING: Fewer matches than GEN electrons!" << std::endl;
    }
  }

//  std::cout << "gen2gsf size: " << gen2gsf.size() << std::endl;
//  size_t idx = 0; 
//  for ( auto iter : gen2gsf ) {
//    reco::CandidatePtr gen = iter.first;
//    reco::GsfTrackRef trk = iter.second;
//    std::cout << "  idx: " << idx
//      	      << " gen key: " << gen.key()
//	      << " trk key: " << trk.key()
//	      << " dr2: " << gen2gsf_dr2[iter.first]
//	      << std::endl
//	      << "  gen pt: " << gen->pt()
//	      << " gen eta: " << gen->eta()
//	      << " gen phi: " << gen->phi()
//	      << std::endl
//	      << "  trk pt: " << trk->ptMode()
//	      << " trk eta: " << trk->etaMode()
//	      << " trk phi: " << trk->phiMode()
//	      << std::endl;
//    idx++;
//  }
  
//  gen2gsf.clear();
//  gen2gsf_dr.clear();
//  other_gsf.clear();
//  std::set<reco::CandidatePtr> matched_gen;
//  std::set<reco::GsfTrackRef> matched_gsf;
//
//  // Match GsfTracks to GEN ele from B decays
//  for ( const auto& gen : electrons_from_B ) {
//    size_t best_idx = 666;
//    //double min_dr2 = dr_max_*dr_max_;
//    double min_dr2 = 0.01;
//    for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
//      reco::GsfTrackRef gsf(gsfTracks,idx);
//      double dr2 = reco::deltaR2(gsf->etaMode(),
//				 gsf->phiMode(),
//				 gen->eta(), 
//				 gen->phi() );
//      if ( dr2 < min_dr2 ) {
//	min_dr2 = dr2;
//	best_idx = idx;
//      }
//    }
//    //if ( min_dr2 < dr_max_*dr_max_ ) {
//    if ( best_idx < 666 && min_dr2 < 0.01 ) {
//      reco::GsfTrackRef gsf(gsfTracks,best_idx);
//      gen2gsf.insert( std::pair<reco::CandidatePtr, reco::GsfTrackRef>( gen, gsf ) );
//      gen2gsf_dr.insert( std::pair<reco::CandidatePtr, float>( gen, min_dr2 ) );
//      matched_gen.insert(gen);
//      matched_gsf.insert(gsf);
//    }
//  }

  // Add null reco::GsfTrack Ref for unmatched GEN particles
  for ( const auto& gen : electrons_from_B ) {
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) != matched_gen.end() ) { continue; }
    gen2gsf[gen] = reco::GsfTrackRef();
    gen2gsf_dr2_mode[gen] = -1.;
    gen2gsf_dr2[gen] = -1.;
  }

  // Identify GsfTracks from background, i.e. not matched to a GEN electron
  other_gsf.reserve( gsfTracks->size() - matched_trk.size() );
  for ( size_t idx = 0; idx < gsfTracks->size(); idx++ ) {
    reco::GsfTrackRef gsf(gsfTracks,idx);
    if ( matched_trk.find(gsf) == matched_trk.end() ) {
      other_gsf.push_back(gsf);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchGenToEle( const std::set<reco::CandidatePtr>& electrons_from_B,
				const edm::Handle< edm::View<reco::GsfElectron> >& gsfElectrons,
				std::map<reco::CandidatePtr, reco::GsfElectronPtr>& gen2ele,
				std::map<reco::CandidatePtr, float>& gen2ele_dr2,
				std::vector<reco::GsfElectronPtr>& other_ele ) {

  gen2ele.clear();
  gen2ele_dr2.clear();
  
  // Identify all possible matches b/w GEN ele and tracks, store indices and dr2
  std::vector<Match> matches;
  matches.reserve( electrons_from_B.size()*gsfElectrons->size() );
  for ( auto gen : electrons_from_B ) {
    for ( size_t iele = 0; iele < gsfElectrons->size(); ++iele ) {
      reco::GsfElectronPtr ele(gsfElectrons,iele);
      matches.emplace_back(gen,iele,reco::deltaR2(ele->eta(),
						  ele->phi(),
						  gen->eta(), 
						  gen->phi()) );
    }
  }
  
  // Sort all matches by dr2
  std::sort(matches.begin(), matches.end(), Match::compare_by_dr2);

  // Store best matches by dr2
  std::set<reco::CandidatePtr> matched_gen;
  std::set<reco::GsfElectronPtr> matched_ele;
  for ( auto iter : matches ) {
    reco::CandidatePtr gen = iter.gen_;
    reco::GsfElectronPtr ele(gsfElectrons,iter.idx_);
    float dr2 = iter.dr2_;
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) == matched_gen.end() &&
	 std::find( matched_ele.begin(), matched_ele.end(), ele ) == matched_ele.end() && 
	 dr2 < dr_max_*dr_max_ ) { //@@ max deltaR considered
      gen2ele.insert( std::pair<reco::CandidatePtr, reco::GsfElectronPtr>( gen, ele ) );
      gen2ele_dr2.insert( std::pair<reco::CandidatePtr, float>( gen, dr2 ) );
      matched_gen.insert(gen);
      matched_ele.insert(ele);
      // Check if match found for all gen particles
      if ( matched_gen.size() >= electrons_from_B.size() &&
	   matched_ele.size() >= electrons_from_B.size() ) {
	break;
      }
    }
  }

  if ( verbose_ > 0 ) {
    if ( gen2ele.size() <  electrons_from_B.size() ) {
      std::cout << "WARNING: Fewer matches than GEN electrons!" << std::endl;
    }
  }

  // Add null reco::GsfElectron Ref for unmatched GEN particles
  for ( const auto& gen : electrons_from_B ) {
    if ( std::find( matched_gen.begin(), matched_gen.end(), gen ) != matched_gen.end() ) { continue; }
    gen2ele[gen] = reco::GsfElectronPtr();
    gen2ele_dr2[gen] = -1.;
  }

  // Identify GsfTracks from background, i.e. not matched to a GEN electron
  other_ele.reserve( gsfElectrons->size() - matched_ele.size() );
  for ( size_t idx = 0; idx < gsfElectrons->size(); idx++ ) {
    reco::GsfElectronPtr ele(gsfElectrons,idx);
    if ( matched_ele.find(ele) == matched_ele.end() ) {
      other_ele.push_back(ele);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////////
//// 
//void IDFeatures::matchTrkToGsf( const edm::Handle< std::vector<reco::Track> >& ctfTracks,	
//				const edm::Handle< std::vector<reco::ElectronSeed> >& eleSeeds,
//				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
//				std::map<reco::TrackRef, reco::ElectronSeedRef>& trk2seed,
//				std::map<reco::ElectronSeedRef, reco::GsfTrackRef>& seed2gsf ) {
//
//
//  // Match recoTracks to gsfTracks
//  std::map<reco::TrackRef, reco::GsfTrackRef> trk2gsf;
//  matchTrkToGsf( ctfTracks_, gsfTracks, trk2gsf );
//  
//  // Match electronSeeds to gsfTracks
//  seed2gsf.clear();
//  matchSeedToGsf( eleSeeds, gsfTracks, seed2gsf );
//
//  // 
//  typedef std::map<reco::TrackRef, reco::GsfTrackRef>::const_iterator Iter;
//  trk2seed.clear();
//  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
//    reco::GsfTrackRef gsf(gsfTracks, idx);
//    if ( gsf.isNull() ) { 
//      std::cout << "ERROR! Null GsfTrackRef!" << std::endl;
//      continue; 
//    }
//    Iter iter = trk2gsf.find(trk);
//    if ( iter == trk2gsf.end() ) { 
//      std::cout << "ERROR! Null TrackRef!" << std::endl;
//      continue;
//    }
//    if ( iter->second.isNull() ) {
//      if ( verbose_ > 2 ) { std::cout << "INFO: Null GsfTrackRef in trk2gsf!" << std::endl; }
//      continue;
//    }
//
//
//
//    trk2gsf[iter->first];
//
//
//    int val = 20;
//    
//    auto result = std::find_if(seed2gsf.begin(),
//			       seed2gsf.end(),
//			       [val](const auto& ) {return mo.second == val; });
//    
//    //RETURN VARIABLE IF FOUND
//    if(result != mapObject.end())
//      int foundkey = result->first;
//    
//    
//      matched.insert(trk);
//    
//  
//}

////////////////////////////////////////////////////////////////////////////////
// 
void IDFeatures::matchTrkToGsf( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				std::map<reco::TrackRef, reco::GsfTrackRef>& trk2gsf ) {
  
  trk2gsf.clear();
  std::set<reco::TrackRef> matched;
  typedef std::map<reco::TrackRef, reco::GsfTrackRef>::const_iterator Iter;

  // Store ctfTrack -> gsfTrack link
  for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
    reco::GsfTrackRef gsf(gsfTracks, idx);
    if ( gsf.isNull() ) { 
      std::cout << "ERROR! Null GsfTrackRef!" << std::endl;
      continue; 
    }
    edm::RefToBase<TrajectorySeed> seed = gsf->seedRef();
    if ( seed.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO! Null TrajectorySeedRef!" << std::endl; }
      continue; 
    }
    reco::ElectronSeedRef ele_seed = seed.castTo<reco::ElectronSeedRef>();
    if ( ele_seed.isNull() ) { 
      std::cout << "ERROR! Null ElectronSeedRef!" << std::endl;
      continue; 
    }
    reco::TrackRef trk = ele_seed->ctfTrack();
    if ( trk.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO! Null TrackRef!" << std::endl; }
      continue;
    }
    if ( trk.id() != ctfTracks.id() ) { 
      std::cout << "ERROR! Collections do not match!" << std::endl;
      continue;
    }
    Iter iter = trk2gsf.find(trk);
    if ( iter == trk2gsf.end() ) { 
      trk2gsf[trk] = gsf; 
      matched.insert(trk);
    } else {
      std::cout << "ERROR! Single TrackRef mapped to many GsfTrackRefs!" << std::endl;
      if ( trk.isNonnull() && iter->first.isNonnull() ) {
	std::cout << " old trk key, pT, eta, phi: " 
		  << std::setw(3) 
		  << iter->first.key() << " "
		  << iter->first->pt() << " "
		  << iter->first->eta() << " "
		  << iter->first->phi() << " " 
		  << std::endl 
		  << " old gsf key, pT, eta, phi: " 
		  << std::setw(3) 
		  << iter->second.key() << " "
		  << iter->second->pt() << " "        //@@ Mode?
		  << iter->second->eta() << " "       //@@ Mode?
		  << iter->second->phi() << ", dR: "  //@@ Mode?
		  << sqrt( reco::deltaR2( iter->first->eta(),
					  iter->first->phi(),
					  iter->second->eta(),     //@@ Mode?
					  iter->second->phi() ) )  //@@ Mode?
		  << std::endl 
		  << " new trk key, pT, eta, phi: " 
		  << std::setw(3) 
		  << trk.key() << " "
		  << trk->pt() << " "
		  << trk->eta() << " "
		  << trk->phi() << " " 
		  << std::endl 
		  << " new gsf key, pT, eta, phi: " 
		  << std::setw(3) 
		  << gsf.key() << " "
		  << gsf->pt() << " "        //@@ Mode?
		  << gsf->eta() << " "       //@@ Mode?
		  << gsf->phi() << ", dR: "  //@@ Mode?
		  << sqrt( reco::deltaR2( trk->eta(), 
					  trk->phi(), 
					  gsf->etaMode(),
					  gsf->phiMode() ) )
		  << std::endl; 
	if ( reco::deltaR2( trk->eta(), 
			    trk->phi(),
			    gsf->etaMode(),
			    gsf->phiMode() ) <
	     reco::deltaR2( iter->first->eta(),
			    iter->first->phi(),
			    iter->second->eta(), 
			    iter->second->phi() )  ) {
	  std::cout << "INFO: Overwriting ctfTrackRef->gsfTrackRef entry!" << std::endl;
	  trk2gsf[trk] = gsf; 
	  matched.insert(trk);
	}
      }
    }
    
  } // GsfTracks loop
    
  // Add null gsfTrack Ref for unmatched ctfTracks
  for ( size_t idx = 0; idx < ctfTracks->size(); ++idx ) {
    reco::TrackRef trk(ctfTracks, idx);
    if ( std::find( matched.begin(), matched.end(), trk ) != matched.end() ) { continue; }
    trk2gsf[trk] = reco::GsfTrackRef();
  }

//  int cntr = 0;
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
void IDFeatures::checkGsfToTrk( const edm::Handle< std::vector<reco::Track> >& ctfTracks,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				const edm::Handle< edm::Association< std::vector<reco::Track> > >& gsf2trk,
				const std::map<reco::TrackRef, reco::GsfTrackRef>& trk2gsf 
				) {
  

  if ( trk2gsf.empty() ) {
    std::cout << "ERROR! (checkGsfToTrk) Empty trk2gsf map!" << std::endl;
  }
  if ( gsfTracks->size() != gsf2trk->size() ) {
    std::cout << "ERROR! (checkGsfToTrk) Different sizes for gsfTracksH_ and gsf2trk!" 
	      << " # gsfTracks->size(): " << gsfTracks->size()
	      << " # gsf2trk entries: " << gsf2trk->size()
	      << std::endl;
  }
  unsigned int cntr = 0;
  typedef std::map<reco::TrackRef, reco::GsfTrackRef>::const_iterator Iter;
  for ( Iter iter = trk2gsf.begin(); iter != trk2gsf.end(); ++iter ) {
    if ( iter->second.isNonnull() ) { cntr++; }
  }
  if ( cntr != gsf2trk->size() ) {
    std::cout << "ERROR! (checkGsfToTrk) Different number of GsfTrackRefs in trk2gsf and gsf2trk!" 
	      << " # nonnull Refs: " << cntr
	      << " # gsf2trk entries: " << gsf2trk->size()
	      << std::endl;

    for ( size_t idx = 0; idx < gsfTracks->size(); ++idx ) {
      reco::GsfTrackRef gsf(gsfTracks,idx);
      if ( gsf.isNull() ) {
	std::cout << "ERROR! (checkGsfToTrk) Null GsfTrackRef!" << std::endl;
	continue; 
      }
      reco::TrackRef trk = (*gsf2trk)[gsf];
      if ( trk.isNull() ) {
	std::cout << "ERROR! (checkGsfToTrk) Null TrackRef in gsf2trk!" << std::endl;
	continue; 
      }
      if ( trk.id() != ctfTracks.id() ) {
	std::cout << "ERROR! (checkGsfToTrk) Collections do not match!" << std::endl;
	continue; 
      }
    }

  }
  
  typedef std::map<reco::TrackRef, reco::GsfTrackRef>::const_iterator Iter;
  for ( Iter iter = trk2gsf.begin(); iter != trk2gsf.end(); ++iter ) {
    if ( iter->first.isNull() ) { 
      std::cout << "ERROR! (checkGsfToTrk) Null TrackRef in trk2gsf!" << std::endl;
      continue;
    }
    if ( iter->second.isNull() ) { 
      if ( verbose_ > 2 ) { std::cout << "INFO: Null GsfTrackRef in trk2gsf!" << std::endl; }
      continue;
    }
    reco::TrackRef trk = (*gsf2trk)[iter->second];
    if ( trk.isNull() ) {
      std::cout << "ERROR! (checkGsfToTrk) Null TrackRef in gsf2trk!" << std::endl;
      continue;
    }
    if ( trk.id() != iter->first.id() ) {
      std::cout << "ERROR! (checkGsfToTrk) Collections do not match!" << std::endl;
      continue;
    }
    if ( trk.key() != iter->first.key() ) {
      std::cout << "ERROR! (checkGsfToTrk) Keys do not match!" << std::endl;
      continue;
    }
//    double dr2 = reco::deltaR2( iter->first->eta(),
//				iter->first->phi(),
//				iter->second->eta(), 
//				iter->second->phi() );
//    double dr2_mode = reco::deltaR2( iter->first->eta(),
//				     iter->first->phi(),
//				     iter->second->etaMode(), 
//				     iter->second->phiMode() );
//    static double threshold = 0.02*0.02;
//    if ( dr2_mode > threshold ) {
//      std::cout << "WARNING! Large deltaR between reco::Tracks and gsfTracks!" 
//		<< " dR: " << sqrt(dr2) 
//		<< " dR(mode): " << sqrt(dr2_mode) 
//		<< " ratio: " << ( dr2_mode > 0. ? sqrt(dr2)/sqrt(dr2_mode) : -1 )
//		<< std::endl;
//    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// 
reco::GsfTrackRef IDFeatures::matchGsfToEgammaGsf( const reco::GsfTrackRef gsfTrack,
						   const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks_egamma ) {
  size_t best_idx = -1;
  double min_dr2 = dr_threshold_*dr_threshold_;
  for ( size_t idx = 0; idx < gsfTracks_egamma->size(); ++idx ) {
    reco::GsfTrackRef gsf_egamma(gsfTracks_egamma,idx);
    double dr2 = reco::deltaR2(gsf_egamma->etaMode(), //@@ use mode value for eta ?!
			       gsf_egamma->phiMode(), //@@ use mode value for phi ?!
			       gsfTrack->etaMode(),  //@@ use mode value for eta ?!
			       gsfTrack->phiMode() ); //@@ use mode value for phi ?!
    if ( dr2 < min_dr2 ) {
      min_dr2 = dr2;
      best_idx = idx;
    }
  }
  if ( min_dr2 < dr_threshold_*dr_threshold_ ) {
    return reco::GsfTrackRef(gsfTracks_egamma,best_idx);
  } else {
    return reco::GsfTrackRef(gsfTracks_egamma.id()); // null Ref
  }  
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDFeatures);








//  for ( auto gen_electron : electrons_from_B ) {
//
//    // Init and set truth label
//    ntuple_.reset();
//    ntuple_.is_e(true);
//
//    // Fill Rho, Event, and GEN branches
//    ntuple_.set_rho(*rho);
//    ntuple_.fill_evt(event.id());
//    ntuple_.fill_gen(gen_electron);
//
//    // Check if GEN electron is matched to a reco::Track
//    const auto& matched_trk = gen2trk.find(gen_electron);
//    if ( matched_trk != gen2trk.end() && matched_trk->second.isNonnull() ) { 
//      ntuple_.has_trk(true);
//      std::cout << "Matched GenParticle #" << matched_trk->first.key()
//		<< " to reco::Track #" << matched_trk->second.key()
//		<< std::endl;
//      
//      // Check if reco::Track is matched to an ElectronSeed
//      const auto& matched_seed = trk2seed.find(matched_trk->second);
//      if ( matched_seed != trk2seed.end() && matched_seed->second.isNonnull() ) { 
//	ntuple_.has_seed(true);
//	std::cout << "Matched reco::Track #" << matched_seed->first.key()
//		  << " to ElectronSeed #" << matched_seed->second.key()
//		  << std::endl;
//	
//	// Check if ElectronSeed is matched to a reco::GsfTrack
//	const auto& matched_gsf = seed2gsf.find(matched_seed->second);
//	if ( matched_gsf != seed2gsf.end() && matched_gsf->second.isNonnull() ) { 
//	  ntuple_.has_gsf(true);
//	  std::cout << "Matched ElectronSeed #" << matched_gsf->first.key()
//		    << " to GsfTrack #" << matched_gsf->second.key()
//		    << std::endl;
//	  
//	  ntuple_.fill_gsf(matched_gsf->second, *beamspot);
//
//	  // Store Seed BDT discrimator values
//	  float unbiased = (*mvaUnbiased)[matched_gsf->second];
//	  float ptbiased = (*mvaPtbiased)[matched_gsf->second];
//	  ntuple_.fill_seed( unbiased, ptbiased );
//	  
//	  // Check if GsfTrack is matched to an electron 
//	  const auto& matched_ele = gsf2ele.find(matched_gsf->second);
//	  if ( matched_ele != gsf2ele.end() ) {
//	    ntuple_.has_ele(true);
//	    std::cout << "Matched GsfTrack #" << matched_ele->first.key()
//		      << " to GsfElectron #" << matched_ele->second.key()
//		      << std::endl;
//	    
//	    //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	    float id_lowpt = -999.;
//	    if ( mvaValueLowPt.isValid() && mvaValueLowPt->size() == gsfElectrons->size() ) {
//	      id_lowpt = mvaValueLowPt->get( matched_ele->second.key() );
//	    }
//
//	    //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	    float mva_value = -999.;
//	    if ( mvaValue.isValid() && mvaValue->size() == gsfElectrons->size() ) {
//	      mva_value = mvaValue->get( matched_ele->second.key() );
//	    }
//	    
//	    //@@ dirty hack as is not in Event nor embedded in pat::Electron
//	    float conv_vtx_fit_prob = -999.;
//	    //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == gsfElectrons->size() ) {
//	    //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
//	    //}
//
//	    ntuple_.fill_ele( matched_ele->second, id_lowpt, mva_value, conv_vtx_fit_prob, *rho );
//
//	    //@@ Add SuperCluster vars?
//	    //ntuple_.fill_supercluster(matched_ele->second);
//
//	  } // Find reco::Track
//	} // Find ElectronSeed
//      } // Find GsfTrack
//    } // Find GenParticle
//
//    // Check if GEN electron is matched to a EGamma GsfTrack
//    const auto& matched_gsf_egamma = gen2gsf_egamma.find(gen_electron);
//    if ( matched_gsf_egamma != gen2gsf_egamma.end() ) { 
//
//      // Record presence of EGAMMA EGamma GsfTrack
//      ntuple_.is_egamma(true);
//      ntuple_.has_gsf(true);
//
//      // Check if GsfTrack is matched to an electron 
//      const auto& matched_egammaele = gsf2egammaele.find(matched_gsf_egamma->second);
//      if ( matched_egammaele != gsf2egammaele.end() ) {
//
//	ntuple_.has_ele(true);
//
//	//@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//	//float mva_value = -999.;
//	//if ( mvaValue.isValid() && mvaValue->size() == gsfElectrons->size() ) {
//	//mva_value = mvaValue->get( matched_ele->second.key() );
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

    
//  void fillOtherGSF( const edm::EventID& id,
//		     const std::vector<reco::GsfTrackRef>& other_gsf, 
//		     const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
//		     const edm::Handle< edm::ValueMap<float> >& mvaValueLowPt,
//		     const edm::Handle< edm::ValueMap<float> >& mvaValue,
//		     const edm::Handle< edm::ValueMap<float> >& mvaUnbiased,
//		     const edm::Handle< edm::ValueMap<float> >& mvaPtbiased,
//		     const edm::Handle<reco::BeamSpot>& beamspot,
//		     float rho, bool is_egamma, size_t nelectrons, bool store_id );

//  //////////
//  // Fill ntuple with background pT electrons (prescaled by 'fakesMultiplier' configurable)
//  //////////
//
//  fillOtherGSF( event.id(), 
//		other_gsf, 
//		gsf2ele, 
//		mvaValueLowPt, 
//		mvaValue, 
//		mvaUnbiased, 
//		mvaPtbiased, 
//		beamspot, 
//		*rho, 
//		false, 
//		electrons->size(), 
//		true );
//
//  fillOtherGSF( event.id(), 
//		other_gsf_egamma, 
//		gsf2egammaele, 
//		mvaValueLowPt, 
//		mvaValue, 
//		mvaUnbiased, 
//		mvaPtbiased, 
//		beamspot, 
//		*rho, 
//		true, 
//		electrons->size(), 
//		false );


//////////////////////////////////////////////////////////////////////////////////
//// Fill fake electrons (either PF or LowPt)
//void IDFeatures::fillOtherGSF( const edm::EventID& id,
//			       const std::vector<reco::GsfTrackRef>& other_gsf, 
//			       const std::map<reco::GsfTrackRef, reco::GsfElectronPtr>& gsf2ele,
//			       const edm::Handle< edm::ValueMap<float> >& mvaValueLowPt,
//			       const edm::Handle< edm::ValueMap<float> >& mvaValue,
//			       const edm::Handle< edm::ValueMap<float> >& mvaUnbiased,
//			       const edm::Handle< edm::ValueMap<float> >& mvaPtbiased,
//			       const edm::Handle<reco::BeamSpot>& beamspot,
//			       float rho, bool is_egamma, size_t nelectrons, bool store_id ) {
//  
//  for ( const auto& other : other_gsf ) {
//
//    if ( gRandom->Rndm() > prescale_ ) { continue; }
//
//    // Init and set truth label
//    ntuple_.reset();
//    ntuple_.is_other(true);
//    ntuple_.is_egamma(is_egamma);
//		
//    // Fill Rho, Event, and GEN branches
//    ntuple_.set_rho(rho);
//    ntuple_.fill_evt(id);
//    
//    ntuple_.fill_gsf(other, *beamspot);
//    
//    // Store Seed BDT discrimator values
//    if ( !is_egamma ) {
//      float unbiased = (*mvaUnbiased)[other];
//      float ptbiased = (*mvaPtbiased)[other];
//      ntuple_.fill_seed( unbiased, ptbiased );
//    }
//    
//    // Check if GsfTrack is matched to an electron 
//    const auto& matched_ele = gsf2ele.find(other);
//    if ( matched_ele != gsf2ele.end() ) {
//      
//      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//      float id_lowpt = -999.;
//      if ( store_id && mvaValueLowPt.isValid() && mvaValueLowPt->size() == nelectrons ) {
//	id_lowpt = mvaValueLowPt->get( matched_ele->second.key() );
//      }
//      
//      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//      float mva_value = -999.;
//      if ( store_id && mvaValue.isValid() && mvaValue->size() == nelectrons ) {
//        mva_value = mvaValue->get( matched_ele->second.key() );
//      }
//      
//      //@@ dirty hack as is not in Event nor embedded in pat::Electron
//      float conv_vtx_fit_prob = -999.;
//      //if ( convVtxFitProb.isValid() && convVtxFitProb->size() == electrons->size() ) {
//      //  conv_vtx_fit_prob = convVtxFitProb->get( matched_ele->second.key() );
//      //}
//      
//      ntuple_.fill_ele( matched_ele->second, id_lowpt, mva_value, conv_vtx_fit_prob, rho );
//      
//      //@@ Add SuperCluster vars?
//      //ntuple_.fill_supercluster(matched_ele->second);
//      
//    } // Find electron
//    
//    // Fill tree
//    tree_->Fill();
//    
//  } // Fill ntuple with background low pT electrons
//  
//}
