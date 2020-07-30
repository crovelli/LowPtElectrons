#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "LowPtElectrons/LowPtElectrons/interface/IDSlimNtuple.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TTree.h"
#include <set>
#include <vector>
#include <math.h>
#include <iostream> 
#include <boost/core/demangle.hpp>

namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<ElectronSeed> ElectronSeedPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
typedef std::map<unsigned long,int> PdgIds;

class IDSlimNtuplizer : public edm::EDAnalyzer {
  
public:
  
  explicit IDSlimNtuplizer( const edm::ParameterSet& );
  ~IDSlimNtuplizer() {}


  // ---------------------------------------
  // Main methods
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );      

  // Delete collections at the end                                        
  void deleteCollections(); 

  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectronsFromZ( std::set<reco::GenParticlePtr>& signal_electrons ); 
  void signalElectronsFromB( std::set<reco::GenParticlePtr>& signal_electrons ); 
  void signalElectronsFromGun( std::set<reco::GenParticlePtr>& signal_electrons ); 

  // GEN-based method to provide a sample of "signal" electrons            
  void genElectronsFromZ( std::set<reco::GenParticlePtr>& electrons_from_Z);  

  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,  
			  float muon_pt = 7., float muon_eta = 1.5 );

  void genElectronsFromGun( std::set<reco::GenParticlePtr>& electrons_from_gun);


  // ---------------------------------------
  // Utility methods
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );

  bool extrapolate_to_ECAL(reco::TrackPtr kfTrackRef, float& eta_ECAL, float& phi_ECAL);
  
private:

  // Regression stuff
  //std::unique_ptr<ModifyObjectValueBase> regression_;     // Low pt 
  //std::unique_ptr<ModifyObjectValueBase> regressionGsf_;  // Gsf
  
  // Misc  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDSlimNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double prescale_; 
  bool isMC_;
  double minTrackPt_;   
  double maxTrackPt_;   
  double maxTrackEta_;   
  bool tag_side_muon;
  
  // Generic collections
  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; 
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_;
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > packedCands_; 
  edm::Handle< edm::View<pat::PackedCandidate> > packedCandsH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > lostTracks_; 
  edm::Handle< edm::View<pat::PackedCandidate> > lostTracksH_;

  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHits_;
  edm::Handle<EcalRecHitCollection> ebRecHitsH_;
  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHits_;
  edm::Handle<EcalRecHitCollection> eeRecHitsH_;
  const edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsEGM_;
  edm::Handle<EcalRecHitCollection> ebRecHitsEGMH_;
  const edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsEGM_;
  edm::Handle<EcalRecHitCollection> eeRecHitsEGMH_;

  noZS::EcalClusterLazyTools *ecalTools_;       

  // Low pT collections
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; 
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;
  // const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; 
  // edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;
  const edm::EDGetTokenT<pat::ElectronCollection> patElectrons_;   
  edm::Handle<pat::ElectronCollection> gsfElectronsH_;

  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; 
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;
  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; 
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_;
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;

  PdgIds pdgids_;  
};

////////////////////////////////////////////////////////////////////////////////
IDSlimNtuplizer::IDSlimNtuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    prescale_(cfg.getParameter<double>("prescale")),  
    isMC_(true),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackPt_(cfg.getParameter<double>("maxTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
    beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
    beamspotH_(),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    ctfTracks_(consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("ctfTracks"))),
    ctfTracksH_(),
    packedCands_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCands"))),
    packedCandsH_(),
    lostTracks_(consumes< edm::View<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("lostTracks"))),
    lostTracksH_(),
    ebRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHits"))),
    ebRecHitsH_(),
    eeRecHits_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHits"))),
    eeRecHitsH_(),
    ebRecHitsEGM_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ebRecHitsEGM"))),
    ebRecHitsEGMH_(),
    eeRecHitsEGM_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("eeRecHitsEGM"))),
    eeRecHitsEGMH_(),
    // Low pT collections
    gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
    gsfTracksH_(),
    // gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
    // patElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
    patElectrons_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("patElectrons") )},
    gsfElectronsH_(),
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
    // 
    pdgids_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "Verbosity level: "<< verbose_ << std::endl;

    /*
    // Regression stuff - lowPtElectrons
    if( cfg.existsAs<edm::ParameterSet>("lowPtRegressionConfig") ) {
      const edm::ParameterSet& iconf = cfg.getParameterSet("lowPtRegressionConfig");
      const std::string& mname = iconf.getParameter<std::string>("modifierName");
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf);
      regression_.reset(plugin);
      edm::ConsumesCollector sumes = consumesCollector();
      regression_->setConsumes(sumes);
    } else {
      regression_.reset(nullptr);
    }

    // Regression stuff - GSF electrons
    if( cfg.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
      const edm::ParameterSet& iconf = cfg.getParameterSet("gsfRegressionConfig");
      const std::string& mname = iconf.getParameter<std::string>("modifierName");
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf);
      regressionGsf_.reset(plugin);
      edm::ConsumesCollector sumes = consumesCollector();
      regressionGsf_->setConsumes(sumes);
    } else {
      regressionGsf_.reset(nullptr);
    }
    */

  }

////////////////////////////////////////////////////////////////////////////////
// Initialise the weights LUT to filter fake tracks
void IDSlimNtuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) { }

////////////////////////////////////////////////////////////////////////////////
//
void IDSlimNtuplizer::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  // ----------------------------------------
  // To be set by hand 
  // ----------------------------------------
  //
  // Slim or Large size
  bool useEleGun=0;
  bool useZ=1;
  // bool useEnergyRegression=1;
  // ----------------------------------------

  // Reset ntuple
  ntuple_.reset();

  // Update all handles - MUST be called every event! 
  readCollections(event,setup);

  /*
  // Setup energy regressions for event
  if (useEnergyRegression) {
    regression_->setEvent(event);
    regression_->setEventContent(setup);
    regressionGsf_->setEvent(event);
    regressionGsf_->setEventContent(setup);
  }
  */

  // Gen level electrons from B                        
  std::set<reco::GenParticlePtr> signal_electrons;
  if (!useEleGun && !useZ) {
    if (isMC_) signalElectronsFromB(signal_electrons);          
    if (!tag_side_muon) return;
  } 
  // Gen level electrons from Z                        
  if (useZ) {
    if (isMC_) signalElectronsFromZ(signal_electrons);          
  } 
  // Gen level electrons from particle gun
  if (useEleGun) {
    if (isMC_) signalElectronsFromGun(signal_electrons);          
  }

  // Loop over low-pT electrons                  
  for( size_t electronlooper = 0; electronlooper < gsfElectronsH_->size(); electronlooper++ ) {
    
    // ---------------------------------
    // General event info - we save 1 entry per electron
    ntuple_.fill_evt( event.id() );
    ntuple_.is_egamma_ = false;
    ntuple_.weight_ = 1.;
  
    // Low pT electrons
    const reco::GsfElectronPtr ele(gsfElectronsH_, electronlooper);
    edm::Ref<pat::ElectronCollection> eleRef(gsfElectronsH_,electronlooper);   
    // std::cout << "E/P Ptr = " << ele->eSuperClusterOverP() << std::endl;
    // std::cout << "E/P Ref = " << eleRef->eSuperClusterOverP() << std::endl;

    // filter candidates: there must be a gsf with mode-pT>0.5
    reco::GsfTrackPtr gsf = edm::refToPtr(ele->gsfTrack());    
    if ( !validPtr(gsf) ) continue;     
    reco::TrackRef trk = ele->closestCtfTrackRef();  
    if ( gsf->ptMode() < minTrackPt_ ) continue;
    if ( fabs(gsf->etaMode()) > maxTrackEta_ ) continue;
    
    // Work on the electron candidate
    TVector3 eleTV3(0,0,0);
    eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());

    
    // ---------------------------------
    // Signal or fake electron, using gen-level info (-999 means nothing found with dR<0.1 )
    float dRGenMin=999.;
    reco::GenParticlePtr theGenParticle;
    TVector3 genTV3(0,0,0);
    for ( auto sig : signal_electrons ) {      
      genTV3.SetPtEtaPhi(sig->pt(), sig->eta(), sig->phi());
      float dR = eleTV3.DeltaR(genTV3);
      if (dR<dRGenMin) { 
	theGenParticle = sig;
	dRGenMin=dR;
      }
    }

    // Keep only electrons with dRGenMin<0.03 (signal) or >0.1 (fakes)
    if (dRGenMin>=0.03 && dRGenMin<0.1) continue;

    genTV3.SetPtEtaPhi(theGenParticle->pt(), theGenParticle->eta(), theGenParticle->phi());
    if (dRGenMin<0.1) {
      ntuple_.fill_gen( theGenParticle ); 
      ntuple_.gen_dR_ = dRGenMin;
      ntuple_.is_e_ = true;
      ntuple_.is_other_ = false;
      ntuple_.gen_tag_side_=tag_side_muon;
    } else { 
      ntuple_.fill_gen_default(); 
      ntuple_.gen_dR_ = dRGenMin;
      ntuple_.is_e_ = false;
      ntuple_.is_other_ = true;
      ntuple_.gen_tag_side_=tag_side_muon;
    }

    // prescale fake electrons 
    if (dRGenMin>=0.1) {
      if ( gRandom->Rndm() > prescale_  ) continue;
      ntuple_.weight_ = 1./prescale_;
    }  


    // ---------------------------------
    // Electron ID: dirty hack as ID is not in Event nor embedded in pat::Electron
    float mva_value = -999.;
    // std::cout << "Dumper: mvaValueLowPtH_->size() = " << mvaValueLowPtH_->size() << ", gsfElectronsH_->size() = " << gsfElectronsH_->size() << std::endl;
    // for (size_t ii=0; ii<mvaValueLowPtH_->size(); ii++ ) std::cout << ii << " " << (*mvaValueLowPtH_).get(ii) << std::endl;

    if ( mvaValueLowPtH_.isValid()) {
      mva_value = float((*mvaValueLowPtH_)[eleRef]);
      // std::cout << "Dumper: mva_value = " << mva_value << std::endl;
    } else {
      std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
    }

    //if ( mvaValueLowPtH_.isValid() && 
    // mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
    // mva_value = mvaValueLowPtH_->get( ele.key() );
    // std::cout << "Dumper: mva_value = " << mva_value << std::endl;
    // } else {
    // std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
    // }
      

      // ---------------------------------
    // GSF track linked to electron
    ntuple_.fill_gsf( gsf, *beamspotH_ );
    // TVector3 gsfTV3(0,0,0);
    // gsfTV3.SetPtEtaPhi(gsf->ptMode(), gsf->etaMode(), gsf->phiMode()); 
    // ntuple_.gsf_dr_ = eleTV3.DeltaR(gsfTV3);  
    // std::cout << "HAND: ntuple_.gsf_dr_ = " << eleTV3.DeltaR(gsfTV3) << std::endl;
    float unbiasedSeedBdt_ = (*mvaUnbiasedH_)[gsf];


    // ---------------------------------
    // Fill ele info 
    ntuple_.fill_ele( ele, mva_value, -999, *rhoH_, unbiasedSeedBdt_ );

    /*
    // Regression stuff
    float pre_ecal     = -999.;
    float pre_ecaltrk  = -999.;
    float post_ecal    = -999.;
    float post_ecaltrk = -999.;
    if (useEnergyRegression) {
      reco::GsfElectron newElectron(*ele);
      pre_ecal = newElectron.correctedEcalEnergy();
      pre_ecaltrk = newElectron.energy();
      regression_->modifyObject(newElectron);
      post_ecal = newElectron.correctedEcalEnergy();
      post_ecaltrk = newElectron.energy();
    }
    */
    
    // ---------------------------------
    // Supercluster linked to electron
    ntuple_.fill_supercluster_miniAOD(ele);  

    // ---------------------------------
    // KTF track linked to electron
    if ( trk.isNonnull() ) {
      // TVector3 trkTV3(0,0,0);
      // trkTV3.SetPtEtaPhi(trk->pt(), trk->eta(), trk->phi());  
      // ntuple_.trk_dr_ = eleTV3.DeltaR(trkTV3);  
      // std::cout << "HAND: ntuple_.trk_dr_ = " << eleTV3.DeltaR(trkTV3) << std::endl;
      ntuple_.trk_pt_  = trk->pt();
      ntuple_.trk_eta_ = trk->eta();
      ntuple_.trk_phi_ = trk->phi();
      ntuple_.trk_p_   = trk->p();
    } else { 
      ntuple_.fill_trk_default();
      // ntuple_.trk_dr_ = -999.;
      // std::cout << "HAND: ntuple_.trk_dr_ = -999" << std::endl;
    }

    reco::SuperClusterRef scp = ele->superCluster(); 

    /*
    // correct variables thanks to regression
    if (useEnergyRegression) {
      ntuple_.eid_match_SC_EoverP_=ntuple_.eid_match_SC_EoverP_*post_ecal/pre_ecal;
      ntuple_.eid_match_eclu_EoverP_=1/post_ecal-1/post_ecaltrk;
    }
    ntuple_.pre_ecal_=pre_ecal;
    ntuple_.pre_ecaltrk_=pre_ecaltrk;
    ntuple_.post_ecal_=post_ecal;
    ntuple_.post_ecaltrk_=post_ecaltrk;
    ntuple_.sc_raw_energy_=scp->rawEnergy();
    ntuple_.sc_energy_=scp->energy();
    */

    tree_->Fill();

    /*
    std::cout << "Dumper" << std::endl;
    std::cout << ntuple_.eid_rho_ << " " 
	      << ntuple_.eid_sc_eta_ << " " 
	      << ntuple_.eid_shape_full5x5_r9_ << " " 
	      << ntuple_.eid_sc_etaWidth_ << " " 
	      << ntuple_.eid_sc_phiWidth_ << " " 
	      << ntuple_.eid_shape_full5x5_HoverE_ << " " 
	      << ntuple_.eid_trk_nhits_ << " " 
	      << ntuple_.eid_trk_chi2red_ << " " 
	      << ntuple_.eid_gsf_chi2red_ << " " 
	      << ntuple_.eid_brem_frac_ << " " 
	      << ntuple_.eid_gsf_nhits_ << " " 
	      << ntuple_.eid_match_SC_EoverP_ << " " 
	      << ntuple_.eid_match_eclu_EoverP_ << " " 
	      << ntuple_.eid_match_SC_dEta_ << " " 
	      << ntuple_.eid_match_SC_dPhi_ << " " 
	      << ntuple_.eid_match_seed_dEta_ << " " 
	      << ntuple_.eid_sc_E_ << " " 
	      << ntuple_.eid_trk_p_ << " " 
	      << ntuple_.gsf_mode_p_ << " " 
	      << ntuple_.core_shFracHits_ << " " 
	      << ntuple_.seed_unbiased_ << " " 
	      << ntuple_.gsf_dr_ << " " 
	      << ntuple_.trk_dr_ << " " 
	      << ntuple_.sc_Nclus_ << " " 
	      << ntuple_.sc_clus1_nxtal_ << " " 
	      << ntuple_.sc_clus1_dphi_ << " " 
	      << ntuple_.sc_clus2_dphi_ << " " 
	      << ntuple_.sc_clus1_deta_ << " " 
	      << ntuple_.sc_clus2_deta_ << " " 
	      << ntuple_.sc_clus1_E_ << " " 
	      << ntuple_.sc_clus2_E_ << " " 
	      << ntuple_.sc_clus1_E_ov_p_ << " " 
	      << ntuple_.sc_clus2_E_ov_p_ << std::endl;
    std::cout << std::endl;
    */

  } // electron looper

  // Delete
  deleteCollections();
}

////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons
  event.getByToken(patElectrons_,gsfElectronsH_);
  if ( gsfElectronsH_.isValid() ) { 
    std::cout << "File contains MINIAOD data tier!" << std::endl;
  } else {
    throw cms::Exception(" Collection not found: ") 
      << " failed to find a standard miniAOD particle collection " << std::endl;
  }

  // Generic collections 
  event.getByToken(rho_, rhoH_);
  event.getByToken(beamspot_, beamspotH_);
  
  // GEN particles
  if ( isMC_ ) {
    event.getByToken(prunedGenParticles_, genParticlesH_);
    if ( !(genParticlesH_.isValid()) ) { 
      isMC_ = false;
      std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
    }
  }

  // KF tracks
  event.getByToken(packedCands_,packedCandsH_);
  event.getByToken(lostTracks_,lostTracksH_);

  // RecHits 
  event.getByToken(ebRecHitsEGM_, ebRecHitsEGMH_);
  event.getByToken(eeRecHitsEGM_, eeRecHitsEGMH_);
  if (!ebRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EB not valid" << std::endl;
  if (!eeRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EE not valid" << std::endl;
  if (ebRecHitsEGMH_.isValid() && eeRecHitsEGMH_.isValid()) ecalTools_ = new noZS::EcalClusterLazyTools(event, setup, ebRecHitsEGM_, eeRecHitsEGM_);

  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  event.getByToken(packedCandLinks_, packedCandLinksH_); 
  event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  
  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
}

void IDSlimNtuplizer::deleteCollections( ) {

  delete ecalTools_;
}

// Gen-level electons from B
void IDSlimNtuplizer::signalElectronsFromB( std::set<reco::GenParticlePtr>& signal_electrons ) {

  signal_electrons.clear();
  std::set<reco::GenParticlePtr> electrons_from_B;
  genElectronsFromB(electrons_from_B);
  for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
}

// Gen-level electons from B 
void IDSlimNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
					 float muon_pt, float muon_eta ) {   
  
  electrons_from_B.clear();
  tag_side_muon = false;

  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
    
    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "ERROR! GenParticlePtr:"
		<< " gen.isNull(): " << gen.isNull()
		<< " gen.isAvailable(): " << gen.isAvailable()
		<< std::endl;
      continue;
    }
    
    // Last copy of GEN electron 
    bool is_ele = std::abs(gen->pdgId()) == 11 && gen->isLastCopy(); //@@ not a method of Candidate

    //if (gen->mother())
    //std::cout << "GEN: " << idx << " << id = " << gen->pdgId() << ", motherID = " << gen->mother()->pdgId() << std::endl;
    
    // Does GEN ele comes from B decay?
    bool non_resonant = gen->numberOfMothers() >= 1 && gen->mother() &&   // has mother
      std::abs(gen->mother()->pdgId()) > 510 &&                           // mother is B
      std::abs(gen->mother()->pdgId()) < 546;                             // mother is B
    bool resonant = gen->numberOfMothers() >= 1 && gen->mother() &&       // has mother
      std::abs(gen->mother()->pdgId()) == 443 &&                          // mother is J/psi
      gen->mother()->numberOfMothers() >= 1 && gen->mother()->mother() && // has grandmother
      std::abs(gen->mother()->mother()->pdgId()) > 510 &&                 // grandmother is B
      std::abs(gen->mother()->mother()->pdgId()) < 546;                   // grandmother is B
    
    //  Check for tag side muon
    tag_side_muon |= ( std::abs(gen->pdgId()) == 13 && gen->isLastCopy() 
		       && 
		       ( ( gen->numberOfMothers() >= 1 &&
			   gen->mother() &&
			   std::abs(gen->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->pdgId()) < 546 && 
			   gen->mother()->pt() > muon_pt && 
			   std::abs(gen->mother()->eta()) < muon_eta ) 
			 ||
			 ( gen->numberOfMothers() >= 1 && 
			   gen->mother() &&
			   gen->mother()->numberOfMothers() >= 1 && 
			   gen->mother()->mother() &&
			   std::abs(gen->mother()->mother()->pdgId()) > 510 &&
			   std::abs(gen->mother()->mother()->pdgId()) < 546 && 
			   gen->mother()->mother()->pt() > muon_pt &&
			   std::abs(gen->mother()->mother()->eta()) < muon_eta ) ) );
    
    // is coming from a B
    if ( is_ele && ( ( resonant || non_resonant ) || !check_from_B_ ) ) {
      electrons_from_B.insert(gen);
      if ( verbose_ > 0 ) {
	std::cout << "electronsFromB: "
		  << " #signal_electrons: " << electrons_from_B.size()
		  << " resonant? " << resonant
		  << " non resonant? " << non_resonant
		  << " tag-side muon? " << tag_side_muon
		  << std::endl;
      }
    }
    
  } // genParticles loop

  // We don't want this for now
  //// if ( !tag_side_muon ) { electrons_from_B.clear(); }
}

// Gen-level electons from Z 
void IDSlimNtuplizer::signalElectronsFromZ( std::set<reco::GenParticlePtr>& signal_electrons ) {

  signal_electrons.clear();
  std::set<reco::GenParticlePtr> electrons_from_Z;
  genElectronsFromZ(electrons_from_Z);
  for ( auto gen : electrons_from_Z ) { signal_electrons.insert(gen); }
}

void IDSlimNtuplizer::genElectronsFromZ( std::set<reco::GenParticlePtr>& electrons_from_Z) {
  
  electrons_from_Z.clear();

  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {

    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "ERROR! GenParticlePtr:"
		<< " gen.isNull(): " << gen.isNull()
		<< " gen.isAvailable(): " << gen.isAvailable()
		<< std::endl;
      continue;
    }

    // Last copy of GEN electron 
    bool is_ele = std::abs(gen->pdgId()) == 11; 
    if (!is_ele) continue;

    // Does GEN ele comes from Z decay?
    if (gen->numberOfMothers()<1) continue;
    if (gen->mother()) {
      bool fromZ = gen->numberOfMothers() >= 1 && gen->mother() &&  std::abs(gen->mother()->pdgId())==23;
      if ( fromZ ) electrons_from_Z.insert(gen);
    }

  } // genParticles loop
}

// Gen-level electons from particle gun
void IDSlimNtuplizer::signalElectronsFromGun( std::set<reco::GenParticlePtr>& signal_electrons ) {

  signal_electrons.clear();
  std::set<reco::GenParticlePtr> electrons_from_gun;
  genElectronsFromGun(electrons_from_gun);
  for ( auto gen : electrons_from_gun ) { signal_electrons.insert(gen); }
}

// Gen-level electons from B 
void IDSlimNtuplizer::genElectronsFromGun( std::set<reco::GenParticlePtr>& electrons_from_gun) {
  
  electrons_from_gun.clear();

  for ( size_t idx = 0; idx < genParticlesH_->size(); idx++ ) {
    
    reco::GenParticlePtr gen(genParticlesH_, idx);
    if ( !validPtr(gen) ) {
      std::cout << "ERROR! GenParticlePtr:"
		<< " gen.isNull(): " << gen.isNull()
		<< " gen.isAvailable(): " << gen.isAvailable()
		<< std::endl;
      continue;
    }

    // Last copy of GEN electron 
    bool is_ele = std::abs(gen->pdgId()) == 11 && gen->isLastCopy(); //@@ not a method of Candidate
    if ( is_ele ) electrons_from_gun.insert(gen);
    
  } // genParticles loop
  
}

////////////////////////////////////////////////////////////////////////////////
template <typename T> 
bool IDSlimNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDSlimNtuplizer);
