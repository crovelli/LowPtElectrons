#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
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
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
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
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }

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
  void signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ); 
  void signalElectronsFromGun( std::set<reco::GenParticlePtr>& signal_electrons ); 

  // GEN-based method to provide a sample of "signal" electrons            
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,  
			  float muon_pt = 7., float muon_eta = 1.5 );
  void genElectronsFromGun( std::set<reco::GenParticlePtr>& electrons_from_gun);


  // ---------------------------------------
  // Utility methods
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );

private:

  // Regression stuff
  std::unique_ptr<ModifyObjectValueBase> regression_;     
  
  // Misc  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDSlimNtuple ntuple_;
  int verbose_;
  bool check_from_B_;
  double prescale_; 
  int isAOD_;
  bool isMC_;
  double minTrackPt_;   
  double maxTrackPt_;   
  double maxTrackEta_;   
  bool tag_side_muon;
  
  // Generic collections
  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_;       // AOD
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  // Low pT collections
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaUnbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaUnbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaPtbiased_; 
  edm::Handle< edm::ValueMap<float> > mvaPtbiasedH_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueLowPt_;
  edm::Handle< edm::ValueMap<float> > mvaValueLowPtH_;
};

////////////////////////////////////////////////////////////////////////////////
IDSlimNtuplizer::IDSlimNtuplizer( const edm::ParameterSet& cfg ) 
  : tree_(nullptr),
    ntuple_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    prescale_(cfg.getParameter<double>("prescale")),  
    isAOD_(-1),
    isMC_(true),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    maxTrackPt_(cfg.getParameter<double>("maxTrackPt")),
    maxTrackEta_(cfg.getParameter<double>("maxTrackEta")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
    genParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genParticles"))),
    prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
    genParticlesH_(),
    // Low pT collections
    gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
    patElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
    gsfElectronsH_(),
    mvaUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaUnbiased"))),
    mvaUnbiasedH_(),
    mvaPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaPtbiased"))),
    mvaPtbiasedH_(),
    mvaValueLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("mvaValueLowPt"))),
    mvaValueLowPtH_()
  {
    tree_ = fs_->make<TTree>("tree","tree");
    ntuple_.link_tree(tree_);
    std::cout << "Verbosity level: "<< verbose_ << std::endl;

    // Regression stuff - lowPtElectrons
    if( cfg.existsAs<edm::ParameterSet>("lowPtRegressionConfig") ) {
      const edm::ParameterSet& iconf = cfg.getParameterSet("lowPtRegressionConfig");
      const std::string& mname = iconf.getParameter<std::string>("modifierName");
      auto cc = consumesCollector();                  
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf, cc);
      regression_.reset(plugin);
      // edm::ConsumesCollector sumes = consumesCollector();
      // regression_->setConsumes(sumes);
    } else {
      regression_.reset(nullptr);
    }
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
  bool useEnergyRegression=1;
  // ----------------------------------------

  // Reset ntuple
  ntuple_.reset();

  // Update all handles - MUST be called every event! 
  readCollections(event,setup);

  // B-field
  edm::ESHandle<MagneticField> field;
  setup.get<IdealMagneticFieldRecord>().get(field);
  math::XYZVector zfield(field->inTesla(GlobalPoint(0, 0, 0)));
  
  // Setup energy regressions for event
  if (useEnergyRegression) {
    regression_->setEvent(event);
    regression_->setEventContent(setup);
  }

  // Gen level electrons from B                        
  std::set<reco::GenParticlePtr> signal_electrons;
  if (!useEleGun) {
    if (isMC_) signalElectrons(signal_electrons);          
    if (!tag_side_muon) return;
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
    int mva_id = -999;
    if ( mvaValueLowPtH_.isValid() && 
	 mvaValueLowPtH_->size() == gsfElectronsH_->size() ) {
      mva_value = mvaValueLowPtH_->get( ele.key() );
    } else {
      std::cout << "ERROR! Issue matching MVA output to GsfElectrons!" << std::endl;
    }

    // ---------------------------------
    // Fill ele info 
    float unbiasedSeedBdt_ = (*mvaUnbiasedH_)[gsf];
    float ptbiasedSeedBdt_ = (*mvaPtbiasedH_)[gsf];
    ntuple_.fill_bdt( unbiasedSeedBdt_, ptbiasedSeedBdt_ );
    ntuple_.fill_ele( ele, mva_value, mva_id, *rhoH_, unbiasedSeedBdt_, zfield.z() );

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

    // ---------------------------------
    // GSF track linked to electron
    ntuple_.fill_gsf( gsf );

    // ---------------------------------
    // Supercluster linked to electron
    ntuple_.fill_supercluster_miniAOD(ele);  

    // fill how many tracks there are around first second and third supercluster within dR<0.1 
    reco::SuperClusterRef scp = ele->superCluster(); 

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

    tree_->Fill();

  } // electron looper

  // Delete
  deleteCollections();
}

////////////////////////////////////////////////////////////////////////////////
void IDSlimNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {

  // Low pT electrons (and identify if data or MC and RECO/AOD or MINIAOD)
  if ( isAOD_ == -1 ) {
    event.getByToken(gsfElectrons_, gsfElectronsH_);
    if ( gsfElectronsH_.isValid() ) {
      isAOD_ = 1;
      std::cout << "File contains AOD data tier!" << std::endl;
    } else {
      event.getByToken(patElectrons_,gsfElectronsH_);
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
    event.getByToken(gsfElectrons_, gsfElectronsH_);
  } else if ( isAOD_ == 0 ) {
    event.getByToken(patElectrons_,gsfElectronsH_);
  } else {
    throw cms::Exception(" Invalid value for isAOD: ") 
      << isAOD_ 
      << std::endl;
  }

  // Generic collections 
  event.getByToken(rho_, rhoH_);
  
  // GEN particles
  if ( isMC_ ) {
    if ( isAOD_ == 1 ) { 
      event.getByToken(genParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "No GEN info found in AOD data tier!" << std::endl;
      }
    } else if ( isAOD_ == 0 ) { 
      event.getByToken(prunedGenParticles_, genParticlesH_);
      if ( !(genParticlesH_.isValid()) ) { 
	isMC_ = false;
	std::cout << "No GEN info found in MINIAOD data tier!" << std::endl;
      }
    }
  }

  // IDs
  event.getByToken(mvaUnbiased_, mvaUnbiasedH_);
  event.getByToken(mvaPtbiased_, mvaPtbiasedH_);
  event.getByToken(mvaValueLowPt_, mvaValueLowPtH_);
}

void IDSlimNtuplizer::deleteCollections( ) {
}

// Gen-level electons from B
void IDSlimNtuplizer::signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ) {

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
