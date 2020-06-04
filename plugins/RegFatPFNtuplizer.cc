#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
#include "DataFormats/TrackReco/interface/DeDxData.h"
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
#include "LowPtElectrons/LowPtElectrons/interface/RegFatPFNtuple.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"



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
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

namespace reco{ 
  class SuperCluster;
  class GenParticle;
}



typedef std::map<unsigned long,int> PdgIds;

class RegFatPFNtuplizer : public edm::EDAnalyzer {
  
public:
  
  explicit RegFatPFNtuplizer( const edm::ParameterSet& );
  ~RegFatPFNtuplizer() {}


  // ---------------------------------------
  // Main methods
  virtual void beginJob(); // from Sam

  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const&){}
  //virtual void endJob();

  
  // Reads all collections from the Event
  void readCollections( const edm::Event&, const edm::EventSetup& );      

  // Delete collections at the end                                        
  void deleteCollections(); 

  // Wraps other methods to provide a sample of "signal" electrons
  void signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ); 

  // GEN-based method to provide a sample of "signal" electrons            
  void genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,  
			  float muon_pt = 7., float muon_eta = 1.5 );

  // ---------------------------------------
  // Utility methods
  template <typename T> bool validPtr( edm::Ptr<T>& ptr );

  bool gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk );
  bool gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );
  //
  bool egmToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk );
  bool egmToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed );
  bool egmSeedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk );

  bool extrapolate_to_ECAL(reco::TrackPtr kfTrackRef, float& eta_ECAL, float& phi_ECAL);

  template<typename T>
  void setToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& tagName){
    token = consumes<T>(iPara.getParameter<edm::InputTag>(tagName));
  }
  template<typename T>
  void setToken(std::vector<edm::EDGetTokenT<T> > & tokens,const edm::ParameterSet& iPara,const std::string& tagName){
    for(auto& tag: iPara.getParameter<std::vector<edm::InputTag> >(tagName)){
      tokens.push_back(consumes<T>(tag));
    }
  }




  static const reco::GenParticle* matchGenPart(float eta,float phi,const std::vector<reco::GenParticle>& genParts);
  const reco::SuperCluster*  matchSC(const reco::SuperCluster* scToMatch,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles);
  static const reco::SuperCluster*  matchSC(float eta,float phi,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles);

  
private:

  // from Sam
  TTree* egRegTree_;
  RegFatPFNtuple egRegTreeData_;
  // TTree* tree_;	
  // IDSlimNtuple ntuple_;


  
  // Misc  
  edm::Service<TFileService> fs_;
  int verbose_;
  bool check_from_B_;
  double prescale_; 
  int isAOD_;
  bool isMC_;
  double minTrackPt_;   
  bool tag_side_muon;


  edm::ESHandle<CaloTopology> caloTopoHandle;
  edm::ESHandle<EcalChannelStatus> chanStatusHandle;

  
  // Generic collections
  const edm::EDGetTokenT<double> rho_;
  edm::Handle<double> rhoH_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  edm::Handle<reco::BeamSpot> beamspotH_;

  const edm::EDGetTokenT< edm::View<reco::GenParticle> > genParticles_;       // AOD
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< edm::View<reco::Track> > ctfTracks_; // AOD
  edm::Handle< edm::View<reco::Track> > ctfTracksH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > packedCands_; // MINIAOD
  edm::Handle< edm::View<pat::PackedCandidate> > packedCandsH_;
  const edm::EDGetTokenT< edm::View<pat::PackedCandidate> > lostTracks_; // MINIAOD
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

  const edm::EDGetTokenT<edm::ValueMap< reco::DeDxData > > dEdx1Tag_;  // AOD only
  edm::Handle< edm::ValueMap< reco::DeDxData > > dEdx1H_;
  std::vector<const edm::ValueMap<reco::DeDxData>*> v_dEdx_;

  // Low pT collections
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_; // AOD and MINIAOD
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksH_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; // AOD
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;
  const edm::EDGetTokenT< edm::Association<reco::TrackCollection> > gsfTrackLinks_; // AOD
  edm::Handle<edm::Association<reco::TrackCollection> > gsfTrackLinksH_;
  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > packedCandLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > packedCandLinksH_;
  const edm::EDGetTokenT< edm::Association<pat::PackedCandidateCollection> > lostTrackLinks_; // MINIAOD
  edm::Handle<edm::Association<pat::PackedCandidateCollection> > lostTrackLinksH_;

  // Sam
  edm::EDGetTokenT<reco::VertexCollection>  verticesToken_;
  edm::Handle<reco::VertexCollection> verticesH_;
  std::vector<edm::EDGetTokenT<std::vector<reco::SuperCluster > > > scTokens_;
  std::vector<edm::Handle<reco::SuperClusterCollection > >   scH_;


  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puSumToken_;
  edm::Handle<std::vector<PileupSummaryInfo> > puSumH_;

  //  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron> > > eleAltTokens_;
  // edm::Handle<std::vector<reco::GsfElectron > > eleAH_;


  // EGamma collections                                                                           

  const edm::EDGetTokenT< std::vector<reco::ElectronSeed> > eleSeedsEGamma_; // AOD           
  edm::Handle< std::vector<reco::ElectronSeed> > eleSeedsEGammaH_; // AOD                     

  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_; // AOD               
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracksEGamma_MAOD_; // MINIAOD      
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracksEGammaH_;
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectronsEGamma_; // AOD              
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectronsEGamma_; // MINIAOD          
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsEGammaH_;


  PdgIds pdgids_;  
};

////////////////////////////////////////////////////////////////////////////////
RegFatPFNtuplizer::RegFatPFNtuplizer( const edm::ParameterSet& cfg ) 
  : egRegTree_(nullptr),
    egRegTreeData_(),
    verbose_(cfg.getParameter<int>("verbose")),
    check_from_B_(cfg.getParameter<bool>("checkFromB")),
    prescale_(cfg.getParameter<double>("prescale")),  
    isAOD_(-1),
    isMC_(true),
    minTrackPt_(cfg.getParameter<double>("minTrackPt")),
    // Generic collections
    rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
    rhoH_(),
    beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
    beamspotH_(),
    genParticles_(consumes<edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag >("genParticles"))),
  prunedGenParticles_(consumes<edm::View<reco::GenParticle> >(cfg.getParameter<edm::InputTag >("prunedGenParticles"))),
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
  dEdx1Tag_(consumes< edm::ValueMap<reco::DeDxData> >(cfg.getParameter<edm::InputTag>("dEdx1Tag"))),
// Low pT collections
  gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  gsfTracksH_(),
  gsfElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("gsfElectrons"))),
  patElectrons_(consumes< edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("patElectrons"))),
  gsfElectronsH_(),
  gsfTrackLinks_(consumes< edm::Association<reco::TrackCollection> >(cfg.getParameter<edm::InputTag>("gsfTrackLinks"))),
  gsfTrackLinksH_(),
  packedCandLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("packedCandLinks"))),
  packedCandLinksH_(),
  lostTrackLinks_(consumes< edm::Association<pat::PackedCandidateCollection> >(cfg.getParameter<edm::InputTag>("lostTrackLinks"))),
  lostTrackLinksH_(),
// EGamma collections                                                                                                             
  eleSeedsEGamma_(consumes< std::vector<reco::ElectronSeed> >(cfg.getParameter<edm::InputTag>("eleSeedsEGamma"))),
  eleSeedsEGammaH_(),
  gsfTracksEGamma_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma"))),
  gsfTracksEGamma_MAOD_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracksEGamma_MAOD"))),
  gsfTracksEGammaH_(),
  gsfElectronsEGamma_(consumes<edm::View<reco::GsfElectron > >(cfg.getParameter<edm::InputTag>("gsfElectronsEGamma"))),
  patElectronsEGamma_(consumes<edm::View<reco::GsfElectron > >(cfg.getParameter<edm::InputTag>("patElectronsEGamma"))),
  gsfElectronsEGammaH_(),
  pdgids_()
{
    //    tree_ = fs_->make<TTree>("tree","tree");
    //    ntuple_.link_tree(tree_);
    std::cout << "Verbosity level: "<< verbose_ << std::endl;
    setToken(scTokens_, cfg , "scTag"); 
    setToken(verticesToken_, cfg, "verticesCollection" );

    //  verticesToken_(consumes<edm::View<reco::VertexCollection> >(cfg.getParameter<edm::InputTag>("verticesCollection"))),
    //verticesH_();
    setToken(puSumToken_, cfg, "puSum" );
 
    //  puSumToken_(consumes< std::vector<PileupSummaryInfo> >(cfg.getParameter<edm::InputTag>("puSum"))),
    // puSumH_();
    //scH_();

  std::cout<<"RegFatPFNtuplizer"<<std::endl;

  edm::Service<TFileService> fs;
  fs->file().cd();
  egRegTree_ = new TTree("egRegTree","");
  //egRegTreeData_.setNrEnergies(eleAltTokens_.size());
  egRegTreeData_.createBranches(egRegTree_);


  }

void RegFatPFNtuplizer::beginJob()
{
}

namespace {
  template<typename T> 
  edm::Handle<T> getHandle(const edm::Event& iEvent,const edm::EDGetTokenT<T>& token)
  {
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
  template<typename T> 
  std::vector<edm::Handle<T> > getHandle(const edm::Event& iEvent,const std::vector<edm::EDGetTokenT<T> >& tokens)
  {
    std::vector<edm::Handle<T> > handles;
    for(auto& token : tokens){
      edm::Handle<T> handle;
      iEvent.getByToken(token,handle);
      handles.emplace_back(std::move(handle));
    }
    return handles;
  }
}



////////////////////////////////////////////////////////////////////////////////

void RegFatPFNtuplizer::beginRun( const edm::Run& run, const edm::EventSetup& es ) { 
  std::cout<<"RegFatPFNtuplizer::beginRun"<<std::endl;



}

////////////////////////////////////////////////////////////////////////////////
//
void RegFatPFNtuplizer::analyze( const edm::Event& event, const edm::EventSetup& setup ) {

  bool debug=false; 
  // Reset ntuple
  egRegTreeData_.reset();
  if(debug)std::cout<<"RegFatPFNtuplizer::analyze"<<std::endl;

  // Update all handles - MUST be called every event! 
  if(debug)std::cout<<"RegFatPFNtuplizer::analyze>> reading collections..."<<std::endl;
  readCollections(event,setup);
  if(debug)std::cout<<"RegFatPFNtuplizer::analyze>> read collections done"<<std::endl;

  // Gen level electrons from B                       
  std::set<reco::GenParticlePtr> signal_electrons;
  if (isMC_) signalElectrons(signal_electrons);          
  if (!tag_side_muon && check_from_B_   ) return;

  int nrVert = verticesH_->size();
  float nrPUInt = -1;
  float nrPUIntTrue = -1;
  if(puSumH_.isValid()){
    for(auto& puInfo : *puSumH_){
      if(puInfo.getBunchCrossing()==0){
        nrPUInt = puInfo.getPU_NumInteractions();
        nrPUIntTrue = puInfo.getTrueNumInteractions();
        break;
      }
    }
  }

  if(debug) std::cout<<"RegFatPFNtuplizer::analyze>> loop on electrons ..."<<std::endl; 


  // Loop over GSF electrons
  for( size_t electronlooper = 0; electronlooper < gsfElectronsEGammaH_->size(); electronlooper++ ) {

    const reco::GsfElectronPtr ele(gsfElectronsEGammaH_, electronlooper);

    // filter candidates: there must be a trk, a gsf track and a SC linked to the electron (should be always the case). 
	// Ele, trk and gsf_trk must have pT>0.5
	// trk must be high purity
	reco::GsfTrackPtr pfgsf = edm::refToPtr(ele->gsfTrack());    
	reco::TrackPtr pftrk;
	if ( !validPtr(pfgsf) )     continue;     
	if ( !egmToTrk(pfgsf,pftrk) ) continue;
	if ( !validPtr(pftrk) )     continue;     
	if ( pftrk->pt() < minTrackPt_ ) continue;
	

	if ( pftrk->pt() < minTrackPt_ ) continue;

	// Work on the electron candidate
	TVector3 eleTV3(0,0,0);
	eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());

    
	// ---------------------------------
	// Signal or fake electron, using gen-level info (-999 means nothing found with dR<0.05 )
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
	genTV3.SetPtEtaPhi(theGenParticle->pt(), theGenParticle->eta(), theGenParticle->phi());
	if (dRGenMin<0.05) {
	  // fill entry  
	  
	  // ---------------------------------
	  // Supercluster linked to electron
	  const reco::SuperClusterRef& scref =  ele->superCluster();

	  if(debug) std::cout<<"RegFatPFNtuplizer::analyze>> found one electron calling fill() ene,SCene="<<ele->energy()<<", "<<scref->energy()<< std::endl; 

      

 

      egRegTreeData_.fill(event,nrVert,*rhoH_,nrPUInt,nrPUIntTrue,
			  *ebRecHitsH_,*eeRecHitsH_,
			  *caloTopoHandle,
			  *chanStatusHandle,
			  &(*scref),&(*theGenParticle),&(*ele), &(*pfgsf) );
      
 
      egRegTree_->Fill();

    }




  } // LowPt electron looper
  if(debug) std::cout<<"RegFatPFNtuplizer::analyze>> end loop on Low Pt electrons"<<std::endl; 

  // Delete
  deleteCollections();
}

////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuplizer::readCollections( const edm::Event& event, const edm::EventSetup& setup ) {


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
  event.getByToken(beamspot_, beamspotH_);
  
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

  // KF tracks
  if ( isAOD_ == 1 ) { 
    event.getByToken(ctfTracks_, ctfTracksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCands_,packedCandsH_);
    event.getByToken(lostTracks_,lostTracksH_);
  }

  // RecHits 
  if ( isAOD_ == 1 ) {   
    event.getByToken(ebRecHits_, ebRecHitsH_);
    event.getByToken(eeRecHits_, eeRecHitsH_);
    if (!ebRecHitsH_.isValid()) std::cout << "rechits EB not valid" << std::endl;
    if (!eeRecHitsH_.isValid()) std::cout << "rechits EE not valid" << std::endl;
    if (ebRecHitsH_.isValid() && eeRecHitsH_.isValid()) ecalTools_ = new noZS::EcalClusterLazyTools(event, setup, ebRecHits_, eeRecHits_);
  }
  if ( isAOD_ == 0 ) {   
    event.getByToken(ebRecHitsEGM_, ebRecHitsEGMH_);
    event.getByToken(eeRecHitsEGM_, eeRecHitsEGMH_);
    if (!ebRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EB not valid" << std::endl;
    if (!eeRecHitsEGMH_.isValid()) std::cout << "rechitsEGM EE not valid" << std::endl;
    if (ebRecHitsEGMH_.isValid() && eeRecHitsEGMH_.isValid()) ecalTools_ = new noZS::EcalClusterLazyTools(event, setup, ebRecHitsEGM_, eeRecHitsEGM_);
  }


  // GsfTracks
  event.getByToken(gsfTracks_, gsfTracksH_);

  // Links
  if ( isAOD_ == 1 ) { 
    event.getByToken(gsfTrackLinks_, gsfTrackLinksH_);
  } else if ( isAOD_ == 0 ) { 
    event.getByToken(packedCandLinks_, packedCandLinksH_); 
    event.getByToken(lostTrackLinks_, lostTrackLinksH_); 
  }
  

  // dEdx
  if ( isAOD_ == 1 ) event.getByToken(dEdx1Tag_, dEdx1H_);

  //FC
  // EGamma collections                                                                      
  if      ( isAOD_ == 1 ) { event.getByToken(eleSeedsEGamma_, eleSeedsEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfTracksEGamma_, gsfTracksEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(gsfTracksEGamma_MAOD_, gsfTracksEGammaH_); }
  if      ( isAOD_ == 1 ) { event.getByToken(gsfElectronsEGamma_, gsfElectronsEGammaH_); }
  else if ( isAOD_ == 0 ) { event.getByToken(patElectronsEGamma_, gsfElectronsEGammaH_); }
  // end FC

  
  // collections from Sam 

  event.getByToken(verticesToken_, verticesH_);
  event.getByToken(puSumToken_, puSumH_);

  scH_ =getHandle(event,scTokens_);

  //  event.getByToken(scTokens_, scH_);
  // event.getByToken(eleAltTokens_, eleAH_);

  setup.get<CaloTopologyRecord>().get(caloTopoHandle);
  setup.get<EcalChannelStatusRcd>().get(chanStatusHandle);





}

void RegFatPFNtuplizer::deleteCollections( ) {

  delete ecalTools_;
}

// Gen-level electons from B
void RegFatPFNtuplizer::signalElectrons( std::set<reco::GenParticlePtr>& signal_electrons ) {

  signal_electrons.clear();
  std::set<reco::GenParticlePtr> electrons_from_B;
  genElectronsFromB(electrons_from_B);
  for ( auto gen : electrons_from_B ) { signal_electrons.insert(gen); }
}

// Gen-level electons from B 
void RegFatPFNtuplizer::genElectronsFromB( std::set<reco::GenParticlePtr>& electrons_from_B,
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

////////////////////////////////////////////////////////////////////////////////
template <typename T> 
bool RegFatPFNtuplizer::validPtr( edm::Ptr<T>& ptr ) {
  return ( ptr.isNonnull() && ptr.isAvailable() );
}

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::gsfToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  edm::RefToBase<TrajectorySeed> traj;
  if ( gsf->extra().isNonnull() && gsf->extra().isAvailable() ) { 
    traj = gsf->seedRef(); 
  } else {
    if ( verbose_ > 2 ) { // TrackExtra are not stored by default in MINIAOD
      std::cout << "ERROR: TrackExtra:" 
		<< " gsf->extra().isNull(): " << gsf->extra().isNull()
		<< " gsf->extra().isAvailable(): " << gsf->extra().isAvailable()
		<< std::endl; 
    }
    return false;
  }
  if ( traj.isNull() || !traj.isAvailable() ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR: TrajectorySeedRef:" 
		<< " traj.isNull(): " << traj.isNull()
		<< " traj.isAvailable(): " << traj.isAvailable()
		<< std::endl; 
    }
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::seedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  if ( !seed->isTrackerDriven() ) {
    if ( verbose_ > 3 ) {
      std::cout << "INFO! ElectronSeedPtr:"
		<< " seed->isTrackerDriven(): " << seed->isTrackerDriven()
		<< std::endl;
    }
  }
  trk = edm::refToPtr(seed->ctfTrack());
  if ( !validPtr(trk) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "INFO! TrackPtr:"
		<< " trk.isNull(): " << trk.isNull()
		<< " trk.isAvailable(): " << trk.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::gsfToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk ) {   

  // Attempt to navigate via Seed (and TrackExtra) to Track
  reco::ElectronSeedPtr seed;
  if ( gsfToSeed(gsf,seed) && seedToTrk(seed,trk) ) { return true; } // works for AOD 

  // ... if above fails (e.g. TrackExtra missing), attempt to use ...
  if ( isAOD_ == 1 ) {
    // ... track Association in AOD
    reco::GsfTrackRef gsf_ref(gsfTracksH_,(unsigned long)gsf.key());
    reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
    trk = edm::refToPtr(trk_ref);
    if ( validPtr(trk) ) { return true; }
  } else if ( isAOD_ == 0 ) {
    // ... "packedCand" Associations in MINIAOD
    reco::GsfTrackRef gsf_ref(gsfTracksH_,(unsigned long)gsf.key());
    pat::PackedCandidateRef packed_ref = (*packedCandLinksH_)[gsf_ref];
    if ( packed_ref.isAvailable() && 
	 packed_ref.isNonnull() && 
	 packed_ref->bestTrack() != nullptr ) { 
      trk = reco::TrackPtr(packed_ref->bestTrack(),packed_ref.key());
      if ( validPtr(trk) ) { return true; }
    }
    // ... "lostTrack" Associations in MINIAOD
    pat::PackedCandidateRef lost_ref = (*lostTrackLinksH_)[gsf_ref];
    if ( lost_ref.isAvailable() && 
	 lost_ref.isNonnull() && 
	 lost_ref->bestTrack() != nullptr ) { 
      trk = reco::TrackPtr(lost_ref->bestTrack(),lost_ref.key());
      if ( validPtr(trk) ) { return true; }
    }
  }

  return false;
}

// here methods for EGM electrons

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::egmToSeed( reco::GsfTrackPtr& gsf, reco::ElectronSeedPtr& seed ) {
  if ( !validPtr(gsf) ) {
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! GsfTrackPtr:"
		<< " gsf.isNull(): " << gsf.isNull()
		<< " gsf.isAvailable(): " << gsf.isAvailable()
		<< std::endl;
    }
    return false;
  }
  edm::RefToBase<TrajectorySeed> traj;
  if ( gsf->extra().isNonnull() && gsf->extra().isAvailable() ) { 
    traj = gsf->seedRef(); 
  } else {
    if ( verbose_ > 2 ) { // TrackExtra are not stored by default in MINIAOD
      std::cout << "ERROR: TrackExtra:" 
		<< " gsf->extra().isNull(): " << gsf->extra().isNull()
		<< " gsf->extra().isAvailable(): " << gsf->extra().isAvailable()
		<< std::endl; 
    }
    return false;
  }
  if ( traj.isNull() || !traj.isAvailable() ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR: TrajectorySeedRef:" 
		<< " traj.isNull(): " << traj.isNull()
		<< " traj.isAvailable(): " << traj.isAvailable()
		<< std::endl; 
    }
    return false;
  }
  seed = edm::refToPtr(traj.castTo<reco::ElectronSeedRef>());
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::egmSeedToTrk( reco::ElectronSeedPtr& seed, reco::TrackPtr& trk ) {
  if ( !validPtr(seed) ) { 
    if ( verbose_ > 0 ) {
      std::cout << "ERROR! ElectronSeedPtr:"
		<< " seed.isNull(): " << seed.isNull()
		<< " seed.isAvailable(): " << seed.isAvailable()
		<< std::endl;
    }
    return false;
  }
  if ( !seed->isTrackerDriven() ) {
    if ( verbose_ > 3 ) {
      std::cout << "INFO! ElectronSeedPtr:"
		<< " seed->isTrackerDriven(): " << seed->isTrackerDriven()
		<< std::endl;
    }
  }
  trk = edm::refToPtr(seed->ctfTrack());
  if ( !validPtr(trk) ) { 
    if ( verbose_ > 3 ) {
      std::cout << "INFO! TrackPtr:"
		<< " trk.isNull(): " << trk.isNull()
		<< " trk.isAvailable(): " << trk.isAvailable()
		<< std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool RegFatPFNtuplizer::egmToTrk( reco::GsfTrackPtr& gsf, reco::TrackPtr& trk ) {   
  bool debug=false; 
  if(debug)std::cout<<"ele loop - egmToTrk called "<<std::endl; 

  // Attempt to navigate via Seed (and TrackExtra) to Track
  reco::ElectronSeedPtr seed;
  if ( egmToSeed(gsf,seed) && egmSeedToTrk(seed,trk) ) { // works for AOD 
    return true; 
  }

  // ... if above fails (e.g. TrackExtra missing), attempt to use ...
  if ( isAOD_ == 1 ) {
    // ... track Association in AOD
    reco::GsfTrackRef gsf_ref(gsfTracksEGammaH_,(unsigned long)gsf.key());
    reco::TrackRef trk_ref = (*gsfTrackLinksH_)[gsf_ref];
    trk = edm::refToPtr(trk_ref);
    if ( validPtr(trk) ) { return true; }
  } else if ( isAOD_ == 0 ) {
    if(debug)std::cout<<"ele loop - last resort "<<std::endl; 
    // last resort... try DR match 
    float dRmin=99.0;
    size_t iptr=0;
    for ( auto& ptr : *packedCandsH_) {
      if(ptr.bestTrack() == nullptr) { continue;}
      reco::TrackPtr trkx(ptr.bestTrack(), iptr);
      float dRcur=deltaR2(trkx->eta(),trkx->phi(),gsf->eta(),gsf->phi());
      if(dRcur<dRmin){
	dRmin=dRcur;
	trk=trkx; 
      }
      ++iptr;
    }
    if(debug)std::cout<<"ele loop - other lost tracks "<<std::endl; 
    for ( auto& ptr : *lostTracksH_) {
      if(ptr.bestTrack() == nullptr) { continue;}
      reco::TrackPtr trkx(ptr.bestTrack(), iptr);
      float dRcur=deltaR2(trkx->eta(),trkx->phi(),gsf->eta(),gsf->phi());
      if(dRcur<dRmin){
	dRmin=dRcur;
	trk=trkx; 
      }
      ++iptr;
    }
    if(debug)std::cout<<"ele loop dR= "<<dRmin<<std::endl; 

    if(dRmin<99.0 && validPtr(trk)) return true;
  }

  return false;
}

bool RegFatPFNtuplizer::extrapolate_to_ECAL(reco::TrackPtr kfTrackRef, float& eta_ECAL, float& phi_ECAL){

  // Propagate 'electron' to ECAL surface
  double mass_=0.000511*0.000511; // ele mass 
  bool result=false;
  if (! validPtr(kfTrackRef) ) return result; 

  float p2=0;
  float px=0;
  float py=0;
  float pz=0;
  float vx=0;
  float vy=0;
  float vz=0;
  if ( kfTrackRef->extra().isAvailable() && kfTrackRef->extra().isNonnull() ) {
    p2=kfTrackRef->outerMomentum().Mag2();
    px=kfTrackRef->outerMomentum().x();
    py=kfTrackRef->outerMomentum().y();
    pz=kfTrackRef->outerMomentum().z();
    vx=kfTrackRef->outerPosition().x();
    vy=kfTrackRef->outerPosition().y();
    vz=kfTrackRef->outerPosition().z();
  } else {
    p2=pow( kfTrackRef->p() ,2 );
    px=kfTrackRef->px();
    py=kfTrackRef->py();
    pz=kfTrackRef->pz();
    vx=kfTrackRef->vx(); // must be in cm
    vy=kfTrackRef->vy();
    vz=kfTrackRef->vz();
  }


  float energy = sqrt(mass_ + p2);
  XYZTLorentzVector mom = XYZTLorentzVector(px,py,pz, energy);
  XYZTLorentzVector pos = XYZTLorentzVector(vx,vy,vz, 0.);

  float field_z=3.8;

  BaseParticlePropagator mypart(RawParticle(mom,pos), 0, 0, field_z);
  mypart.setCharge(kfTrackRef->charge());
  mypart.propagateToEcalEntrance(true); // true only first half loop , false more than one loop
  bool reach_ECAL=mypart.getSuccess(); // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps 

  // ECAL entry point for track
  GlobalPoint ecal_pos(mypart.x(), mypart.y(), mypart.z());

  eta_ECAL=ecal_pos.eta();
  phi_ECAL=ecal_pos.phi();

  return reach_ECAL; 

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RegFatPFNtuplizer);
