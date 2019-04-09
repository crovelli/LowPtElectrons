#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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

////////////////////////////////////////////////////////////////////////////////
//
class IDFeatures : public edm::EDAnalyzer {
  
public:
  
  explicit IDFeatures( const edm::ParameterSet& );
  ~IDFeatures() {}
  
  virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;
  virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;
  
private:
  
  void electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& prunedGenParticles,
		       const edm::Handle< std::vector<pat::PackedGenParticle> >& packedGenParticles,
		       std::set<pat::PackedGenParticleRef>& electrons_from_B );
  
  void electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& genParticles, 
		       std::set<reco::GenParticleRef>& electrons_from_B );
  
  bool isAncestor( const reco::Candidate* ancestor, 
		   const reco::Candidate* particle );

  void matchGenToGsf( const std::set<pat::PackedGenParticleRef>& electrons_from_B,
		      const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      std::map<pat::PackedGenParticleRef, reco::GsfTrackRef>& gen2gsf,
		      std::vector<reco::GsfTrackRef>& other_gsf );

  void matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
		      const edm::Handle< std::vector<pat::Electron> >& electrons,
		      std::map<reco::GsfTrackRef, pat::ElectronRef>& gsf2ele );

  reco::GsfTrackRef matchGsfToGsf( const reco::GsfTrackRef gsfTrack,
				   const edm::Handle< std::vector<reco::GsfTrack> >& pfGsfTracks );
    
private:
  
  edm::Service<TFileService> fs_;
  TTree* tree_;	
  IDNtuple ntuple_;
  bool check_from_B_;
  double dr_max_; // DeltaR matching
  double fakes_multiplier_;
  const edm::EDGetTokenT<double> rho_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  const edm::EDGetTokenT< std::vector<reco::GenParticle> > prunedGenParticles_;
  const edm::EDGetTokenT< std::vector<pat::PackedGenParticle> > packedGenParticles_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > gsfTracks_;
  const edm::EDGetTokenT< std::vector<reco::GsfTrack> > pfGsfTracks_;
  const edm::EDGetTokenT< std::vector<pat::Electron> > electrons_;
  const edm::EDGetTokenT< std::vector<pat::Electron> > pfElectrons_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaSeedUnbiased_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaSeedPtbiased_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaIDLowPt_;
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaIDv2_;
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
  rho_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
  beamspot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
  prunedGenParticles_(consumes< std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("prunedGenParticles"))),
  packedGenParticles_(consumes< std::vector<pat::PackedGenParticle> >(cfg.getParameter<edm::InputTag>("packedGenParticles"))),
  gsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("gsfTracks"))),
  pfGsfTracks_(consumes< std::vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("pfGsfTracks"))),
  electrons_(consumes< std::vector<pat::Electron> >(cfg.getParameter<edm::InputTag>("electrons"))),
  pfElectrons_(consumes< std::vector<pat::Electron> >(cfg.getParameter<edm::InputTag>("pfElectrons"))),
  mvaSeedUnbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVASeedUnbiased"))),
  mvaSeedPtbiased_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVASeedPtbiased"))),
  mvaIDLowPt_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDLowPt"))),
  mvaIDv2_(consumes< edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("MVAIDV2")))//,
  //convVtxFitProb_(consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("convVtxFitProb")))
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

  // Event collections 
  edm::Handle<double> rho;
  event.getByToken(rho_, rho);
  
  edm::Handle<reco::BeamSpot> beamspot;
  event.getByToken(beamspot_, beamspot);

  edm::Handle< std::vector<reco::GenParticle> > prunedGenParticles;
  event.getByToken(prunedGenParticles_, prunedGenParticles);
  
  edm::Handle< std::vector<pat::PackedGenParticle> > packedGenParticles;
  event.getByToken(packedGenParticles_, packedGenParticles);
  
  edm::Handle< std::vector<reco::GsfTrack> > gsfTracks;
  event.getByToken(gsfTracks_, gsfTracks);

  edm::Handle< std::vector<pat::Electron> > electrons;
  event.getByToken(electrons_, electrons);

  edm::Handle< std::vector<reco::GsfTrack> > pfGsfTracks;
  event.getByToken(pfGsfTracks_, pfGsfTracks);

  edm::Handle< std::vector<pat::Electron> > pfElectrons;
  event.getByToken(pfElectrons_, pfElectrons);

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
  std::set<pat::PackedGenParticleRef> electrons_from_B;
  electronsFromB( prunedGenParticles, packedGenParticles, electrons_from_B );

  // Match GEN electrons to low pT GsfTracks
  std::map<pat::PackedGenParticleRef, reco::GsfTrackRef> gen2gsf;
  std::vector<reco::GsfTrackRef> other_gsf;
  matchGenToGsf( electrons_from_B, gsfTracks, gen2gsf, other_gsf );

  // Match GEN electrons to EGamma GsfTracks
  std::map<pat::PackedGenParticleRef, reco::GsfTrackRef> gen2pfgsf;
  std::vector<reco::GsfTrackRef> other_pfgsf;
  matchGenToGsf( electrons_from_B, pfGsfTracks, gen2pfgsf, other_pfgsf );
  
  // Match GsfTracks to low pT electrons
  std::map<reco::GsfTrackRef, pat::ElectronRef> gsf2ele;
  matchGsfToEle( gsfTracks, electrons, gsf2ele );

  // Match GsfTracks to EGamma PF electrons
  std::map<reco::GsfTrackRef, pat::ElectronRef> gsf2pfele;
  matchGsfToEle( pfGsfTracks, pfElectrons, gsf2pfele );
  
  //////////
  // Fill ntuple with signal electrons
  //////////

  for ( auto gen : electrons_from_B ) {

    // Init and set truth label
    ntuple_.reset();
    ntuple_.is_e(true);
    ntuple_.is_other(false);

    // Fill Rho, Event, and GEN branches
    ntuple_.set_rho(*rho);
    ntuple_.fill_evt(event.id());
    ntuple_.fill_gen(gen);

    // Check if GEN electron is matched to a GsfTrack
    const auto& matched_gsf = gen2gsf.find(gen);
    if ( matched_gsf != gen2gsf.end() ) { 

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

      } // Find GsfTrack

    } // Find GenParticle

    // Check if GEN electron is matched to a PF GsfTrack
    const auto& matched_pfgsf = gen2pfgsf.find(gen);
    if ( matched_pfgsf != gen2pfgsf.end() ) { 

      // Record presence of CMS EGamma GsfTrack
      ntuple_.has_pfgsf(true);

      // Check if GsfTrack is matched to an electron 
      const auto& matched_pfele = gsf2pfele.find(matched_pfgsf->second);
      if ( matched_pfele != gsf2pfele.end() ) {

	ntuple_.has_pfele(true);

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
    reco::GsfTrackRef matched_pfgsf = matchGsfToGsf( other, pfGsfTracks );
    if ( matched_pfgsf.isNonnull() ) {

      ntuple_.has_pfgsf(true);

      // Check if GsfTrack is matched to a PF electron 
      const auto& matched_pfele = gsf2pfele.find(matched_pfgsf);
      if ( matched_pfele != gsf2pfele.end() ) {
	ntuple_.has_pfele(true);
      }

    } // Find GsfTrack

    // Fill tree
    tree_->Fill();

  } // Fill ntuple with background low pT electrons
  
}

////////////////////////////////////////////////////////////////////////////////
// Assumes pruned and packed GenParticles from MINIAOD
void IDFeatures::electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& prunedGenParticles,
				 const edm::Handle< std::vector<pat::PackedGenParticle> >& packedGenParticles,
				 std::set<pat::PackedGenParticleRef>& electrons_from_B ) {

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Examples
  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/packedGenParticles_cfi.py
  // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py

  electrons_from_B.clear();
  
  // Iterate through ("status 1") packed candidates
  for ( size_t ipa = 0; ipa < packedGenParticles->size(); ipa++ ) {
    pat::PackedGenParticleRef packed(packedGenParticles, ipa);
    
    // Find GEN electrons
    bool is_ele = std::abs(packed->pdgId()) == 11; // && packed->isLastCopy() 

    // Search for ancestors only for GEN electrons ... 
    if ( !is_ele ) { continue; }
    
    // Does GEN ele come from B decay?
    bool comes_from_B = false;
    for ( size_t ipr = 0; ipr < prunedGenParticles->size(); ipr++ ) {
      reco::GenParticleRef pruned(prunedGenParticles, ipr);
      comes_from_B = 
	std::abs(pruned->pdgId()) > 510 &&
	std::abs(pruned->pdgId()) < 546 &&
	packed->mother(0) != nullptr && 
	isAncestor( &*pruned, &*packed );

      // Is coming from a B
      if ( is_ele && ( comes_from_B || !check_from_B_ ) ) {
	electrons_from_B.insert(packed);
      }

    }

    //std::cout << "is_ele: " << is_ele << " comes_from_B: " << comes_from_B << std::endl;
    
  } // packedGenParticles loop
  
}

////////////////////////////////////////////////////////////////////////////////
// Assumes GenParticles from RECO/AOD
void IDFeatures::electronsFromB( const edm::Handle< std::vector<reco::GenParticle> >& genParticles,
				 std::set<reco::GenParticleRef>& electrons_from_B ) {
  
  electrons_from_B.clear();
  for ( size_t idx = 0; idx < genParticles->size(); idx++ ) {
    reco::GenParticleRef genp(genParticles, idx);
    
    // Last copy of GEN electron 
    bool is_ele = genp->isLastCopy() && std::abs(genp->pdgId()) == 11;
    
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
void IDFeatures::matchGenToGsf( const std::set<pat::PackedGenParticleRef>& electrons_from_B,
				const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				std::map<pat::PackedGenParticleRef, reco::GsfTrackRef>& gen2gsf,
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
      gen2gsf.insert( std::pair<pat::PackedGenParticleRef, reco::GsfTrackRef>( gen, gsf ) );
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
void IDFeatures::matchGsfToEle( const edm::Handle< std::vector<reco::GsfTrack> >& gsfTracks,
				const edm::Handle< std::vector<pat::Electron> >& electrons,
				std::map<reco::GsfTrackRef, pat::ElectronRef>& gsf2ele ) {
  gsf2ele.clear();
  for ( size_t idx = 0; idx < electrons->size(); ++idx ) {
    pat::ElectronRef ele(electrons, idx);
    reco::GsfTrackRef gsf = ele->gsfTrack();
    if ( gsf2ele.find(gsf) != gsf2ele.end() ) {
      std::cout << "THIS SHOULD NEVER HAPPEN! Multiple low pT electrons matched to the same GSFTrack?!"
		<< std::endl;
    } else {
      gsf2ele.insert( std::pair<reco::GsfTrackRef, pat::ElectronRef>(gsf, ele) );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// 
reco::GsfTrackRef IDFeatures::matchGsfToGsf( const reco::GsfTrackRef gsfTrack,
					     const edm::Handle< std::vector<reco::GsfTrack> >& pfGsfTracks ) {
  size_t best_idx = -1;
  double min_dr2 = dr_max_*dr_max_;
  for ( size_t idx = 0; idx < pfGsfTracks->size(); ++idx ) {
    reco::GsfTrackRef pfgsf(pfGsfTracks,idx);
    double dr2 = reco::deltaR2(pfgsf->etaMode(), //@@ use mode value for eta ?!
			       pfgsf->phiMode(), //@@ use mode value for phi ?!
			       gsfTrack->etaMode(),  //@@ use mode value for eta ?!
			       gsfTrack->phiMode() ); //@@ use mode value for phi ?!
    if ( dr2 < min_dr2 ) {
      min_dr2 = dr2;
      best_idx = idx;
    }
  }
  if ( min_dr2 < dr_max_*dr_max_ ) {
    return reco::GsfTrackRef(pfGsfTracks,best_idx);
  } else {
    return reco::GsfTrackRef(pfGsfTracks.id()); // null Ref
  }  
}

////////////////////////////////////////////////////////////////////////////////
//
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDFeatures);








  //////////
  // Fill ntuple with background PF electrons (prescaled by 'fakesMultiplier' configurable)
  //////////

//  indices.clear(); // reset
//  nfakes = 0; // reset
//  while ( indices.size() < other_pfgsf.size() && // stop when all tracks considered
//	  ( nfakes < fakes_multiplier_ || fakes_multiplier_ < 0 ) ) { // stop when fakesMultiplier is satisfied
//    int index = int( gRandom->Rndm() * other_pfgsf.size() ); // pick a track at random
//    if ( std::find( indices.begin(),
//		    indices.end(),
//		    index ) != indices.end() ) { continue; } // consider each track only once
//    indices.push_back(index); // record tracks used
//    nfakes++;
//    reco::GsfTrackRef other = other_pfgsf.at(index);
//
//    ntuple_.has_pfgsf(true);
//
//    // Check if GsfTrack is matched to a PF electron 
//    const auto& matched_pfele = gsf2pfele.find(other);
//    if ( matched_pfele != gsf2pfele.end() ) {
//
//      ntuple_.has_pfele(true);
//
//      //@@ dirty hack as ID is not in Event nor embedded in pat::Electron
//      float id_v2 = -999.;
//      //if ( mvaIDv2.isValid() && mvaIDv2->size() == electrons->size() ) {
//      //  id_v2 = mvaIDv2->get( matched_ele->second.key() );
//      //}
//      
//    } // Find GsfTrack
//    
//    // Fill tree
//    tree_->Fill();
//    
//  } // Fill ntuple with background PF electrons

