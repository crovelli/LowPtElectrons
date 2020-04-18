#ifndef LowPtElectrons_LowPtElectrons_RegFatNtuple
#define LowPtElectrons_LowPtElectrons_RegFatNtuple

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "FWCore/Framework/interface/Event.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
// nuove 
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "TTree.h"



#include <vector>

class TTree;

namespace reco { typedef edm::Ptr<GenParticle> GenParticlePtr; }
namespace reco { typedef edm::Ptr<Track> TrackPtr; }
namespace reco { typedef edm::Ptr<GsfTrack> GsfTrackPtr; }
namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

namespace reco{
  class SuperCluster;
  class GenParticle;
}

namespace edm{
  class Event;
}
class CaloTopology;

struct ClustStruct {
  float clusterRawEnergy,clusterDEtaToSeed,clusterDPhiToSeed;
  static std::string contents(){return "clusterRawEnergy/F:clusterDEtaToSeed:clusterDPhiToSeed";}
  void clear(){clusterRawEnergy=clusterDEtaToSeed=clusterDPhiToSeed=0.;}
  void fill(float iClusterRawEnergy,float iClusterDEtaToSeed,float iClusterDPhiToSeed){
    clusterRawEnergy = iClusterRawEnergy;
    clusterDEtaToSeed = iClusterDEtaToSeed;
    clusterDPhiToSeed = iClusterDPhiToSeed;
  }
};
struct EleStruct {
  float et,energy,energyErr,ecalEnergy,ecalEnergyErr,eta,phi,trkEtaMode,trkPhiMode,trkPMode,trkPModeErr,fbrem,corrMean,corrSigma,hademTow,hademCone,trkPInn,trkPtInn,trkPVtx,trkPOut,trkChi2,trkNDof,ecalDrivenSeed,nrSatCrys,scRawEnergy,scRawESEnergy;
  static std::string contents(){return "et/F:energy:energyErr:ecalEnergy:ecalEnergyErr:eta:phi:trkEtaMode:trkPhiMode:trkPMode:trkPModeErr:fbrem:corrMean:corrSigma:hademTow:hademCone:trkPInn:trkPtInn:trkPVtx:trkPOut:trkChi2:trkNDof:ecalDrivenSeed:nrSatCrys:scRawEnergy:scRawESEnergy";}
  void clear(){et=energy=energyErr=ecalEnergy=ecalEnergyErr=eta=phi=trkEtaMode=trkPhiMode=trkPMode=trkPModeErr=fbrem=corrMean=corrSigma=hademTow=hademCone=trkPInn=trkPtInn=trkPVtx=trkPOut=trkChi2=trkNDof=ecalDrivenSeed=nrSatCrys=scRawEnergy=scRawESEnergy=0.;}
  void fill(const reco::GsfElectron& ele, const  reco::GsfTrack& gsf );
};
struct PFEleStruct {
  float et,energy,ecalEnergy,trkPMode,trkPtMode,trkP,scRawEnergy,scRawESEnergy;
  static std::string contents(){return "et/F:energy:ecalEnergy:trkPMode:trkPtMode:trkP:scRawEnergy:scRawESEnergy";}
  void clear(){et=energy=ecalEnergy=trkPMode=trkPtMode=trkP=scRawEnergy=scRawESEnergy=0.;}
  void fill(const reco::GsfElectron& ipfele, const  reco::GsfTrack& ipfgsf,  const  reco::Track& ipftrk );
};

struct EleEnergyStruct {
  float ecalTrk,ecalTrkErr,ecal,ecalErr;
  static std::string contents(){return "ecalTrk/F:ecalTrkErr:ecal:ecalErr";}
  void clear(){ecalTrk=ecalTrkErr=ecal=ecalErr=0.;}
  void fill(const reco::GsfElectron& ele){ecalTrk=ele.energy();ecalTrkErr=ele.p4Error(reco::GsfElectron::P4_COMBINATION);ecal=ele.ecalEnergy();ecalErr=ele.ecalEnergyError();}
};





struct SuperClustStruct {
  float rawEnergy,rawESEnergy,etaWidth,phiWidth,seedClusEnergy,numberOfClusters,numberOfSubClusters,clusterMaxDR,clusterMaxDRDPhi,clusterMaxDRDEta,clusterMaxDRRawEnergy,corrEnergy,scEta,scPhi,seedEta,seedPhi,dEtaSeedSC,dPhiSeedSC,isEB,iEtaOrX,iPhiOrY,iEtaMod5,iPhiMod2,iEtaMod20,iPhiMod20,etaGapCode,phiGapCode,nearbyChanStatus,corrEnergyAlt,rawEnergyAlt,nrClusAlt,scSinTheta,seedSinTheta;
  static std::string contents(){return "rawEnergy/F:rawESEnergy:etaWidth:phiWidth:seedClusEnergy:numberOfClusters:numberOfSubClusters:clusterMaxDR:clusterMaxDRDPhi:clusterMaxDRDEta:clusterMaxDRRawEnergy:corrEnergy:scEta:scPhi:seedEta:seedPhi:dEtaSeedSC:dPhiSeedSC:isEB:iEtaOrX:iPhiOrY:iEtaMod5:iPhiMod2:iEtaMod20:iPhiMod20:etaGapCode:phiGapCode:nearbyChanStatus:corrEnergyAlt:rawEnergyAlt:nrClusAlt:scSinTheta:seedSinTheta";}
  void clear(){
    rawEnergy=rawESEnergy=etaWidth=phiWidth=seedClusEnergy=numberOfClusters=numberOfSubClusters=clusterMaxDR=clusterMaxDRDPhi=clusterMaxDRDEta=clusterMaxDRRawEnergy=corrEnergy=scEta=scPhi=seedEta=seedPhi=dEtaSeedSC=dPhiSeedSC=isEB=iEtaOrX=iPhiOrY=iEtaMod5=iPhiMod2=iEtaMod20=iPhiMod20=etaGapCode=phiGapCode=nearbyChanStatus=corrEnergyAlt=rawEnergyAlt=nrClusAlt=scSinTheta=seedSinTheta=0.;
  }

  void fill(const reco::SuperCluster& sc,const EcalChannelStatus& ecalChanStatus);
};

struct ShowerShapeStruct {
  float e3x3,e5x5,seedClusEnergy,eMax,e2nd,eLeftRightDiffSumRatio,eTopBottomDiffSumRatio,sigmaIEtaIEta,sigmaIEtaIPhi,sigmaIPhiIPhi,e2x5Max,e2x5Top,e2x5Bottom,e2x5Left,e2x5Right,eTop,eBottom,eLeft,eRight;
  static std::string contents(){return "e3x3:e5x5:seedClusEnergy:eMax:e2nd:eLeftRightDiffSumRatio:eTopBottomDiffSumRatio:sigmaIEtaIEta:sigmaIEtaIPhi:sigmaIPhiIPhi:e2x5Max:e2x5Top:e2x5Bottom:e2x5Left:e2x5Right:eTop:eBottom:eLeft:eRight";}
  void clear(){
    e3x3=e5x5=seedClusEnergy=eMax=e2nd=eLeftRightDiffSumRatio=eTopBottomDiffSumRatio=sigmaIEtaIEta=sigmaIEtaIPhi=sigmaIPhiIPhi=e2x5Max=e2x5Top=e2x5Bottom=e2x5Left=e2x5Right=eTop=eBottom=eLeft=eRight=0.;
  }
  template<bool full5x5>
  void fill(const reco::CaloCluster& clus,const EcalRecHitCollection& ecalHitsEB,const EcalRecHitCollection& ecalHitsEE,const CaloTopology& topo);

  void fill(const reco::GsfElectron::ShowerShape& eleSS,const reco::GsfElectron& ele);

  void fillShowerShape(const reco::GsfElectron* ele);

};

struct EvtStruct {
  int runnr,lumiSec,eventnr;
  static std::string contents(){return "runnr/I:lumiSec:eventnr";}
  void clear(){runnr=lumiSec=eventnr=0;}
  void fill(const edm::Event& event);
};

struct GenInfoStruct {
  float energy,pt,eta,phi,pdgId,status,dR;
  static std::string contents(){return "energy/F:pt:eta:phi:pdgId:status:dR";}
  void clear(){energy=pt=eta=phi=pdgId=status=dR=0;}
  void fill(const reco::GenParticle& genPart, float iDR);
};



template<bool full5x5>
void ShowerShapeStruct::fill(const reco::CaloCluster& clus,const EcalRecHitCollection& ecalHitsEB,const EcalRecHitCollection& ecalHitsEE,const CaloTopology& topo)
{
  bool debug=false;
  if(debug)std::cout<<"entering fill shower shape"<<std::endl;

  const bool isEB = clus.seed().subdetId()==EcalBarrel;
  const EcalRecHitCollection& ecalHits = isEB ? ecalHitsEB : ecalHitsEE;
  if(debug)std::cout<<"shower shape is eb"<<std::endl;

  e3x3 = EcalClusterToolsT<full5x5>::e3x3(clus,&ecalHits,&topo);
  e5x5 = EcalClusterToolsT<full5x5>::e5x5(clus,&ecalHits,&topo);
  eMax = EcalClusterToolsT<full5x5>::eMax(clus,&ecalHits);
  e2nd = EcalClusterToolsT<full5x5>::e2nd(clus,&ecalHits);
  eLeft = EcalClusterToolsT<full5x5>::eLeft(clus,&ecalHits,&topo);
  eRight = EcalClusterToolsT<full5x5>::eRight(clus,&ecalHits,&topo);
  const float eLeftRightSum  = eLeft + eRight;
  const float eLeftRightDiff  = eLeft - eRight;
  eLeftRightDiffSumRatio  = eLeftRightSum !=0.f ? eLeftRightDiff/eLeftRightSum : 0.f;
  eTop = EcalClusterToolsT<full5x5>::eTop(clus,&ecalHits,&topo);
  eBottom = EcalClusterToolsT<full5x5>::eBottom(clus,&ecalHits,&topo);
  const float eTopBottomSum  = eTop + eBottom;
  const float eTopBottomDiff  = eTop - eBottom;
  eTopBottomDiffSumRatio  = eTopBottomSum !=0.f ? eTopBottomDiff/eTopBottomSum : 0.f;
  if(debug)std::cout<<"shower shape mid"<<std::endl;

  e2x5Bottom =  EcalClusterToolsT<full5x5>::e2x5Bottom(clus,&ecalHits,&topo);
  e2x5Top =  EcalClusterToolsT<full5x5>::e2x5Top(clus,&ecalHits,&topo);
  e2x5Left =  EcalClusterToolsT<full5x5>::e2x5Left(clus,&ecalHits,&topo);
  e2x5Right =  EcalClusterToolsT<full5x5>::e2x5Right(clus,&ecalHits,&topo);
  e2x5Max = EcalClusterToolsT<full5x5>::e2x5Max(clus,&ecalHits,&topo);
  if(debug)std::cout<<"shower shape dopo cluster tools"<<std::endl;

  const auto localCovs =  EcalClusterToolsT<full5x5>::localCovariances(clus,&ecalHits,&topo);
  if(debug)std::cout<<"shower shape dopo cov"<<std::endl;

  sigmaIEtaIEta =std::sqrt(localCovs[0]);
  sigmaIEtaIPhi = std::numeric_limits<float>::max();
  sigmaIPhiIPhi = std::numeric_limits<float>::max();
  if (!edm::isNotFinite(localCovs[2])) sigmaIPhiIPhi = std::sqrt(localCovs[2]) ;
  if(debug)std::cout<<"shower shape dopo numeric limit"<<std::endl;

  const bool applySPPBug = false;
  const float seeBySpp = applySPPBug ? sigmaIEtaIEta*std::numeric_limits<float>::max() : sigmaIEtaIEta*sigmaIPhiIPhi;

  if(  seeBySpp > 0 ) {
    sigmaIEtaIPhi = localCovs[1] / seeBySpp;
  } else if ( localCovs[1] > 0 ) {
    sigmaIEtaIPhi = 1.f;
  } else {
    sigmaIEtaIPhi = -1.f;
  }
  if(debug)std::cout<<"shower shape end"<<std::endl;

}


// 



void ShowerShapeStruct::fillShowerShape(const reco::GsfElectron* ele){

  bool debug=false;
  if(debug)std::cout<<"entering fill shower shape"<<std::endl;

  if ( ele->superCluster().isNull() ) { return; }
  const reco::SuperClusterRef& sc = ele->superCluster();

  //  const bool isEB = sc->seed()->subdetId()==EcalBarrel;

  // if(debug)std::cout<<"shower shape is eb"<<std::endl;



  e3x3 = ele->full5x5_r9()/sc->rawEnergy();
  
  float maxEne1=-1;
  float maxEne2=-1;
  int i1=-1;
  //  int i2=-1;
  int clusNum=0;

  if(sc->clustersSize()>0&& sc->clustersBegin()!=sc->clustersEnd()){
    for(auto& cluster : sc->clusters()) {
      if (cluster->energy() > maxEne1){
	maxEne1=cluster->energy();
	i1=clusNum;
      }
      clusNum++;
    }
    if(debug)std::cout<<"fill_supercl mini 3 - i1="<<i1<<"clusNum="<<clusNum<<std::endl;
    
    if(sc->clustersSize()>1){
      clusNum=0;
      for(auto& cluster : sc->clusters()) {
	if (clusNum!=i1) {
	  if (cluster->energy() > maxEne2){
	    maxEne2=cluster->energy();
	    //    i2=clusNum;
	  }
	}
	clusNum++;
      }
    }
  }
  
  eMax=maxEne1;
  e2nd = maxEne2;

  eLeft = ele->full5x5_eLeft();
  eRight =ele->full5x5_eRight();
  const float eLeftRightSum  = eLeft + eRight;
  const float eLeftRightDiff  = eLeft - eRight;
  eLeftRightDiffSumRatio  = eLeftRightSum !=0.f ? eLeftRightDiff/eLeftRightSum : 0.f;
  eTop = ele->full5x5_eTop();
  eBottom = ele->full5x5_eBottom();
  const float eTopBottomSum  = eTop + eBottom;
  const float eTopBottomDiff  = eTop - eBottom;
  eTopBottomDiffSumRatio  = eTopBottomSum !=0.f ? eTopBottomDiff/eTopBottomSum : 0.f;
  if(debug)std::cout<<"shower shape mid"<<std::endl;

  e2x5Bottom =  ele->full5x5_e2x5Bottom();
  e2x5Top =  ele->full5x5_e2x5Top();
  e2x5Left =  ele->full5x5_e2x5Left();
  e2x5Right =  ele->full5x5_e2x5Right();
  e2x5Max = ele->full5x5_e2x5Max();
  if(debug)std::cout<<"shower shape dopo cluster tools"<<std::endl;


  sigmaIEtaIEta = ele->full5x5_sigmaIetaIeta();
  sigmaIEtaIPhi = ele->full5x5_showerShape().sigmaIetaIphi;
  sigmaIPhiIPhi = ele->full5x5_sigmaIphiIphi();

  if(debug)std::cout<<"shower shape end"<<std::endl;

}




// Small class to provide fillers and hide tree I/O
class RegFatNtuple {

 public:

  static constexpr size_t NHITS_MAX = 30;
  static constexpr int NCLUS_MAX = 50;
  static constexpr int NEG_INT = -999;
  static constexpr float NEG_FLOAT = -999.;
  
  RegFatNtuple() {}
  
  void reset() {
    RegFatNtuple dummy; // create a new object 
    *this = dummy; // use assignment to reset
  }

  int nrVert=0;
  float rho=0;
  float nrPUInt=0;
  float nrPUIntTrue=0;
  EvtStruct evt;
  SuperClustStruct sc;
  ShowerShapeStruct ssFull;
  ShowerShapeStruct ssFrac;
  EleStruct ele;
  ShowerShapeStruct eleSSFull;
  GenInfoStruct mc;
  PFEleStruct pfele;
  ClustStruct clus1;
  ClustStruct clus2;
  ClustStruct clus3;
  std::vector<EleEnergyStruct> eleEnergies;
  void setNrEnergies(unsigned int nrEleEnergies){
    eleEnergies.resize(nrEleEnergies);
  }
  void createBranches(TTree* tree);
  void setBranchAddresses(TTree* tree);
  void fill(const edm::Event& event,int iNrVert,float iRho,float nrPUInt,float nrTruePUInt,
	    const EcalRecHitCollection& ecalHitsEB,const EcalRecHitCollection& ecalHitsEE,
	    const CaloTopology& topo,const EcalChannelStatus& ecalChanStatus,const reco::SuperCluster* iSC,
	    const reco::GenParticle* iMC,const reco::GsfElectron* iEle,
	    const  reco::GsfTrack* gsf ,
	    const reco::GsfElectron* iPFEle, const  reco::GsfTrack* ipfGsf, const  reco::Track* ipfTrk);

  void clear(){
    nrVert=0;
    rho=0.;
    nrPUInt=0.;
    nrPUIntTrue=0.;
    evt.clear();
    sc.clear();
    ssFull.clear();
    ssFrac.clear();
    ele.clear();
    eleSSFull.clear();
    mc.clear();
    pfele.clear();
    clus1.clear();
    clus2.clear();
    clus3.clear();
    for(auto& x : eleEnergies) x.clear();
  }

  
  void link_tree( TTree* tree );

  void set_rho( float r ) { rho_ = r; }

  void fill_evt( const edm::EventID& id );
  void fill_gen( const pat::PackedGenParticleRef );   
  void fill_gen( const reco::GenParticlePtr );
  void fill_gen_default( );

  void fill_trk( const reco::TrackPtr& trk,
		 const reco::BeamSpot& spot );

  void fill_trk_dEdx( const reco::TrackPtr& trk,
                      std::vector< const edm::ValueMap< reco::DeDxData > * > & x );
  void fill_trk_dEdx_default();

  void fill_bdt( double seed_unbiased,
		 double seed_ptbiased );

  void fill_gsf( const reco::GsfTrackPtr trk,
		 const reco::BeamSpot& spot );

  void fill_ele( const reco::GsfElectronPtr ele,
		 float mva_value,
		 int mva_id,
		 float ele_conv_vtx_fit_prob,
		 const double rho );

  void fill_supercluster(const reco::GsfElectronPtr ele, noZS::EcalClusterLazyTools *ecalTools_);
  void fill_supercluster_miniAOD(const reco::GsfElectronPtr ele); 

  template < typename T> bool validPtr( edm::Ptr<T>& ptr);

 public:

  // Event
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long evt_ = 0;
  float weight_ = 1.;   
  float rho_ = RegFatNtuple::NEG_FLOAT;

  // Data sample
  int is_aod_ = -1;
  int is_mc_ = -1;

  // Labels
  bool is_e_ = false;
  bool is_e_not_matched_ = false;  
  bool is_other_ = false;             
  bool is_egamma_ = false;
  bool has_trk_  = true;
  bool has_seed_ = true;
  bool has_gsf_  = true;
  bool has_ele_  = true;

  // GEN electrons
  float gen_dR_  = RegFatNtuple::NEG_FLOAT;
  float gen_pt_  = RegFatNtuple::NEG_FLOAT;
  float gen_eta_ = RegFatNtuple::NEG_FLOAT;
  float gen_phi_ = RegFatNtuple::NEG_FLOAT;
  float gen_p_   = RegFatNtuple::NEG_FLOAT;
  int gen_charge_ = RegFatNtuple::NEG_INT;
  int gen_pdgid_  = RegFatNtuple::NEG_INT;
  int gen_mom_pdgid_ = RegFatNtuple::NEG_INT;
  int gen_gran_pdgid_ = RegFatNtuple::NEG_INT;
  int gen_tag_side_ = RegFatNtuple::NEG_INT;  
  float gen_trk_dr_  = RegFatNtuple::NEG_FLOAT;
  float gen_gsf_dr_  = RegFatNtuple::NEG_FLOAT;

  // RECO steps
  float gsf_dr_ = RegFatNtuple::NEG_FLOAT;
  float trk_dr_ = RegFatNtuple::NEG_FLOAT;

  // GSF tracks: kine
  float gsf_pt_   = RegFatNtuple::NEG_FLOAT;
  float gsf_eta_  = RegFatNtuple::NEG_FLOAT;
  float gsf_phi_  = RegFatNtuple::NEG_FLOAT;
  float gsf_p_    = RegFatNtuple::NEG_FLOAT;
  int gsf_charge_ = RegFatNtuple::NEG_INT;
  float gsf_inp_  = RegFatNtuple::NEG_FLOAT;
  float gsf_outp_ = RegFatNtuple::NEG_FLOAT;

  // GSF tracks: kine (mode)
  float gsf_mode_pt_  = RegFatNtuple::NEG_FLOAT;
  float gsf_mode_eta_ = RegFatNtuple::NEG_FLOAT;
  float gsf_mode_phi_ = RegFatNtuple::NEG_FLOAT;
  float gsf_mode_p_   = RegFatNtuple::NEG_FLOAT;

  // GSF tracks: quality
  int gsf_missing_inner_hits_ = RegFatNtuple::NEG_INT;

  // GSF tracks: displacement
  float gsf_dxy_ = RegFatNtuple::NEG_FLOAT;
  float gsf_dxy_err_ = RegFatNtuple::NEG_FLOAT;
  float gsf_dz_ = RegFatNtuple::NEG_FLOAT;
  float gsf_dz_err_ = RegFatNtuple::NEG_FLOAT;

  // GSF pos at ECAL
  float gsf_x_ = RegFatNtuple::NEG_FLOAT;     
  float gsf_y_ = RegFatNtuple::NEG_FLOAT;      
  float gsf_z_ = RegFatNtuple::NEG_FLOAT;     

  // GSF tracks: tangents
  int gsf_ntangents_ = 0; //@@ RegFatNtuple::NEG_INT;
  //float gsf_hit_dpt_[NHITS_MAX] = {0}; //@@ {RegFatNtuple::NEG_FLOAT};
  //float gsf_hit_dpt_unc_[NHITS_MAX] = {0}; //@@ {RegFatNtuple::NEG_FLOAT};
  //std::vector<float> gsf_extapolated_eta_;
  //std::vector<float> gsf_extapolated_phi_;

  // Seed BDT discriminator values at GsfTrack level
  float seed_unbiased_ = RegFatNtuple::NEG_FLOAT;
  float seed_ptbiased_ = RegFatNtuple::NEG_FLOAT;

  // KF tracks: kine
  float trk_pt_  = RegFatNtuple::NEG_FLOAT;
  float trk_eta_  = RegFatNtuple::NEG_FLOAT;
  float trk_phi_  = RegFatNtuple::NEG_FLOAT;
  float trk_p_    = RegFatNtuple::NEG_INT;
  int trk_charge_ = RegFatNtuple::NEG_INT;
  float trk_inp_  = RegFatNtuple::NEG_FLOAT;
  float trk_outp_ = RegFatNtuple::NEG_FLOAT;
  int pdg_id_     = RegFatNtuple::NEG_INT;

  // KF tracks: quality
  int trk_nhits_ = RegFatNtuple::NEG_INT;
  int trk_missing_inner_hits_ = RegFatNtuple::NEG_INT;
  float trk_chi2red_   = RegFatNtuple::NEG_FLOAT;
  int trk_high_purity_ = RegFatNtuple::NEG_INT;

  // KF tracks: displ
  float trk_dxy_     = RegFatNtuple::NEG_FLOAT;
  float trk_dxy_err_ = RegFatNtuple::NEG_FLOAT;
  float trk_dz_      = RegFatNtuple::NEG_FLOAT;
  float trk_dz_err_  = RegFatNtuple::NEG_FLOAT;

  // KF tracks: dE/dx
  float trk_dEdx1_   = RegFatNtuple::NEG_FLOAT;
  int trk_dEdx1_Nm_  = RegFatNtuple::NEG_INT;
  int trk_dEdx1_NSm_ = RegFatNtuple::NEG_INT;

  // GSF electrons: kinematics
  float ele_p_   = RegFatNtuple::NEG_FLOAT;
  float ele_pt_  = RegFatNtuple::NEG_FLOAT;
  float ele_eta_ = RegFatNtuple::NEG_FLOAT;
  float ele_phi_ = RegFatNtuple::NEG_FLOAT;
  int p4kind_    = RegFatNtuple::NEG_INT;
  float core_shFracHits_ = RegFatNtuple::NEG_FLOAT;
  float ele_p_atvtx_  = RegFatNtuple::NEG_FLOAT;
  float ele_p_atcalo_ = RegFatNtuple::NEG_FLOAT;
  int fiducial_isEB_      = RegFatNtuple::NEG_INT;
  int fiducial_isEE_      = RegFatNtuple::NEG_INT; 
  int fiducial_isEBEEGap_ = RegFatNtuple::NEG_INT;
    
  // Electron: Charge
  int chPix_ = RegFatNtuple::NEG_INT;
  int chGCP_ = RegFatNtuple::NEG_INT;
  int chGP_  = RegFatNtuple::NEG_INT;
  int chGC_  = RegFatNtuple::NEG_INT;  

  // Electrons: IDs
  float ele_mva_value_ = RegFatNtuple::NEG_FLOAT;
  int ele_mva_id_      = RegFatNtuple::NEG_INT;
  float ele_conv_vtx_fit_prob_ = RegFatNtuple::NEG_FLOAT;

  float eid_rho_    = RegFatNtuple::NEG_FLOAT;
  float eid_ele_pt_ = RegFatNtuple::NEG_FLOAT;
  float eid_sc_eta_ = RegFatNtuple::NEG_FLOAT;
  float eid_shape_full5x5_sigmaIetaIeta_ = RegFatNtuple::NEG_FLOAT;
  float eid_shape_full5x5_sigmaIphiIphi_ = RegFatNtuple::NEG_FLOAT;
  float eid_shape_full5x5_circularity_   = RegFatNtuple::NEG_FLOAT;
  float eid_shape_full5x5_r9_            = RegFatNtuple::NEG_FLOAT;
  float eid_sc_etaWidth_ = RegFatNtuple::NEG_FLOAT;
  float eid_sc_phiWidth_ = RegFatNtuple::NEG_FLOAT;
  float eid_shape_full5x5_HoverE_ = RegFatNtuple::NEG_FLOAT;
  float eid_trk_nhits_   = RegFatNtuple::NEG_FLOAT;
  float eid_trk_chi2red_ = RegFatNtuple::NEG_FLOAT;
  float eid_gsf_chi2red_ = RegFatNtuple::NEG_FLOAT;
  float eid_brem_frac_   = RegFatNtuple::NEG_FLOAT;
  float eid_gsf_nhits_   = RegFatNtuple::NEG_FLOAT;
  float eid_match_SC_EoverP_   = RegFatNtuple::NEG_FLOAT;
  float eid_match_eclu_EoverP_ = RegFatNtuple::NEG_FLOAT;
  float eid_match_SC_dEta_   = RegFatNtuple::NEG_FLOAT;
  float eid_match_SC_dPhi_   = RegFatNtuple::NEG_FLOAT;
  float eid_match_seed_dEta_ = RegFatNtuple::NEG_FLOAT;
  float eid_sc_E_  = RegFatNtuple::NEG_FLOAT;
  float eid_trk_p_ = RegFatNtuple::NEG_FLOAT;

  // Electron: firther track-Cluster matching
  float match_seed_EoverP_    = RegFatNtuple::NEG_FLOAT;
  float match_seed_EoverPout_ = RegFatNtuple::NEG_FLOAT;
  float match_seed_dPhi_      = RegFatNtuple::NEG_FLOAT;
  float match_seed_dEta_vtx_  = RegFatNtuple::NEG_FLOAT;
  float match_eclu_EoverPout_ = RegFatNtuple::NEG_FLOAT;
  float match_eclu_dEta_      = RegFatNtuple::NEG_FLOAT;
  float match_eclu_dPhi_      = RegFatNtuple::NEG_FLOAT;

  // Further full 5x5 shower shape 
  float shape_full5x5_e1x5_    = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_e2x5Max_ = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_e5x5_    = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_HoverEBc_     = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth1_   = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth2_   = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth1Bc_ = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_hcalDepth2Bc_ = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_eLeft_   = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_eRight_  = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_eTop_    = RegFatNtuple::NEG_FLOAT;
  float shape_full5x5_eBottom_ = RegFatNtuple::NEG_FLOAT;
  
  // Electron, brem fractions
  float brem_fracTrk_ = RegFatNtuple::NEG_FLOAT;
  float brem_fracSC_  = RegFatNtuple::NEG_FLOAT;
  int brem_N_         = RegFatNtuple::NEG_INT;

  // SuperClusters 
  float sc_Et_  = RegFatNtuple::NEG_FLOAT;
  int sc_Nclus_ = RegFatNtuple::NEG_INT;
  int sc_Nclus_deta01_ = RegFatNtuple::NEG_INT;
  int sc_Nclus_deta02_ = RegFatNtuple::NEG_INT;
  int sc_Nclus_deta03_ = RegFatNtuple::NEG_INT;
  bool sc_goodSeed_ = false;

  float sc_E_ps_= RegFatNtuple::NEG_FLOAT;
  float sc_E_ps1_= RegFatNtuple::NEG_FLOAT;
  float sc_E_ps2_= RegFatNtuple::NEG_FLOAT;


  // Clusters 
  float sc_cluster_et_[NCLUS_MAX]    = {RegFatNtuple::NEG_FLOAT}; 
  float sc_cluster_E_[NCLUS_MAX]     = {RegFatNtuple::NEG_FLOAT}; 
  float sc_cluster_eta_[NCLUS_MAX]   = {RegFatNtuple::NEG_FLOAT}; 
  float sc_cluster_phi_[NCLUS_MAX]   = {RegFatNtuple::NEG_FLOAT}; 
  int   sc_cluster_nxtal_[NCLUS_MAX] = {RegFatNtuple::NEG_INT}; 
  float sc_cluster_e1x3_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e1x5_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2x2_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e3x3_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e5x5_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eMax_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2nd_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Right_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Left_[NCLUS_MAX]   = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Top_[NCLUS_MAX]    = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_e2x5Bottom_[NCLUS_MAX] = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eRight_[NCLUS_MAX]  = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eLeft_[NCLUS_MAX]   = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eTop_[NCLUS_MAX]    = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eBottom_[NCLUS_MAX] = {RegFatNtuple::NEG_FLOAT};  
  float sc_cluster_eMaxOver2x2_[NCLUS_MAX] = {RegFatNtuple::NEG_FLOAT};
  float sc_cluster_eMaxOver3x3_[NCLUS_MAX] = {RegFatNtuple::NEG_FLOAT};
  float sc_cluster_eMaxOver1x3_[NCLUS_MAX] = {RegFatNtuple::NEG_FLOAT};

  // When running on miniaod, only 3 highest ET clusters only
  float sc_clus1_et_=0;
  float sc_clus1_E_=0;
  float sc_clus1_E_ov_p_=0;
  float sc_clus1_E_ov_E_=0;
  float sc_clus1_eta_=0;
  float sc_clus1_phi_=0;
  int sc_clus1_nxtal_=0;
  float sc_clus1_dphi_=0;
  float sc_clus1_deta_=0;
  float sc_clus1_ntrk_deta01_=0;
  //
  float sc_clus2_et_=0;
  float sc_clus2_E_=0;
  float sc_clus2_E_ov_p_=0;
  float sc_clus2_E_ov_E_=0;
  float sc_clus2_eta_=0;
  float sc_clus2_phi_=0;
  int sc_clus2_nxtal_=0;
  float sc_clus2_dphi_=0;
  float sc_clus2_deta_=0;
  float sc_clus2_ntrk_deta01_=0;
  //
  float sc_clus3_et_=0;
  float sc_clus3_E_=0;
  float sc_clus3_E_ov_p_=0;
  float sc_clus3_E_ov_E_=0;
  float sc_clus3_eta_=0;
  float sc_clus3_phi_=0;
  int sc_clus3_nxtal_=0;
  float sc_clus3_dphi_=0;
  float sc_clus3_deta_=0;
  float sc_clus3_ntrk_deta01_=0;

};

#endif // LowPtElectrons_LowPtElectrons_RegFatNtuple
