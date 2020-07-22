#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "LowPtElectrons/LowPtElectrons/interface/RegFatPFNtuple.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronFeatures.h"
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"


#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"



#include "TMath.h"
#include "TTree.h"




#include <iostream>

using namespace std;


///////////////////////////////////////

void RegFatPFNtuple::createBranches(TTree* tree)
{
  tree->Branch("nrVert",&nrVert,"nrVert/I");
  tree->Branch("rho",&rho,"rho/F");
  tree->Branch("nrPUInt",&nrPUInt,"nrPUInt/F");
  tree->Branch("nrPUIntTrue",&nrPUIntTrue,"nrPUIntTrue/F");
  tree->Branch("evt",&evt,evt.contents().c_str());
  tree->Branch("sc",&sc,sc.contents().c_str());
  tree->Branch("ssFull",&ssFull,ssFull.contents().c_str());
  tree->Branch("ssFrac",&ssFrac,ssFrac.contents().c_str());
  tree->Branch("ele",&ele,ele.contents().c_str());
  tree->Branch("eleSSFull",&eleSSFull,eleSSFull.contents().c_str());
  tree->Branch("mc",&mc,mc.contents().c_str());
  tree->Branch("pfele",&pfele,pfele.contents().c_str());
  tree->Branch("clus1",&clus1,clus1.contents().c_str());
  tree->Branch("clus2",&clus2,clus2.contents().c_str());
  tree->Branch("clus3",&clus3,clus3.contents().c_str());
  for(size_t eleNr=0;eleNr<eleEnergies.size();eleNr++){
    std::string name = "eleAltEnergy"+std::to_string(eleNr+1);
    tree->Branch(name.c_str(),&eleEnergies[eleNr],eleEnergies[eleNr].contents().c_str());
  }
}

void RegFatPFNtuple::setBranchAddresses(TTree* tree)
{
  tree->SetBranchAddress("nrVert",&nrVert);
  tree->SetBranchAddress("rho",&rho);
  tree->SetBranchAddress("nrPUInt",&nrPUInt);
  tree->SetBranchAddress("nrPUIntTrue",&nrPUIntTrue);
  tree->SetBranchAddress("evt",&evt);
  tree->SetBranchAddress("sc",&sc);
  tree->SetBranchAddress("ssFull",&ssFull);
  tree->SetBranchAddress("ssFrac",&ssFrac);
  tree->SetBranchAddress("ele",&ele);
  tree->SetBranchAddress("eleSSFull",&eleSSFull);
  tree->SetBranchAddress("mc",&mc);
  tree->SetBranchAddress("clus1",&clus1);
  tree->SetBranchAddress("clus2",&clus2);
  tree->SetBranchAddress("clus3",&clus3);
  for(size_t eleNr=0;eleNr<eleEnergies.size();eleNr++){
    std::string name = "eleAltEnergy"+std::to_string(eleNr+1);
    tree->SetBranchAddress(name.c_str(),&eleEnergies[eleNr]);
  }
  
}

void EvtStructPF::fill(const edm::Event& event)
{
  runnr = event.id().run();
  lumiSec = event.luminosityBlock();
  eventnr = event.id().event();
}

void GenInfoStructPF::fill(const reco::GenParticle& genPart,float iDR)
{
  energy = genPart.energy();
  pt = genPart.pt();
  eta = genPart.eta();
  phi = genPart.phi();
  pdgId = genPart.pdgId();
  status = genPart.status();
  dR = iDR;
}

void RegFatPFNtuple::fill( const edm::Event& event,int iNrVert,float iRho,float iNrPUInt,float iNrPUIntTrue,
	    const EcalRecHitCollection& ecalHitsEB,const EcalRecHitCollection& ecalHitsEE,
	    const CaloTopology& topo,const EcalChannelStatus& ecalChanStatus,const reco::SuperCluster* iSC,
			   const reco::GenParticle* iMC,const reco::GsfElectron* iEle, const  reco::GsfTrack* gsf)
{
  if(debug)std::cout<<"RegFatPFNtuple::fill calling clear()"<<std::endl;
  clear();
  if(debug)std::cout<<"RegFatPFNtuple::fill done clear()"<<std::endl;

  nrVert = iNrVert;
  rho = iRho;
  nrPUInt = iNrPUInt;
  nrPUIntTrue = iNrPUIntTrue;
  evt.fill(event);
  if(debug)std::cout<<"RegFatPFNtuple::fill done fill event"<<std::endl;
  if(debug)std::cout<<"RegFatPFNtuple::fill super-cluster energy is "<<iSC->energy()<<std::endl;
  if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC "<<iSC->energy()<<std::endl;

  if(iSC){
    if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -1 "<<iSC->energy()<<std::endl;
    sc.fill(*iSC,ecalChanStatus);
    if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -2 "<<iSC->energy()<<std::endl;
    ssFull.fillShowerShape(iEle);
    // if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -3 "<<iSC->energy()<<std::endl;
    // ssFrac.fillShowerShape<false>(iEle,ecalHitsEB,ecalHitsEE);
    if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -4"<<iSC->energy()<<std::endl;
    auto fillClus = [&iSC](ClustStructPF& clus,size_t index){
      if(index<static_cast<size_t>(iSC->clusters().size())){
	auto iClus = iSC->clusters()[index];
	clus.fill( iClus->energy(),
		   iClus->eta()-iSC->seed()->eta(),
		   reco::deltaPhi(iClus->phi(),iSC->seed()->phi()));
      }else{
	clus.fill( 0, 0, 0 );
      }
    };  
    if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -5"<<iSC->energy()<<std::endl;
    try{ 
      fillClus(clus1,1);  
      fillClus(clus2,2);  
      fillClus(clus3,3);
    } catch (...) {
      if(debug) std::cout<<"there are no clusters"<<std::endl; 
    }  
    if(debug)std::cout<<"RegFatPFNtuple::fill going to fill SC -6"<<iSC->energy()<<std::endl;
  }
  if(debug)std::cout<<"RegFatPFNtuple::fill done fill clus"<<std::endl;

  if(iMC) mc.fill(*iMC, iSC ? std::sqrt(reco::deltaR2(iSC->eta(),iSC->phi(),iMC->eta(),iMC->phi())) : 999);
  if(iEle){
    if(debug)std::cout<<"RegFatPFNtuple::fill filling ele and gsf trk"<<std::endl;
    ele.fill(*iEle,*gsf);
    if(debug)std::cout<<"RegFatPFNtuple::fill filling ele and gsf trk done... now shower shapes " <<std::endl;
    eleSSFull.fill(iEle->full5x5_showerShape(),*iEle);
    if(debug)std::cout<<"RegFatPFNtuple::fill filling shower shapes done "<<std::endl;
  }
  if(debug)std::cout<<"RegFatPFNtuple::fill done all"<<std::endl;

  //  if(altEles.size() != eleEnergies.size()){
  //  throw cms::Exception("LogicError") <<" alt electrons && eleEnergies are not equal in size "<<altEles.size()<<" "<<eleEnergies.size();
  // }

  //for(size_t eleNr=0;eleNr<altEles.size();eleNr++){
  //  if(altEles[eleNr]) eleEnergies[eleNr].fill(*altEles[eleNr]);
  // }

}

template<typename DetIdType>
int getDeadCrysCode(const DetIdType& id,const EcalChannelStatus& ecalChanStatus)
{
  //constexpr int deadChanMask = EcalChannelStatusCode::kDeadFE | EcalChannelStatusCode::kDeadVFE | EcalChannelStatusCode::kNoDataNoTP | EcalChannelStatusCode::kNonRespondingIsolated;
  constexpr int deadChanMask = EcalChannelStatusCode::kNoDataNoTP;
  
  int crysCode =0 ;

  //meh I could derive this from iEtaOrXNr/iPhiOrYNr indices but sometimes simpliest is best
  int bitNr=0;
  for(int iEtaOrXNr=-1;iEtaOrXNr<=1;iEtaOrXNr++){
    for(int iPhiOrYNr=-1;iPhiOrYNr<=1;iPhiOrYNr++){
      if(iEtaOrXNr==0 && iPhiOrYNr==0) continue;
      DetIdType currId = id.offsetBy(iEtaOrXNr,iPhiOrYNr);
      if(currId.rawId()!=0){
	int statusCode = ecalChanStatus[currId.rawId()].getEncodedStatusCode();
	int bit = 0x1 << bitNr;
	if((statusCode&deadChanMask)!=0) crysCode |=bit;
      }
      bitNr++;
    }
  }
  return crysCode;
}

void SuperClustStructPF::fill(const reco::SuperCluster& sc,const EcalChannelStatus& ecalChanStatus)
{
  if(debug)std::cout<<"SuperClustStructPF::fill >> SC raw ene="<<sc.rawEnergy()<<std::endl;
  auto& seedClus = *sc.seed(); 
  isEB = seedClus.seed().subdetId()==EcalBarrel;
  
  rawEnergy = sc.rawEnergy();
  rawESEnergy = sc.preshowerEnergy();
  etaWidth = sc.etaWidth();
  phiWidth = sc.phiWidth();
  seedClusEnergy = seedClus.energy();
  corrEnergy = sc.energy();
  scEta = sc.eta();
  scPhi = sc.phi();
  auto divideWithZeroCheck= [](float numer,float denom)->float{
    return denom!=0 ? numer/denom : 0.;
  };
  scSinTheta = divideWithZeroCheck(sc.position().rho(),sc.position().r());
  seedEta = seedClus.eta();
  seedPhi = seedClus.phi();
  seedSinTheta = divideWithZeroCheck(seedClus.position().rho(),seedClus.position().r());
  dPhiSeedSC = reco::deltaPhi(seedClus.phi(),sc.position().Phi()); //needs this way due to rounding errors
  dEtaSeedSC = seedClus.eta() - sc.position().Eta(); //needs this way due to rounding errors
  numberOfClusters = sc.clusters().size();
  numberOfSubClusters = std::max(0,static_cast<int>(sc.clusters().size())-1);

  if(debug)std::cout<<"SuperClustStructPF::fill >> SC basic var done SC raw ene="<<sc.rawEnergy()<<std::endl;

  if(isEB){
    EBDetId ebDetId(seedClus.seed());
    iEtaOrX = ebDetId.ieta();
    iPhiOrY = ebDetId.iphi();
    const int iEtaCorr = ebDetId.ieta() - (ebDetId.ieta() > 0 ? +1 : -1);
    const int iEtaCorr26 = ebDetId.ieta() - (ebDetId.ieta() > 0 ? +26 : -26);
    iEtaMod5 = iEtaCorr%5;
    iEtaMod20 = std::abs(ebDetId.ieta())<=25 ? iEtaCorr%20 : iEtaCorr26%20;
    iPhiMod2 = (ebDetId.iphi()-1)%2;
    iPhiMod20 = (ebDetId.iphi()-1)%20;
    auto gapCode=[](int iEtaAbs){
      if(iEtaAbs==25 || iEtaAbs==45 || iEtaAbs==65 || iEtaAbs==85) return -1;//before gap
      else if(iEtaAbs==1 || iEtaAbs==26 || iEtaAbs==46 || iEtaAbs==66) return 1;//after gap
      else return 0; //no gap
    };
    etaGapCode = gapCode(ebDetId.ietaAbs());
    phiGapCode = ebDetId.iphi()%20 == 0 ? -1 : ebDetId.iphi()%20 ==1 ? 1 : 0;
    nearbyChanStatus = getDeadCrysCode(ebDetId,ecalChanStatus);
    
  }else{
    EEDetId eeDetId(seedClus.seed());
    iEtaOrX = eeDetId.ix();
    iPhiOrY = eeDetId.iy();
    auto gapCode=[](int iAbs){
      if(iAbs==45 || iAbs==50 || iAbs==55) return -1;//before gap
      else if(iAbs==46 || iAbs==51 || iAbs==56) return 1;//before gap
      else return 0; //no gap
    };
    etaGapCode = gapCode(eeDetId.ix());
    phiGapCode = gapCode(eeDetId.iy());
    nearbyChanStatus = getDeadCrysCode(eeDetId,ecalChanStatus);
  }
  clusterMaxDR     = 999.;
  clusterMaxDRDPhi = 999.;
  clusterMaxDRDEta = 999.;
  clusterMaxDRRawEnergy = 0.;

  if(debug)std::cout<<"SuperClustStructPF::fill >> SC basic var done rec hits SC raw ene="<<sc.rawEnergy()<<std::endl;
  if(sc.clusters().size()>0 ){
    try{
      float maxDR2 = 0;
      for(auto& clus : sc.clusters()){
	if(clus == sc.seed()) continue;
	float dR2 = reco::deltaR2(seedEta,seedPhi,clus->eta(),clus->phi());
	if(dR2 > maxDR2 ){
	  maxDR2 = dR2;
	  clusterMaxDR = std::sqrt(dR2);
	  clusterMaxDRDPhi = reco::deltaPhi(clus->phi(),seedPhi);
	  clusterMaxDRDEta = clus->eta()-seedEta;
	  clusterMaxDRRawEnergy = clus->energy();
	}
      }
    } catch (...){
      if(debug) std::cout<<"SuperClustStructPF::fill >> there were no clusters attached to this SC "<<sc.rawEnergy()<<std::endl;
    }
 
  }
  if(debug)std::cout<<"SuperClustStructPF::fill >> SC basic all done SC raw ene="<<sc.rawEnergy()<<std::endl;


  /*  if(altSC){
    corrEnergyAlt = altSC->energy();
    rawEnergyAlt = altSC->rawEnergy();
    nrClusAlt = altSC->clusters().size();
    } */

}


void EleStructPF::fill(const reco::GsfElectron& ele, const  reco::GsfTrack& gsf )
{
  et = ele.et();
  energy = ele.energy();
  energyErr = ele.p4Error(reco::GsfElectron::P4_COMBINATION);
  ecalEnergy = ele.ecalEnergy();
  ecalEnergyErr = ele.ecalEnergyError();
  eta = ele.eta();
  phi = ele.phi();
 
  trkEtaMode = gsf.etaMode();
  trkPhiMode = gsf.phiMode();
  trkPMode = gsf.pMode();
  trkPModeErr = std::abs(gsf.qoverpModeError())*gsf.pMode()*gsf.pMode();
  trkPInn = gsf.p();
  trkPtInn = gsf.pt();
  trkPVtx = std::sqrt(ele.trackMomentumAtVtx().Mag2());
  trkPOut = std::sqrt(ele.trackMomentumOut().Mag2());
  trkChi2 = gsf.chi2();
  trkNDof = gsf.ndof();
  fbrem = ele.fbrem();
  corrMean = 1.;
  corrSigma = 0.;
  hademTow = ele.hcalOverEcalBc();
  hademCone = ele.hcalOverEcal();
  ecalDrivenSeed = ele.ecalDrivenSeed();
  nrSatCrys = ele.nSaturatedXtals();
  scRawEnergy = ele.superCluster()->rawEnergy();
  scRawESEnergy = ele.superCluster()->preshowerEnergy();
}

void PFEleStructPF::fill(const reco::GsfElectron& ele, const  reco::GsfTrack& gsf,const  reco::Track& trk )
{


    et = ele.et();
    energy = ele.energy();
    ecalEnergy = ele.ecalEnergy();
    scRawEnergy = ele.superCluster()->rawEnergy();
    scRawESEnergy = ele.superCluster()->preshowerEnergy();
    trkPMode = gsf.pMode();
    trkPtMode = gsf.ptMode();
    trkP = trk.p();

}


void ShowerShapeStructPF::fill(const reco::GsfElectron::ShowerShape& eleSS,const reco::GsfElectron& ele)
{
  e3x3 = eleSS.r9*ele.superCluster()->rawEnergy();
  e5x5 = eleSS.e5x5;
  eMax = eleSS.eMax;
  e2nd = eleSS.e2nd;
  eLeft = eleSS.eLeft;
  eRight = eleSS.eRight;
  const float eLeftRightSum  = eLeft + eRight;
  const float eLeftRightDiff  = eLeft - eRight;
  eLeftRightDiffSumRatio  = eLeftRightSum !=0.f ? eLeftRightDiff/eLeftRightSum : 0.f;
  eTop = eleSS.eTop;
  eBottom = eleSS.eBottom;
  const float eTopBottomSum  = eTop + eBottom;
  const float eTopBottomDiff  = eTop - eBottom;
  eTopBottomDiffSumRatio  = eTopBottomSum !=0.f ? eTopBottomDiff/eTopBottomSum : 0.f;

  e2x5Bottom = eleSS.e2x5Bottom;
  e2x5Top = eleSS.e2x5Top;
  e2x5Left = eleSS.e2x5Left;
  e2x5Right = eleSS.e2x5Right;
  sigmaIEtaIEta = eleSS.sigmaIetaIeta;
  sigmaIEtaIPhi = eleSS.sigmaIetaIphi;
  sigmaIPhiIPhi = eleSS.sigmaIphiIphi;
}




/////////////////////////////////////////////////////////////////////////////////
//
void RegFatPFNtuple::link_tree( TTree *tree ) {

  std::cout<<"I am running RegFatPFNtuple::link_tree"<<std::endl; 
  // general

  createBranches(tree);

}

/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_evt( const edm::EventID& id ) {  
  run_  = id.run();
  lumi_ = id.luminosityBlock();
  evt_  = id.event();
}

/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_gen( const reco::GenParticlePtr genp ) {

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = genp->pdgId();
  if ( genp->mother(0) )
    gen_mom_pdgid_ = (*genp->mother(0)).pdgId();
  else 
    gen_mom_pdgid_ = 0;
  if ( (genp->mother(0))->mother(0) )
    gen_gran_pdgid_ = ( *(*genp->mother(0)).mother(0)).pdgId();
  else 
    gen_gran_pdgid_ = 0;
}

void RegFatPFNtuple::fill_gen( const pat::PackedGenParticleRef genp ) {  

  gen_pt_  = genp->pt();
  gen_eta_ = genp->eta();
  gen_phi_ = genp->phi();
  gen_p_ = genp->p();
  gen_charge_ = genp->charge();
  gen_pdgid_ = genp->pdgId();
  if ( genp->mother(0) )
    gen_mom_pdgid_ = (*genp->mother(0)).pdgId();
  else 
    gen_mom_pdgid_ = 0;
  if ( (genp->mother(0))->mother(0) )
    gen_gran_pdgid_ = ( *(*genp->mother(0)).mother(0)).pdgId();
  else 
    gen_gran_pdgid_ = 0;
}

void RegFatPFNtuple::fill_gen_default() {

  gen_pt_  = -999.;
  gen_eta_ = -999.;
  gen_phi_ = -999.;
  gen_p_   = -999.;
  gen_charge_ = -999; 
  gen_pdgid_  = -999;
  gen_mom_pdgid_  = -999; 
  gen_gran_pdgid_ = -999;
}

/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_trk( const reco::TrackPtr& trk,
			 const reco::BeamSpot& spot ) {       
  
  if ( trk.isNonnull() ) {
    // kine
    trk_pt_ = trk->pt();
    trk_eta_ = trk->eta();
    trk_phi_ = trk->phi();
    trk_p_ = trk->p();
    trk_charge_ = trk->charge();
    if ( trk->extra().isAvailable() && trk->extra().isNonnull() ) {
      trk_inp_ = sqrt( trk->innerMomentum().mag2() );
      trk_outp_ = sqrt( trk->outerMomentum().mag2() );
    }
    // quality
    trk_nhits_ = trk->found();
    trk_missing_inner_hits_ = trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    trk_chi2red_ = trk->normalizedChi2();
    // displ
    trk_dxy_ = trk->dxy(spot);
    trk_dxy_err_ = trk->dxyError();
    trk_dz_ = trk->dz(spot.position());
    trk_dz_err_ = trk->dzError();
    trk_high_purity_=trk->quality( reco::TrackBase::qualityByName("highPurity") ) ;

  }
  
}

void RegFatPFNtuple::fill_trk_dEdx( const reco::TrackPtr& trk,
				  std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx ) {

  if ( trk.isNonnull() ) {
    const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[0]);
    const reco::DeDxData& dedx = dEdxTrack[trk];
    trk_dEdx1_=dedx.dEdx();
    trk_dEdx1_Nm_=dedx.numberOfMeasurements();
    trk_dEdx1_NSm_=dedx.numberOfSaturatedMeasurements();
  }

}

void RegFatPFNtuple::fill_trk_dEdx_default( ) {

  trk_dEdx1_     = -999.;
  trk_dEdx1_Nm_  = -999;
  trk_dEdx1_NSm_ = -999;
}

/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_bdt( double seed_unbiased, 
			 double seed_ptbiased ) {          
  seed_unbiased_ = seed_unbiased;
  seed_ptbiased_ = seed_ptbiased;
}

/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_gsf( const reco::GsfTrackPtr gsf, 
			 const reco::BeamSpot& spot ) {        

  if ( gsf.isNull() ) {
    //@@ Shouldn't happen, but do we just take dummy values...? 
  } else {

    // Kinematics
    gsf_pt_ = gsf->pt();
    gsf_eta_ = gsf->eta();
    gsf_phi_ = gsf->phi();
    gsf_p_ = gsf->p();
    gsf_charge_ = gsf->charge();
    if ( gsf->extra().isAvailable() && gsf->extra().isNonnull() ) {
      gsf_inp_ = sqrt(gsf->innerMomentum().mag2());
      gsf_outp_ = sqrt(gsf->outerMomentum().mag2());
    }
    
    // Kinematics (MODE)
    gsf_mode_pt_ = gsf->ptMode();
    gsf_mode_eta_ = gsf->etaMode();
    gsf_mode_phi_ = gsf->phiMode();
    gsf_mode_p_ = gsf->pMode();

    // Quality
    gsf_missing_inner_hits_ = gsf->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

    // Displacement
    gsf_dxy_ = gsf->dxy(spot);
    gsf_dxy_err_ = gsf->dxyError();
    gsf_dz_ = gsf->dz(spot.position());
    gsf_dz_err_ = gsf->dzError();

    // Tangents (requires TrackExtra)
    //    const auto& extra = gsf->gsfExtra(); //@@ Collection does not exist?! 
    //    if ( extra.isNonnull() ) {
    //      gsf_ntangents_ = (extra->tangentsSize() > NHITS_MAX) ? NHITS_MAX : extra->tangentsSize();
    //      for (int idx = 0; idx < gsf_ntangents_; idx++ ) {
    //	gsf_hit_dpt_[idx] = extra->tangents().at(idx).deltaP().value();
    //	gsf_hit_dpt_unc_[idx] = extra->tangents().at(idx).deltaP().error();
    //      }
    //    }
    
  } 
}



/////////////////////////////////////////////////////////////////////////////////
void RegFatPFNtuple::fill_ele( const reco::GsfElectronPtr ele,
			 float mva_value,
			 int mva_id,
			 float ele_conv_vtx_fit_prob,
			 const double rho ) {       


  if ( ele.isNonnull() ) {  

    // Kinematics 
    ele_p_   = ele->p();
    ele_pt_  = ele->pt();
    ele_eta_ = ele->eta();
    ele_phi_ = ele->phi();
    
    // Momentum
    p4kind_  = ele->candidateP4Kind();
    core_shFracHits_ = ele->shFracInnerHits();
    ele_p_atvtx_  = sqrt(ele->trackMomentumAtVtx().mag2());
    ele_p_atcalo_ = sqrt(ele->trackMomentumAtCalo().mag2());
    // Fiducial flags 
    fiducial_isEB_ = ele->isEB();
    fiducial_isEE_ = ele->isEE();
    fiducial_isEBEEGap_ = ele->isEBEEGap();

    // Charge 
    chPix_ = ele->scPixCharge();
    chGCP_ = ele->isGsfCtfScPixChargeConsistent();
    chGP_  = ele->isGsfScPixChargeConsistent();
    chGC_  = ele->isGsfCtfChargeConsistent();    

    // MVA IDs: only filled if 'ValueMap->size() == electrons->size()' in IDFeatures::analyze()
    if ( mva_value > -666. ) { ele_mva_value_ = mva_value; }
    if ( mva_id > -666 ) { ele_mva_id_ = mva_id; }
    if ( ele_conv_vtx_fit_prob > -666. ) { ele_conv_vtx_fit_prob_ = ele_conv_vtx_fit_prob; }
    
    // ElectronID variables
    std::vector<float> vfeatures;
    //@@ ORDER IS IMPORTANT!
    size_t idx = 0;
    eid_rho_ = vfeatures[idx++];
    eid_ele_pt_ = vfeatures[idx++];
    eid_sc_eta_ = vfeatures[idx++];
    eid_shape_full5x5_sigmaIetaIeta_ = vfeatures[idx++];
    eid_shape_full5x5_sigmaIphiIphi_ = vfeatures[idx++];
    eid_shape_full5x5_circularity_ = vfeatures[idx++];
    eid_shape_full5x5_r9_ = vfeatures[idx++];
    eid_sc_etaWidth_ = vfeatures[idx++];
    eid_sc_phiWidth_ = vfeatures[idx++];
    eid_shape_full5x5_HoverE_ = vfeatures[idx++];
    eid_trk_nhits_ = vfeatures[idx++];
    eid_trk_chi2red_ = vfeatures[idx++];
    eid_gsf_chi2red_ = vfeatures[idx++];
    eid_brem_frac_ = vfeatures[idx++];
    eid_gsf_nhits_ = vfeatures[idx++];
    eid_match_SC_EoverP_ = vfeatures[idx++];
    eid_match_eclu_EoverP_ = vfeatures[idx++];
    eid_match_SC_dEta_ = vfeatures[idx++];
    eid_match_SC_dPhi_ = vfeatures[idx++];
    eid_match_seed_dEta_ = vfeatures[idx++];
    eid_sc_E_ = vfeatures[idx++];
    eid_trk_p_ = vfeatures[idx++];

    // Further track-Cluster matching 
    match_seed_EoverP_    = ele->eSeedClusterOverP();
    match_seed_EoverPout_ = ele->eSeedClusterOverPout();
    match_seed_dPhi_      = ele->deltaPhiSeedClusterTrackAtCalo();
    match_seed_dEta_vtx_  = ele->deltaEtaSeedClusterTrackAtVtx();
    match_eclu_EoverPout_ = ele->eEleClusterOverPout();
    match_eclu_dEta_ = ele->deltaEtaEleClusterTrackAtCalo();
    match_eclu_dPhi_ = ele->deltaPhiEleClusterTrackAtCalo();

    // Further full 5x5 shower shape 
    shape_full5x5_e1x5_     = ele->full5x5_e1x5();
    shape_full5x5_e2x5Max_  = ele->full5x5_e2x5Max();
    shape_full5x5_e5x5_     = ele->full5x5_e5x5();
    shape_full5x5_HoverEBc_ = ele->full5x5_hcalOverEcalBc();
    shape_full5x5_hcalDepth1_    = ele->full5x5_hcalDepth1OverEcal();
    shape_full5x5_hcalDepth2_    = ele->full5x5_hcalDepth2OverEcal();
    shape_full5x5_hcalDepth1Bc_  = ele->full5x5_hcalDepth1OverEcalBc();
    shape_full5x5_hcalDepth2Bc_  = ele->full5x5_hcalDepth2OverEcalBc();
    shape_full5x5_eLeft_   = ele->full5x5_eLeft();
    shape_full5x5_eRight_  = ele->full5x5_eRight();
    shape_full5x5_eTop_    = ele->full5x5_eTop();
    shape_full5x5_eBottom_ = ele->full5x5_eBottom();

    // Brem fractions and classification 
    brem_fracTrk_ = ele->trackFbrem();
    brem_fracSC_  = ele->superClusterFbrem();
    brem_N_ = ele->numberOfBrems();



  }    
}

void RegFatPFNtuple::fill_supercluster(const reco::GsfElectronPtr ele, noZS::EcalClusterLazyTools *ecalTools_ ) {

  if ( ele.isNull() ) { return; }

  if ( ele->superCluster().isNull() ) { return; }
  const reco::SuperClusterRef& sc = ele->superCluster();

  int clusNum=0;
  float maxEne=-1;
  for(auto& cluster : sc->clusters()) {
    if (clusNum<NCLUS_MAX) {
      float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) + pow(cluster->z(),2) );
      sc_cluster_et_[clusNum]    = clusterEt;
      sc_cluster_E_[clusNum]     = cluster->energy();
      sc_cluster_eta_[clusNum]   = cluster->eta();
      sc_cluster_phi_[clusNum]   = cluster->phi();
      sc_cluster_nxtal_[clusNum] = cluster->size(); 
      sc_cluster_e1x3_[clusNum] = ecalTools_->e1x3(*cluster);
      sc_cluster_e1x5_[clusNum] = ecalTools_->e1x5(*cluster);
      sc_cluster_e2x2_[clusNum] = ecalTools_->e2x2(*cluster);
      sc_cluster_e3x3_[clusNum] = ecalTools_->e3x3(*cluster);
      sc_cluster_e5x5_[clusNum] = ecalTools_->e5x5(*cluster);
      sc_cluster_eMax_[clusNum] = ecalTools_->eMax(*cluster);
      sc_cluster_e2nd_[clusNum] = ecalTools_->e2nd(*cluster);
      sc_cluster_e2x5Right_[clusNum]  = ecalTools_->e2x5Right(*cluster);
      sc_cluster_e2x5Left_[clusNum]   = ecalTools_->e2x5Left(*cluster);
      sc_cluster_e2x5Top_[clusNum]    = ecalTools_->e2x5Top(*cluster);
      sc_cluster_e2x5Bottom_[clusNum] = ecalTools_->e2x5Bottom(*cluster);
      sc_cluster_eRight_[clusNum]  = ecalTools_->eRight(*cluster);
      sc_cluster_eLeft_[clusNum]   = ecalTools_->eLeft(*cluster);
      sc_cluster_eTop_[clusNum]    = ecalTools_->eTop(*cluster);
      sc_cluster_eBottom_[clusNum] = ecalTools_->eBottom(*cluster);
      sc_cluster_eMaxOver2x2_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e2x2(*cluster);
      sc_cluster_eMaxOver3x3_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e3x3(*cluster);
      sc_cluster_eMaxOver1x3_[clusNum] = ecalTools_->eMax(*cluster)/ecalTools_->e1x3(*cluster);
      if (cluster->energy() > maxEne) maxEne=cluster->energy();
      clusNum++;
    }
  }
  sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );
  sc_Nclus_ = sc->clustersSize();

  float seedEne = sc->seed()->energy();
  if ( fabs(seedEne-maxEne)<0.001 ) sc_goodSeed_ = true;
}

// FC new method 
void RegFatPFNtuple::fill_supercluster_miniAOD(const reco::GsfElectronPtr ele  ) {
  bool debug=false; 
  if(debug)std::cout<<"fill_supercl mini"<<std::endl; 
  if ( ele.isNull() ) { return; }

  if ( ele->superCluster().isNull() ) { return; }
  const reco::SuperClusterRef& sc = ele->superCluster();

  // Propagate 'electron' to ECAL surface
  double mass_=0.000511*0.000511; // ele mass 

  reco::GsfTrackPtr kfTrackRef = edm::refToPtr(ele->gsfTrack());
  if (! validPtr(kfTrackRef) ) return; 
  if(debug)std::cout<<"fill_supercl mini 1"<<std::endl; 

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

  if(debug)std::cout<<"fill_supercl mini 2"<<std::endl; 

  //  std::cout<<"ele part: pt/eta/phi"<<kfTrackRef->pt()<<"/"<<kfTrackRef->eta()<<"/"<<kfTrackRef->phi()<<" vtx"<<vx<<"/"<<vy<<"/"<<vz<<" energy="<<energy<<" charge="<< kfTrackRef->charge()<<std::endl;  

  //  math::XYZVector field(field_->inTesla(GlobalPoint(0, 0, 0)));
  float field_z=3.8;

  BaseParticlePropagator mypart(RawParticle(mom,pos,kfTrackRef->charge()), 0, 0, field_z);
  mypart.propagateToEcalEntrance(true); // true only first half loop , false more than one loop
  bool reach_ECAL=mypart.getSuccess(); // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps 

  // ECAL entry point for track
  // GlobalPoint ecal_pos(mypart.propagated().vertex().x(), mypart.propagated().vertex().y(), mypart.propagated().vertex().z());
  GlobalPoint ecal_pos(mypart.particle().vertex().x(), mypart.particle().vertex().y(), mypart.particle().vertex().z());
  // Preshower limit
  //  bool below_ps = pow(ecal_pos.z(), 2.) > boundary_ * ecal_pos.perp2();
  // Iterate through ECAL clusters

  gsf_x_=mypart.particle().vertex().x();
  gsf_y_=mypart.particle().vertex().y();
  gsf_z_=mypart.particle().vertex().z();

  if(debug)std::cout<<"fill_supercl mini 3"<<std::endl; 
  if(debug)std::cout<<"fill_supercl mini 3 - size="<<sc->clustersSize()<<std::endl; 


  int clusNum=0;

  float maxEne1=-1;
  float maxEne2=-1;
  float maxEne3=-1;
  int i1=-1;
  int i2=-1;
  int i3=-1;

  try{
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
	      i2=clusNum;
	    }
	  }
	  clusNum++;
	}
      }


      if(debug)std::cout<<"fill_supercl mini 3 - i2="<<i2<<"clusNum="<<clusNum<<std::endl; 
    
      if(sc->clustersSize()>2){
	clusNum=0;
	for(auto& cluster : sc->clusters()) {
	  if (clusNum!=i1 && clusNum!=i2) {
	    if (cluster->energy() > maxEne3){
	      maxEne3=cluster->energy();
	      i3=clusNum;
	    }
	  }
	  clusNum++;
	}
      }
      if(debug)std::cout<<"fill_supercl mini 3 - i3="<<i3<<"clusNum="<<clusNum<<std::endl; 
    }
  }catch(...){
    std::cout<<"exception caught clusNum="<<clusNum<<" clus size"<<sc->clustersSize()<<" energy="<< sc->energy()<<std::endl;
  }
  if(debug)std::cout<<"fill_supercl mini 3 end clusnum="<<clusNum<<std::endl; 
  
  sc_clus1_et_    = 0;
  sc_clus1_E_     = 0;
  sc_clus1_E_ov_p_     = 0;
  sc_clus1_E_ov_E_     = 0;
  sc_clus1_eta_   = 0;
  sc_clus1_phi_   = 0;
  sc_clus1_nxtal_ = 0;
  sc_clus1_deta_ = -100;
  sc_clus1_dphi_ = -100;
  //
  sc_clus2_et_    = 0;
  sc_clus2_E_     = 0;
  sc_clus2_E_ov_p_     = 0;
  sc_clus2_E_ov_E_     = 0;
  sc_clus2_eta_   = 0;
  sc_clus2_phi_   = 0;
  sc_clus2_nxtal_ = 0;
  sc_clus2_deta_ = -100;
  sc_clus2_dphi_ = -100;
  //
  sc_clus3_et_    = 0;
  sc_clus3_E_     = 0;
  sc_clus3_E_ov_p_     = 0;
  sc_clus3_E_ov_E_     = 0;
  sc_clus3_eta_   = 0;
  sc_clus3_phi_   = 0;
  sc_clus3_nxtal_ = 0;
  sc_clus3_deta_ = -100;
  sc_clus3_dphi_ = -100;

  if(debug)std::cout<<"fill_supercl mini 4"<<std::endl; 


  // trovati i 3 cluster piu` energetici 
  // riempio i primi 3 cluster in E 
  // i1 i2 i3 
  clusNum=0;

  try{
  
    if(sc->clustersSize()>0&& sc->clustersBegin()!=sc->clustersEnd()){
      
      for(auto& cluster : sc->clusters()) {
	

    // Correct ecal_pos for shower depth
    //    double shower_depth = reco::PFCluster::getDepthCorrection(cluster->correctedEnergy(), below_ps, false);
    //GlobalPoint showerPos = ecal_pos + GlobalVector(mypart.propagated().momentum().x(),
    //						  mypart.propagated().momentum().y(),
    //						  mypart.propagated().momentum().z()).unit() * shower_depth;


    double pi_=3.1415926535;
    // Determine dR squared
    //    float dr2 = reco::deltaR2(cluRef->positionREP(), showerPos);
    // determine deta and dphi 
    float deta = std::abs(ecal_pos.eta()-cluster->eta()) ;
    float dphi = std::abs(ecal_pos.phi()-cluster->phi());
    if (dphi > pi_)  dphi -= 2 * pi_;
    if(ecal_pos.phi()-cluster->phi()<0) dphi=-dphi;
    if(ecal_pos.eta()-cluster->eta()<0) deta=-deta;

    if (clusNum==i1) {
      float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
      sc_clus1_et_    = clusterEt;
      sc_clus1_E_     = cluster->energy();
      if(gsf_mode_p_>0)sc_clus1_E_ov_p_     = cluster->energy()/gsf_mode_p_;
      if(eid_sc_E_>0)sc_clus1_E_ov_E_     = cluster->energy()/eid_sc_E_;
      sc_clus1_eta_   = cluster->eta();
      sc_clus1_phi_   = cluster->phi();
      sc_clus1_nxtal_ =(int) cluster->size();
      if(reach_ECAL<=0){
	//	std::cout<<"did not reach ECAL"<<reach_ECAL<<std::endl; 
      } else {
	sc_clus1_deta_ = deta;
	sc_clus1_dphi_ = dphi;
	//	std::cout<<"cluster 1 ene="<<cluster->energy()<<" eta="<<cluster->eta()<<" extra_eta="<<ecal_pos.eta()<<" phi="<<cluster->phi()<<" extra_phi="<<ecal_pos.phi()<<std::endl;
      }

    } else if(clusNum==i2){
      float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
      sc_clus2_et_    = clusterEt;
      sc_clus2_E_     = cluster->energy();
      if(gsf_mode_p_>0)sc_clus2_E_ov_p_     = cluster->energy()/gsf_mode_p_;
      if(eid_sc_E_>0)sc_clus2_E_ov_E_     = cluster->energy()/eid_sc_E_;
      sc_clus2_eta_   = cluster->eta();
      sc_clus2_phi_   = cluster->phi();
      sc_clus2_nxtal_ = (int) cluster->size();
  if(reach_ECAL<=0){
    //std::cout<<"did not reach ECAL"<<reach_ECAL<<std::endl; 
  } else {
    sc_clus2_deta_ = deta;
    sc_clus2_dphi_ = dphi;
    //std::cout<<"cluster 2 ene="<<cluster->energy()<<" eta="<<cluster->eta()<<" extra_eta="<<ecal_pos.eta()<<" phi="<<cluster->phi()<<" extra_phi="<<ecal_pos.phi()<<std::endl;
  }
    } else if(clusNum==i3){
      float clusterEt = cluster->energy() * sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) ) / sqrt( pow(cluster->x(),2) + pow(cluster->y(),2) +pow(cluster->z(),2) );
      sc_clus3_et_    = clusterEt;
      sc_clus3_E_     = cluster->energy();
      if(gsf_mode_p_>0)sc_clus3_E_ov_p_     = cluster->energy()/gsf_mode_p_;
      if(eid_sc_E_>0)sc_clus3_E_ov_E_     = cluster->energy()/eid_sc_E_;
      sc_clus3_eta_   = cluster->eta();
      sc_clus3_phi_   = cluster->phi();
      sc_clus3_nxtal_ = (int) cluster->size();
      if(reach_ECAL<=0){
	//std::cout<<"did not reach ECAL"<<reach_ECAL<<std::endl; 
      } else {
      sc_clus3_deta_ = deta;
      sc_clus3_dphi_ = dphi;
      //std::cout<<"cluster 3 ene="<<cluster->energy()<<" eta="<<cluster->eta()<<" extra_eta="<<ecal_pos.eta()<<" phi="<<cluster->phi()<<" extra_phi="<<ecal_pos.phi()<<std::endl;
      }
    }
    clusNum++;
      }
    }
  }catch(...){
    //    std::cout<<"caught an exception"<<std::endl;
  }

  sc_Et_ = sc->energy() * sqrt( pow(sc->x(),2) + pow(sc->y(),2) ) / 
    sqrt( pow(sc->x(),2) + pow(sc->y(),2) + pow(sc->z(),2) );
  sc_Nclus_ = sc->clustersSize();
  if(debug)std::cout<<"fill_supercl mini almost end"<<std::endl; 

  float seedEne = sc->seed()->energy();
  if ( fabs(seedEne-maxEne1)<0.001 ) sc_goodSeed_ = true;

  sc_E_ps_=sc->preshowerEnergy();
  sc_E_ps1_=sc->preshowerEnergyPlane1();
  sc_E_ps2_=sc->preshowerEnergyPlane2();

  sc_Nclus_deta01_=0;
  sc_Nclus_deta02_=0;
  sc_Nclus_deta03_=0;
  try{
    if(sc->clustersSize()>0&& sc->clustersBegin()!=sc->clustersEnd()){
      for(auto& cluster : sc->clusters()) {
	float deta = std::abs(ecal_pos.eta()-cluster->eta()) ;
	if(deta<0.1) sc_Nclus_deta01_=sc_Nclus_deta01_+1; 
	if(deta<0.2) sc_Nclus_deta02_=sc_Nclus_deta02_+1; 
	if(deta<0.3) sc_Nclus_deta03_=sc_Nclus_deta03_+1; 
      }
    }
  }catch(...){
    //    std::cout<<"caught an exception"<<std::endl;
  }

  if(debug)std::cout<<"fill_supercl mini end"<<std::endl; 

}

// end FC

template < typename T> 
bool RegFatPFNtuple::validPtr(edm::Ptr<T>& ptr){
  return (ptr.isNonnull() && ptr.isAvailable() );
}

