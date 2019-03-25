# Low pT electrons

## Recipe to add low pT electrons to RECO data tier

Details can be found [here](LowPtElectrons.md). 

## Recipe to produce ntuples used for BDT training 

Ntuple production:
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git cms-init
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun mc_features.py 
```
