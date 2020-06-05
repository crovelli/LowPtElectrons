```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

### Get latest code and model for electron ID (Aug07) - and energy regression
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
git cms-addpkg RecoEgamma/EgammaTools
git cms-addpkg RecoEgamma/ElectronIdentification
git cms-merge-topic crovelli:from-CMSSW_10_2_15_Aug07Id_and_regr
```

### Install ntuplizer code
```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 
cd LowPtElectrons/LowPtElectrons
git remote add crovelli git@github.com:crovelli/LowPtElectrons.git
git fetch crovelli
git checkout -b regression crovelli/regression
scram b
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun ntuplizer.py 
```
