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

### This is required if running on CRAB!
```
git cms-addpkg RecoEgamma/ElectronIdentification
mv $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data # this is required if running on CRAB
 ```

### Install ntuplizer code
```
git clone git@github.com:crovelli/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 
git checkout -b regression
scram b
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun ntuplizer.py 
```
