```
cmsrel CMSSW_10_6_X_2020-11-24-2300 
cd CMSSW_10_6_X_2020-11-24-2300/src
cmsenv
git cms-init
```

### Get latest code and model for electron ID (2020Sept15, trained for 2018 b-park) - and energy regression
git cms-merge-topic bainbrid:LowPtElectrons_userFloats_106X

git clone --single-branch --branch from-CMSSW_10_2_15_2020Sept15 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data


### Install ntuplizer code
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 

cd LowPtElectrons/LowPtElectrons

git remote add crovelli git@github.com:crovelli/LowPtElectrons.git

git fetch crovelli

git checkout -b for106x crovelli/for106x

scram b

cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run

voms-proxy-init --voms cms

cmsRun ntuplizer.py 
```
