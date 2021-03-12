```
cmsrel CMSSW_10_6_20 
cd CMSSW_10_6_20/src
cmsenv
git cms-init
```

### Install ntuplizer code
```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons 
cd LowPtElectrons/LowPtElectrons
git remote add crovelli git@github.com:crovelli/LowPtElectrons.git
git fetch crovelli
git checkout -b for106x_bkgFix crovelli/for106x_bkgFix
```

### Compile
```
scram b
```

### Run
```
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun ntuplizer.py 
```
