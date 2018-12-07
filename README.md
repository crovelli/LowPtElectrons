# LowPtElectrons (production branch)

First, setup your env. Then, create release area.
```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
```

Some initialisation.
```
git cms-init
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git remote -v
```

Merge topic for new electron reco (production branch).
```
git cms-merge-topic bainbrid:LowPtElectrons_prod
```

Checkout ntuplizer code.
```
cd $CMSSW_BASE/src
mkdir LowPtElectrons
cd LowPtElectrons
git clone git@github.com:bainbrid/LowPtElectrons.git
cd LowPtElectrons
git remote add bainbrid git@github.com:bainbrid/LowPtElectrons.git
git fetch bainbrid LowPtElectrons_prod
git checkout -b LowPtElectrons_prod bainbrid/LowPtElectrons_prod
```

Build and run.
``` 
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
export INPUT_FILES=root://cms-xrd-global.cern.ch//store/cmst3/group/bpark/BToKee_Pythia_PUMix_18_03_18_180318_112206_0000/BToKee_PUMix_10.root
cmsRun mc_features.py inputFiles=$INPUT_FILES maxEvents=1
```

