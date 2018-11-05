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
git clone https://github.com/bainbrid/LowPtElectrons.git
git fetch bainbrid LowPtElectons_prod
```

Build and run.
``` 
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
cmsRun mc_features.py maxEvents=1
```

