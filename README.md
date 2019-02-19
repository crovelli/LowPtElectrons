# LowPtElectrons (10X branches)

## 10.5.X development

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_5_X_2019-02-18-2300
cd CMSSW_10_5_X_2019-02-18-2300/src
cmsenv
```

Some initialisation.
```
git cms-init
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git remote -v
```

Merge topics.
```
git cms-merge-topic bainbrid: LowPtElectronsFull_105X_SCbugfix
```

Build.
``` 
cd $CMSSW_BASE/src
scram b -j8
```

Add models from cms-data.
```
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

Optional: checkout and run simple ntuplizer.
```
git clone git@github.com:bainbrid/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
. test_10X.sh
```

## 10.2.X development

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_X_2019-02-18-2300
cd CMSSW_10_2_X_2019-02-18-2300/src
cmsenv
```

Some initialisation.
```
git cms-init
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git remote -v
```

Merge topics.
```
git cms-merge-topic bainbrid: LowPtElectronsFull_102X
git cms-merge-topic bainbrid: LowPtElectronsFull_102X_update
```

Build.
``` 
cd $CMSSW_BASE/src
scram b -j8
```

Add models from cms-data.
```
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

Optional: checkout and run simple ntuplizer.
```
git clone git@github.com:bainbrid/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
. test_10X.sh
```
