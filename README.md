# LowPtElectrons (production branch)

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_4_X_2018-12-12-1100 # i.e. latest IB
cd CMSSW_10_4_X_2018-12-12-1100/src
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
git cms-merge-topic bainbrid:LowPtElectrons_10X
```

Checkout ntuplizer code.
```
mkdir $CMSSW_BASE/src/LowPtElectrons
cd $CMSSW_BASE/src/LowPtElectrons
git clone git@github.com:bainbrid/LowPtElectrons.git
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
git fetch origin LowPtElectrons_10X
git checkout LowPtElectrons_10X
```

Build.
``` 
cd $CMSSW_BASE/src
scram b -j8
```

Add models from cms-data.
```
cd $CMSSW_BASE/external/$SCRAM_ARCH
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
ls -l $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons
```

Run.
``` 
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
. test_10X.sh
```
