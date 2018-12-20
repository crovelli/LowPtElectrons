# LowPtElectrons (10X branches)

Choose option ```#1``` or ```#2``` (in two places, be consistent)

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_4_0_pre4; cd CMSSW_10_4_0_pre4/src.  #1
cmsrel CMSSW_10_2_9; cd CMSSW_10_2_9/src             #2
cmsenv
```

Some initialisation.
```
git cms-initproduction
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git remote -v
```

Merge topic for new electron reco (10X branches).
```
git cms-merge-topic bainbrid:LowPtElectrons_10X      #1
git cms-merge-topic bainbrid:LowPtElectrons_10_2_9   #2
```

Checkout ntuplizer code.
```
mkdir $CMSSW_BASE/src/LowPtElectrons
cd $CMSSW_BASE/src/LowPtElectrons
git clone git@github.com:bainbrid/LowPtElectrons.git
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
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
