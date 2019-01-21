# LowPtElectrons (10X branches)

Link to ```RecoEgamma/EgammaElectronProducers``` package on ```bainbrid/cmssw:LowPtElectrons_10X branch``` [here](https://github.com/bainbrid/cmssw/tree/LowPtElectrons_10X/RecoEgamma/EgammaElectronProducers).

Recommended:
- Choose ```#4``` for full chain in ```CMSSW_10_2_X```
- Choose ```#5``` for full chain in ```CMSSW_10_5_X```

Obsolete:
- Choose option ```#1``` for the branch used for the "up to GsfTracks" PR.
- Or choose ```#2``` or ```#3``` for the full chain (up to GsfElectrons) in 10.4 or 10.2.

The options should be chosen consistently in two places.

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_4_0_pre4; cd CMSSW_10_4_0_pre4/src   #1
cmsrel CMSSW_10_4_0_pre4; cd CMSSW_10_4_0_pre4/src   #2
cmsrel CMSSW_10_2_9; cd CMSSW_10_2_9/src             #3
cmsrel CMSSW_10_2_10; cd CMSSW_10_2_10/src           #4
cmsrel CMSSW_10_5_X_2019-01-17-2300; cd CMSSW_10_5_X_2019-01-17-2300/src #5
cmsenv
```

Some initialisation.
```
git cms-init
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git remote -v
```

Merge topic for new electron reco (10X branches).
```
git cms-merge-topic bainbrid:LowPtElectrons_10X      #1
git cms-merge-topic bainbrid:LowPtElectronsFull_10X  #2
git cms-merge-topic bainbrid:LowPtElectrons_10_2_9   #3
git cms-merge-topic bainbrid:LowPtElectronsFull_102X #4
git cms-merge-topic bainbrid:LowPtElectronsFull_105X #5
```

Extra instructions for ```#4``` only.
```
cd $CMSSW_BASE/src
rm -r Validation/ SimTracker/ SimG4CMS/ L1Trigger/ DQM/ CalibTracker/ HLTrigger/ 
enw .git/info/sparse-checkout # remove from this list all packages consistent with above rm command
git read-tree -mu HEAD
git status # there should be no add/rm files
```

Build.
``` 
cd $CMSSW_BASE/src
scram b -j8
```

Add models from cms-data.
```
cd $CMSSW_BASE/external/$SCRAM_ARCH
git clone git@github.com:bainbrid/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
ls -l $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons
```

The following instructions are optional only. 

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
Run.
``` 
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
. test_10X.sh
```
