# Recipes to add low pT electrons to RECO data tier

## Overview 

The new reconstruction needs to produce GsfTracks from scratch and so the RAW data tier is needed as input. 

The low pT reconstruction is fully enabled by specifying the ```bParking``` era, either via a ```cmsRun --era=bParking ...``` command or in the cfg python configuration file directly, e.g. ```process = cms.Process('RECO',eras.Run2_2018,eras.bParking)```.

__bParking era__ 

In the __10.2.X__ release cycle, the ```bParking``` era *__must__* be specified, otherwise low pT electrons will not be produced. The region pT > 0.5 GeV is covered, and a logical OR of the Loose working points for the two BDT models employed in the Seeding module is used. All relevant collections (listed below) are stored in the RECO, AOD, and MINIAOD data tiers. 

In the __10.5.X__ release cycle (and beyond), the same behaviour is observed as above if the ```bParking``` era is specified. If not, then low pT electrons are still produced but with tighter requirements: only the range pT > 1.0 GeV is available by default, and the Tight working points in the Seeding module are employed. All collections are still stored in RECO and AOD but *__not__* MINIAOD. 

The table below summarises the behaviour for different release cycles. 

| Release | bParking era? | pT threshold | Seed WP | Data tiers         |
| ---     | ---           | ---          | ---     | ---                |
| 10.5.X  | Yes           | 0.5          | Loose   | RECO, AOD, MINIAOD |
| 10.5.X  | No            | 1.0          | Tight   | RECO, AOD          |
| 10.2.X  | Yes           | 0.5          | Loose   | RECO, AOD, MINIAOD |
| 10.2.X  | No            | -            | -       | N/A                |

__Collections__

The low pT reconstruction uses the EGamma data formats. The GsfElectrons and their corresponding GsfElectronCores, the GsfTracks, the SuperClusters and constituent CaloClusters, are all stored. In MINIAOD, the "slimmed" PatElectrons replace the GsfElectrons. There are also three ValueMaps that stored the discriminant values of the BDTs used in the chain: the "unbiased" and "ptbiased" models in the Seeding module, plus the model defining the electron ID. The former two ValueMaps are keyed by a GsfTrackRef, while the latter is keyed by a GsfElectronRef. 

RECO and AOD:
```
recoGsfElectrons_lowPtGsfElectrons__RECO
recoGsfElectronCores_lowPtGsfElectronCores__RECO
recoGsfTracks_lowPtGsfEleGsfTracks__RECO
recoSuperClusters_lowPtGsfElectronSuperClusters__RECO
recoCaloClusters_lowPtGsfElectronSuperClusters__RECO
floatedmValueMap_lowPtGsfElectronSeedValueMaps_unbiased_RECO
floatedmValueMap_lowPtGsfElectronSeedValueMaps_ptbiased_RECO
floatedmValueMap_lowPtGsfElectronID__RECO
```

 MINIAOD:
```
patElectrons_slimmedLowPtElectrons__RECO
recoCaloClusters_lowPtGsfElectronSuperClusters__RECO
recoGsfElectronCores_lowPtGsfElectronCores__RECO
floatedmValueMap_lowPtGsfElectronSeedValueMaps_ptbiased_RECO
floatedmValueMap_lowPtGsfElectronID__RECO
floatedmValueMap_lowPtGsfElectronSeedValueMaps_unbiased_RECO
recoGsfTracks_lowPtGsfEleGsfTracks__RECO
recoSuperClusters_lowPtGsfElectronSuperClusters__RECO
```
## CMSSW_10_2_13 (recommended release)

```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git cms-init
```

Optional: check out and build.
```
git cms-addpkg RecoEgamma/EgammaElectronProducers
cd $CMSSW_BASE/src
scram b -j8 # repeat if necesssary
```

Optional: check out the BDT models:
```
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

Optional: checkout and run simple analyzer.
```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun LowPtGsfElectronsAnalyzer_cfg.py
```

## CMSSW_10_6_X (SLC7, development cycle)

```
ssh -X username@lxplus7.cern.ch
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_X_2019-02-27-1100
cd CMSSW_10_6_X_2019-02-27-1100/src
cmsenv
git cms-init
```

Optional: check out and build.
``` 
git cms-addpkg RecoEgamma/EgammaElectronProducers
cd $CMSSW_BASE/src
scram b -j8
```

Optional: earlier IBs would require the following:
```
git cms-merge-topic bainbrid:LowPtElectronsFull_105X_Autumn18
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

Optional: checkout and run simple analyzer.
```
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun LowPtGsfElectronsAnalyzer_cfg.py
```

## CMSSW_10_5_0_pre2 (closed release cycle)

Setup env and create release area.
```
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_5_0_pre2
cd CMSSW_10_5_0_pre2/src
cmsenv
git cms-init
```

Check out and build.
```
git cms-merge-topic bainbrid:LowPtElectronsFull_105X_Autumn18
cd $CMSSW_BASE/src
scram b -j8
```

Add models from cms-data.
```
git clone git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

Optional: checkout and run simple analyzer.
```
git remote add bainbrid git@github.com:bainbrid/cmssw.git
git clone git@github.com:CMSBParking/LowPtElectrons.git $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons
cd $CMSSW_BASE/src
scram b -j8
cd $CMSSW_BASE/src/LowPtElectrons/LowPtElectrons/run
voms-proxy-init --voms cms
cmsRun LowPtGsfElectronsAnalyzer_cfg.py
```
