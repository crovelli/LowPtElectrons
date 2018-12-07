#INPUT_FILES=root://cms-xrd-global.cern.ch//store/cmst3/group/bpark/BToKee_Pythia_PUMix_18_03_18_180318_112206_0000/BToKee_PUMix_10.root # 94X
#INPUT_FILES=root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_pre3/RelValZEE_13/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2018_realistic_v8-v1/20000/03156BB2-7F39-5B4E-B636-AB7C58CFF19D.root # 10X (Zee)
#INPUT_FILES=file:/eos/cms/store/group/phys_egamma/lowPtEle/doubleElePt0p5To10_1040pre2_NoPU_GEN_SIM_RAW.root # 10X (ele gun, PU=0)
INPUT_FILES=file:/eos/cms/store/group/phys_egamma/lowPtEle/doubleElePt0p5To10_1040pre2_PUPMX_GEN_SIM_RAW.root # 10X (ele gun, PU=nominal)
#INPUT_FILES=file:/eos/cms/store/group/phys_egamma/lowPtEle/doubleElePt0p5To10_1040pre2_PU20_GEN_SIM_RAW.root # 10X (ele gun, PU=20)

cmsRun test_10X.py inputFiles=$INPUT_FILES maxEvents=1
