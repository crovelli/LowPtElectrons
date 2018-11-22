OUTPUT=track_features.root
FAKES=-1
EDOUT=0
EVENTS=1
INPUT_FILES=root://cms-xrd-global.cern.ch//store/cmst3/group/bpark/BToKee_Pythia_PUMix_18_03_18_180318_112206_0000/BToKee_PUMix_10.root

cmsRun mc_features.py \
    inputFiles=$INPUT_FILES \
    outname=$OUTPUT \
    fakesMultiplier=$FAKES \
    edout=$EDOUT \
    maxEvents=$EVENTS
