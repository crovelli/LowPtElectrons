cmsenv
OLD=`pwd`
cd $CMSSW_BASE/src/
git cms-addpkg RecoEcal/EgammaCoreTools
git apply LowPtElectrons/LowPtElectrons/scripts/energy_matrix.patch
sed -i -e 's|//#define|#define|g' LowPtElectrons/LowPtElectrons/src/ElectronNtuple.cc
scram b
cd $OLD
