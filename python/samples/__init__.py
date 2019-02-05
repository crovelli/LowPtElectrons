import fill6371 as fill6371AOD 
import RAWBlocks
import BtoKee as BtoKeeRAW
import os

def file2list(fname):
   fullname = os.path.join(
      os.environ['CMSSW_BASE'], 
      'src/LowPtElectrons/LowPtElectrons/python/samples', 
      fname
      )   
   return [i.strip() for i in open(fullname)]

all_samples = {
   'fill6371AOD' : fill6371AOD.input_files,
   'RAWBlocks' : RAWBlocks.input_files,
   'RAWTest' : ['/store/data/Run2017F/SingleElectron/RAW/v1/000/306/459/00000/B6410DA6-C8C5-E711-947F-02163E01A45E.root'],
   'EleGunTest' : ['file:/afs/cern.ch/work/m/mverzett/public/ELEGUN_GEN-SIM-RAW.root'],
   'RAWMC' : BtoKeeRAW.all_files,
   'RAWMCTest' : BtoKeeRAW.test_file,
   'BToJPsieeK_0' : file2list('BToJPsieeK_0.list'),
   'BToJPsieeK_1' : file2list('BToJPsieeK_1.list'),
   'BToJPsieeK' : file2list('BToJPsieeK.list'),
   'BToKStaree' : file2list('BToKStaree.list'),
   'BsToJPsieePhi' : file2list('BsToJPsieePhi.list'),
}
