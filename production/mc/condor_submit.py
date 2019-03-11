import os
from argparse import ArgumentParser
import datetime
import stat
from LowPtElectrons.LowPtElectrons.samples import all_samples

parser = ArgumentParser()
parser.add_argument('what', choices=all_samples.keys())
parser.add_argument('--group', type=int, default=1)
args = parser.parse_args()

date = datetime.date.today().strftime('%Y%b%d')

condir = 'condor_%s_%s' % (date, args.what)
outdir = '/eos/cms/store/cmst3/group/bpark/%s_%s' % (date, args.what)
try:
   os.makedirs(condir)
   os.makedirs(outdir)
except OSError:
   pass

dotsh = '''#! /bin/bash

sbox=$PWD
cd {CMSSW}
eval `scramv1 runtime -sh`
cd $sbox

cmsRun {CFG} $@
mv *.root {OUTDIR}/.
'''.format(
   CMSSW = os.environ['CMSSW_BASE'],
   CFG = os.path.join(
      os.environ['CMSSW_BASE'],
      'src/LowPtElectrons/LowPtElectrons/run/mc_features.py'
      ),
   OUTDIR = outdir,
   )

batch = '%s/batch.sh' % condir
with open(batch, 'w') as out:
   out.write(dotsh)

st = os.stat(batch)
os.chmod(batch, st.st_mode | stat.S_IEXEC)

nchuncks = len(all_samples[args.what])//args.group

jdl='''Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = True
executable = {CMSSW}/src/LowPtElectrons/LowPtElectrons/production/mc/{CONDIR}/batch.sh
use_x509userproxy = True
+MaxRuntime = 160600
+AccountingGroup = "group_u_CMST3.all"

Output = con_df_$(ProcId).out
Error = con_df_$(ProcId).err
Log = con_df_$(ProcId).log
Arguments = ichunk=$(ProcId) nchunks={NCHUNKS} outname={SAMPLE}_$(ProcId).root data={SAMPLE}
Queue {NCHUNKS}
'''.format(
   CMSSW = os.environ['CMSSW_BASE'],
   CONDIR = condir,
   NCHUNKS = nchuncks,
   SAMPLE = args.what,
)
with open('%s/condor.jdl' % condir, 'w') as out:
   out.write(jdl)

