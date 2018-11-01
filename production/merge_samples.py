#! /bin/env python

import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('jobid')
parser.add_argument('--se', default='/eos/cms/store/cmst3/user/mverzett/')
parser.add_argument('--group', default=120, type=int)
args = parser.parse_args()

import ROOT as R

def syscall(cmd):
   print 'Executing: %s' % cmd
   retval = os.system(cmd)
   if retval != 0:
      raise RuntimeError('Command failed!')

import sys
toolbar_width=50
def progbar(frac):
   done = int(frac*toolbar_width)
   msg = "[%s%s] %d%%" % (
      '#'*done, " "*(toolbar_width-done),
      int(frac*100)
      )
   sys.stdout.write(msg)
   sys.stdout.flush()
   sys.stdout.write("\b" * len(msg)) # return to start of line, after '['

def hadd(ins, out, tree_path='features/tree'):
   print 'producing', out,' ',
   tc = R.TChain(tree_path)
   for i in ins:
      tc.Add(i)
   
   tf = R.TFile(out, 'RECREATE')
   td = tf.mkdir(os.path.dirname(tree_path))
   td.cd()
   otree = tc.CloneTree(0)
   entries = tc.GetEntries()
   print '(', entries, 'entries)'
   for idx in xrange(entries):
      if idx%500 == 0: progbar(idx/float(entries))
      tc.GetEntry(idx)
      otree.Fill()
   #otree.AutoSave()
   otree.Write()
   #tf.Write()
   print ''
   tf.Close()   

def chunks(l, n):
   """Yield successive n-sized chunks from l."""
   for i in range(0, len(l), n):
      yield l[i:i + n]

from glob import glob
from pdb import set_trace
indirs = glob('%s/*/crab_%s_*' % (args.se, args.jobid))

for sample in indirs:
   ins = glob('%s/*/*/*.root' % sample)
   print 'sample: ', sample
   print 'found', len(ins), 'input files'
   base_name = sample.split('_%s_' % args.jobid)[1]
   ins_chuncks = [i for i in chunks(ins, args.group)]
   for idx, chunk in enumerate(ins_chuncks):
      hadd(
         chunk,
         '%s_%s_%d.root' % (base_name, args.jobid, idx)         
         )
   
