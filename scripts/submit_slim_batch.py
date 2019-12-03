#! /usr/bin/env python
#
# example: python submit_slim_batch.py -c -N -1 -n 1 -p testslim --cfg=slimntuplizerbatch_cfg.py BuToKJpsiToee.dat
# -c = just create and do not submit (remove it to submit)
# this is to write 1 job per file in the dataset (-N = run on all the events in a file; -n = 1job/file)
#
# to copy the output on EOS: python submit_slim_batch.py -c -N -1 -n 1 -p testslim --cfg=slimntuplizerbatch_cfg.py --eos=ok BuToKJpsiToee.dat 
#
# if files are accessed via Grid, before you need:
# voms-proxy-init --voms cms --valid 168:00
# mv /tmp/x509up_u7291 ~/X509_USER_PROXY
# setenv X509_USER_PROXY ~/X509_USER_PROXY
# echo $X509_USER_PROXY 
#
# BuToKJpsiToee.dat contains the dataset files. To be made with 
# das_client --query 'file dataset=/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM' --limit 0 > BuToKJpsiToee.dat

import os
import sys
import re
import time
import commands
import optparse
import datetime

def makeCondorFile(jobdir, srcFiles, options):
    dummy_exec = open(jobdir+'/dummy_exec.sh','w')
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('bash $*\n')
    dummy_exec.close()
     
    condor_file_name = jobdir+'/condor_submit.condor'
    condor_file = open(condor_file_name,'w')
    condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = true
Log        = {jd}/$(ProcId).log
Output     = {jd}/$(ProcId).out
Error      = {jd}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
request_memory = 2000
+MaxRuntime = {rt}\n
'''.format(de=os.path.abspath(dummy_exec.name), jd=os.path.abspath(jobdir), rt=int(options.runtime*3600), here=os.environ['PWD'] ) )
    if os.environ['USER'] in ['mdunser', 'psilva']:
        condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n\n\n')
    for sf in srcFiles:
        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(sf)))
    condor_file.close()
    return condor_file_name

def main():
#######################################
### usage  
#######################################
    usage = '''usage: %prog [opts] --cfg cmssw.py dataset'''
    parser = optparse.OptionParser(usage=usage)
    now = datetime.datetime.now()
    defaultoutputdir='job_'+str(now.year)+str(now.month)+str(now.day)+"_"+str(now.hour)+str(now.minute)+str(now.second)
        
    parser.add_option('-q', '--queue',       action='store',     dest='queue',       help='run in batch in queue specified as option (default -q 8nh)', default='8nh')
    parser.add_option('-n', '--nfileperjob', action='store',     dest='nfileperjob', help='split the jobs with n files read/batch job'                , default=1,   type='int')
    parser.add_option('-p', '--prefix',      action='store',     dest='prefix',      help='the prefix to be added to the output'                      , default=defaultoutputdir)
    parser.add_option('-a', '--application', action='store',     dest='application', help='the executable to be run'                                  , default='cmsRun')
    parser.add_option('-d', '--download',    action='store',     dest='download',    help='download the output on a local computer'                   , default='')
    parser.add_option('-c', '--create',      action='store_true',dest='create',      help='create only the jobs, do not submit them'                  , default=False)
    parser.add_option('-t', '--testnjobs',   action='store',     dest='testnjobs',   help='submit only the first n jobs'                              , default=1000000, type='int')
    parser.add_option('-N', '--neventsjob', action='store',      dest='neventsjob',  help='split the jobs with n events  / batch job'                 , default=200,   type='int')
    parser.add_option('-T', '--eventsperfile', action='store',   dest='eventsperfile',  help='number of events per input file'                        , default=-1,   type='int')
    parser.add_option('-r', '--runtime',     action='store',     dest='runtime',     help='New runtime for condor resubmission in hours. default None: will take the original one.', default=4        , type=int);
    parser.add_option('--eos',               action='store',     dest='eos',         help='copy the output in the specified EOS path'                 , default='')
    parser.add_option('--cfg',               action='store',     dest='cfg',         help='the cfg to be run'                                         , default='pippo_cfg.py')
    parser.add_option('--scheduler',         action='store',     dest='scheduler',   help='select the batch scheduler (lsf,condor). Default=condor'   , default='condor')
    (opt, args) = parser.parse_args()
    print args

    if len(args) != 1:
        print len(args)
        print usage
        sys.exit(1)
    inputlist = args[0]

    output = os.path.splitext(os.path.basename(inputlist))[0]

    print "the outputs will be in the directory: "+opt.prefix

    if opt.download=='pccmsrm':
        diskoutputdir = "/cmsrm/pc24_2/emanuele/data/EcalReco7.1.X/"
    else: diskoutputdir = ''
    diskoutputmain = diskoutputdir+"/"+opt.prefix+"/"+output

    jobdir = opt.prefix+"/"+output
    logdir = jobdir+"/log/"
    os.system("mkdir -p "+jobdir)
    os.system("mkdir -p "+logdir)
    os.system("mkdir -p "+jobdir+"/src/")
    os.system("mkdir -p "+jobdir+"/cfg/")

    outputroot = diskoutputmain+"/root/"

    if (diskoutputdir != "none" and opt.download=='pccmsrm'): 
        os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm24 mkdir -p "+diskoutputmain)


    #look for the current directory
    #######################################
    pwd = os.environ['PWD']
    scramarch = os.environ['SCRAM_ARCH']
    #######################################
    inputListfile=open(inputlist)
    inputfiles = inputListfile.readlines()
    ijob=0
    jobdir = opt.prefix+"/"+output

    srcfiles = []
    while (len(inputfiles) > 0):
        L = []
        for line in range(min(opt.nfileperjob,len(inputfiles))):
            ntpfile = inputfiles.pop()
            ntpfile = ntpfile.rstrip('\n')
            ntpfile = re.sub(r'/eos/cms','',ntpfile.rstrip())
            if ntpfile != '':
                L.append("\'"+ntpfile+"\',\n")

        firstEvent = 1
        while (firstEvent < opt.eventsperfile or opt.eventsperfile == -1):
            lastEvent = firstEvent+opt.neventsjob
                
            # prepare the cfg
            icfgfilename = pwd+"/"+opt.prefix+"/"+output+"/cfg/cmssw"+str(ijob)+"_cfg.py"
            icfgfile = open(icfgfilename,'w')
            icfgfile.write('import sys\n')
            cfgfile=open(opt.cfg,'r')
            stringtoreplace = ''.join(L)
            stringtoreplace = stringtoreplace[:-2] # remove the "," and the end of line for the last input
            stringtoreplace = 'fileNames = cms.untracked.vstring('+stringtoreplace+')\n#'
            if (opt.eventsperfile == -1): maxEventsString = 'process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )#'
            else: maxEventsString = 'process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32('+str(opt.neventsjob)+') )#'

            for line in cfgfile:
                line = re.sub(r'fileNames = cms.untracked.vstring',stringtoreplace, line.rstrip())
                line = re.sub(r'fileName = cms.untracked.string','fileName = cms.untracked.string("/tmp/'+output+'_'+str(ijob)+'.root"),#', line.rstrip())
                line = re.sub(r'fileName = cms.string','fileName = cms.string("/tmp/'+output+'_'+str(ijob)+'.root") #', line.rstrip())
                line = re.sub(r'process.maxEvents = cms.untracked.PSet', maxEventsString, line.rstrip())
                icfgfile.write(line+'\n')

            if (opt.eventsperfile > -1): icfgfile.write('process.source.skipEvents=cms.untracked.uint32('+str(firstEvent-1)+')\n')
            icfgfile.close()

            # prepare the script to run
            rootoutputfile = output+'_'+str(ijob)+'.root'
            outputname = jobdir+"/src/submit_"+str(ijob)+".src"
            outputfile = open(outputname,'w')
            outputfile.write('#!/bin/bash\n')
            outputfile.write('export SCRAM_ARCH='+scramarch+'\n')
            outputfile.write('cd '+pwd+'\n')
            outputfile.write('eval `scramv1 runtime -sh`\n')
            if opt.scheduler=='lsf':
                outputfile.write('cd $WORKDIR\n')
            elif opt.scheduler=='condor':
                #outputfile.write('cd /tmp/'+os.environ['USER']+'\n')
                rootoutputfile = '/tmp/'+rootoutputfile
            outputfile.write('echo $PWD\n')
            #outputfile.write(opt.application+' '+icfgfilename+' "'+rootoutputfile+'" \n')
            outputfile.write(opt.application+' '+icfgfilename+' \n')
            if(opt.download=='pccmsrm'): outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm24:'+diskoutputmain+'/{}\n') 
            if(opt.eos!=''): 
                #outputfile.write('cp '+rootoutputfile+' /eos/cms/store/user/crovelli/LowPtEle/Batch'+opt.eos+'/'+output+'_'+str(ijob)+'.root\n')
                outputfile.write('cp '+rootoutputfile+' /eos/cms/store/user/crovelli/LowPtEle/Batch2/'+output+'_'+str(ijob)+'.root\n')
                outputfile.write('rm '+rootoutputfile)
            outputfile.close()
            logfile = logdir+output+"_"+str(ijob)+".log"
            scriptfile = pwd+"/"+outputname
            if opt.scheduler=='lsf':
                cmd = "bsub -q "+opt.queue+" -o "+logfile+" source "+scriptfile
                print cmd
                if not opt.create:
                    os.system(cmd)
            elif opt.scheduler=='condor':
                srcfiles.append(outputname)
            else:
                print "ERROR. Scheduler ",opt.scheduler," not implemented. Choose either 'lsf' or 'condor'."
                sys.exit(1)
            ijob = ijob+1
            if(ijob==opt.testnjobs): break
            if (opt.eventsperfile == -1): break
            else: firstEvent = lastEvent

    if opt.scheduler=='condor':
        cf = makeCondorFile(jobdir,srcfiles,opt)
        subcmd = 'condor_submit {rf} '.format(rf = cf)
        if opt.create:
            print 'running dry, printing the commands...'
            print subcmd
        else:
            print 'submitting for real...'
            os.system(subcmd)

if __name__ == "__main__":
        main()

