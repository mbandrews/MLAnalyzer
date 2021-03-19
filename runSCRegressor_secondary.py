import os

cfg='RecHitAnalyzer/python/SCRegressor_cfg_secondary.py'
#cfg='RecHitAnalyzer/python/SCRegressor_cfg_data.py'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/04B84384-FF41-E811-ACA3-7845C4FC3B18.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePhotonPt10To100_pythia8_ReAOD_PU2017_MINIAODSIM/190611_010300/0000/step_MINIAODSIMext2_1.root'
#inputFiles_='file:/uscms/home/mba2012/nobackup/GUN/CMSSW_9_4_13/src/step_MINIAODSIMext2.root'
#inputFiles_='file:step_full_filtered.root'
#inputFiles_='file:QCD_Pt-40toInf.root'
#inputFiles_='root://cmseos.fnal.gov/%s'%inputFiles_
#inputFiles_='file:/eos/uscms%s'%inputFiles_
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0C57DB95-FA41-E811-B58D-008CFAFBEBF2.root'
#root://cmsxrootd-site.fnal.gov

maxEvents_=-1
#maxEvents_=1000
#maxEvents_=100
skipEvents_=0
#outputFile_='output.root'
#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
#inputTag='TEST'
#evtList_ = '../flashgg/eventsToProc_h24g_ma1GeV_2photons_EBonly.txt'
evtList_ = 'DoubleEG_2017B_2photons_ggskim_event_list.txt'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
#for ievt in range(1):
#if not os.path.isdir(inputTag):
#    os.system('mkdir %s'%(inputTag))
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
cmd="cmsRun %s maxEvents=%d skipEvents=%d"%(cfg,maxEvents_,skipEvents_)
#cmd="cmsRun %s maxEvents=%d skipEvents=%d eventsToProcess_load=%s"%(cfg,maxEvents_,skipEvents_,evtList_)
print '%s'%cmd
os.system(cmd)
#os.system('mv c*.eps %s/'%(inputTag))
#    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
