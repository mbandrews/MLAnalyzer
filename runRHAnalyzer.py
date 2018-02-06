import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:SinglePhotonPt50_FEVTDEBUG.root'
#inputFiles_='file:pickeventsRECO.root'
#inputFiles_='file:SingleElectronPt50_FEVTDEBUG.root'
#inputFiles_='file:../step_full.root'
#inputFiles_='file:../step_full_H2GG_MoriondPU.root'
#inputFiles_='file:../step_full_H2GG_noPU.root'
#inputFiles_='file:../step_full_H2GG_HLPU.root'
inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/H125GGgluonfusion_13TeV_TuneCUETP8M1_FEVTDEBUG/170529_103138/0000/step_full_1.root'
#inputFiles='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/DoublePhotonGaussPt55_StdDev20_FEVTDEBUG/170606_115953/0002/step_full_2000.root'
#inputFiles_='file:/uscms/home/mba2012/work/MLHEP/GUN/CMSSW_8_0_24_patch1/src/step_full.root'
#inputFiles_='root://cmseos.fnal.gov:1094//store/user/mba2012/FEVTDEBUG/step_full_17.root'
#inputFiles_='file:../step3.root'
maxEvents_=-1
skipEvents_=0
#outputFile_='output.root'
inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
inputTag='TEST'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
#for ievt in range(1):
if not os.path.isdir(inputTag):
    os.system('mkdir %s'%(inputTag))
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
#os.system('mv c*.eps %s/'%(inputTag))
#    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
