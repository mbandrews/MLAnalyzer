import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents',
    default=0,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
options.register('eventsToProcess',
     default='',
     mult=VarParsing.VarParsing.multiplicity.list,
     mytype=VarParsing.VarParsing.varType.string,
     info = "Events to process")
options.parseArguments()

process = cms.Process("FEVTAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("DQM.Integration.config.FrontierCondition_GT_Offline_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1)
    input = cms.untracked.int32(options.maxEvents)
    )

print " >> Loaded",len(options.inputFiles),"input files from list."
#evtsToProc = open('list_pi0_ptrecoOgen1p2To1p6_eventsToProcess_.txt').read().splitlines()
#evtsToProc = open('pi0_evtsToProc.txt').read().splitlines()
#evtsToProc = open('eta_evtsToProc.txt').read().splitlines()
#print evtsToProc
slim_files = ['root://cmseos.fnal.gov//store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_1GeV_PU2017_AODSIM_slim/190719_005502/0000/step_aodsim_slim_%d.root'%i for i in range(13,15)]
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'file:myfile.root'
      #'/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_1GeV_TEST_PU2017_MINIAODSIM/190712_044728/0000/step_miniaodsim_ext3_1.root'
      "root://cms-xrd-global.cern.ch//store/group/phys_higgs/cmshgg/mandrews/flashgg/h24g_26Sep2019/RunIIFall18-4_0_0-119-g2d54185d/MINIAODSIM/h24g_26Sep2019-RunIIFall18-4_0_0-119-g2d54185d-v0-mandrews-h24gamma_1j_1M_1GeV_PU2017_MINIAODSIM-919c80a76a70185609d372d13ecbc645/190926_214616/0000/myMicroAODOutputFile_22.root"
      ),
    secondaryFileNames = cms.untracked.vstring(
      #'file:myfile.root'
      #'/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_1GeV_TEST_PU2017_AODSIM_trunc/190712_005252/0000/step_aodsim_trunc_1.root'
      #'/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_1GeV_TEST_PU2017_AODSIM_slim/190713_203315/0000/step_aodsim_slim_1.root'
      *slim_files
      )
    #, skipEvents = cms.untracked.uint32(options.skipEvents)
    #, eventsToProcess = cms.untracked.VEventRange('1:6931:1723687928','1:6932:1723895372')
    #, eventsToProcess = cms.untracked.VEventRange(*evtsToProc)
    #, lumisToProcess = cms.untracked.VLuminosityBlockRange('1:2133-1:2133')
    #, lumisToProcess = cms.untracked.VLuminosityBlockRange('1:3393-1:3393')
    )

if options.eventsToProcess:
    process.source.eventsToProcess = \
           cms.untracked.VEventRange (options.eventsToProcess)
#process.options = cms.untracked.PSet(
#)
#process.options.numberOfThreads=cms.untracked.uint32(4)

#process.GlobalTag.globaltag = cms.string('80X_dataRun2_HLT_v12')
process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v17')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.fevt = cms.EDAnalyzer('SCRegressor'
    #, EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , gsfElectronCollection = cms.InputTag('gedGsfElectrons')
    #, photonCollection = cms.InputTag('gedPhotons')
    , photonCollection = cms.InputTag('slimmedPhotons')
    , EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , EERecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEE')
    , ESRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsES')
    , reducedAODEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    , reducedAODEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    , reducedAODESRecHitCollection = cms.InputTag('reducedEcalRecHitsES')
    #, reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    #, reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    #, reducedESRecHitCollection = cms.InputTag('reducedEcalRecHitsES')
    , reducedEBRecHitCollection = cms.InputTag('reducedEgamma:reducedEBRecHits')
    , reducedEERecHitCollection = cms.InputTag('reducedEgamma:reducedEERecHits')
    , reducedESRecHitCollection = cms.InputTag('reducedEgamma:reducedESRecHits')
    , genParticleCollection = cms.InputTag('genParticles')
    #, genParticleCollection = cms.InputTag('prunedGenParticles')
    , genJetCollection = cms.InputTag('ak4GenJets')
    , trackCollection = cms.InputTag("generalTracks")
    , rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
    , trgResults = cms.InputTag("TriggerResults","","HLT")
    )

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('histo.root')
    fileName = cms.string(options.outputFile)
    )

process.hltFilter = cms.EDFilter("HLTHighLevel",
                                          eventSetupPathsKey = cms.string(''),
                                          TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                          HLTPaths = cms.vstring('HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_*_Mass55_v*'),
                                          andOr = cms.bool(True),
                                          throw = cms.bool(False)
                                          )

### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       era='2017-Nov17ReReco',
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
                                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
                       )

### reduce effect of high eta EE noise on the PF MET measurement
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
        process,
        isData = False, # false for MC
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
)

process.p = cms.Path(
  #process.hltFilter*
  #process.fullPatMetSequenceModifiedMET*
  #process.egammaPostRecoSeq*
  process.fevt
)
