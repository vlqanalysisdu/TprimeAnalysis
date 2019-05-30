import os
import FWCore.ParameterSet.Config as cms


process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.Tracer = cms.Service("Tracer")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

Summer16_23_dbFile = "sqlite:" + os.getcwd() + "/Summer16_23Sep2016V4_MC.db"

#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
			connect = cms.string('sqlite:Summer16_23Sep2016V4_MC.db'),
			#connect = cms.string(Summer16_23_dbFile),
			toGet = cms.VPSet(
				cms.PSet(
					record = cms.string('JetCorrectionsRecord'),
					tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
					label = cms.untracked.string('AK4PFchs')
					),
				cms.PSet(
					record = cms.string('JetCorrectionsRecord'),
					tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK8PFchs'),
					label = cms.untracked.string('AK8PFchs')
					)))
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

import FWCore.Utilities.FileUtils as FileUtils 
#readFiles = cms.untracked.vstring(FileUtils.loadListFromFile('TprimeBToTH_M-1500_LH_TuneCUETP8M1_13TeV.txt'))
#readFiles = cms.untracked.vstring(FileUtils.loadListFromFile('TEMP_NAME.txt'))

process.source = cms.Source("PoolSource",
			fileNames = cms.untracked.vstring(
			'/store/mc/RunIISummer16MiniAODv2/WW_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/00E95BAA-C3D7-E611-A416-0025905A60BC.root'
				))
		   #     fileNames = readFiles,  
	           #     skipBadFiles = cms.untracked.bool(True) 
		   #		)

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

### EGM 80X regression
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

### EGM scale and smearing correction         
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
			calibratedPatElectrons = cms.PSet(
				initialSeed = cms.untracked.uint32(12345),
				engineName = cms.untracked.string('TRandom3')
				),
			calibratedPatPhotons = cms.PSet(
				initialSeed = cms.untracked.uint32(12345),
				engineName = cms.untracked.string('TRandom3')
				),
			ggNtuplizer  = cms.PSet(
				initialSeed = cms.untracked.uint32(67890),
				engineName = cms.untracked.string('TRandom3')
				)
			)

process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
process.calibratedPatElectrons.isMC = cms.bool(True)
process.calibratedPatPhotons.isMC = cms.bool(True)

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.TFileService = cms.Service("TFileService", fileName = cms.string('QCD_Pt_MUEnrichedPt_bkg.root'))

jecLevels = [
		'Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt',
		'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt'
	]

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
			process,
			jetSource = cms.InputTag('slimmedJets'),
			labelName = 'UpdatedJEC',
			jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			btagDiscriminators = ['deepFlavourJetTags:probudsg', 'deepFlavourJetTags:probb', 'deepFlavourJetTags:probc', 'deepFlavourJetTags:probbb', 'deepFlavourJetTags:probcc']
			)
updateJetCollection(
			process,
			jetSource = cms.InputTag('slimmedJetsAK8'),
			labelName = 'UpdatedJECAK8',
			jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			btagDiscriminators = ['pfBoostedDoubleSecondaryVertexAK8BJetTags'],
			btagPrefix = 'newV4' # optional, in case interested in accessing both the old and new discriminator values
			)

# jettool box
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
			PUMethod='CHS',
			updateCollection='slimmedJetsAK8', 
			JETCorrPayload='AK8PFchs', 
			addNsub=True,
			maxTau=4,
		  )

jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
			PUMethod='Puppi',
			addPruning=True,
			addSoftDrop=True ,           # add basic grooming
			addTrimming=True,
			addFiltering=True, 
			addSoftDropSubjets=True,
			addNsub=True,
			maxTau=4,                       # add Nsubjettiness tau1, tau2, tau3, tau4
			JETCorrPayload = 'AK8PFPuppi',
			JETCorrLevels = ['L2Relative', 'L3Absolute']
			#PUMethod='Puppi',
			#JETCorrPayload = 'AK8PFPuppi',
			#addNsub=True,
			#maxTau=4,
		  )

# MET correction and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
		isData=False
		)

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.load("ggAnalysis.ggNtuplizer.ggPhotonIso_CITK_PUPPI_cff")
process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.jecAK8PayloadNames=cms.vstring(jecLevels)
process.ggNtuplizer.runHFElectrons=cms.bool(True)
process.ggNtuplizer.isAOD=cms.bool(False)
process.ggNtuplizer.doGenParticles=cms.bool(True)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.patTriggerResults=cms.InputTag("TriggerResults", "", "PAT")

#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
	'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
	'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
	'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
	'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
	'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedPhotons')
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedPhotons')

#process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
#                                   src = cms.InputTag("slimmedMuons"),
#                                   preselection = cms.string("track.isNonnull"),
#                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
#                                   fractionOfSharedSegments = cms.double(0.499))

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.p = cms.Path(
			process.regressionApplication*
			process.ggMETFiltersSequence*
			process.calibratedPatElectrons*
			process.calibratedPatPhotons* 
			process.egmGsfElectronIDSequence*
			process.egmPhotonIDSequence*
			#process.cleanedMu*
			process.ggNtuplizer
			)

#print process.dumpPython()
