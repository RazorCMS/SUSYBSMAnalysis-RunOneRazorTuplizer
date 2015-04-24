import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("razorTuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12_DR53X/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V19-v1/00000/00903D2F-3E44-E311-8AAB-00266CF9B274.root' #MC example file        
        #'/store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/30002/FEDB1808-B48A-E211-96A8-20CF3027A59B.root' #Data example file
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#For Data
#process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
#For MC
process.GlobalTag.globaltag = 'START53_V27::All'


#For Transient Track Builder
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

#For BTagging Discriminators
process.load('SUSYBSMAnalysis.RunOneRazorTuplizer.btagging_cff')

#For MC jet matching
process.load('SUSYBSMAnalysis.RunOneRazorTuplizer.jetflavorMatching_cff')

#For MET Filters
process.load("SUSYBSMAnalysis.RunOneRazorTuplizer.MetOptionalFilters_cff")

#For Pileup Jet ID
from CMGTools.External.pujetidsequence_cff import puJetId
from CMGTools.External.pujetidsequence_cff import puJetMva

process.recoPuJetId = puJetId.clone(
   jets = cms.InputTag("ak5PFJets"),
   applyJec = cms.bool(True),
   inputIsCorrected = cms.bool(False),                
)

process.recoPuJetMva = puJetMva.clone(
   jets = cms.InputTag("ak5PFJets"),
   jetids = cms.InputTag("recoPuJetId"),
   produceJetIds = cms.bool(True),
   applyJec = cms.bool(True),
   inputIsCorrected = cms.bool(False),                
)



#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorTuplizer', 
    isData = cms.bool(False),    
    useGen = cms.bool(True),
    enableTriggerInfo = cms.bool(True),                                 
    triggerPathNamesFile = cms.string("SUSYBSMAnalysis/RunOneRazorTuplizer/data/RunOneRazorHLTPathnames.dat"),
    vertices = cms.string("offlinePrimaryVertices"),    
    muons = cms.string("muons"),
    electrons = cms.string("gsfElectrons"),
    #taus = cms.string("slimmedTaus"),
    photons = cms.string("photons"),
    jets = cms.string("ak5PFJets"),
    jetsAK8 = cms.string("ak8PFJets"),
    genMets = cms.string("genMetTrue"),
    mets = cms.string("pfMet"),
    pfCands = cms.string("particleFlow"),
    genParticles = cms.string("genParticles"),
    genJets = cms.string("ak5GenJets"),
    triggerBits = cms.string("TriggerResults"),
    lheInfo = cms.string("externalLHEProducer"),
    genInfo = cms.string("generator"),
    puInfo = cms.string("addPileupInfo"), #uncomment if no pre-mixing
    beamSpot = cms.string("offlineBeamSpot"),
    ebRecHits = cms.string("reducedEcalRecHitsEB"),
    eeRecHits = cms.string("reducedEcalRecHitsEE"),
    conversions = cms.string("allConversions"),
    csvBJetTags = cms.string("newCombinedSecondaryVertexBJetTags"),              
    cisvBJetTags = cms.string("myCombinedInclusiveSecondaryVertexBJetTags"),
    jetFlavorMatch = cms.string("myAK5PFJetFlavourAssociation"),

)

#run
process.p = cms.Path( process.metOptionalFilterSequence*
                      process.myJetFlavourId*
                      process.newJetBtagging*
                      process.recoPuJetId*
                      process.recoPuJetMva*
                      process.ntuples)

