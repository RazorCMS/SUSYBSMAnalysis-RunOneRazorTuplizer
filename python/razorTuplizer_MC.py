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

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#For Data
#process.GlobalTag.globaltag = 'FT_53_V6_AN1::All'
#For MC
process.GlobalTag.globaltag = 'START53_V10::All'


#For Transient Track Builder
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

#For BTagging Discriminators
process.load('SUSYBSMAnalysis.RunOneRazorTuplizer.btagging_cff')

#For MC jet matching
process.load('SUSYBSMAnalysis.RunOneRazorTuplizer.jetflavorMatching_cff')

# --- good vertex filter ---
process.load("SUSYBSMAnalysis.RunOneRazorTuplizer.vertexFiltering_cff")

#For MET Filters
process.load("SUSYBSMAnalysis.RunOneRazorTuplizer.MetOptionalFilters_cff")

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorTuplizer', 
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
    #triggerPrescales = cms.string("patTrigger"),
    #triggerObjects = cms.string("selectedPatTrigger"),
    #metFilterBits = cms.string("TriggerResults", "", "PAT"),

    lheInfo = cms.string("externalLHEProducer"),
    genInfo = cms.string("generator"),
    puInfo = cms.string("addPileupInfo"), #uncomment if no pre-mixing
    #puInfo = cms.string("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    #hcalNoiseInfo = cms.string("hcalnoise", "", "RECO"),

    #secondaryVertices = cms.string("slimmedSecondaryVertices", "", "PAT"),

    #rhoAll = cms.string("kt6PFJets", "rho"),
    #rhoFastjetAll = cms.string("fixedGridRhoFastjetAll", "", "RECO"),
    #rhoFastjetAllCalo = cms.string("fixedGridRhoFastjetAllCalo", "", "RECO"),
    #rhoFastjetCentralCalo = cms.string("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    #rhoFastjetCentralChargedPileUp = cms.string("kt6PFJetsCentralChargedPileUp", "rho"),
    #rhoFastjetCentralNeutral = cms.string("kt6PFJetsCentralNeutral", "rho"),

    beamSpot = cms.string("offlineBeamSpot"),

    ebRecHits = cms.string("reducedEcalRecHitsEB"),
    eeRecHits = cms.string("reducedEcalRecHitsEE"),
    #esRecHits = cms.string("reducedEgamma", "reducedESRecHits", "PAT"),
    #ebeeClusters = cms.string("reducedEgamma", "reducedEBEEClusters", "PAT"),
    #esClusters = cms.string("reducedEgamma", "reducedESClusters", "PAT"),
    conversions = cms.string("allConversions"),
    csvBJetTags = cms.string("newCombinedSecondaryVertexBJetTags"),              
    cisvBJetTags = cms.string("myCombinedInclusiveSecondaryVertexBJetTags"),
    jetFlavorMatch = cms.string("myAK5PFJetFlavourAssociation"),
   #singleLegConversions = cms.string("reducedEgamma", "reducedSingleLegConversions", "PAT"),
    #gedGsfElectronCores = cms.string("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    #gedPhotonCores = cms.string("reducedEgamma", "reducedGedPhotonCores", "PAT"),
    #superClusters = cms.string("reducedEgamma", "reducedSuperClusters", "PAT"),

    #lostTracks = cms.string("lostTracks", "", "PAT")
)

#run
process.p = cms.Path( process.goodPrimaryVertices*
                      process.metOptionalFilterSequence*
                      process.myJetFlavourId*
                      process.newJetBtagging*
                      process.ntuples)
