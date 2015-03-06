import FWCore.ParameterSet.Config as cms

myJetPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(False),
    src = cms.InputTag("genParticles")                            
)

myAK5PFJetPartonAssociation = cms.EDProducer("JetPartonMatcher",
    jets    = cms.InputTag("ak5PFJets"),
    partons = cms.InputTag("myJetPartons"),
    coneSizeToAssociate = cms.double(0.3),
)

myAK5PFJetFlavourAssociation = cms.EDProducer("JetFlavourIdentifier",
    srcByReference    = cms.InputTag("myAK5PFJetPartonAssociation"),
    physicsDefinition = cms.bool(False)
)

myAK5GenJetPartonAssociation = cms.EDProducer("JetPartonMatcher",
    jets    = cms.InputTag("ak5GenJets"),
    partons = cms.InputTag("myJetPartons"),
    coneSizeToAssociate = cms.double(0.3),
)

myAK5GenJetFlavourAssociation = cms.EDProducer("JetFlavourIdentifier",
    srcByReference    = cms.InputTag("myAK5GenJetPartonAssociation"),
    physicsDefinition = cms.bool(False)
)

myJetPartonMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = cms.InputTag("ak5PFJets"),         # RECO objects to match
    matched     = cms.InputTag("genParticles"),      # mc-truth particle collection
    mcPdgId     = cms.vint32(1, 2, 3, 4, 5, 21),     # one or more PDG ID (quarks except top; gluons)
    mcStatus    = cms.vint32(3),                     # PYTHIA status code (3 = hard scattering)
    checkCharge = cms.bool(False),                   # False = any value of the charge of MC and RECO is ok
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                   # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)

myJetGenJetMatch = cms.EDProducer("GenJetMatcher",   # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = cms.InputTag("ak5PFJets"),         # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("ak5GenJets"),        # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                   # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)

# default PAT sequence for jet flavour identification
myJetFlavourId = cms.Sequence(myJetPartons * myAK5PFJetPartonAssociation * myAK5PFJetFlavourAssociation * myAK5GenJetPartonAssociation * myAK5GenJetFlavourAssociation)
myGenJetFlavourId = cms.Sequence(myJetPartons * myAK5GenJetPartonAssociation * myAK5GenJetFlavourAssociation)
