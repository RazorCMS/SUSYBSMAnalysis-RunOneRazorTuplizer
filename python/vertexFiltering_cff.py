import FWCore.ParameterSet.Config as cms

VERTEX_SEL=("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")

goodPrimaryVertices = cms.EDFilter("VertexSelector",
                                   src = cms.InputTag("offlinePrimaryVertices"),
                                   cut = cms.string(VERTEX_SEL),
                                   filter = cms.bool(True),
                                   )
