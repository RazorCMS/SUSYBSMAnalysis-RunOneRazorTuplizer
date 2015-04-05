from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *

# in 42X, it has to be customized to use run/lumi/event filter, otherwise it needs AOD
from RecoMET.METFilters.hcalLaserEventFilter_cfi import *
hcalLaserEventFilter.taggingMode = True
hcalLaserEventFilter.vetoByRunEventNumber = False
hcalLaserEventFilter.vetoByHBHEOccupancy= True

from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi import *
EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)

from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = True

from RecoMET.METFilters.ecalLaserCorrFilter_cfi import *
ecalLaserCorrFilter.taggingMode = True


goodVertices = cms.EDFilter(
   "VertexSelector",
   filter = cms.bool(False),
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)


## tracking failure filter
from RecoMET.METFilters.trackingFailureFilter_cfi import *
trackingFailureFilter.taggingMode = cms.bool(True)

## The tracking POG filters 
from RecoMET.METFilters.trackingPOGFilters_cff import *
## NOTE: to make tagging mode of the tracking POG filters (three of them), please do:
manystripclus53X.taggedMode = cms.untracked.bool(True)
manystripclus53X.forcedValue = cms.untracked.bool(False)
toomanystripclus53X.taggedMode = cms.untracked.bool(True)
toomanystripclus53X.forcedValue = cms.untracked.bool(False)
logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
## Also the stored boolean for the three filters is opposite to what we usually
## have for other filters, i.e., true means rejected bad events while false means 
## good events.


############################################
# Extra Sequences ported from Vecbos
############################################
from SUSYBSMAnalysis.RunOneRazorTuplizer.EcalDeadCellEventFlagProducer_cfi import *
from SUSYBSMAnalysis.RunOneRazorTuplizer.simpleDRFlagProducer_cfi import *

from SUSYBSMAnalysis.RunOneRazorTuplizer.ecalanomalouseventfilter_cfi import *
EcalAnomalousEventFilter.recHitsEB = 'reducedEcalRecHitsEB'
EcalAnomalousEventFilter.recHitsEE = 'reducedEcalRecHitsEE'


metOptionalFilterSequence = cms.Sequence( HBHENoiseFilterResultProducer
* hcalLaserEventFilter
* EcalDeadCellTriggerPrimitiveFilter 
* EcalDeadCellBoundaryEnergyFilter
* goodVertices * trackingFailureFilter
* eeBadScFilter
* ecalLaserCorrFilter
* trkPOGFilters
* EcalDeadCellEventFlagProducer
* simpleDRFlagProducer
* EcalAnomalousEventFilter
)
