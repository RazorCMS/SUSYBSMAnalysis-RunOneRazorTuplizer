from SUSYBSMAnalysis.RunOneRazorTuplizer.EcalDeadCellEventFlagProducer_cfi import *
from SUSYBSMAnalysis.RunOneRazorTuplizer.simpleDRFlagProducer_cfi import *

# the CSC Halo filter only works in the RECO
#from SUSYBSMAnalysis.RunOneRazorTuplizer.CSCHaloFlagProducer_cfi import *

from RecoMET.METFilters.trackingFailureFilter_cfi import *
trackingFailureFilter.taggingMode = True
trackingFailureFilter.VertexSource = "goodPrimaryVertices"

from SUSYBSMAnalysis.RunOneRazorTuplizer.ecalanomalouseventfilter_cfi import *
EcalAnomalousEventFilter.recHitsEB = 'reducedEcalRecHitsEB'
EcalAnomalousEventFilter.recHitsEE = 'reducedEcalRecHitsEE'

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *

# in 42X, it has to be customized to use run/lumi/event filter, otherwise it needs AOD
from RecoMET.METFilters.hcalLaserEventFilter_cfi import *
hcalLaserEventFilter.taggingMode = True
hcalLaserEventFilter.vetoByRunEventNumber = False
hcalLaserEventFilter.vetoByHBHEOccupancy= True

from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = True

from RecoMET.METFilters.ecalLaserCorrFilter_cfi import *
ecalLaserCorrFilter.taggingMode = True

metOptionalFilterSequence = cms.Sequence( EcalDeadCellEventFlagProducer
* HBHENoiseFilterResultProducer
* hcalLaserEventFilter
* eeBadScFilter
* simpleDRFlagProducer
* trackingFailureFilter
* EcalAnomalousEventFilter
* ecalLaserCorrFilter)
