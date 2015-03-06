from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola__Summer12_DR53X-PU_S10_START53_V19-v1__06March2015__V1' #change the request name for each new task
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/razorTuplizer_MC.py'
config.JobType.outputFiles = ['razorNtuple.root']

config.section_("Data")

#MC example
config.Data.inputDataset = '/TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM'

config.Data.dbsUrl = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
#config.Data.dbsUrl = 'phys03' 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03' #enable for publishing
#config.Data.publishDataName = 'razorNtuple'
#config.Data.ignoreLocality = False #disable AAA
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
