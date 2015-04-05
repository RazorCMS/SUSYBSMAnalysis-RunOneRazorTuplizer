// -*- C++ -*-
// Class:      RazorTuplizer
/*
Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#ifndef RAZORTUPLIZER_H
#define RAZORTUPLIZER_H

// system include files
#include <memory>
#include <string>
#include <vector>

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//CMSSW package includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
//#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//------ Class declaration ------//

class RazorTuplizer : public edm::EDAnalyzer {
public:
  //analyzer constructor and destructor
  explicit RazorTuplizer(const edm::ParameterSet&);
  ~RazorTuplizer();
  
  void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config
  virtual void resetBranches();
  
  //enable desired output variables
  virtual void setBranches();
  virtual void enableEventInfoBranches();
  virtual void enablePileUpBranches();
  virtual void enableMuonBranches();
  virtual void enableElectronBranches();
  virtual void enableTauBranches();
  virtual void enableIsoPFCandidateBranches();
  virtual void enablePhotonBranches();
  virtual void enableJetBranches();
  virtual void enableJetAK8Branches();
  virtual void enableMetBranches();
  virtual void enableRazorBranches();
  virtual void enableTriggerBranches();
  virtual void enableMCBranches();
  virtual void enableGenParticleBranches();
  
  //select objects and fill tree branches
  virtual bool fillEventInfo(const edm::Event& iEvent);
  virtual bool fillPileUp();//Fill summary PU info
  virtual bool fillMuons();//Fills looseID muon 4-momentum only. PT > 5GeV
  virtual bool fillElectrons();//Fills Ele 4-momentum only. PT > 5GeV
  virtual bool fillTaus();//Fills Tau 4-momentum only. PT > 20GeV
  virtual bool fillIsoPFCandidates();//Fills Isolated PF Candidates, PT > 5 GeV
  virtual bool fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup);//Fills photon 4-momentum only. PT > 20GeV && ISO < 0.3
  virtual bool fillJets(const edm::Event& iEvent );//Fills AK5 Jet 4-momentum, CSV, and CISV. PT > 20GeV 
  virtual bool fillJetsAK8();//Fills AK8 Jet 4-momentum.
  virtual bool fillMet(const edm::Event& iEvent);//Fills MET(mag, phi)
  virtual bool fillRazor();//Fills MR and RSQ
  virtual bool fillTrigger(const edm::Event& iEvent);//Fills trigger information
  virtual bool fillMC();
  virtual bool fillGenParticles();
  
  //------ HELPER FUNCTIONS ------//
  
  //splits jets into two hemisperes for razor variable calculation
  //(minimizes sum of mass^2's of hemispheres)
  vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
  //compute M_R using two hemispheres
  float computeMR(TLorentzVector hem1, TLorentzVector hem2);
  //compute R^2 using two hemispheres and MET vector
  float computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet);
  //returns true if particle 1 is an ancestor of particle 2, false otherwise
  //(takes two members of prunedGenParticles)
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
  //follows the particle's ancestry back until finding a particle of different type
  const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
  //follows the particle's ancestry back and finds the "oldest" particle with the same ID
  const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);
  //electron veto for photons (for use until an official recipe exists)
  bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<reco::GsfElectron> > &eleCol,
				const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot, 
				float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);
  bool isGoodPV( const reco::Vertex *v);
  bool isPFNoPU( const reco::PFCandidate candidate, const reco::Vertex *PV, edm::Handle<reco::VertexCollection> vertices);
  
  double getLeptonPtRel(edm::Handle<reco::PFJetCollection> jets, const reco::Candidate* lepton);

  double getPFMiniIsolation(edm::Handle<reco::PFCandidateCollection> pfcands,
			    const reco::Candidate* ptcl,
			    double r_iso_min = 0.05, double r_iso_max = 0.2 , double kt_scale = 10.0,
			    bool use_pfweight = false, bool charged_only = false);  
  
protected:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //MVAs for triggering and non-triggering electron ID
  EGammaMvaEleEstimator* myMVATrig;
  EGammaMvaEleEstimator* myMVANonTrig;
  
  //Photon energy regression helpers
  EGEnergyCorrector photonEnergyCorrector;

  //----- Member data ------//

  // Control Switches
  bool    useGen_;
  bool enableTriggerInfo_;
  
  // Input file containing the mapping of the HLT Triggers
  string triggerPathNamesFile_;

  // //EDM tokens for each miniAOD input object
  // edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  // edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  // edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
  // edm::EDGetTokenT<reco::TauCollection> tausToken_;
  // edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
  // edm::EDGetTokenT<reco::PFJetCollection> jetsToken_;
  // edm::EDGetTokenT<reco::PFJetCollection> jetsAK8Token_;
  // edm::EDGetTokenT<reco::PFCandidateCollection> packedPFCandsToken_;
  // edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  // edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  // edm::EDGetTokenT<reco::PFMETCollection> metToken_;
  // edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
  // edm::EDGetTokenT<LHEEventProduct> lheInfoToken_;
  // edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  // edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
  // edm::EDGetTokenT<HcalNoiseSummary> hcalNoiseInfoToken_;
  // //  edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > secondaryVerticesToken_;
  // edm::EDGetTokenT<double> rhoAllToken_;
  // edm::EDGetTokenT<double> rhoFastjetAllToken_;
  // edm::EDGetTokenT<double> rhoFastjetAllCaloToken_;
  // edm::EDGetTokenT<double> rhoFastjetCentralCaloToken_;
  // edm::EDGetTokenT<double> rhoFastjetCentralChargedPileUpToken_;
  // edm::EDGetTokenT<double> rhoFastjetCentralNeutralToken_;
  // edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  // edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
  // edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;
  // edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHitsToken_;
  // edm::EDGetTokenT<vector<reco::CaloCluster> > ebeeClustersToken_;
  // edm::EDGetTokenT<vector<reco::CaloCluster> > esClustersToken_;
  // edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
  // edm::EDGetTokenT<vector<reco::Conversion> > singleLegConversionsToken_;
  // edm::EDGetTokenT<vector<reco::GsfElectronCore> > gedGsfElectronCoresToken_;
  // edm::EDGetTokenT<vector<reco::PhotonCore> > gedPhotonCoresToken_;
  // edm::EDGetTokenT<vector<reco::SuperCluster> > superClustersToken_;
    //EDM tokens for each miniAOD input object


  string  verticesSrcName;
  string  muonsSrcName;
  string  electronsSrcName;
  string  tausSrcName;
  string  photonsSrcName;
  string  jetsSrcName;
  string  csvBJetTagsSrcName;
  string  cisvBJetTagsSrcName;
  string  jetFlavorMatchSrcName;
  string  jetsAK8SrcName;
  string  pfCandsSrcName;
  string  genParticlesSrcName;
  string  genJetsSrcName;
  string  triggerBitsSrcName;
  string  genMetSrcName;
  string  metSrcName;
  string  metFilterBitsSrcName;
  string  lheInfoSrcName;
  string  genInfoSrcName;
  string  puInfoSrcName;
  string  hcalNoiseInfoSrcName;
  //  string  secondaryVerticesSrcName;
  string  rhoAllSrcName;
  string  rhoFastjetAllSrcName;
  string  rhoFastjetAllCaloSrcName;
  string  rhoFastjetCentralCaloSrcName;
  string  rhoFastjetCentralChargedPileUpSrcName;
  string  rhoFastjetCentralNeutralSrcName;
  string  beamSpotSrcName;
  string  ebRecHitsSrcName;
  string  eeRecHitsSrcName;
  string  esRecHitsSrcName;
  string  ebeeClustersSrcName;
  string  esClustersSrcName;
  string  conversionsSrcName;
  string  singleLegConversionsSrcName;
  string  gedGsfElectronCoresSrcName;
  string  gedPhotonCoresSrcName;
  string  superClustersSrcName;
  

  //EDM handles for each miniAOD input object
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<edm::TriggerResults> metFilterBits;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::PFCandidateCollection> PFCands;
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::PhotonCollection> photons;
  //edm::Handle<reco::TauCollection> taus;
  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::PFJetCollection> jetsAK8;
  edm::Handle<reco::GenMETCollection> genMets;
  edm::Handle<reco::PFMETCollection> mets;
  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<LHEEventProduct> lheInfo;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  edm::Handle<HcalNoiseSummary> hcalNoiseInfo;
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > secondaryVertices;
  edm::Handle<double> rhoAll;
  edm::Handle<double> rhoFastjetCentralChargedPileUp;
  edm::Handle<double> rhoFastjetCentralNeutral;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHits;
  edm::Handle<vector<reco::CaloCluster> > ebeeClusters;
  edm::Handle<vector<reco::CaloCluster> > esClusters;
  edm::Handle<vector<reco::Conversion> > conversions;
  edm::Handle<vector<reco::Conversion>> singleLegConversions;
  edm::Handle<vector<reco::GsfElectronCore> > gedGsfElectronCores;
  edm::Handle<vector<reco::PhotonCore> > gedPhotonCores;
  edm::Handle<vector<reco::SuperCluster> > superClusters;
  edm::Handle<reco::JetTagCollection> csvBJetTags;
  edm::Handle<reco::JetTagCollection> cisvBJetTags;
  edm::Handle<reco::JetFlavourMatchingCollection> jetFlavorMatch;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilderHandle;
  EcalClusterLazyTools *ecalLazyTools;
  EcalClusterLocal *ecalLocal;  
  const reco::Vertex *myPV;

  //output tree
  TTree *RazorEvents;
  TH1F *NEvents;

  //------ Variables for tree ------//

  //PU
  int nBunchXing;
  int BunchXing[99];
  int nPU[99];
  float nPUmean[99];

  //Muons
  int nMuons;
  float muonE[99];
  float muonPt[99];
  float muonEta[99];
  float muonPhi[99];
  int muonCharge[99];//muon charge
  bool muonIsLoose[99];
  bool muonIsTight[99];
  float muon_d0[99];//transverse impact paramenter
  float muon_dZ[99];//impact parameter
  float muon_ip3d[99];//3d impact paramenter
  float muon_ip3dSignificance[99];//3d impact paramenter/error
  unsigned int muonType[99];//muonTypeBit: global, tracker, standalone 
  unsigned int muonQuality[99];//muonID Quality Bits
  float muon_pileupIso[99];
  float muon_chargedIso[99];
  float muon_photonIso[99];
  float muon_neutralHadIso[99];
  float muon_ptrel[99];
  float muon_miniiso[99];

  //Electrons
  int nElectrons;
  float eleE[99];
  float elePt[99];
  float eleEta[99];
  float elePhi[99];
  float eleCharge[99];
  float eleE_SC[99];
  //float SC_ElePt[99]; 
  float eleEta_SC[99];
  float elePhi_SC[99];
  float eleSigmaIetaIeta[99];
  float eleFull5x5SigmaIetaIeta[99];
  float eleR9[99];
  float ele_dEta[99];
  float ele_dPhi[99];
  float ele_HoverE[99];
  float ele_d0[99];
  float ele_dZ[99];
  float ele_pileupIso[99];
  float ele_chargedIso[99];
  float ele_photonIso[99];
  float ele_neutralHadIso[99];
  int ele_MissHits[99];
  bool ele_PassConvVeto[99];
  float ele_OneOverEminusOneOverP[99];
  float ele_IDMVATrig[99];
  float ele_IDMVANonTrig[99];
  float ele_RegressionE[99];
  float ele_CombineP4[99];
  float ele_ptrel[99];
  float ele_miniiso[99];

  //Taus
  int nTaus;
  float tauE[99];
  float tauPt[99];
  float tauEta[99];
  float tauPhi[99];
  bool tau_IsLoose[99];
  bool tau_IsMedium[99];
  bool tau_IsTight[99];
  bool tau_passEleVetoLoose[99];
  bool tau_passEleVetoMedium[99];
  bool tau_passEleVetoTight[99];
  bool tau_passMuVetoLoose[99];
  bool tau_passMuVetoMedium[99];
  bool tau_passMuVetoTight[99];  
  UInt_t tau_ID[99];//tauID Bits
  float tau_combinedIsoDeltaBetaCorr3Hits[99];
  float tau_eleVetoMVA[99];
  int tau_eleVetoCategory[99];
  float tau_muonVetoMVA[99];
  float tau_isoMVAnewDMwLT[99];
  float tau_isoMVAnewDMwoLT[99]; 
  float tau_leadCandPt[99];
  int tau_leadCandID[99];
  float tau_leadChargedHadrCandPt[99];
  int tau_leadChargedHadrCandID[99];

  //IsolatedChargedPFCandidates
  int nIsoPFCandidates;
  float isoPFCandidatePt[99];
  float isoPFCandidateEta[99];
  float isoPFCandidatePhi[99];
  float isoPFCandidateIso04[99];
  float isoPFCandidateD0[99];
  int   isoPFCandidatePdgId[99];

  //Photons
  int nPhotons;
  float phoE[99];
  float phoPt[99];
  float phoEta[99];
  float phoPhi[99];
  float phoSigmaIetaIeta[99];
  float phoFull5x5SigmaIetaIeta[99];
  float phoR9[99];
  float pho_HoverE[99];
  float pho_sumChargedHadronPt[99];
  float pho_sumNeutralHadronEt[99];
  float pho_sumPhotonEt[99];
  float pho_sumWorstVertexChargedHadronPt[99];
  bool  pho_isConversion[99];
  bool  pho_passEleVeto[99];
  float pho_RegressionE[99];
  float pho_RegressionEUncertainty[99];
  float pho_IDMVA[99];
  float pho_superClusterEta[99];
  float pho_superClusterPhi[99];
  bool pho_hasPixelSeed[99];

  //AK4 Jets
  int nJets;
  float jetE[99];
  float jetPt[99];
  float jetEta[99];
  float jetPhi[99];
  float jetCSV[99];
  float jetCISV[99];
  float jetMass[99];
  float jetJetArea[99];
  float jetPileupE[99];  
  float jetPileupId[99];
  int   jetPileupIdFlag[99];
  bool  jetPassIDLoose[99];
  bool  jetPassIDTight[99];
  int   jetPartonFlavor[99];
  int   jetHadronFlavor[99];

  //AK8 Jets
  int nFatJets;
  float fatJetE[99];
  float fatJetPt[99];
  float fatJetEta[99];
  float fatJetPhi[99];
  float fatJetTrimmedM[99];
  float fatJetPrunedM[99];
  float fatJetFilteredM[99];
  float fatJetTau1[99];
  float fatJetTau2[99];
  float fatJetTau3[99];

  //MET 
  float metPt;
  float metPhi;
  float sumMET;
  float UncMETdpx;
  float UncMETdpy;
  float UncMETdSumEt;
  bool Flag_HBHENoiseFilter;
  bool Flag_CSCTightHaloFilter;
  bool Flag_hcalLaserEventFilter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_EcalDeadCellBoundaryEnergyFilter;
  bool Flag_goodVertices;
  bool Flag_trackingFailureFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalLaserCorrFilter;
  bool Flag_trkPOGFilters;  
  bool Flag_trkPOG_manystripclus53X;
  bool Flag_trkPOG_toomanystripclus53X;
  bool Flag_trkPOG_logErrorTooManyClusters;
  bool Flag_METFilters;

  bool Flag_EcalDeadCellEvent;
  bool Flag_IsNotDeadEcalCluster;
  bool Flag_EcalDeadDR;
  bool Flag_EcalBoundaryDR;

  //MC
  int nGenJets;
  float genJetE[99];
  float genJetPt[99];
  float genJetEta[99];
  float genJetPhi[99];
  float genMetPt;
  float genMetPhi;
  float genVertexX;
  float genVertexY;
  float genVertexZ;
  float genWeight;
  unsigned int genSignalProcessID;
  float genQScale;
  float genAlphaQCD;
  float genAlphaQED;

  //gen info
  int nGenParticle;
  int gParticleMotherId[500];
  int gParticleMotherIndex[500];
  int gParticleId[500];
  int gParticleStatus[500];
  float gParticleE[500];
  float gParticlePt[500];
  float gParticleEta[500];
  float gParticlePhi[500];

  //razor variables
  float MR, RSQ;
  float MR_AK8, RSQ_AK8;
  
  //event info
  int nPV;
  int runNum;
  int lumiNum;
  int eventNum;
  float pvX;
  float pvY;
  float pvZ;
  float fixedGridRhoAll;
  float fixedGridRhoFastjetAll;
  float fixedGridRhoFastjetAllCalo;
  float fixedGridRhoFastjetCentralCalo;
  float fixedGridRhoFastjetCentralChargedPileUp;
  float fixedGridRhoFastjetCentralNeutral;

  //trigger info
  vector<string>  *nameHLT;
  static const int NTriggersMAX = 100;
  string triggerPathNames[NTriggersMAX];
  bool triggerDecision[NTriggersMAX];

};

#endif
