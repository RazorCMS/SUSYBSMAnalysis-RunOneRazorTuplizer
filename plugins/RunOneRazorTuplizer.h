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

//------ Array Size Constants ------//
#define OBJECTARRAYSIZE 99
#define GENPARTICLEARRAYSIZE 500

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
			    const reco::Vertex *PV, 
			    edm::Handle<reco::VertexCollection> vertices,
			    double r_iso_min = 0.05, double r_iso_max = 0.2 , double kt_scale = 10.0,
			    bool use_pfweight = false, bool charged_only = false);  
  
  //Get photon 4-momentum after primary vertex correction
  TLorentzVector photonP4FromVtx( TVector3 vtx, TVector3 phoPos, double E );
  
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
  bool    isData_;
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
  edm::Handle<reco::PFMETCollection> type0Mets;
  edm::Handle<reco::PFMETCollection> type0Plus1Mets;
  edm::Handle<reco::PFMETCollection> type1Mets;
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
  int BunchXing[OBJECTARRAYSIZE];
  int nPU[OBJECTARRAYSIZE];
  float nPUmean[OBJECTARRAYSIZE];

  //Muons
  int nMuons;
  float muonE[OBJECTARRAYSIZE];
  float muonPt[OBJECTARRAYSIZE];
  float muonEta[OBJECTARRAYSIZE];
  float muonPhi[OBJECTARRAYSIZE];
  int muonCharge[OBJECTARRAYSIZE];//muon charge
  bool muonIsLoose[OBJECTARRAYSIZE];
  bool muonIsTight[OBJECTARRAYSIZE];
  float muon_d0[OBJECTARRAYSIZE];//transverse impact paramenter
  float muon_dZ[OBJECTARRAYSIZE];//impact parameter
  float muon_ip3d[OBJECTARRAYSIZE];//3d impact paramenter
  float muon_ip3dSignificance[OBJECTARRAYSIZE];//3d impact paramenter/error
  unsigned int muonType[OBJECTARRAYSIZE];//muonTypeBit: global, tracker, standalone 
  unsigned int muonQuality[OBJECTARRAYSIZE];//muonID Quality Bits
  float muon_pileupIso[OBJECTARRAYSIZE];
  float muon_chargedIso[OBJECTARRAYSIZE];
  float muon_photonIso[OBJECTARRAYSIZE];
  float muon_neutralHadIso[OBJECTARRAYSIZE];
  float muon_ptrel[OBJECTARRAYSIZE];
  float muon_miniiso[OBJECTARRAYSIZE];

  //Electrons
  int nElectrons;
  float eleE[OBJECTARRAYSIZE];
  float elePt[OBJECTARRAYSIZE];
  float eleEta[OBJECTARRAYSIZE];
  float elePhi[OBJECTARRAYSIZE];
  float eleCharge[OBJECTARRAYSIZE];
  float eleE_SC[OBJECTARRAYSIZE];
  //float SC_ElePt[OBJECTARRAYSIZE]; 
  float eleEta_SC[OBJECTARRAYSIZE];
  float elePhi_SC[OBJECTARRAYSIZE];
  float eleSigmaIetaIeta[OBJECTARRAYSIZE];
  float eleFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
  float eleR9[OBJECTARRAYSIZE];
  float ele_dEta[OBJECTARRAYSIZE];
  float ele_dPhi[OBJECTARRAYSIZE];
  float ele_HoverE[OBJECTARRAYSIZE];
  float ele_d0[OBJECTARRAYSIZE];
  float ele_dZ[OBJECTARRAYSIZE];
  float ele_pileupIso[OBJECTARRAYSIZE];
  float ele_chargedIso[OBJECTARRAYSIZE];
  float ele_photonIso[OBJECTARRAYSIZE];
  float ele_neutralHadIso[OBJECTARRAYSIZE];
  int ele_MissHits[OBJECTARRAYSIZE];
  bool ele_PassConvVeto[OBJECTARRAYSIZE];
  float ele_OneOverEminusOneOverP[OBJECTARRAYSIZE];
  float ele_IDMVATrig[OBJECTARRAYSIZE];
  float ele_IDMVANonTrig[OBJECTARRAYSIZE];
  float ele_RegressionE[OBJECTARRAYSIZE];
  float ele_CombineP4[OBJECTARRAYSIZE];
  float ele_ptrel[OBJECTARRAYSIZE];
  float ele_miniiso[OBJECTARRAYSIZE];

  //Taus
  int nTaus;
  float tauE[OBJECTARRAYSIZE];
  float tauPt[OBJECTARRAYSIZE];
  float tauEta[OBJECTARRAYSIZE];
  float tauPhi[OBJECTARRAYSIZE];
  bool tau_IsLoose[OBJECTARRAYSIZE];
  bool tau_IsMedium[OBJECTARRAYSIZE];
  bool tau_IsTight[OBJECTARRAYSIZE];
  bool tau_passEleVetoLoose[OBJECTARRAYSIZE];
  bool tau_passEleVetoMedium[OBJECTARRAYSIZE];
  bool tau_passEleVetoTight[OBJECTARRAYSIZE];
  bool tau_passMuVetoLoose[OBJECTARRAYSIZE];
  bool tau_passMuVetoMedium[OBJECTARRAYSIZE];
  bool tau_passMuVetoTight[OBJECTARRAYSIZE];  
  UInt_t tau_ID[OBJECTARRAYSIZE];//tauID Bits
  float tau_combinedIsoDeltaBetaCorr3Hits[OBJECTARRAYSIZE];
  float tau_eleVetoMVA[OBJECTARRAYSIZE];
  int tau_eleVetoCategory[OBJECTARRAYSIZE];
  float tau_muonVetoMVA[OBJECTARRAYSIZE];
  float tau_isoMVAnewDMwLT[OBJECTARRAYSIZE];
  float tau_isoMVAnewDMwoLT[OBJECTARRAYSIZE]; 
  float tau_leadCandPt[OBJECTARRAYSIZE];
  int tau_leadCandID[OBJECTARRAYSIZE];
  float tau_leadChargedHadrCandPt[OBJECTARRAYSIZE];
  int tau_leadChargedHadrCandID[OBJECTARRAYSIZE];

  //IsolatedChargedPFCandidates
  int nIsoPFCandidates;
  float isoPFCandidatePt[OBJECTARRAYSIZE];
  float isoPFCandidateEta[OBJECTARRAYSIZE];
  float isoPFCandidatePhi[OBJECTARRAYSIZE];
  float isoPFCandidateIso04[OBJECTARRAYSIZE];
  float isoPFCandidateD0[OBJECTARRAYSIZE];
  int   isoPFCandidatePdgId[OBJECTARRAYSIZE];

  //Photons
  int nPhotons;
  float phoE[OBJECTARRAYSIZE];
  float phoPt[OBJECTARRAYSIZE];
  float phoEta[OBJECTARRAYSIZE];
  float phoPhi[OBJECTARRAYSIZE];
  float phoSigmaIetaIeta[OBJECTARRAYSIZE];
  float phoFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
  float phoR9[OBJECTARRAYSIZE];
  float pho_HoverE[OBJECTARRAYSIZE];
  float pho_sumChargedHadronPt[OBJECTARRAYSIZE];
  float pho_sumNeutralHadronEt[OBJECTARRAYSIZE];
  float pho_sumPhotonEt[OBJECTARRAYSIZE];
  float pho_sumWorstVertexChargedHadronPt[OBJECTARRAYSIZE];
  bool  pho_isConversion[OBJECTARRAYSIZE];
  bool  pho_passEleVeto[OBJECTARRAYSIZE];
  float pho_RegressionE[OBJECTARRAYSIZE];
  float pho_RegressionEUncertainty[OBJECTARRAYSIZE];
  float pho_IDMVA[OBJECTARRAYSIZE];
  float pho_superClusterEta[OBJECTARRAYSIZE];
  float pho_superClusterPhi[OBJECTARRAYSIZE];
  bool pho_hasPixelSeed[OBJECTARRAYSIZE];

  //AK4 Jets
  int nJets;
  float jetE[OBJECTARRAYSIZE];
  float jetPt[OBJECTARRAYSIZE];
  float jetEta[OBJECTARRAYSIZE];
  float jetPhi[OBJECTARRAYSIZE];
  float jetCSV[OBJECTARRAYSIZE];
  float jetCISV[OBJECTARRAYSIZE];
  float jetMass[OBJECTARRAYSIZE];
  float jetJetArea[OBJECTARRAYSIZE];
  float jetPileupE[OBJECTARRAYSIZE];  
  float jetPileupId[OBJECTARRAYSIZE];
  int   jetPileupIdFlag[OBJECTARRAYSIZE];
  bool  jetPassIDLoose[OBJECTARRAYSIZE];
  bool  jetPassIDTight[OBJECTARRAYSIZE];
  bool  jetPassMuFrac[OBJECTARRAYSIZE];
  bool  jetPassEleFrac[OBJECTARRAYSIZE];
  int   jetPartonFlavor[OBJECTARRAYSIZE];
  int   jetHadronFlavor[OBJECTARRAYSIZE];

  //AK8 Jets
  int nFatJets;
  float fatJetE[OBJECTARRAYSIZE];
  float fatJetPt[OBJECTARRAYSIZE];
  float fatJetEta[OBJECTARRAYSIZE];
  float fatJetPhi[OBJECTARRAYSIZE];
  float fatJetTrimmedM[OBJECTARRAYSIZE];
  float fatJetPrunedM[OBJECTARRAYSIZE];
  float fatJetFilteredM[OBJECTARRAYSIZE];
  float fatJetTau1[OBJECTARRAYSIZE];
  float fatJetTau2[OBJECTARRAYSIZE];
  float fatJetTau3[OBJECTARRAYSIZE];

  //MET 
  float metPt;
  float metPhi;
  float sumMET;
  float metType0Pt;
  float metType0Phi;
  float metType1Pt;
  float metType1Phi;
  float metType0Plus1Pt;
  float metType0Plus1Phi;
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
  float genJetE[OBJECTARRAYSIZE];
  float genJetPt[OBJECTARRAYSIZE];
  float genJetEta[OBJECTARRAYSIZE];
  float genJetPhi[OBJECTARRAYSIZE];
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
  int gParticleMotherId[GENPARTICLEARRAYSIZE];
  int gParticleMotherIndex[GENPARTICLEARRAYSIZE];
  int gParticleId[GENPARTICLEARRAYSIZE];
  int gParticleStatus[GENPARTICLEARRAYSIZE];
  float gParticleE[GENPARTICLEARRAYSIZE];
  float gParticlePt[GENPARTICLEARRAYSIZE];
  float gParticleEta[GENPARTICLEARRAYSIZE];
  float gParticlePhi[GENPARTICLEARRAYSIZE];

  //razor variables
  float MR, RSQ;
  float MR_AK8, RSQ_AK8;
  
  //event info
  bool isData;
  int nPV;
  uint runNum;
  uint lumiNum;
  uint eventNum;
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
