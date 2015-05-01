// -*- C++ -*-
// Class:      RazorTuplizer
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RunOneRazorTuplizer.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include <PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h>
#include <PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h>
#include "DataFormats/METReco/interface/AnomalousECALVariables.h"

//------ Constructors and destructor ------//
RazorTuplizer::RazorTuplizer(const edm::ParameterSet& iConfig): 
    //get inputs from config file  
    isData_(iConfig.getParameter<bool> ("isData")),
    useGen_(iConfig.getParameter<bool> ("useGen")),
    enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
    triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),  
    verticesSrcName(iConfig.getParameter<string> ("vertices")),
    muonsSrcName(iConfig.getParameter<string> ("muons")),
    electronsSrcName(iConfig.getParameter<string> ("electrons")),
    //tausSrcName(iConfig.getParameter<string> ("taus")),
    photonsSrcName(iConfig.getParameter<string> ("photons")),
    jetsSrcName(iConfig.getParameter<string> ("jets")),
    csvBJetTagsSrcName(iConfig.getParameter<string> ("csvBJetTags")),
    cisvBJetTagsSrcName(iConfig.getParameter<string> ("cisvBJetTags")),
    jetFlavorMatchSrcName(iConfig.getParameter<string> ("jetFlavorMatch")),
    jetsAK8SrcName(iConfig.getParameter<string> ("jetsAK8")),
    pfCandsSrcName(iConfig.getParameter<string> ("pfCands")),
    genParticlesSrcName(iConfig.getParameter<string> ("genParticles")),
    genJetsSrcName(iConfig.getParameter<string> ("genJets")),
    triggerBitsSrcName(iConfig.getParameter<string> ("triggerBits")),
    genMetSrcName(iConfig.getParameter<string> ("genMets")),
    metSrcName(iConfig.getParameter<string> ("mets")),
    lheInfoSrcName(iConfig.getParameter<string> ("lheInfo")),
    genInfoSrcName(iConfig.getParameter<string> ("genInfo")),
    puInfoSrcName(iConfig.getParameter<string> ("puInfo")),
    beamSpotSrcName(iConfig.getParameter<string> ("beamSpot")),
    ebRecHitsSrcName(iConfig.getParameter<string> ("ebRecHits")),
    eeRecHitsSrcName(iConfig.getParameter<string> ("eeRecHits")),
    conversionsSrcName(iConfig.getParameter<string> ("conversions"))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;
  
  //set up output tree
  RazorEvents = fs->make<TTree>("RazorEvents", "selected miniAOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);

  //set up electron MVA ID
  std::vector<std::string> myTrigWeights;
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat1.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat2.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat3.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat4.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat5.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_TrigV0_Cat6.weights.xml").fullPath().c_str());

  myMVATrig = new EGammaMvaEleEstimator();
  myMVATrig->initialize("BDT",
  			EGammaMvaEleEstimator::kTrig,
  			true,
  			myTrigWeights);

  std::vector<std::string> myNonTrigWeights;
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml").fullPath().c_str());
  
  myMVANonTrig = new EGammaMvaEleEstimator();
  myMVANonTrig->initialize("BDT",
			   EGammaMvaEleEstimator::kNonTrig,
			   true,
			   myNonTrigWeights);



  //Read in HLT Trigger Path List from config file
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open()) {
    string line;
    int index;
    string hltpathname;
    
    while(myfile>>index>>hltpathname) {
      if (index < NTriggersMAX) {
  	triggerPathNames[index] = hltpathname;
      }    
    }    
    myfile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
  
  if(enableTriggerInfo_) {
    cout << "\n";
    cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";    
    for (int i = 0; i<NTriggersMAX; ++i) {
      if (triggerPathNames[i] != "") cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
    }
    cout << "****************************************************************************\n";    
    cout << "\n";
  }

}

RazorTuplizer::~RazorTuplizer()
{
}

//------ Enable the desired set of branches ------//
void RazorTuplizer::setBranches(){
  enableEventInfoBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enableIsoPFCandidateBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  
  if (enableTriggerInfo_) enableTriggerBranches();
  enableMCBranches();
  enableGenParticleBranches();
}

void RazorTuplizer::enableEventInfoBranches(){
  RazorEvents->Branch("isData", &isData, "isData/O");
  RazorEvents->Branch("nPV", &nPV, "nPV/I");
  RazorEvents->Branch("runNum", &runNum, "runNum/I");
  RazorEvents->Branch("lumiNum", &lumiNum, "lumiNum/I");
  RazorEvents->Branch("eventNum", &eventNum, "eventNum/I");
  RazorEvents->Branch("pvX", &pvX, "pvX/F");
  RazorEvents->Branch("pvY", &pvY, "pvY/F");
  RazorEvents->Branch("pvZ", &pvZ, "pvZ/F");
  RazorEvents->Branch("fixedGridRhoAll", &fixedGridRhoAll, "fixedGridRhoAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral/F");
}

void RazorTuplizer::enablePileUpBranches(){
  RazorEvents->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  RazorEvents->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  RazorEvents->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  RazorEvents->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void RazorTuplizer::enableMuonBranches(){
  RazorEvents->Branch("nMuons", &nMuons,"nMuons/I");
  RazorEvents->Branch("muonE", muonE,"muonE[nMuons]/F");
  RazorEvents->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  RazorEvents->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  RazorEvents->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  RazorEvents->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  RazorEvents->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  RazorEvents->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  RazorEvents->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  RazorEvents->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  RazorEvents->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  RazorEvents->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  RazorEvents->Branch("muonType", muonType, "muonType[nMuons]/i");
  RazorEvents->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  RazorEvents->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  RazorEvents->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  RazorEvents->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  RazorEvents->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");
  RazorEvents->Branch("muon_ptrel", muon_ptrel, "muon_ptrel[nMuons]/F");
  RazorEvents->Branch("muon_miniiso", muon_miniiso, "muon_miniiso[nMuons]/F");
}

void RazorTuplizer::enableElectronBranches(){
  RazorEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
  RazorEvents->Branch("eleE", eleE,"eleE[nElectrons]/F");
  RazorEvents->Branch("elePt", elePt,"elePt[nElectrons]/F");
  RazorEvents->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  RazorEvents->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
  RazorEvents->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  //RazorEvents->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");
  RazorEvents->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  //RazorEvents->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");
  RazorEvents->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  RazorEvents->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  RazorEvents->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  RazorEvents->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  RazorEvents->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  RazorEvents->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  RazorEvents->Branch("ele_pileupIso", ele_pileupIso, "ele_pileupIso[nElectrons]/F");
  RazorEvents->Branch("ele_chargedIso", ele_chargedIso, "ele_chargedIso[nElectrons]/F");
  RazorEvents->Branch("ele_photonIso", ele_photonIso, "ele_photonIso[nElectrons]/F");
  RazorEvents->Branch("ele_neutralHadIso", ele_neutralHadIso, "ele_neutralHadIso[nElectrons]/F");
  RazorEvents->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  RazorEvents->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  RazorEvents->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVATrig", ele_IDMVATrig, "ele_IDMVATrig[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVANonTrig", ele_IDMVANonTrig, "ele_IDMVANonTrig[nElectrons]/F");
  RazorEvents->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  RazorEvents->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
  RazorEvents->Branch("ele_ptrel", ele_ptrel, "ele_ptrel[nElectrons]/F");
  RazorEvents->Branch("ele_miniiso", ele_miniiso, "ele_miniiso[nElectrons]/F");
}

void RazorTuplizer::enableTauBranches(){
  RazorEvents->Branch("nTaus", &nTaus,"nTaus/I");
  RazorEvents->Branch("tauE", tauE,"tauE[nTaus]/F");
  RazorEvents->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  RazorEvents->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  RazorEvents->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
  RazorEvents->Branch("tau_IsLoose", tau_IsLoose, "tau_IsLoose[nTaus]/O");
  RazorEvents->Branch("tau_IsMedium", tau_IsMedium, "tau_IsMedium[nTaus]/O");
  RazorEvents->Branch("tau_IsTight", tau_IsTight, "tau_IsTight[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoLoose", tau_passEleVetoLoose, "tau_passEleVetoLoose[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoMedium", tau_passEleVetoMedium, "tau_passEleVetoMedium[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoTight", tau_passEleVetoTight, "tau_passEleVetoTight[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoLoose", tau_passMuVetoLoose, "tau_passMuVetoLoose[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoMedium", tau_passMuVetoMedium, "tau_passMuVetoMedium[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoTight", tau_passMuVetoTight, "tau_passMuVetoTight[nTaus]/O");
  RazorEvents->Branch("tau_ID", tau_ID, "tau_ID[nTaus]/i");
  RazorEvents->Branch("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, "tau_combinedIsoDeltaBetaCorr3Hits[nTaus]/F");
  RazorEvents->Branch("tau_eleVetoMVA", tau_eleVetoMVA, "tau_eleVetoMVA[nTaus]/F");
  RazorEvents->Branch("tau_eleVetoCategory", tau_eleVetoCategory, "tau_eleVetoCategory[nTaus]/I");
  RazorEvents->Branch("tau_muonVetoMVA", tau_muonVetoMVA, "tau_muonVetoMVA[nTaus]/F");
  RazorEvents->Branch("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, "tau_isoMVAnewDMwLT[nTaus]/F");
  RazorEvents->Branch("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, "tau_isoMVAnewDMwoLT[nTaus]/F");
  RazorEvents->Branch("tau_leadCandPt", tau_leadCandPt, "tau_leadCandPt[nTaus]/F");
  RazorEvents->Branch("tau_leadCandID", tau_leadCandID, "tau_leadCandID[nTaus]/I");
  RazorEvents->Branch("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, "tau_leadChargedHadrCandPt[nTaus]/F");
  RazorEvents->Branch("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, "tau_leadChargedHadrCandID[nTaus]/I"); 
}

void RazorTuplizer::enableIsoPFCandidateBranches(){
  RazorEvents->Branch("nIsoPFCandidates", &nIsoPFCandidates, "nIsoPFCandidates/i");
  RazorEvents->Branch("isoPFCandidatePt", isoPFCandidatePt, "isoPFCandidatePt[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateEta", isoPFCandidateEta, "isoPFCandidateEta[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidatePhi", isoPFCandidatePhi, "isoPFCandidatePhi[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateIso04", isoPFCandidateIso04, "isoPFCandidateIso04[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateD0", isoPFCandidateD0, "isoPFCandidateD0[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidatePdgId", isoPFCandidatePdgId, "isoPFCandidatePdgId[nIsoPFCandidates]/I");  
}

void RazorTuplizer::enablePhotonBranches(){
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("phoE", phoE,"phoE[nPhotons]/F");
  RazorEvents->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  RazorEvents->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  RazorEvents->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  RazorEvents->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  RazorEvents->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  RazorEvents->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, "pho_sumWorstVertexChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  RazorEvents->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  RazorEvents->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  RazorEvents->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  RazorEvents->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  RazorEvents->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
}

void RazorTuplizer::enableJetBranches(){
  RazorEvents->Branch("nJets", &nJets,"nJets/I");
  RazorEvents->Branch("jetE", jetE,"jetE[nJets]/F");
  RazorEvents->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  RazorEvents->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  RazorEvents->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  RazorEvents->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  RazorEvents->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  RazorEvents->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  RazorEvents->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  RazorEvents->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  RazorEvents->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  RazorEvents->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  RazorEvents->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  RazorEvents->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  RazorEvents->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  RazorEvents->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  RazorEvents->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  RazorEvents->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
}

void RazorTuplizer::enableJetAK8Branches(){
  RazorEvents->Branch("nFatJets", &nFatJets,"nFatJets/i");
  RazorEvents->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  RazorEvents->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  RazorEvents->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  RazorEvents->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
  RazorEvents->Branch("fatJetPrunedM", fatJetPrunedM,"fatJetPrunedM[nFatJets]/F");
  RazorEvents->Branch("fatJetTrimmedM", fatJetTrimmedM,"fatJetTrimmedM[nFatJets]/F");
  RazorEvents->Branch("fatJetFilteredM", fatJetFilteredM,"fatJetFilteredM[nFatJets]/F");
  RazorEvents->Branch("fatJetTau1", fatJetTau1,"fatJetTau1[nFatJets]/F");
  RazorEvents->Branch("fatJetTau2", fatJetTau2,"fatJetTau2[nFatJets]/F");
  RazorEvents->Branch("fatJetTau3", fatJetTau3,"fatJetTau3[nFatJets]/F");
}

void RazorTuplizer::enableMetBranches(){
  RazorEvents->Branch("metPt", &metPt, "metPt/F");
  RazorEvents->Branch("metPhi", &metPhi, "metPhi/F");
  RazorEvents->Branch("sumMET", &sumMET, "sumMET/F");
  RazorEvents->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  RazorEvents->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  RazorEvents->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  RazorEvents->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  RazorEvents->Branch("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, "Flag_EcalDeadCellBoundaryEnergyFilter/O");
  RazorEvents->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  RazorEvents->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  RazorEvents->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  RazorEvents->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  RazorEvents->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  RazorEvents->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  RazorEvents->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  RazorEvents->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  RazorEvents->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");  
  RazorEvents->Branch("Flag_EcalDeadCellEvent", &Flag_EcalDeadCellEvent, "Flag_EcalDeadCellEvent/O");  
  RazorEvents->Branch("Flag_IsNotDeadEcalCluster", &Flag_IsNotDeadEcalCluster, "Flag_IsNotDeadEcalCluster/O");  
  RazorEvents->Branch("Flag_EcalDeadDR", &Flag_EcalDeadDR, "Flag_EcalDeadDR/O");  
  RazorEvents->Branch("Flag_EcalBoundaryDR", &Flag_EcalBoundaryDR, "Flag_EcalBoundaryDR/O");  
}

void RazorTuplizer::enableRazorBranches(){
  RazorEvents->Branch("MR", &MR, "MR/F");
  RazorEvents->Branch("RSQ", &RSQ, "RSQ/F");
  RazorEvents->Branch("MR_AK8", &MR_AK8, "MR_AK8/F");
  RazorEvents->Branch("RSQ_AK8", &RSQ_AK8, "RSQ_AK8/F");
}

void RazorTuplizer::enableTriggerBranches(){
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  //RazorEvents->Branch("nameHLT", "std::vector<std::string>",&nameHLT);
  RazorEvents->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
}

void RazorTuplizer::enableMCBranches(){
  RazorEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
  RazorEvents->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  RazorEvents->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  RazorEvents->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  RazorEvents->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  RazorEvents->Branch("genMetPt", &genMetPt, "genMetPt/F");
  RazorEvents->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
  RazorEvents->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  RazorEvents->Branch("genWeight", &genWeight, "genWeight/F");
  RazorEvents->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  RazorEvents->Branch("genQScale", &genQScale, "genQScale/F");
  RazorEvents->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  RazorEvents->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
}

void RazorTuplizer::enableGenParticleBranches(){
  RazorEvents->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  RazorEvents->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  RazorEvents->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  RazorEvents->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  RazorEvents->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  RazorEvents->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  RazorEvents->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  RazorEvents->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  RazorEvents->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
}

//------ Load the miniAOD objects and reset tree variables for each event ------//
void RazorTuplizer::loadEvent(const edm::Event& iEvent){
  //load all miniAOD objects for the current event
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),triggerBits);
  iEvent.getByLabel(verticesSrcName, vertices);
  iEvent.getByLabel(pfCandsSrcName, PFCands);
  iEvent.getByLabel(muonsSrcName, muons);
  iEvent.getByLabel(electronsSrcName, electrons);
  iEvent.getByLabel(photonsSrcName, photons);
  //iEvent.getByLabel(tausSrcName, taus);
  iEvent.getByLabel(jetsSrcName, jets);
  iEvent.getByLabel(csvBJetTagsSrcName, csvBJetTags);
  //iEvent.getByLabel(cisvBJetTagsSrcName, cisvBJetTags);
  iEvent.getByLabel(jetsAK8SrcName, jetsAK8);
  iEvent.getByLabel(metSrcName, mets);
  iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoAll);
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralChargedPileUp","rho"),rhoFastjetCentralChargedPileUp);
  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutral","rho"),rhoFastjetCentralNeutral);
  iEvent.getByLabel(beamSpotSrcName,beamSpot);
  iEvent.getByLabel(ebRecHitsSrcName,ebRecHits);
  iEvent.getByLabel(eeRecHitsSrcName,eeRecHits);
  iEvent.getByLabel(conversionsSrcName,conversions);
  
  if (useGen_) {
    iEvent.getByLabel(jetFlavorMatchSrcName, jetFlavorMatch);  
    iEvent.getByLabel(genParticlesSrcName,genParticles);
    iEvent.getByLabel(genJetsSrcName,genJets);
    iEvent.getByLabel(lheInfoSrcName, lheInfo);
    iEvent.getByLabel(genInfoSrcName,genInfo);
    iEvent.getByLabel(puInfoSrcName,puInfo);
    iEvent.getByLabel(genMetSrcName, genMets);
  }

}   

//called by the loadEvent() method
void RazorTuplizer::resetBranches(){
    //reset tree variables
    nBunchXing = 0;
    nMuons = 0;
    nElectrons = 0;
    nTaus = 0;
    nPhotons = 0;
    nJets = 0;
    nFatJets = 0;
    nGenJets = 0;
    nGenParticle = 0;

    //nameHLT->clear();

    for(int i = 0; i < 99; i++){
        //PU
        BunchXing[i] = -99;
        nPU[i] = -99;
        nPUmean[i] = -99.0;

        //Muon
        muonE[i] = 0.0;
        muonPt[i] = 0.0;
        muonEta[i] = 0.0;
        muonPhi[i] = 0.0;
        muonCharge[i] = -99;
        muonIsLoose[i] = false;
        muonIsTight[i] = false;
        muon_d0[i] = -99.0;
        muon_dZ[i] = -99.0;
        muon_ip3d[i] = -99.0;
        muon_ip3dSignificance[i] = -99.0;
        muonType[i] = 0;
        muonQuality[i] = 0;
        muon_pileupIso[i] = -99.0;
        muon_chargedIso[i] = -99.0;
        muon_photonIso[i] = -99.0;
        muon_neutralHadIso[i] = -99.0;
	muon_ptrel[i] = -99.0;
	muon_miniiso[i] = -99.0;

        //Electron
        eleE[i] = 0.0;
        elePt[i] = 0.0;
        eleEta[i] = 0.0;
        elePhi[i] = 0.0;
        eleE_SC[i] = -99.0;
        eleEta_SC[i] = -99.0;
        elePhi_SC[i] = -99.0;
        eleSigmaIetaIeta[i] = -99.0;
        eleFull5x5SigmaIetaIeta[i] = -99.0;
        eleR9[i] = -99;
        ele_dEta[i] = -99;
        ele_dPhi[i] = -99;
        ele_HoverE[i] = -99;
        ele_d0[i] = -99;
        ele_dZ[i] = -99;
	ele_pileupIso[i] = -99.0;
        ele_chargedIso[i] = -99.0;
        ele_photonIso[i] = -99.0;
        ele_neutralHadIso[i] = -99.0;
	ele_MissHits[i] = -99;
        ele_PassConvVeto[i] = false;
        ele_OneOverEminusOneOverP[i] = -99.0;
        ele_IDMVATrig[i] = -99.0;
        ele_IDMVANonTrig[i] = -99.0;
        ele_RegressionE[i] = -99.0;
        ele_CombineP4[i] = -99.0;
	ele_ptrel[i] = -99.0;
	ele_miniiso[i] = -99.0;

        //Tau
        tauE[i] = 0.0;
        tauPt[i] = 0.0;
        tauEta[i] = 0.0;
        tauPhi[i] = 0.0;
        tau_IsLoose[i] = false;
        tau_IsMedium[i] = false;
        tau_IsTight[i] = false;
        tau_passEleVetoLoose[i] = false;
        tau_passEleVetoMedium[i] = false;
        tau_passEleVetoTight[i] = false;
        tau_passMuVetoLoose[i] = false;
        tau_passMuVetoMedium[i] = false;
        tau_passMuVetoTight[i] = false;
        tau_ID[i] = 0;
        tau_combinedIsoDeltaBetaCorr3Hits[i] = -99.0;
        tau_eleVetoMVA[i] = -99.0;
        tau_eleVetoCategory[i] = -1;
        tau_muonVetoMVA[i] = -99.0;
        tau_isoMVAnewDMwLT[i] = -99.0;
        tau_isoMVAnewDMwoLT[i] = -99.0;
        tau_leadCandPt[i] = -99.0;
        tau_leadCandID[i] = 0;
        tau_leadChargedHadrCandPt[i] = -99.0;
        tau_leadChargedHadrCandID[i] = 0;

        //IsoPFCandidates
        nIsoPFCandidates = 0;
        isoPFCandidatePt[i] = -99.0;
        isoPFCandidateEta[i] = -99.0;
        isoPFCandidatePhi[i] = -99.0;
        isoPFCandidateIso04[i] = -99.0;
        isoPFCandidateD0[i] = -99.0;
        isoPFCandidatePdgId[i] = 0;

        //Photon
        phoE[i] = 0.0;
        phoPt[i] = 0.0;
        phoEta[i] = 0.0;
        phoPhi[i] = 0.0;
        phoSigmaIetaIeta[i] = -99.0;
        phoFull5x5SigmaIetaIeta[i] = -99.0;
        phoR9[i] = -99.0;
        pho_HoverE[i] = -99.0;
        pho_sumChargedHadronPt[i] = -99.0;
        pho_sumNeutralHadronEt[i] = -99.0;
        pho_sumPhotonEt[i] = -99.0;
	pho_sumWorstVertexChargedHadronPt[i] = -99.0;
        pho_isConversion[i] = false;
        pho_passEleVeto[i] = false;    
        pho_RegressionE[i] = -99.0;
        pho_RegressionEUncertainty[i] = -99.0;
        pho_IDMVA[i] = -99.0;
        pho_superClusterEta[i] = -99.0;
        pho_superClusterPhi[i] = -99.0;
        pho_hasPixelSeed[i] = false;

        //Jet
        jetE[i] = 0.0;
        jetPt[i] = 0.0;
        jetEta[i] = 0.0;
        jetPhi[i] = 0.0;
        jetCSV[i] = 0.0;
        jetCISV[i] = 0.0;
        jetMass[i] =  -99.0;
        jetJetArea[i] = -99.0;
        jetPileupE[i] = -99.0;
        jetPileupId[i] = -99.0;
        jetPartonFlavor[i] = 0;
        jetHadronFlavor[i] = 0;
	jetPassIDLoose[i] = false;
	jetPassIDTight[i] = false;
	jetPassMuFrac[i] = false;
	jetPassEleFrac[i] = false;

        //AK8 Jet
        fatJetE[i] = 0.0;
        fatJetPt[i] = 0.0;
        fatJetEta[i] = 0.0;
        fatJetPhi[i] = 0.0;
        fatJetPrunedM[i] = 0.0;
        fatJetTrimmedM[i] = 0.0;
        fatJetFilteredM[i] = 0.0;
        fatJetTau1[i] = 0.0;
        fatJetTau2[i] = 0.0;
        fatJetTau3[i] = 0.0;

        genJetE[i] = 0.0;
        genJetPt[i] = 0.0;
        genJetEta[i] = 0.0;
        genJetPhi[i] = 0.0;
    }

    for(int i = 0; i < 500; i++){
        //Gen Particle
        gParticleMotherId[i] = -99999;
        gParticleMotherIndex[i] = -99999;
        gParticleId[i] = -99999;
        gParticleStatus[i] = -99999;
        gParticleE[i] = -99999.0;
        gParticlePt[i] = -99999.0;
        gParticleEta[i] = -99999.0;
        gParticlePhi[i] = -99999.0;

    }

    //MET
    metPt = -999;
    metPhi = -999;
    sumMET = -99.0;
    UncMETdpx = -99.0;
    UncMETdpy = -99.0;
    UncMETdSumEt = -99.0;
    Flag_HBHENoiseFilter = false;
    Flag_CSCTightHaloFilter = false;
    Flag_hcalLaserEventFilter = false;
    Flag_EcalDeadCellTriggerPrimitiveFilter = false;
    Flag_EcalDeadCellBoundaryEnergyFilter = false;
    Flag_goodVertices = false;
    Flag_trackingFailureFilter = false;
    Flag_eeBadScFilter = false;
    Flag_ecalLaserCorrFilter = false;
    Flag_trkPOGFilters = false;  
    Flag_trkPOG_manystripclus53X = false;
    Flag_trkPOG_toomanystripclus53X = false;
    Flag_trkPOG_logErrorTooManyClusters = false;
    Flag_METFilters = false;
    Flag_EcalDeadCellEvent = false;
    Flag_IsNotDeadEcalCluster = false;
    Flag_EcalDeadDR = false;
    Flag_EcalBoundaryDR = false;

    genMetPt = -999;
    genMetPhi = -999;

    MR = -999;
    RSQ = -999;

    //Event
    nPV = -1;
    eventNum = 0;
    lumiNum = 0;
    runNum = 0;
    pvX = -99.0;
    pvY = -99.0;
    pvZ = -99.0;
    fixedGridRhoAll = -99.0;
    fixedGridRhoFastjetAll = -99.0;
    fixedGridRhoFastjetAllCalo = -99.0;
    fixedGridRhoFastjetCentralCalo = -99.0;
    fixedGridRhoFastjetCentralChargedPileUp = -99.0;
    fixedGridRhoFastjetCentralNeutral = -99.0;

}

//------ Methods to fill tree variables ------//

bool RazorTuplizer::fillEventInfo(const edm::Event& iEvent){
  //store basic event info
  isData = isData_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();
  
  if (vertices->empty()) {
    std::cout << "Warning :  Event has no primary vertices\n";
    return false; // skip the event if no PV found
  }

  bool foundPV = false;
  myPV = 0;

  pvX = -999;
  pvY = -999;
  pvZ = -999;
  nPV = 0;

  //Check for good vertices
  for(unsigned int i = 0; i < vertices->size(); i++){
    if (isGoodPV(&(vertices->at(i)))) {
      if (!foundPV) { 
	myPV = &(vertices->at(i));
	foundPV = true;
      }
      nPV++;
    }
  }

  //save vertex position
  if (foundPV) {
    pvX = myPV->x();
    pvY = myPV->y();
    pvZ = myPV->z();
  } else {
    return false;
  }

  //get rho
  fixedGridRhoAll = *rhoAll;
  fixedGridRhoFastjetCentralChargedPileUp = *rhoFastjetCentralChargedPileUp;
  fixedGridRhoFastjetCentralNeutral = *rhoFastjetCentralNeutral;

  return true;
}

bool RazorTuplizer::fillPileUp(){
  for(const PileupSummaryInfo &pu : *puInfo){
    BunchXing[nBunchXing] = pu.getBunchCrossing();
    nPU[nBunchXing] = pu.getPU_NumInteractions();
    nPUmean[nBunchXing] = pu.getTrueNumInteractions();
    nBunchXing++;
    //std::cout << "BC: " << pu.getBunchCrossing() << std::endl;
  }
  return true;
};

bool RazorTuplizer::fillMuons(){
  //PV required for Tight working point
  for(const reco::Muon &mu : *muons){
    if(mu.pt() < 5) continue;
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    muonCharge[nMuons] = mu.charge();
    muonIsLoose[nMuons] = muon::isLooseMuon(mu);
    muonIsTight[nMuons] = muon::isTightMuon(mu,*myPV);
    muon_d0[nMuons] = -mu.muonBestTrack()->dxy(myPV->position());
    muon_dZ[nMuons] = mu.muonBestTrack()->dz(myPV->position());

    muon_ip3d[nMuons] = 0;
    muon_ip3dSignificance[nMuons] = 0;
    if (mu.track().isNonnull()) {
      const TransientTrackBuilder *transientTrackBuilder = transientTrackBuilderHandle.product();
      const reco::TransientTrack &tt = transientTrackBuilder->build(mu.track());
      const std::pair<bool,Measurement1D> &ip3dpv =  
	IPTools::signedImpactParameter3D(tt, GlobalVector(mu.track()->px(),mu.track()->py(),mu.track()->pz()),*myPV);      
      muon_ip3d[nMuons] = ip3dpv.second.value();
      muon_ip3dSignificance[nMuons] = ip3dpv.second.value() / ip3dpv.second.error();
    }


    muonType[nMuons] = mu.isMuon() + mu.isGlobalMuon() + mu.isTrackerMuon() + mu.isStandAloneMuon()
      + mu.isCaloMuon() + mu.isPFMuon() + mu.isRPCMuon();
    muonQuality[nMuons] = 
      muon::isGoodMuon(mu,muon::All)
    + muon::isGoodMuon(mu,muon::AllGlobalMuons)
    + muon::isGoodMuon(mu,muon::AllStandAloneMuons)
    + muon::isGoodMuon(mu,muon::AllTrackerMuons)
    + muon::isGoodMuon(mu,muon::TrackerMuonArbitrated)
    + muon::isGoodMuon(mu,muon::AllArbitrated)      
    + muon::isGoodMuon(mu,muon::GlobalMuonPromptTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationTight)      
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityLoose)      
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityTight)      
    + muon::isGoodMuon(mu,muon::TMOneStationLoose)      
    + muon::isGoodMuon(mu,muon::TMOneStationTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtTight)      
    + muon::isGoodMuon(mu,muon::GMTkChiCompatibility)      
    + muon::isGoodMuon(mu,muon::GMStaChiCompatibility)      
    + muon::isGoodMuon(mu,muon::GMTkKinkTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationAngLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationAngTight)      
    + muon::isGoodMuon(mu,muon::TMOneStationAngLoose)      
    + muon::isGoodMuon(mu,muon::TMOneStationAngTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtTight);

    muon_pileupIso[nMuons] = mu.pfIsolationR04().sumPUPt;
    muon_chargedIso[nMuons] = mu.pfIsolationR04().sumChargedHadronPt;
    muon_photonIso[nMuons] = mu.pfIsolationR04().sumPhotonEt;
    muon_neutralHadIso[nMuons] = mu.pfIsolationR04().sumNeutralHadronEt;
    muon_ptrel[nMuons] = getLeptonPtRel( jets, &mu );
    muon_miniiso[nMuons] = getPFMiniIsolation(PFCands, dynamic_cast<const reco::Candidate *>(&mu), myPV, vertices, 0.05, 0.2, 10., false, false);

    nMuons++;
  }
  
  assert(nMuons < OBJECTARRAYSIZE);

  return true;
};

bool RazorTuplizer::fillElectrons(){
  for(const reco::GsfElectron &ele : *electrons){
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    eleCharge[nElectrons] = ele.charge();
    eleE_SC[nElectrons] = ele.superCluster()->energy();
    eleEta_SC[nElectrons] = ele.superCluster()->eta();
    elePhi_SC[nElectrons] = ele.superCluster()->phi();
    eleSigmaIetaIeta[nElectrons] = ele.sigmaIetaIeta();
    eleR9[nElectrons] = ecalLazyTools->e3x3(*(ele.superCluster()->seed())) / ele.superCluster()->rawEnergy();
    ele_dEta[nElectrons] = ele.deltaEtaSuperClusterTrackAtVtx();
    ele_dPhi[nElectrons] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_HoverE[nElectrons] = ele.hcalOverEcal();
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(myPV->position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(myPV->position());    
   ele_MissHits[nElectrons] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    //Conversion Veto
    ele_PassConvVeto[nElectrons] = false;
    if( beamSpot.isValid() && conversions.isValid() ) {
      ele_PassConvVeto[nElectrons] = !ConversionTools::hasMatchedConversion(ele,conversions,
									    beamSpot->position());
    } else {
      cout << "\n\nERROR!!! conversions not found!!!\n";
    }
  
    // 1/E - 1/P
    if( ele.ecalEnergy() == 0 ){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else if( !std::isfinite(ele.ecalEnergy())){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else {
    ele_OneOverEminusOneOverP[nElectrons] = 1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy();    
    }

    //ID MVA
    ele_IDMVATrig[nElectrons] = myMVATrig->mvaValue(ele, *myPV, *transientTrackBuilderHandle.product(), *ecalLazyTools, false);
    ele_IDMVANonTrig[nElectrons] = myMVANonTrig->mvaValue(ele,*myPV, *transientTrackBuilderHandle.product(), *ecalLazyTools,false);



    double tmpPUPt = 0;
    double tmpChargedHadronPt = 0;
    double tmpPhotonPt = 0;
    double tmpNeutralHadronPt = 0;
    
    for (const reco::PFCandidate &candidate : *PFCands) {

      //isolation cone is 0.4
      double tmpDR = deltaR(candidate.eta(), candidate.phi(), ele.eta(), ele.phi());
      if (tmpDR > 0.4) continue;
      
      //Determine if particle is PU or not
      bool tmpIsPFNoPU = isPFNoPU(candidate, myPV, vertices);

      //don't include PF ele or PF mu
      if (candidate.particleId() == reco::PFCandidate::e || candidate.particleId() == reco::PFCandidate::mu ) {
	continue;
      }
      
      if (!tmpIsPFNoPU) {

	//in vecbos, the isolation sequence structure is such that 
	//the veto below is also applied to pfPU particles.
	//so we will also apply them.
	if ( candidate.charge() != 0 &&
	     fabs(ele.superCluster()->eta()) > 1.479 && tmpDR < 0.015
	     ) {
	  continue;
	}

	tmpPUPt += candidate.pt();

      } else {
	if ( candidate.charge() != 0) {
	  if (fabs(ele.superCluster()->eta()) > 1.479 && tmpDR < 0.015) {
	    continue;	 
	  }
	  tmpChargedHadronPt += candidate.pt();
	} else if ( candidate.particleId() == reco::PFCandidate::gamma ) {

	  //veto PF photons that have the same supercluster as the electron
	  if (candidate.superClusterRef().isNonnull() && ele.superCluster().isNonnull()) {
	    if ( ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0 
		 && candidate.mva_nothing_gamma() > 0.99 
		 && &(*ele.superCluster()) == &(*candidate.superClusterRef())
		 ) {	      
	      continue;
	    }	  
	  }
	  
	  if (fabs(ele.superCluster()->eta()) > 1.479) {	    
	    if (deltaR(ele.eta(),ele.phi(), candidate.eta(), candidate.phi()) < 0.08) {
	      continue;
	    }
	  }
	  tmpPhotonPt += candidate.pt();
	} else if (candidate.particleId() == reco::PFCandidate::h0) {
	  tmpNeutralHadronPt += candidate.pt();
	}
	  
      }	
    }
    
    ele_pileupIso[nElectrons] = tmpPUPt;
    ele_chargedIso[nElectrons] = tmpChargedHadronPt;
    ele_photonIso[nElectrons] = tmpPhotonPt;
    ele_neutralHadIso[nElectrons] = tmpNeutralHadronPt;        
    //ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    //ele_CombineP4[nElectrons] = ele.ecalTrackRegressionEnergy();
    ele_ptrel[nElectrons] = getLeptonPtRel( jets, &ele );
    ele_miniiso[nElectrons] = getPFMiniIsolation(PFCands, dynamic_cast<const reco::Candidate *>(&ele), myPV, vertices, 0.05, 0.2, 10., false, false);

    nElectrons++;
  }

  assert(nElectrons < OBJECTARRAYSIZE);

  return true;
};

bool RazorTuplizer::fillTaus(){
  // for (const pat::Tau &tau : *taus) {
 //   if (tau.pt() < 20) continue;
  //   tauE[nTaus] = tau.energy();
  //   tauPt[nTaus] = tau.pt();
  //   tauEta[nTaus] = tau.eta();
  //   tauPhi[nTaus] = tau.phi();
    
  //   tau_IsLoose[nTaus] = bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
  //   tau_IsMedium[nTaus] = bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
  //   tau_IsTight[nTaus] = bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
  //   tau_passEleVetoLoose[nTaus] = bool(tau.tauID("againstElectronLooseMVA5"));
  //   tau_passEleVetoMedium[nTaus] = bool(tau.tauID("againstElectronMediumMVA5"));
  //   tau_passEleVetoTight[nTaus] = bool(tau.tauID("againstElectronTightMVA5"));
  //   tau_passMuVetoLoose[nTaus] = bool(tau.tauID("againstMuonLooseMVA"));
  //   tau_passMuVetoMedium[nTaus] = bool(tau.tauID("againstMuonMediumMVA"));
  //   tau_passMuVetoTight[nTaus] = bool(tau.tauID("againstMuonTightMVA") );  
  
  //   tau_combinedIsoDeltaBetaCorr3Hits[nTaus] = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  //   tau_eleVetoMVA[nTaus] = tau.tauID("againstElectronMVA5raw") ;
  //   tau_eleVetoCategory[nTaus] = tau.tauID("againstElectronMVA5category");
  //   tau_muonVetoMVA[nTaus] = tau.tauID("againstMuonMVAraw");
  //   tau_isoMVAnewDMwLT[nTaus] = tau.tauID("byIsolationMVA3newDMwLTraw");
  //   tau_isoMVAnewDMwoLT[nTaus] = tau.tauID("byIsolationMVA3newDMwoLTraw") ; 

  //   tau_ID[nTaus] = 
  //     bool(tau.tauID("decayModeFinding")) +
  //     bool(tau.tauID("decayModeFindingNewDMs")) +
  //     bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) +
  //     bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) +
  //     bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")) +
  //     bool(tau.tauID("againstElectronVLooseMVA5")) +
  //     bool(tau.tauID("againstElectronLooseMVA5")) +
  //     bool(tau.tauID("againstElectronMediumMVA5")) +
  //     bool(tau.tauID("againstElectronTightMVA5")) +
  //     bool(tau.tauID("againstElectronVTightMVA5")) +
  //     bool(tau.tauID("againstMuonLoose3")) +
  //     bool(tau.tauID("againstMuonTight3")) +
  //     bool(tau.tauID("againstMuonLooseMVA")) +
  //     bool(tau.tauID("againstMuonMediumMVA")) +
  //     bool(tau.tauID("againstMuonTightMVA")) +
  //     bool(tau.tauID("byVLooseIsolationMVA3newDMwLT")) +
  //     bool(tau.tauID("byLooseIsolationMVA3newDMwLT")) +
  //     bool(tau.tauID("byMediumIsolationMVA3newDMwLT")) +
  //     bool(tau.tauID("byTightIsolationMVA3newDMwLT")) +
  //     bool(tau.tauID("byVTightIsolationMVA3newDMwLT")) +
  //     bool(tau.tauID("byVVTightIsolationMVA3newDMwLT"));

  //   tau_leadCandPt[nTaus] = 0;
  //   tau_leadCandID[nTaus] = 0;
  //   tau_leadChargedHadrCandPt[nTaus] = 0;
  //   tau_leadChargedHadrCandID[nTaus] = 0;
  //   if (tau.leadCand().isNonnull()) {
  //     tau_leadCandPt[nTaus] = tau.leadCand()->pt();
  //     tau_leadCandID[nTaus] = tau.leadCand()->pdgId();
  //   }
  //   if (tau.leadChargedHadrCand().isNonnull()) { 
  //     tau_leadChargedHadrCandPt[nTaus] = tau.leadChargedHadrCand()->pt();
  //     tau_leadChargedHadrCandID[nTaus] = tau.leadChargedHadrCand()->pdgId();
  //   }
      
  //   nTaus++;
  // 
  //assert(nTaus < OBJECTARRAYSIZE);
  return true;
};

bool RazorTuplizer::fillIsoPFCandidates(){
  for (const reco::PFCandidate &candidate : *PFCands) {

    if (candidate.charge() != 0 && candidate.pt() > 5 && 
	candidate.trackRef().isNonnull() && myPV &&
	myPV->trackWeight(candidate.trackRef()) > 0 ) {

      //************************************************************
      //compute isolation of the PF candidate
      //************************************************************
      double tmpIsoPFNoPU = 0;
      double tmpIsoPFPU = 0;
      for (const reco::PFCandidate &isoCandidate : *PFCands) {	
  	if ( (candidate.pdgId() != 1 && candidate.pdgId() != 2)
	     && deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) < 0.4
  	     && !(candidate.eta() == isoCandidate.eta() && candidate.phi() == isoCandidate.phi())
	     ) {

	  bool tmpIsPFNoPU = isPFNoPU( isoCandidate, myPV, vertices);

 	  if (tmpIsPFNoPU) {
  	    tmpIsoPFNoPU += isoCandidate.pt();
  	  } else { 
  	    tmpIsoPFPU += isoCandidate.pt();
  	  }
  	}
      }

      if ( 
  	  (candidate.pt() > 50 ) ||
  	  (candidate.pt() > 20 && (tmpIsoPFNoPU - 0.5*tmpIsoPFPU)/candidate.pt() < 3.0) ||
  	  (candidate.pt() <= 20 && tmpIsoPFNoPU - 0.5*tmpIsoPFPU < 25)
  	   ) {

  	isoPFCandidatePt[nIsoPFCandidates] = candidate.pt();
  	isoPFCandidateEta[nIsoPFCandidates] = candidate.eta();
  	isoPFCandidatePhi[nIsoPFCandidates] = candidate.phi();
  	isoPFCandidateIso04[nIsoPFCandidates] = max(0.0, tmpIsoPFNoPU - 0.5*tmpIsoPFPU) ;

	if (myPV && candidate.gsfTrackRef().isNonnull()) {
	  isoPFCandidateD0[nIsoPFCandidates] = candidate.gsfTrackRef()->dxy(myPV->position());
	} else if (myPV && candidate.trackRef().isNonnull()) {
	  isoPFCandidateD0[nIsoPFCandidates] = candidate.trackRef()->dxy(myPV->position());
	} else { 
	  isoPFCandidateD0[nIsoPFCandidates] = 0.0;
	}

  	isoPFCandidatePdgId[nIsoPFCandidates] = candidate.pdgId();
	
  	nIsoPFCandidates++;
  	
      } // if candidate passes isolation
    }
  }

  assert(nIsoPFCandidates < OBJECTARRAYSIZE);

  return true;
}

bool RazorTuplizer::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //set up photon regression correction
  photonEnergyCorrector.Initialize(iSetup,edm::FileInPath("SUSYBSMAnalysis/RunOneRazorTuplizer/data/gbrv3ph_52x.root").fullPath().c_str());

  for (const reco::Photon &pho : *photons) {
    //if (pho.pt() < 20) continue;

    phoE[nPhotons] = pho.energy();
    /*
      Use Pt, Eta, and Phi are value corrected from the primary vertex position and energy regression
    */
    //phoPt[nPhotons] = pho.pt();
    //phoEta[nPhotons] = pho.eta();
    //phoPhi[nPhotons] = pho.phi();
    //phoSigmaIetaIeta[nPhotons] = pho.sigmaIetaIeta();
    phoSigmaIetaIeta[nPhotons] = sqrt( ecalLazyTools->localCovariances( *pho.superCluster()->seed() )[0] );;
    phoR9[nPhotons] = pho.r9();

    pho_HoverE[nPhotons] = pho.hadTowOverEm();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();
    pho_passEleVeto[nPhotons] = !ConversionTools::hasMatchedPromptElectron(pho.superCluster(),electrons, 
									   conversions, beamSpot->position());

    //**********************************************************
    //Compute PF isolation
    //**********************************************************
    vector<double> tmpChargedHadronPt;
    int bestVertexIndex = -1;
    for(unsigned int i = 0; i < vertices->size(); i++){
      tmpChargedHadronPt.push_back(0);
      if (&(vertices->at(i)) == myPV) {
    	bestVertexIndex = i;
      }      
    }

    double tmpPhotonPt = 0;
    double tmpNeutralHadronPt = 0;

    // First, find photon direction with respect to the good PV
    /*
      phoPos and vtxPos will be used with the Regresion Energy to get the corrected photon 4-momentum
    */
    TVector3 phoPos( pho.superCluster()->x(), pho.superCluster()->y(), pho.superCluster()->z() );
    TVector3 vtxPos( pvX, pvY, pvZ );
   
    for (const reco::PFCandidate &candidate : *PFCands) {

      //don't include PF ele or PF mu
      if (candidate.particleId() == reco::PFCandidate::e || candidate.particleId() == reco::PFCandidate::mu ) {
    	continue;
      }
      
      if ( candidate.particleId() == reco::PFCandidate::h ) {
	
    	for(unsigned int v = 0; v < vertices->size(); v++){
    	  if (!isGoodPV(&(vertices->at(v)))) continue;

    	  TVector3 vtxPos(vertices->at(v).x(),vertices->at(v).y(),vertices->at(v).z());
    	  TVector3 phoDirFromVtx = (phoPos-vtxPos).Unit();
    	  double tmpDR = 0.0;
    	  tmpDR = deltaR(phoDirFromVtx.Eta(), phoDirFromVtx.Phi(), candidate.eta(), candidate.phi());

    	  //dR cone
    	  if ( tmpDR > 0.3 ) continue;
	 
    	  //dR veto
    	  if ( tmpDR < 0.02 ) continue;	 
	  
    	  //dxy, dz to vertex 
    	  float dz = fabs(candidate.trackRef()->dz(vertices->at(v).position()));
    	  if (dz > 0.2) continue;  //vertex of this cand not compatible with this reco vtx
    	  float dxy = fabs(candidate.trackRef()->dxy(vertices->at(v).position()));
    	  if(fabs(dxy) > 0.1) continue;
	 
	  assert(v < tmpChargedHadronPt.size()); //safe guard
    	  tmpChargedHadronPt[v] += candidate.pt();
    	}
    
      } else if ( candidate.particleId() == reco::PFCandidate::gamma ) {

    	TVector3 candVtx(candidate.vertex().x(),candidate.vertex().y(),candidate.vertex().z());
    	TVector3 phoDirFromCandVtx = phoPos-candVtx;

    	//veto PF photons that have the same supercluster as the electron
    	if (candidate.superClusterRef().isNonnull() && pho.superCluster().isNonnull()) {
    	  if (&(*pho.superCluster()) == &(*candidate.superClusterRef())) {	      
    	    continue;
    	  }	  
    	}

    	double tmpDR = deltaR(phoDirFromCandVtx.Eta(), phoDirFromCandVtx.Phi(), candidate.eta(), candidate.phi());
    	if (tmpDR > 0.3) continue;

    	if (fabs(pho.superCluster()->eta()) <= 1.479) {
    	  if ( fabs(phoDirFromCandVtx.Eta() - candidate.eta()) < 0.015) continue;
    	} else {
    	  if (tmpDR < 0.07 
    	      //< 0.00864*fabs(sinh(pho.superCluster()->eta()))*4  // eta dependent veto
    	      ) continue;
    	}
    	tmpPhotonPt += candidate.pt();
      } else if ( candidate.particleId() == reco::PFCandidate::h0 )  {
    	TVector3 candVtx(candidate.vertex().x(),candidate.vertex().y(),candidate.vertex().z());
    	TVector3 phoDirFromCandVtx = phoPos-candVtx;
    	double tmpDR = deltaR(phoDirFromCandVtx.Eta(), phoDirFromCandVtx.Phi(), candidate.eta(), candidate.phi());
    	if (tmpDR > 0.3) continue;
    	tmpNeutralHadronPt += candidate.pt();
      }	  

    }

    if (bestVertexIndex >= 0) {
      pho_sumChargedHadronPt[nPhotons] = tmpChargedHadronPt[bestVertexIndex];
    } else  {
      pho_sumChargedHadronPt[nPhotons] = -99.0;
    }
    pho_sumNeutralHadronEt[nPhotons] = tmpNeutralHadronPt;
    pho_sumPhotonEt[nPhotons] = tmpPhotonPt;

    pho_sumWorstVertexChargedHadronPt[nPhotons] = 0;
    for (int v=0; v<int(tmpChargedHadronPt.size()); ++v) {
      if (!isGoodPV(&(vertices->at(v)))) continue;
      if (pho_sumWorstVertexChargedHadronPt[nPhotons] < tmpChargedHadronPt[v]) {
    	pho_sumWorstVertexChargedHadronPt[nPhotons] = tmpChargedHadronPt[v];
      }
    }

    std::pair<double,double> photonEnergyCorrections = photonEnergyCorrector.CorrectedEnergyWithErrorV3(pho, *vertices, *rhoAll, *ecalLazyTools, iSetup, false);
    pho_RegressionE[nPhotons] = photonEnergyCorrections.first;
    pho_RegressionEUncertainty[nPhotons] = photonEnergyCorrections.second;
    //compute photon corrected 4-mometum 
    TLorentzVector phoP4 = photonP4FromVtx( vtxPos, phoPos, pho_RegressionE[nPhotons] );
    phoPt[nPhotons]  = phoP4.Pt();
    phoEta[nPhotons] = phoP4.Eta();
    phoPhi[nPhotons] = phoP4.Phi();

    pho_IDMVA[nPhotons] = pho.pfMVA();
    pho_superClusterEta[nPhotons] = pho.superCluster()->eta();
    pho_superClusterPhi[nPhotons] = pho.superCluster()->phi();
    pho_hasPixelSeed[nPhotons] = pho.hasPixelSeed();

    /*
    const reco::Candidate* genPhoton = pho.genPhoton();
    if(genPhoton != NULL)std::cout << "======>gen PT: " << genPhoton->pt() <<
      " recoPT: " << pho.pt() << std::endl;
    */
    nPhotons++;
  }

  assert(nPhotons < OBJECTARRAYSIZE);

  //delete lazyToolnoZS;
  return true;
};

bool RazorTuplizer::fillJets(const edm::Event& iEvent){

  edm::Handle<edm::ValueMap<float> > puJetIdMVA;
  iEvent.getByLabel(edm::InputTag("recoPuJetMva","full53xDiscriminant"),puJetIdMVA);
  edm::Handle<edm::ValueMap<int> > puJetIdFlag;
  iEvent.getByLabel(edm::InputTag("recoPuJetMva","full53xId"),puJetIdFlag);

  edm::Handle<reco::JetIDValueMap> hJetIDMap;
  iEvent.getByLabel( "ak5JetID", hJetIDMap );
  PFJetIDSelectionFunctor jetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
  PFJetIDSelectionFunctor jetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
  pat::strbitset ret_loose = jetIDLoose.getBitTemplate(); 
  pat::strbitset ret_tight = jetIDTight.getBitTemplate(); 

  const reco::PFJetCollection inJets = *(jets.product());  

  for (reco::PFJetCollection::const_iterator j = inJets.begin(); 
       j != inJets.end(); ++j) {
    
    reco::PFJetRef jetRef(jets, j - inJets.begin());    
    reco::JetBaseRef jetBaseRef(jetRef);
    
    if (j->pt() < 20) continue;
    jetE[nJets] = j->energy();
    jetPt[nJets] = j->pt();
    jetEta[nJets] = j->eta();
    jetPhi[nJets] = j->phi(); 
    jetCSV[nJets] = (*(csvBJetTags.product()))[jetBaseRef];
    //jetCISV[nJets] = (*(cisvBJetTags.product()))[jetBaseRef];
    jetMass[nJets] = j->mass();
    jetJetArea[nJets] = j->jetArea();
    jetPileupE[nJets] = j->pileup();
    jetPileupId[nJets] =  (*puJetIdMVA)[jetBaseRef];
    jetPileupIdFlag[nJets] = (*puJetIdFlag)[jetBaseRef];

    ret_loose.set(false);
    ret_tight.set(false);
    jetPassIDLoose[nJets] = jetIDLoose(*j, ret_loose);
    jetPassIDTight[nJets] = jetIDTight(*j, ret_loose);        
    jetPassMuFrac[nJets]  = ( j->muonEnergyFraction() < 0.80 );
    jetPassEleFrac[nJets]  = ( j->electronEnergyFraction() < 0.90 );

    if (useGen_) {
      jetPartonFlavor[nJets] = (*jetFlavorMatch)[edm::RefToBase<reco::Jet>(jetRef)].getFlavour();
    } else {
      jetPartonFlavor[nJets] = -999;
    }    
    nJets++;
  }

  assert(nJets < OBJECTARRAYSIZE);

  return true;
};

bool RazorTuplizer::fillJetsAK8(){
  // for (const pat::Jet &j : *jetsAK8) {
  //   fatJetE[nFatJets] = j.correctedP4(0).E();
  //   fatJetPt[nFatJets] = j.correctedP4(0).Pt();
  //   fatJetEta[nFatJets] = j.correctedP4(0).Eta();
  //   fatJetPhi[nFatJets] = j.correctedP4(0).Phi();
  //   fatJetPrunedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSPrunedLinks");
  //   fatJetTrimmedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSTrimmedLinks");
  //   fatJetFilteredM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
  //   fatJetTau1[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau1");
  //   fatJetTau2[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau2");
  //   fatJetTau3[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau3");
  //   nFatJets++;
  // }
  //assert(nFatJets < OBJECTARRAYSIZE);

  return true;
};

bool RazorTuplizer::fillMet(const edm::Event& iEvent){
  const reco::PFMET &Met = mets->front();
  metPt = Met.pt();
  metPhi = Met.phi();
  sumMET = Met.sumEt();
  
  edm::Handle< bool > trackerFailureFilter;
  iEvent.getByLabel("trackingFailureFilter", trackerFailureFilter);
  Flag_trackingFailureFilter = bool(*trackerFailureFilter);

  Flag_goodVertices = true;

  edm::Handle< reco::BeamHaloSummary > beamHaloH;
  iEvent.getByLabel("BeamHaloSummary", beamHaloH);
  Flag_CSCTightHaloFilter = ! beamHaloH->CSCTightHaloId();

  edm::Handle< bool > ecalDeadCellTriggerPrimitiveFilter;
  iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter", ecalDeadCellTriggerPrimitiveFilter);
  Flag_EcalDeadCellTriggerPrimitiveFilter = bool(*ecalDeadCellTriggerPrimitiveFilter);

  edm::Handle< bool > ecalDeadCellBoundaryEnergyFilter;
  iEvent.getByLabel("EcalDeadCellBoundaryEnergyFilter", ecalDeadCellBoundaryEnergyFilter);
  Flag_EcalDeadCellBoundaryEnergyFilter = bool(*ecalDeadCellBoundaryEnergyFilter);

  edm::Handle< bool > ecalLaserCorrFilter;
  iEvent.getByLabel("ecalLaserCorrFilter", ecalLaserCorrFilter);
  Flag_ecalLaserCorrFilter = bool(*ecalLaserCorrFilter);

  edm::Handle< bool > eeBadScFilter;
  iEvent.getByLabel("eeBadScFilter", eeBadScFilter);  
  Flag_eeBadScFilter = bool(*eeBadScFilter);

  edm::Handle< bool > HBHENoiseFilterResult;
  iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), HBHENoiseFilterResult);  
  Flag_HBHENoiseFilter = bool(*HBHENoiseFilterResult);

  edm::Handle< bool > hcalLaserEventFilter;
  iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventFilter);   
  Flag_hcalLaserEventFilter = bool(*hcalLaserEventFilter);

  edm::Handle< bool > manystripclus53X;
  iEvent.getByLabel("manystripclus53X", manystripclus53X);
  Flag_trkPOG_manystripclus53X = *manystripclus53X;

  edm::Handle< bool > toomanystripclus53X;
  iEvent.getByLabel("toomanystripclus53X", toomanystripclus53X);
  Flag_trkPOG_toomanystripclus53X = *toomanystripclus53X;

  edm::Handle< bool > logErrorTooManyClusters;
  iEvent.getByLabel("logErrorTooManyClusters", logErrorTooManyClusters);
  Flag_trkPOG_logErrorTooManyClusters = *logErrorTooManyClusters;

  Flag_trkPOGFilters = !(Flag_trkPOG_manystripclus53X) && 
    !(Flag_trkPOG_toomanystripclus53X) &&
    !(Flag_trkPOG_logErrorTooManyClusters);

  Flag_METFilters = Flag_trackingFailureFilter && Flag_CSCTightHaloFilter && Flag_EcalDeadCellTriggerPrimitiveFilter
    && Flag_ecalLaserCorrFilter && Flag_eeBadScFilter && Flag_HBHENoiseFilter && Flag_hcalLaserEventFilter
    && Flag_trkPOGFilters;

  //***************************************************************
  //Special recipes ported from vecbos sequence
  //Is not part of MET twiki recommendation now
  //***************************************************************
  edm::Handle< bool > EcalDeadCellEventFlag;
  iEvent.getByLabel("EcalDeadCellEventFlagProducer", EcalDeadCellEventFlag);
  Flag_EcalDeadCellEvent = *EcalDeadCellEventFlag;
  
  edm::Handle<AnomalousECALVariables> anomalousECALvarsHandle;
  edm::InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter","anomalousECALVariables");
  iEvent.getByLabel(ecalAnomalousFilterTag, anomalousECALvarsHandle);
  AnomalousECALVariables anomalousECALvars;
  if (anomalousECALvarsHandle.isValid()) {
  anomalousECALvars = *anomalousECALvarsHandle;
  } else {
  edm::LogWarning("anomalous ECAL Vars not valid/found");
  }
  Flag_IsNotDeadEcalCluster = !(anomalousECALvars.isDeadEcalCluster());

  edm::Handle< int > ECALDeadDRFilter;
  iEvent.getByLabel("simpleDRFlagProducer","deadCellStatus", ECALDeadDRFilter);
  int ECALDeadDRFilterFlag = *ECALDeadDRFilter;
  edm::Handle< int > ECALBoundaryDRFilter;
  iEvent.getByLabel("simpleDRFlagProducer","boundaryStatus", ECALBoundaryDRFilter);
  int ECALBoundaryDRFilterFlag = *ECALBoundaryDRFilter;
  Flag_EcalDeadDR = (ECALDeadDRFilterFlag>0) ? true : false;
  Flag_EcalBoundaryDR = (ECALBoundaryDRFilterFlag>0) ? true : false;

  return true;
};

bool RazorTuplizer::fillMC(){
    for(const reco::GenJet &j : *genJets){
        genJetE[nGenJets] = j.energy();
        genJetPt[nGenJets] = j.pt();
        genJetEta[nGenJets] = j.eta();
        genJetPhi[nGenJets] = j.phi();
        nGenJets++;
    }

    const reco::GenMET &genMet = genMets->front();
    genMetPt = genMet.pt();
    genMetPhi = genMet.phi();

    bool foundGenVertex = false;
    for(size_t i=0; i<genParticles->size();i++){
      if (!foundGenVertex) {
	for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j) {
	  const reco::Candidate *dau = (*genParticles)[i].daughter(j);
	  if (dau) {
	    genVertexX = dau->vx();
	    genVertexY = dau->vy();
	    genVertexZ = dau->vz();
	    foundGenVertex = true;
	    break;
	  }
	}
      }
    }
 
    genWeight = genInfo->weight();
    genSignalProcessID = genInfo->signalProcessID();
    genQScale = genInfo->qScale();
    genAlphaQCD = genInfo->alphaQCD();
    genAlphaQED = genInfo->alphaQED();

    return true;
}

bool RazorTuplizer::fillRazor(){   
  return true;
}

bool RazorTuplizer::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding

  //Fills selected gen particles
  for(size_t i=0; i<genParticles->size();i++){
    if ((*genParticles)[i].status() == 11) continue; //skip sherpa documentation lines

    if(   (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 )
	  || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
	  || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() != 2 )	   
	  || (abs((*genParticles)[i].pdgId()) >= 23 && abs((*genParticles)[i].pdgId()) <= 25
	      )
	  || (abs((*genParticles)[i].pdgId()) == 22 && 
	      findFirstMotherWithDifferentID(&(*genParticles)[i]) && abs(findFirstMotherWithDifferentID(&(*genParticles)[i])->pdgId()) <= 25
	      )
	  || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
	  || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
	  ){
      prunedV.push_back(&(*genParticles)[i]);
    }

    // cout << i << " : " << (*genParticles)[i].pdgId() << " " << (*genParticles)[i].status() << " " 
    // 	 << (*genParticles)[i].pt() << " "  << (*genParticles)[i].eta() << " "  << (*genParticles)[i].phi() << " : " 
    // 	 << (*genParticles)[i].numberOfMothers() << " ";
    // if ((*genParticles)[i].numberOfMothers() > 0 ) cout << (*genParticles)[i].mother(0)->pdgId() << " ";
    // cout << "\n";

    //if (prunedV.size()<99) prunedV.push_back(&(*genParticles)[i]); //keep all pruned particles
  }

  //Total number of gen particles
  nGenParticle = prunedV.size();

  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++){
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;
    if(prunedV[i]->numberOfMothers() > 0){
      
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID) {
      	gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
      }
      
      //find the mother and keep going up the mother chain if the ID's are the same
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++){	
      	if(prunedV[j] == originalMotherWithSameID){
      	  gParticleMotherIndex[i] = j;
      	  break;
      	}
      }
    } else {
      gParticleMotherIndex[i] = -1;
    }
  }

  assert(nGenParticle < GENPARTICLEARRAYSIZE);

  return true;
};

bool RazorTuplizer::fillTrigger(const edm::Event& iEvent){

  //initialize trigger decisions
  for (unsigned int j = 0; int(j) < int(NTriggersMAX); ++j) {
    triggerDecision[j] = false;
  }

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  //********************************************************************
  // Debug : Print Trigger Menu
  //********************************************************************
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  // std::cout << "trigger: " << names.size() << "\n";
  // for (int i=0; i<int(names.size()); ++i) {
  //   std::cout << names.triggerName(i) << "\n";
  // }
  // std::cout << "\n\n";

  //********************************************************************
  // Save trigger decisions in array of booleans
  //********************************************************************
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {    
    string hltPathNameReq = "HLT_";   
    
    if ((names.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((names.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);
    
    for (unsigned int j = 0; int(j) < int(NTriggersMAX); ++j) {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j]) {
  	triggerDecision[j] = triggerBits->accept(i);
      }
    }    
  }
    

  return true;
}


//------ Method called for each event ------//

void RazorTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  myPV = 0;
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches
  
  NEvents->Fill(0); //increment event count

  //set up the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilderHandle);

  //set up ecal tools
  ecalLazyTools = new EcalClusterLazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), 
					   edm::InputTag("reducedEcalRecHitsEE"));
  ecalLocal = new EcalClusterLocal;

  // //filler methods should fill relevant tree variables and return false if the event should be rejected
  //bool isGoodEvent = true;


  bool isGoodEvent =
  fillEventInfo(iEvent)
    && fillMuons() 
    && fillElectrons()
    //   && fillTaus()
    && fillIsoPFCandidates()
    && fillPhotons(iEvent,iSetup)
    && fillJets(iEvent)
    && fillJetsAK8()
    && fillMet(iEvent)
    ;

  bool isGoodMCEvent = true;
  if (useGen_) {
    isGoodMCEvent = fillMC()
      && fillPileUp()
      && fillGenParticles();
  }
  
  isGoodEvent = isGoodEvent&&isGoodMCEvent;
  
  //NOTE: if any of the above functions return false, the event will be rejected immediately with no further processing

  if (enableTriggerInfo_) isGoodEvent = (isGoodEvent && fillTrigger(iEvent));

  //fill the tree if the event wasn't rejected
  if(isGoodEvent) {
    RazorEvents->Fill();
  } else {
    cout << "Warning: Event is not good, and was not saved.\n";
  }

  //memory cleanup
  delete ecalLazyTools;
  delete ecalLocal;

}

//------ Method called once each job just before starting event loop ------//
void RazorTuplizer::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void RazorTuplizer::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(RazorTuplizer);
