#include "RunOneRazorTuplizer.h"

//------ Auxiliary tools for RazorTuplizer class ------//

vector<TLorentzVector> RazorTuplizer::getHemispheres(vector<TLorentzVector> jets){
  int nJets = jets.size();
  vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
  vector<TLorentzVector> possibleHem2s;

  //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
  int nComb = pow(2, nJets);
  
  //step 1: store all possible partitions of the input jets
  int j_count;
  for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = nComb/2;
    int count = 0;
    while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
      if(itemp/j_count == 1){
	j_temp1 += jets[count];
      } else {
	j_temp2 += jets[count];
      }
      itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
      j_count /= 2;
      count++;
    }
    possibleHem1s.push_back(j_temp1);
    possibleHem2s.push_back(j_temp2);
  }
  
  //step 2: choose the partition that minimizes m1^2 + m2^2
  double mMin = -1;
  TLorentzVector myHem1;
  TLorentzVector myHem2;
  for(size_t i=0; i < possibleHem1s.size(); i++){
    double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
    if(mMin < 0 || mTemp < mMin){
      mMin = mTemp;
      myHem1 = possibleHem1s[i];
      myHem2 = possibleHem2s[i];
    }
  }
  
  //return the hemispheres in decreasing order of pt
  vector<TLorentzVector> hemsOut;
  if(myHem1.Pt() > myHem2.Pt()){
    hemsOut.push_back(myHem1);
    hemsOut.push_back(myHem2);
  } else {
    hemsOut.push_back(myHem2);
    hemsOut.push_back(myHem1);
  }
  
  return hemsOut;
}

float RazorTuplizer::computeMR(TLorentzVector hem1, TLorentzVector hem2){
  return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

float RazorTuplizer::computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
  double mR = computeMR(hem1, hem2);
  double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
  double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
  double mTR = sqrt(term1 - term2);
  return (mTR / mR) * (mTR / mR);
}

bool RazorTuplizer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++){
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

const reco::Candidate* RazorTuplizer::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
}

const reco::Candidate* RazorTuplizer::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 || 
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {      
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }
  return 0;
}

// //A copy from RecoEgamma/EgammaTools/src/ConversionTools.cc
// //temporary solution because I'm not sure how to convert a PAT electron handle
// //to a GSF electron handle
// bool RazorTuplizer::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<reco::GsfElectron> > &eleCol,
// 					const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot, 
// 					float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {

//    //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
//    //and not matching any conversion in the collection passing the quality cuts
 
//    if (sc.isNull()) return false;
   
//    for (std::vector<reco::GsfElectron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
//      //match electron to supercluster
//      if (it->superCluster()!=sc) continue;
 
//      //check expected inner hits
//      if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
 
//      //check if electron is matching to a conversion
//      if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
        
//      return true;
//    }
   
//    return false;
  
// }

bool RazorTuplizer::isGoodPV( const reco::Vertex *v) {

  if(v->isValid() && !v->isFake()
     && v->ndof() > 4 && fabs(v->z()) <= 24 && v->position().Rho() <= 2
     ) {
    return true;
  } else {
    return false;
  }

}



// **********************************************************
// Old PFNoPU procedure
// **********************************************************
// bool RazorTuplizer::isPFNoPU( const reco::PFCandidate candidate,  const reco::Vertex *PV, edm::Handle<reco::VertexCollection> vertices) {

//   bool tmpIsPFNoPU = true;
//   if(candidate.particleId() == reco::PFCandidate::h) {
//     if(candidate.trackRef().isNonnull() && PV &&
//        PV->trackWeight(candidate.trackRef()) > 0) {
//       tmpIsPFNoPU = true;
//     } else { 
      
//       bool vertexFound = false;
//       const reco::Vertex *closestVtx = 0;
//       double dzmin = 10000;

//       // loop over vertices
//       for (reco::VertexCollection::const_iterator inV = vertices->begin(); 
// 	   inV != vertices->end(); ++inV) {
// 	if(candidate.trackRef().isNonnull() && 
// 	   inV->trackWeight(candidate.trackRef()) > 0) {
// 	  vertexFound = true;
// 	  closestVtx = &(*inV);
// 	  break;
// 	}
// 	double dz = fabs(candidate.vertex().z() - inV->z());
// 	if(dz < dzmin) {
// 	  closestVtx = &(*inV);
// 	  dzmin = dz;
// 	}            
//       }

//       bool fCheckClosestZVertex = true; //we use option 1
//       if(fCheckClosestZVertex) {
// 	// Fallback: if track is not associated with any vertex,
// 	// associate it with the vertex closest in z
// 	if(vertexFound || closestVtx != PV) {
// 	  tmpIsPFNoPU = kFALSE;
// 	} else {
// 	  tmpIsPFNoPU = kTRUE;
// 	}
//       } else {
// 	if(vertexFound && closestVtx != PV) {
// 	  tmpIsPFNoPU = kFALSE;
// 	} else {
// 	  tmpIsPFNoPU = kTRUE;
// 	}
//       }
//     } //charged hadron & trk stuff
//   } else { // neutrals 
//     //
//     tmpIsPFNoPU = kTRUE;
//   }

//   return tmpIsPFNoPU;
// }
  
// **********************************************************
// PFNoPU procedure synchronized with the Vecbos sequence
// **********************************************************
bool RazorTuplizer::isPFNoPU( const reco::PFCandidate candidate,  const reco::Vertex *PV, edm::Handle<reco::VertexCollection> vertices) {

  bool tmpIsPFNoPU = true;
  if(candidate.particleId() == reco::PFCandidate::h) {
   
    bool vertexFound = false;
    const reco::Vertex *closestVtx = 0;
    double bestweight = -1;
    double dzmin = 10000;

    // loop over vertices and find the vertex with largest weight for given track
    for (reco::VertexCollection::const_iterator inV = vertices->begin(); 
	 inV != vertices->end(); ++inV) {
      if(!isGoodPV(&(*inV))) continue;
      if(candidate.trackRef().isNonnull()) {	  
	double w = inV->trackWeight(candidate.trackRef());
	if ( w > 0 && w > bestweight ) {
	  vertexFound = true;
	  bestweight = w;
	  closestVtx = &(*inV);
	  break;
	}
      }
    }
    
    //cout << "best weight: " << bestweight << " : " << vertexFound << " " << bool(closestVtx == PV) << "\n";
    //fall back option is to find the vertex closest in z to the track
    if (!vertexFound) {
      for (reco::VertexCollection::const_iterator inV = vertices->begin(); 
	   inV != vertices->end(); ++inV) {
	if(!isGoodPV(&(*inV))) continue;
	double dz = fabs( candidate.vertex().z() - inV->z() );
	//cout << "dz " << candidate.vertex().z() << " " << inV->z() << " " << dz << "\n";
	if (dz < dzmin) {
	  dzmin = dz;
	  vertexFound = true;
	  closestVtx = &(*inV);
	}
      }
      //cout << "fall back dz : " << dzmin << " " << vertexFound << " " << bool(closestVtx == PV) << "\n";
    }
    
    if (!vertexFound || closestVtx == PV) {
      tmpIsPFNoPU = true;
    } else {
      tmpIsPFNoPU = false;
    }
  } else { // neutrals 
    //
    tmpIsPFNoPU = true;
  }

  return tmpIsPFNoPU;
}
  




double RazorTuplizer::getPFMiniIsolation(edm::Handle<reco::PFCandidateCollection> pfcands,
					 const reco::Candidate* ptcl,
					 const reco::Vertex *PV, 
					 edm::Handle<reco::VertexCollection> vertices,
					 double r_iso_min, double r_iso_max, double kt_scale,
					 bool use_pfweight, bool charged_only) {
  
  if (ptcl->pt()<5.) return 99999.;
  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if(ptcl->isElectron()) {
    if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  } else if(ptcl->isMuon()) {
    deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
  } else {
    //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }
  double iso_nh(0.); double iso_ch(0.);
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh(0.5);
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
  for (const reco::PFCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, *ptcl);
    if (dr > r_iso) continue;

    //is the candidate PFNoPU?
    bool tmpIsPFNoPU = isPFNoPU( pfc , PV, vertices ); 

    ////////////////// NEUTRALS /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {
  	double wpf(1.);
  	if (use_pfweight){
  	  double wpv(0.), wpu(0.);
  	  for (const reco::PFCandidate &jpfc : *pfcands) {
  	    double jdr = deltaR(pfc, jpfc);
  	    if (pfc.charge()!=0 || jdr<0.00001) continue;
  	    double jpt = jpfc.pt();
  	    if (tmpIsPFNoPU) wpv *= jpt/jdr;
  	    else wpu *= jpt/jdr;
  	  }
  	  wpv = log(wpv);
  	  wpu = log(wpu);
  	  wpf = wpv/(wpv+wpu);
  	}
  	/////////// PHOTONS ////////////
  	if (abs(pfc.pdgId())==22) {
  	  if(dr < deadcone_ph) continue;
  	  iso_ph += wpf*pfc.pt();
  	  /////////// NEUTRAL HADRONS ////////////
  	} else if (abs(pfc.pdgId())==130) {
  	  if(dr < deadcone_nh) continue;
  	  iso_nh += wpf*pfc.pt();
  	}
      }
      ////////////////// CHARGED from PV /////////////////////////
    } else if (tmpIsPFNoPU){
      if (abs(pfc.pdgId())==211) {
  	if(dr < deadcone_ch) continue;
  	iso_ch += pfc.pt();
      }
      ////////////////// CHARGED from PU /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
  	if(dr < deadcone_pu) continue;
  	iso_pu += pfc.pt();
      }
    }
  }
  double iso(0.);
  if (charged_only){
    iso = iso_ch;
  } else {
    iso = iso_ph + iso_nh;
    if (!use_pfweight) iso -= 0.5*iso_pu;
    if (iso>0) iso += iso_ch;
    else iso = iso_ch;
  }
  iso = iso/ptcl->pt();
  return iso;

  return 0;
}


//**************************************************************
//Compute ptRel for leptons
//1) find closest jet
//2) subtract lepton from jet
//3) project lepton momentum perpendicular to closest jet
//**************************************************************
double RazorTuplizer::getLeptonPtRel(edm::Handle<reco::PFJetCollection> jets, const reco::Candidate* lepton) {

    const reco::PFJet *closestJet = 0;
    double minDR = 9999;
    for (const reco::PFJet &j : *jets) {
      if (j.pt() < 20) continue;
      double tmpDR = deltaR(j.eta(),j.phi(),lepton->eta(),lepton->phi());
      if (tmpDR < minDR) {
    	minDR = tmpDR;
    	closestJet = &j;
      }
    }

    //if no jet was found nearby, return some large default value
    if (!closestJet) return 9999;

    TLorentzVector closestJetFourVector(closestJet->px(),closestJet->py(),closestJet->pz(),closestJet->energy());    

    for (unsigned int i = 0, n = closestJet->numberOfSourceCandidatePtrs(); i < n; ++i) {
      
      const reco::PFCandidate &candidate = dynamic_cast<const reco::PFCandidate &>(*(closestJet->sourceCandidatePtr(i)));
      bool isPartOfLepton = false;

      if (lepton->isMuon()) {
	if ( candidate.muonRef().isNonnull() &&   &(*candidate.muonRef()) == lepton) {
	  isPartOfLepton = true;
	}
      }

      if (lepton->isElectron()) {
 	if ( (candidate.superClusterRef().isNonnull() && &(*(((reco::GsfElectron*)lepton)->superCluster())) == &(*candidate.superClusterRef())) || 
	     (candidate.gsfElectronRef().isNonnull() && lepton ==  &(*candidate.gsfElectronRef())) || 
	     (candidate.gsfTrackRef().isNonnull() && &(*(((reco::GsfElectron*)lepton)->gsfTrack()))  ==  &(*candidate.gsfTrackRef())) 
	     ) {
	  isPartOfLepton = true;
	}
      }
      //if the PF candidate is part of the muon, subtract its momentum from the jet momentum
      if (isPartOfLepton) {
    	closestJetFourVector.SetPxPyPzE( closestJetFourVector.Px() - candidate.px(), 
    					 closestJetFourVector.Py() - candidate.py(),
    					 closestJetFourVector.Pz() - candidate.pz(),
    					 closestJetFourVector.E() - candidate.energy());
      }
    }
    TLorentzVector lepFourVector(lepton->px(),lepton->py(),lepton->pz(),lepton->energy());    
    return lepFourVector.Perp(closestJetFourVector.Vect());

  return 0;

}

TLorentzVector RazorTuplizer::photonP4FromVtx( TVector3 vtx, TVector3 phoPos, double E )
{
  TVector3 p_hat;//Corrected photon direction
  p_hat = phoPos - vtx;
  TVector3 phoP3 = p_hat.Unit()*E;
  TLorentzVector phoP4;
  phoP4.SetVectM( phoP3, .0 );
  return phoP4;
};
