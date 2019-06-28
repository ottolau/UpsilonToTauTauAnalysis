#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"

#include "UpsilonToTauTauAnalysis/UpsilonToTauTauAnalyzer/interface/UpsilonToTauTauAnalyzer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
float          gen_upsilon_pt_;
float          gen_upsilon_eta_;
float          gen_upsilon_phi_;
float          gen_upsilon_m_;

bool           found_sv_reco_;
float          reco_sv_chi2_;
float          reco_sv_ndof_;
float          reco_sv_prob_;
float          reco_sv_ctxy_;
float          reco_sv_cosangle_;
float          reco_sv_lxy_;
float          reco_sv_lxyerror_;


void UpsilonToTauTauAnalyzer::branchesMatchGenParticles(TTree* tree) {

  tree->Branch("gen_upsilon_pt",            &gen_upsilon_pt_);
  tree->Branch("gen_upsilon_eta",           &gen_upsilon_eta_);
  tree->Branch("gen_upsilon_phi",           &gen_upsilon_phi_);
  tree->Branch("gen_upsilon_m",             &gen_upsilon_m_);

  tree->Branch("found_sv_reco",             &found_sv_reco_);
  tree->Branch("reco_sv_chi2",              &reco_sv_chi2_);
  tree->Branch("reco_sv_ndof",              &reco_sv_ndof_);
  tree->Branch("reco_sv_prob",              &reco_sv_prob_);
  tree->Branch("reco_sv_ctxy",              &reco_sv_ctxy_);
  tree->Branch("reco_sv_cosangle",          &reco_sv_cosangle_);
  tree->Branch("reco_sv_lxy",               &reco_sv_lxy_);
  tree->Branch("reco_sv_lxyerror",          &reco_sv_lxyerror_);

}

void UpsilonToTauTauAnalyzer::fillMatchGenParticles(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv, reco::Vertex vtx) {
    
  gen_upsilon_pt_       = -99;
  gen_upsilon_eta_      = -99;
  gen_upsilon_phi_      = -99;
  gen_upsilon_m_        = -99;

  found_sv_reco_        = false;
  reco_sv_chi2_         = -99;
  reco_sv_ndof_         = -99;
  reco_sv_prob_         = -99;
  reco_sv_ctxy_         = -99;
  reco_sv_cosangle_     = -99;
  reco_sv_lxy_          = -99;
  reco_sv_lxyerror_     = -99;


  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  edm::Handle<edm::View<pat::PackedGenParticle> > packedGenParticlesHandle;
  e.getByToken(packedGenParticlesCollection_, packedGenParticlesHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::Handle<pat::PackedCandidateCollection> losttracks;
  e.getByToken(lostTracksLabel_, losttracks);

  std::vector<pat::PackedCandidate> alltracks;
  alltracks.reserve(pfcands->size() + losttracks->size());
  alltracks.insert(alltracks.end(), pfcands->begin(), pfcands->end());
  alltracks.insert(alltracks.end(), losttracks->begin(), losttracks->end());

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  VertexDistanceXY vertTool;

  int found_pi_gen = 0;
  vector<int> pi_gen_id;

  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    // get Upsilon (pdgid=553)
    if (abs(ip->pdgId()) != 553) continue;
    found_pi_gen = 0;
    const Candidate* Upilson = &(*ip);

    // loop over all generated final state particles, get 
    for (edm::View<pat::PackedGenParticle>::const_iterator ipp = packedGenParticlesHandle->begin(); ipp != packedGenParticlesHandle->end(); ++ipp) {
      // check if the generated final state particle comes from the Y(1S)
      if (! (ipp->mother(0) && isAncestor(Upilson, ipp->mother(0)))) continue;
      if (abs(ipp->pdgId()) == 211) {
	found_pi_gen++;
	pi_gen_id.push_back(ipp - packedGenParticlesHandle->begin());
      }
    }
    if (found_pi_gen == 6) break;

  }

  if (found_pi_gen != 6) return;


  float pionM = 0.139;
  TLorentzVector pi_gen_lv, upsilon_gen_lv;
  for (int i = 0; i < 6; ++i) {
    pi_gen_lv.SetPtEtaPhiM(packedGenParticlesHandle->at(pi_gen_id[i]).pt(), packedGenParticlesHandle->at(pi_gen_id[i]).eta(), packedGenParticlesHandle->at(pi_gen_id[i]).phi(), pionM);
    upsilon_gen_lv += pi_gen_lv;
  }

  gen_upsilon_pt_ = upsilon_gen_lv.Pt();
  gen_upsilon_eta_ = upsilon_gen_lv.Eta();
  gen_upsilon_phi_ = upsilon_gen_lv.Phi();
  gen_upsilon_m_ = upsilon_gen_lv.M();


  /*
  if (found_pi1_reco && found_pi2_reco && found_pi3_reco) {

    // kalman vertex fit
    KalmanVertexFitter theKalmanFitter(false);
    TransientVertex TauKalmanVertex;
    std::vector<reco::TransientTrack> tempTracks;
    tempTracks.push_back(getTransientTrack( tracksHandle->at(kp_reco_id) ));
    tempTracks.push_back(getTransientTrack( tracksHandle->at(km_reco_id) ));
    tempTracks.push_back(getTransientTrack( *(electronHandle->at(elep_reco_id).gsfTrack()) ));
    tempTracks.push_back(getTransientTrack( *(electronHandle->at(elem_reco_id).gsfTrack()) ));

    TauKalmanVertex = theKalmanFitter.vertex(tempTracks);
    if (TauKalmanVertex.isValid() && TauKalmanVertex.totalChiSquared() > 0.0) {
      elep_reco_lv.SetPtEtaPhiM(electronHandle->at(elep_reco_id).pt(), electronHandle->at(elep_reco_id).eta(), electronHandle->at(elep_reco_id).phi(), eleM);
      elem_reco_lv.SetPtEtaPhiM(electronHandle->at(elem_reco_id).pt(), electronHandle->at(elem_reco_id).eta(), electronHandle->at(elem_reco_id).phi(), eleM);
      kp_reco_lv.SetPtEtaPhiM(tracksHandle->at(kp_reco_id).pt(), tracksHandle->at(kp_reco_id).eta(), tracksHandle->at(kp_reco_id).phi(), kaonM);
      km_reco_lv.SetPtEtaPhiM(tracksHandle->at(km_reco_id).pt(), tracksHandle->at(km_reco_id).eta(), tracksHandle->at(km_reco_id).phi(), kaonM);
      ee_reco_lv = elep_reco_lv + elem_reco_lv;
      phi_reco_lv = kp_reco_lv + km_reco_lv;
      bs_reco_lv = ee_reco_lv + phi_reco_lv; 
   
      float ctxy = ((TauKalmanVertex.position().x() - pv.x())*bs_reco_lv.Px() + (TauKalmanVertex.position().y() - pv.y())*bs_reco_lv.Py())/(pow(bs_reco_lv.Pt(),2))*bs_reco_lv.M();
      
      math::XYZVector perp(bs_reco_lv.Px(), bs_reco_lv.Py(), 0.);
      math::XYZPoint dxybs(-1*(pv.x() - TauKalmanVertex.position().x()), -1*(pv.y() - TauKalmanVertex.position().y()), 0.);
      math::XYZVector vperp(dxybs.x(), dxybs.y(), 0.);
      float cosAngle = vperp.Dot(perp)/(vperp.R()*perp.R());

      found_sv_reco_ = true;
      reco_sv_chi2_ = TauKalmanVertex.totalChiSquared();
      reco_sv_ndof_ = TauKalmanVertex.degreesOfFreedom();
      reco_sv_prob_ = TMath::Prob(TauKalmanVertex.totalChiSquared(), TauKalmanVertex.degreesOfFreedom());
      reco_sv_ctxy_ = ctxy;
      reco_sv_cosangle_ = cosAngle;
      reco_sv_lxy_ = vertTool.distance(vtx, TauKalmanVertex).value();
      reco_sv_lxyerror_ = vertTool.distance(vtx, TauKalmanVertex).error();
      reco_ee_m_ = ee_reco_lv.M();
      reco_phi_m_ = phi_reco_lv.M();
      reco_bs_m_ = bs_reco_lv.M();

    }

  }
  */

}

bool UpsilonToTauTauAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle)
{
//particle is already the ancestor
        if(ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
//if we did not return yet, then particle and ancestor are not relatives
        return false;
}

Double_t UpsilonToTauTauAnalyzer::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

Double_t UpsilonToTauTauAnalyzer::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

