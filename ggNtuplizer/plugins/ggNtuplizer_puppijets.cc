#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"


using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

Int_t          npuppijet_;
vector<float>  puppijetPt_;
vector<float>  puppijetEn_;
vector<float>  puppijetEta_;
vector<float>  puppijetPhi_;
vector<float>  puppijetRawPt_;
vector<float>  puppijetRawEn_;
vector<float>  puppijetMt_;
vector<float>  puppijetArea_;
vector<float>  puppijetLeadTrackPt_;
vector<float>  puppijetLeadTrackEta_;
vector<float>  puppijetLeadTrackPhi_;
vector<int>    puppijetLepTrackPID_;
vector<float>  puppijetLepTrackPt_;
vector<float>  puppijetLepTrackEta_;
vector<float>  puppijetLepTrackPhi_;
vector<float>  puppijetCHF_;
vector<float>  puppijetNHF_;
vector<float>  puppijetCEF_;
vector<float>  puppijetNEF_;
vector<int>    puppijetNCH_;
vector<int>    puppijetNNP_;
vector<float>  puppijetMUF_;
vector<float>  puppijetHFHAE_;
vector<float>  puppijetHFEME_;
vector<int>    puppijetNConstituents_;
vector<float>  puppijetVtxPt_;
vector<float>  puppijetVtxMass_;
vector<float>  puppijetVtxNtrks_;
vector<float>  puppijetVtx3DVal_;
vector<float>  puppijetVtx3DSig_;
vector<float>  puppijetCSV2BjetTags_;
vector<float>  puppijetDeepCSVTags_b_;
vector<float>  puppijetDeepCSVTags_bb_;
vector<float>  puppijetDeepCSVTags_c_;
vector<float>  puppijetDeepCSVTags_udsg_;
vector<int>    puppijetPartonID_;
vector<int>    puppijetHadFlvr_;
vector<bool>   puppijetPFLooseId_;
vector<int>    puppijetID_; 
vector<float>  puppijetPUID_;
vector<int>    puppijetPUFullID_;
vector<float>  puppijetJECUnc_;
vector<float>  puppijetP4Smear_;
vector<float>  puppijetP4SmearUp_;
vector<float>  puppijetP4SmearDo_;
vector<ULong64_t> puppijetFiredTrgs_;
//gen-info for ak4
vector<float>  puppijetGenJetEn_;
vector<float>  puppijetGenJetPt_;
vector<float>  puppijetGenJetEta_;
vector<float>  puppijetGenJetPhi_;
vector<int>    puppijetGenPartonID_;
vector<float>  puppijetGenEn_;
vector<float>  puppijetGenPt_;
vector<float>  puppijetGenEta_;
vector<float>  puppijetGenPhi_;
vector<int>    puppijetGenPartonMomID_;

void ggNtuplizer::branchesPuppiJets(TTree* tree) {
  
  tree->Branch("npuppijet",                &npuppijet_);
  tree->Branch("puppijetPt",               &puppijetPt_);
  tree->Branch("puppijetEn",               &puppijetEn_);
  tree->Branch("puppijetEta",              &puppijetEta_);
  tree->Branch("puppijetPhi",              &puppijetPhi_);
  tree->Branch("puppijetRawPt",            &puppijetRawPt_);
  tree->Branch("puppijetRawEn",            &puppijetRawEn_);
  tree->Branch("puppijetMt",               &puppijetMt_);
  tree->Branch("puppijetArea",             &puppijetArea_);
  tree->Branch("puppijetLeadTrackPt",      &puppijetLeadTrackPt_);
  tree->Branch("puppijetLeadTrackEta",     &puppijetLeadTrackEta_);
  tree->Branch("puppijetLeadTrackPhi",     &puppijetLeadTrackPhi_);
  tree->Branch("puppijetLepTrackPID",      &puppijetLepTrackPID_);
  tree->Branch("puppijetLepTrackPt",       &puppijetLepTrackPt_);
  tree->Branch("puppijetLepTrackEta",      &puppijetLepTrackEta_);
  tree->Branch("puppijetLepTrackPhi",      &puppijetLepTrackPhi_);
  tree->Branch("puppijetCSV2BjetTags",     &puppijetCSV2BjetTags_);
  tree->Branch("puppijetDeepCSVTags_b",    &puppijetDeepCSVTags_b_);
  tree->Branch("puppijetDeepCSVTags_bb",   &puppijetDeepCSVTags_bb_);
  tree->Branch("puppijetDeepCSVTags_c",    &puppijetDeepCSVTags_c_);
  tree->Branch("puppijetDeepCSVTags_udsg", &puppijetDeepCSVTags_udsg_);
  if (doGenParticles_){
    tree->Branch("puppijetPartonID",       &puppijetPartonID_);
    tree->Branch("puppijetHadFlvr",        &puppijetHadFlvr_);
    tree->Branch("puppijetGenJetEn",       &puppijetGenJetEn_);
    tree->Branch("puppijetGenJetPt",       &puppijetGenJetPt_);
    tree->Branch("puppijetGenJetEta",      &puppijetGenJetEta_);
    tree->Branch("puppijetGenJetPhi",      &puppijetGenJetPhi_);
    tree->Branch("puppijetGenPartonID",    &puppijetGenPartonID_);
    tree->Branch("puppijetGenEn",          &puppijetGenEn_);
    tree->Branch("puppijetGenPt",          &puppijetGenPt_);
    tree->Branch("puppijetGenEta",         &puppijetGenEta_);
    tree->Branch("puppijetGenPhi",         &puppijetGenPhi_);
    tree->Branch("puppijetGenPartonMomID", &puppijetGenPartonMomID_);
    tree->Branch("puppijetP4Smear",        &puppijetP4Smear_);
    tree->Branch("puppijetP4SmearUp",      &puppijetP4SmearUp_);
    tree->Branch("puppijetP4SmearDo",      &puppijetP4SmearDo_);
  }  
  tree->Branch("puppijetPFLooseId", &puppijetPFLooseId_);
  tree->Branch("puppijetID",        &puppijetID_);
  tree->Branch("puppijetPUID",      &puppijetPUID_);
  tree->Branch("puppijetPUFullID",  &puppijetPUFullID_);
  tree->Branch("puppijetJECUnc",    &puppijetJECUnc_);
  tree->Branch("puppijetFiredTrgs", &puppijetFiredTrgs_);
  tree->Branch("puppijetCHF",       &puppijetCHF_);
  tree->Branch("puppijetNHF",       &puppijetNHF_);
  tree->Branch("puppijetCEF",       &puppijetCEF_);
  tree->Branch("puppijetNEF",       &puppijetNEF_);
  tree->Branch("puppijetNCH",       &puppijetNCH_);
  tree->Branch("puppijetNNP",       &puppijetNNP_);
  tree->Branch("puppijetMUF",       &puppijetMUF_);
  tree->Branch("puppijetVtxPt",     &puppijetVtxPt_);
  tree->Branch("puppijetVtxMass",   &puppijetVtxMass_);
  tree->Branch("puppijetVtxNtrks",  &puppijetVtxNtrks_);
  tree->Branch("puppijetVtx3DVal",  &puppijetVtx3DVal_);
  tree->Branch("puppijetVtx3DSig",  &puppijetVtx3DSig_);
  if (development_) {
    tree->Branch("puppijetHFHAE",         &puppijetHFHAE_);
    tree->Branch("puppijetHFEME",         &puppijetHFEME_);
    tree->Branch("puppijetNConstituents", &puppijetNConstituents_);
  }

}

void ggNtuplizer::fillPuppiJets(const edm::Event& e, const edm::EventSetup& es) {

  puppijetPt_                                  .clear();
  puppijetEn_                                  .clear();
  puppijetEta_                                 .clear();
  puppijetPhi_                                 .clear();
  puppijetRawPt_                               .clear();
  puppijetRawEn_                               .clear();
  puppijetMt_                                  .clear();
  puppijetArea_                                .clear();
  puppijetLeadTrackPt_                         .clear();
  puppijetLeadTrackEta_                        .clear();
  puppijetLeadTrackPhi_                        .clear();
  puppijetLepTrackPt_                          .clear();
  puppijetLepTrackPID_                         .clear();
  puppijetLepTrackEta_                         .clear();
  puppijetLepTrackPhi_                         .clear();
  puppijetCSV2BjetTags_                        .clear();
  puppijetDeepCSVTags_b_                       .clear();
  puppijetDeepCSVTags_bb_                      .clear();
  puppijetDeepCSVTags_c_                       .clear();
  puppijetDeepCSVTags_udsg_                    .clear();
  puppijetPartonID_                            .clear();
  puppijetHadFlvr_                             .clear();
  puppijetPFLooseId_                           .clear();
  puppijetID_                                  .clear();
  puppijetPUID_                                .clear();
  puppijetPUFullID_                            .clear();
  puppijetJECUnc_                              .clear();
  puppijetP4Smear_                             .clear();
  puppijetP4SmearUp_                           .clear();
  puppijetP4SmearDo_                           .clear();
  puppijetFiredTrgs_                           .clear();
  puppijetCHF_                                 .clear();
  puppijetNHF_                                 .clear();
  puppijetCEF_                                 .clear();
  puppijetNEF_                                 .clear();
  puppijetNCH_                                 .clear();
  puppijetNNP_                                 .clear();
  puppijetMUF_                                 .clear();
  puppijetVtxPt_                               .clear();
  puppijetVtxMass_                             .clear();
  puppijetVtxNtrks_                            .clear();
  puppijetVtx3DVal_                            .clear();
  puppijetVtx3DSig_                            .clear();
  if (development_) {
    puppijetHFHAE_                               .clear();
    puppijetHFEME_                               .clear();
    puppijetNConstituents_                       .clear();
  }
  puppijetGenJetEn_.clear();
  puppijetGenJetPt_.clear();
  puppijetGenJetEta_.clear();
  puppijetGenJetPhi_.clear();
  puppijetGenPartonID_.clear();
  puppijetGenEn_.clear();
  puppijetGenPt_.clear();
  puppijetGenEta_.clear();
  puppijetGenPhi_.clear();
  puppijetGenPartonMomID_.clear();


  edm::Handle<edm::View<pat::Jet> > puppijetHandle;
  e.getByToken(puppijetsAK4Label_, puppijetHandle);

  if (!puppijetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4 puppi) in event";
    return;
  }


  edm::Handle<vector<reco::GenParticle> > genParticlesPuppiHandle;
  if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesPuppiHandle);
  
  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  float rho = *(rhoHandle.product());
  
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  if (!vtxHandle.isValid()) edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";
  
  // Accessing the JEC uncertainties 
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *puppijecUnc=0;
  puppijecUnc = new JetCorrectionUncertainty(JetCorPar);
  npuppijet_ = 0;

  for (edm::View<pat::Jet>::const_iterator iPuppiJet = puppijetHandle->begin(); iPuppiJet != puppijetHandle->end(); ++iPuppiJet) {

    if (iPuppiJet->pt() < 20) continue;
    puppijetPt_.push_back(    iPuppiJet->pt());
    puppijetEn_.push_back(    iPuppiJet->energy());
    puppijetEta_.push_back(   iPuppiJet->eta());
    puppijetPhi_.push_back(   iPuppiJet->phi());
    puppijetRawPt_.push_back( (*iPuppiJet).correctedJet("Uncorrected").pt());
    puppijetRawEn_.push_back( (*iPuppiJet).correctedJet("Uncorrected").energy());
    puppijetMt_.push_back(    iPuppiJet->mt());
    puppijetArea_.push_back(  iPuppiJet->jetArea());

    puppijetCEF_.push_back(   iPuppiJet->chargedEmEnergyFraction());

    puppijetNEF_.push_back(   iPuppiJet->neutralEmEnergyFraction());
    puppijetCHF_.push_back(   iPuppiJet->chargedHadronEnergyFraction());
    puppijetNHF_.push_back(   iPuppiJet->neutralHadronEnergyFraction());
    puppijetNCH_.push_back(   iPuppiJet->chargedMultiplicity());
    puppijetNNP_.push_back(   iPuppiJet->neutralMultiplicity());
    puppijetMUF_.push_back(   iPuppiJet->muonEnergyFraction());

    if (development_) {
      puppijetHFHAE_.push_back( iPuppiJet->HFHadronEnergy());
      puppijetHFEME_.push_back( iPuppiJet->HFEMEnergy());
      puppijetNConstituents_.push_back(iPuppiJet->numberOfDaughters());
    }

    if (fabs(iPuppiJet->eta()) < 5.2) {
      puppijecUnc->setJetEta(iPuppiJet->eta());
      puppijecUnc->setJetPt(iPuppiJet->pt()); // here you must use the CORRECTED jet pt
      puppijetJECUnc_.push_back(puppijecUnc->getUncertainty(true));
    } else {
      puppijetJECUnc_.push_back(-1.);
    }

    puppijetFiredTrgs_.push_back(matchJetTriggerFilters(iPuppiJet->pt(), iPuppiJet->eta(), iPuppiJet->phi()));    

    //Searching for leading track and lepton
    float leadTrkPt  = -99;
    float leadTrkEta = -99;
    float leadTrkPhi = -99;
    int   lepTrkPID  = -99;
    float lepTrkPt   = -99;
    float lepTrkEta  = -99;
    float lepTrkPhi  = -99;

    for (unsigned id = 0; id < iPuppiJet->getJetConstituents().size(); id++) {

      const edm::Ptr<reco::Candidate> daughter = iPuppiJet->getJetConstituents().at(id);

      if (daughter.isNonnull() && daughter.isAvailable()) {
	if (daughter->charge() != 0 && daughter->pt() > leadTrkPt) {
	  leadTrkPt  = daughter->pt();
	  leadTrkEta = daughter->eta();
	  leadTrkPhi = daughter->phi();
	}

	if (abs(daughter->pdgId()) == 11 || abs(daughter->pdgId()) == 13) {
	  if (daughter->pt() > lepTrkPt) {
	    lepTrkPID = daughter->pdgId();
	    lepTrkPt  = daughter->pt();
	    lepTrkEta = daughter->eta();
	    lepTrkPhi = daughter->phi();
	  }
	}
      }
    }

    puppijetLeadTrackPt_ .push_back(leadTrkPt);
    puppijetLeadTrackEta_.push_back(leadTrkEta);
    puppijetLeadTrackPhi_.push_back(leadTrkPhi);
    puppijetLepTrackPID_ .push_back(lepTrkPID);
    puppijetLepTrackPt_  .push_back(lepTrkPt);
    puppijetLepTrackEta_ .push_back(lepTrkEta);
    puppijetLepTrackPhi_ .push_back(lepTrkPhi);    
    //jetVtxPt_       .push_back(sqrt(pow(iPuppiJet->userFloat("vtxPx"),2)+pow(iPuppiJet->userFloat("vtxPy"),2)));
    //jetVtxMass_     .push_back(iPuppiJet->userFloat("vtxMass"));
    //jetVtxNtrks_    .push_back(iPuppiJet->userFloat("vtxNtracks"));
    //jetVtx3DVal_    .push_back(iPuppiJet->userFloat("vtx3DVal"));
    //jetVtx3DSig_    .push_back(iPuppiJet->userFloat("vtx3DSig"));
    
    //b/c-tagging
    puppijetCSV2BjetTags_    .push_back(iPuppiJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    puppijetDeepCSVTags_b_   .push_back(iPuppiJet->bDiscriminator("pfDeepCSVJetTags:probb"));
    puppijetDeepCSVTags_bb_  .push_back(iPuppiJet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    puppijetDeepCSVTags_c_   .push_back(iPuppiJet->bDiscriminator("pfDeepCSVJetTags:probc"));
    puppijetDeepCSVTags_udsg_.push_back(iPuppiJet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
  
    //parton id
    puppijetPartonID_.push_back(iPuppiJet->partonFlavour());
    puppijetHadFlvr_.push_back(iPuppiJet->hadronFlavour());

    //jet PF Loose ID
    double NHF      = iPuppiJet->neutralHadronEnergyFraction();
    double NEMF     = iPuppiJet->neutralEmEnergyFraction();
    double NumConst = iPuppiJet->chargedMultiplicity()+iPuppiJet->neutralMultiplicity();
    double CHF      = iPuppiJet->chargedHadronEnergyFraction();
    double CHM      = iPuppiJet->chargedMultiplicity();
    double CEMF     = iPuppiJet->chargedEmEnergyFraction();
    double NNP      = iPuppiJet->neutralMultiplicity();

    bool loosepuppijetID = false;    
    bool tightpuppijetID = false;
    if (fabs(iPuppiJet->eta()) <= 2.7) {
      loosepuppijetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(iPuppiJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(iPuppiJet->eta())>2.4);
      tightpuppijetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(iPuppiJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(iPuppiJet->eta())>2.4);
    } else if (fabs(iPuppiJet->eta()) <= 3.0) {
      loosepuppijetID = (NEMF>0.01 && NHF<0.98 && NNP>2);
      tightpuppijetID = (NEMF>0.01 && NHF<0.98 && NNP>2);
    } else {
      loosepuppijetID = (NEMF<0.90 && NNP>10); 
      tightpuppijetID = (NEMF<0.90 && NNP>10);
    }
    puppijetPFLooseId_.push_back(loosepuppijetID);
    Int_t puppijetIDdecision = 0;
    if (loosepuppijetID) puppijetIDdecision += pow(2, 1);
    if (tightpuppijetID) puppijetIDdecision += pow(2, 2);
    puppijetID_.push_back(puppijetIDdecision);    

    // PUpuppijet ID from slimmedpuppijets
    //puppijetPUID_.push_back(iPuppiJet->userFloat("pileupJetId:fullDiscriminant")); //Not available in 94X
    //puppijetPUFullID_.push_back(iPuppiJet->userInt("pileupJetId:fullId"));

    // gen jet and parton
    if (doGenParticles_ && genParticlesPuppiHandle.isValid()) {
      int puppijetGenPartonID    = -99;
      int puppijetGenPartonMomID = -99;
      float puppijetGenEn        = -999.;
      float puppijetGenPt        = -999.;
      float puppijetGenEta       = -999.;
      float puppijetGenPhi       = -999.;      
      if ((*iPuppiJet).genParton()) {
	puppijetGenPartonID = (*iPuppiJet).genParton()->pdgId();
	puppijetGenEn = (*iPuppiJet).genParton()->energy();
	puppijetGenPt = (*iPuppiJet).genParton()->pt();
	puppijetGenEta = (*iPuppiJet).genParton()->eta();
	puppijetGenPhi = (*iPuppiJet).genParton()->phi();
  //cout<<"pdgID"<<puppijetGenPartonID<<endl;

	// if ((*iPuppiJet).genParton()->mother()) {
	//   puppijetGenPartonMomID = (*iPuppiJet).genParton()->mother()->pdgId();
	// }
      }
      
      puppijetGenPartonID_.push_back(puppijetGenPartonID);
      puppijetGenPartonMomID_.push_back(puppijetGenPartonMomID);
      puppijetGenEn_ .push_back(puppijetGenEn);
      puppijetGenPt_ .push_back(puppijetGenPt);
      puppijetGenEta_ .push_back(puppijetGenEta);
      puppijetGenPhi_ .push_back(puppijetGenPhi);
      
      float puppijetGenJetEn  = -999.;
      float puppijetGenJetPt  = -999.;
      float puppijetGenJetEta = -999.;
      float puppijetGenJetPhi = -999.;
      if ((*iPuppiJet).genJet()) {
	puppijetGenJetEn = (*iPuppiJet).genJet()->energy();
	puppijetGenJetPt = (*iPuppiJet).genJet()->pt();
	puppijetGenJetEta = (*iPuppiJet).genJet()->eta();
	puppijetGenJetPhi = (*iPuppiJet).genJet()->phi();
      }
      puppijetGenJetEn_.push_back(puppijetGenJetEn);
      puppijetGenJetPt_.push_back(puppijetGenJetPt);
      puppijetGenJetEta_.push_back(puppijetGenJetEta);
      puppijetGenJetPhi_.push_back(puppijetGenJetPhi);
      
      // access jet resolution             
      JME::JetParameters parameters;
      parameters.setJetPt(iPuppiJet->pt()).setJetEta(iPuppiJet->eta()).setRho(rho);
      float puppijetResolution = PuppijetResolution_.getResolution(parameters);

      edm::Service<edm::RandomNumberGenerator> rng;
      if (!rng.isAvailable()) edm::LogError("JET : random number generator is missing !");
      CLHEP::HepRandomEngine & engine = rng->getEngine( e.streamID() );
      float rnd = CLHEP::RandGauss::shoot(&engine, 0., puppijetResolution);

      float puppijetResolutionSF   = PuppijetResolutionSF_.getScaleFactor(parameters);
      float puppijetResolutionSFUp = PuppijetResolutionSF_.getScaleFactor(parameters, Variation::UP);
      float puppijetResolutionSFDo = PuppijetResolutionSF_.getScaleFactor(parameters, Variation::DOWN);

      float puppijetP4Smear   = -1.;
      float puppijetP4SmearUp = -1.;
      float puppijetP4SmearDo = -1.;
      if (puppijetGenJetPt > 0 && deltaR(iPuppiJet->eta(), iPuppiJet->phi(), puppijetGenJetEta, puppijetGenJetPhi) < 0.2 && fabs(iPuppiJet->pt()-puppijetGenJetPt) < 3*puppijetResolution*iPuppiJet->pt()) {
	puppijetP4Smear   = 1. + (puppijetResolutionSF   - 1.)*(iPuppiJet->pt() - puppijetGenJetPt)/iPuppiJet->pt();
	puppijetP4SmearUp = 1. + (puppijetResolutionSFUp - 1.)*(iPuppiJet->pt() - puppijetGenJetPt)/iPuppiJet->pt();
	puppijetP4SmearDo = 1. + (puppijetResolutionSFDo - 1.)*(iPuppiJet->pt() - puppijetGenJetPt)/iPuppiJet->pt();
      } else {
	puppijetP4Smear   = 1. + rnd*sqrt(max(pow(puppijetResolutionSF,   2)-1, 0.));
        puppijetP4SmearUp = 1. + rnd*sqrt(max(pow(puppijetResolutionSFUp, 2)-1, 0.));
	puppijetP4SmearDo = 1. + rnd*sqrt(max(pow(puppijetResolutionSFDo, 2)-1, 0.));
      }
      puppijetP4Smear_  .push_back(puppijetP4Smear);
      puppijetP4SmearUp_.push_back(puppijetP4SmearUp);
      puppijetP4SmearDo_.push_back(puppijetP4SmearDo);      
    }
    
    npuppijet_++;
    //cout<<"count no of puppi jets="<<npuppijet_<<endl;
  }


  delete puppijecUnc;
}
