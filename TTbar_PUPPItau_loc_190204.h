
#ifndef TTbar_PUPPItau_loc_190204_h
#define TTbar_PUPPItau_loc_190204_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//Header file for the classes stored in the TTree if any.
#include <TVector.h>
#include "TString.h"
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>

#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
//#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include <map>
//#include "TRFIOFile.h"
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <TLatex.h>
#include <Riostream.h>
#include <set>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TKDE.h"

#include <iostream>
#include <array>

#ifdef __MAKECINT__ 
#pragma link C++ class vector<bool>+;
#pragma link C++ class map<TString,TH1D*>+;
#pragma link C++ class vector<unsigned long long>+;
#pragma link C++ class vector<vector<unsigned long long> >+;
#pragma link C++ class vector<ULong64_t>+;
#pragma link C++ class vector<vector<ULong64_t> >+;
#pragma link C++ class vector<long long>+;
#pragma link C++ class vector<vector<long long> >+;
#pragma link C++ class vector<Long64_t>+;
#pragma link C++ class vector<vector<Long64_t> >+;
#endif

//---------------------------------------------


using namespace std;
using namespace TMath; // To include TMath class functions
using namespace ROOT;

class TTbar_PUPPItau_loc_190204 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           v_event;
   Int_t           n_Vtx;
   Int_t           n_GoodVtx;
   Int_t           n_TrksPV;
   Bool_t          is_PVGood;
   Float_t         v_vtx;
   Float_t         v_vty;
   Float_t         v_vtz;
   Float_t         rho_Central;
   Int_t           n_MC;
   vector<int>     *mc_PID;
   vector<float>   *mc_Vtx;
   vector<float>   *mc_Vty;
   vector<float>   *mc_Vtz;
   vector<float>   *mc_Pt;
   vector<float>   *mc_Mass;
   vector<float>   *mc_Eta;
   vector<float>   *mc_Phi;
   vector<float>   *mc_E;
   vector<float>   *mc_Et;
   vector<int>     *mc_GMomPID;
   vector<int>     *mc_MomPID;
   vector<float>   *mc_MomPt;
   vector<float>   *mc_MomMass;
   vector<float>   *mc_MomEta;
   vector<float>   *mc_MomPhi;
   vector<unsigned short> *mc_StatusFlag;
   vector<int>     *mc_Parentage;
   vector<int>     *mc_Status;
   vector<float>   *mc_CalIsoDR03;
   vector<float>   *mc_TrkIsoDR03;
   vector<float>   *mc_CalIsoDR04;
   vector<float>   *mc_TrkIsoDR04;
   Float_t         gen_MET;
   Float_t         gen_METPhi;
   Float_t         pf_MET;
   Float_t         pf_METPhi;
   Float_t         pf_METsumEt;
   Float_t         pf_METmEtSig;
   Float_t         pf_METSig;
   Int_t           n_Ele;
   vector<int>     *ele_Charge;
   vector<float>   *ele_En;
   vector<float>   *ele_D0;
   vector<float>   *ele_Dz;
   vector<float>   *ele_Pt;
   vector<float>   *ele_Eta;
   vector<float>   *ele_Phi;
   vector<float>   *ele_R9;
   vector<float>   *ele_SCEta;
   vector<float>   *ele_SCPhi;
   vector<float>   *ele_HoverE;
   vector<float>   *ele_EoverP;
   vector<float>   *ele_EoverPout;
   vector<float>   *ele_EoverPInv;
   vector<float>   *ele_dEtaAtVtx;
   vector<float>   *ele_dPhiAtVtx;
   vector<float>   *ele_SigmaIEtaIEtaFull5x5;
   vector<float>   *ele_SigmaIPhiIPhiFull5x5;
   vector<int>     *ele_ConvVeto;
   vector<int>     *ele_MissHits;
   vector<float>   *ele_PFChIso;
   vector<float>   *ele_PFPhoIso;
   vector<float>   *ele_PFNeuIso;
   vector<float>   *ele_PFMiniIso;
   vector<float>   *ele_dEtaseedAtVtx;
   Int_t           N_Mu;
   vector<float>   *mu_Pt;
   vector<float>   *mu_En;
   vector<float>   *mu_Eta;
   vector<float>   *mu_Phi;
   vector<int>     *mu_Charge;
   vector<unsigned short> *mu_IDbit;
   vector<float>   *mu_D0;
   vector<float>   *mu_Dz;
   vector<float>   *mu_Chi2NDF;
   vector<float>   *mu_InnerD0;
   vector<float>   *mu_InnerDz;
   vector<float>   *mu_InnervalidFraction;
   vector<float>   *mu_segmentCompatibility;
   vector<float>   *mu_chi2LocalPosition;
   vector<float>   *mu_trkKink;
   vector<float>   *mu_PFChIso;
   vector<float>   *mu_PFPhoIso;
   vector<float>   *mu_PFNeuIso;
   vector<float>   *mu_PFMiniIso;
   vector<int>     *mu_TrkLayers;
   vector<int>     *mu_PixelLayers;
   vector<int>     *mu_PixelHits;
   vector<int>     *mu_MuonHits;
   vector<int>     *mu_Stations;
   vector<int>     *mu_Matches;
   vector<int>     *mu_TrkQuality;
   vector<float>   *mu_IsoTrk;
   Int_t           n_Jet;
   vector<float>   *jet_Pt;
   vector<float>   *jet_En;
   vector<float>   *jet_Eta;
   vector<float>   *jet_Phi;
   vector<float>   *jet_Area;
   vector<float>   *jet_Mt;
   vector<float>   *jet_CSV2BJetTags;
   vector<float>   *jet_JetProbabilityBJetTags;
   vector<float>   *jet_pfCombinedMVAV2BJetTags;
   vector<float>   *jet_DeepCSVTags_b;
   vector<float>   *jet_DeepCSVTags_bb;
   vector<float>   *jet_DeepCSVTags_c;
   vector<float>   *jet_DeepCSVTags_cc;
   vector<float>   *jet_DeepCSVTags_udsg;
   vector<int>     *jet_PartonID;
   vector<int>     *jet_HadFlvr;
   vector<float>   *jet_GenJetEn;
   vector<float>   *jet_GenJetPt;
   vector<float>   *jet_GenJetEta;
   vector<float>   *jet_GenJetPhi;
   vector<int>     *jet_GenPartonID;
   vector<float>   *jet_GenEn;
   vector<float>   *jet_GenPt;
   vector<float>   *jet_GenEta;
   vector<float>   *jet_GenPhi;
   vector<int>     *jet_GenPartonMomID;
   vector<bool>    *jet_PFLooseId;
   vector<int>     *jet_ID;
   vector<float>   *jet_PUID;
   vector<int>     *jet_PUFullID;
   vector<float>   *jet_CHF;
   vector<float>   *jet_NHF;
   vector<float>   *jet_CEF;
   vector<float>   *jet_NEF;
   vector<int>     *jet_NCH;
   vector<int>     *jet_NNP;
   vector<float>   *jet_MUF;
   Int_t           N_AK8Jet;
   Int_t           n_AK8Jetpuppi;
   vector<float>   *AK8_JetPt;
   vector<float>   *AK8_JetEn;
   vector<float>   *AK8_JetEta;
   vector<float>   *AK8_JetPhi;
   vector<float>   *AK8_JetMass;
   vector<float>   *AK8_Jet_tau1;
   vector<float>   *AK8_Jet_tau2;
   vector<float>   *AK8_Jet_tau3;
   vector<float>   *AK8_Jet_CHStau1;
   vector<float>   *AK8_Jet_CHStau2;
   vector<float>   *AK8_Jet_CHStau3;
   vector<float>   *AK8_Jet_CHStau4;
   vector<float>   *AK8_JetCHF;
   vector<float>   *AK8_JetNHF;
   vector<float>   *AK8_JetCEF;
   vector<float>   *AK8_JetNEF;
   vector<int>     *AK8_JetNCH;
   vector<int>     *AK8_JetNNP;
   vector<float>   *AK8_JetMUF;
   vector<int>     *AK8_Jetnconstituents;
   vector<bool>    *AK8_JetPFLooseId;
   vector<bool>    *AK8_JetPFTightLepVetoId;
   vector<float>   *AK8_JetSoftDropMass;
   vector<float>   *AK8_JetSoftDropMassCorr;
   vector<float>   *AK8_JetPrunedMassCorr;
   vector<float>   *AK8_JetL2L3corr;
   vector<float>   *AK8_puppiSDMassL2L3Corr;
   vector<float>   *AK8_puppiSDL2L3Corr;
   vector<float>   *AK8_JetPrunedMass;
   vector<float>   *AK8_JetpfBoostedDSVBTag;
   vector<float>   *AK8_JetCSV;
   vector<float>   *AK8_JetDSVnewV4;
   vector<float>   *AK8_puppiPt;
   vector<float>   *AK8_puppiMass;
   vector<float>   *AK8_puppiEta;
   vector<float>   *AK8_puppiPhi;
   vector<float>   *AK8_puppiTau1;
   vector<float>   *AK8_puppiTau2;
   vector<float>   *AK8_puppiTau3;
   vector<float>   *AK8_puppiTau4;
   vector<float>   *AK8_puppiSDMass;
   vector<int>     *AK8_JetPartonID;
   vector<int>     *AK8_JetHadFlvr;
   vector<int>     *AK8_JetGenJetIndex;
   vector<float>   *AK8_JetGenJetEn;
   vector<float>   *AK8_JetGenJetPt;
   vector<float>   *AK8_JetGenJetEta;
   vector<float>   *AK8_JetGenJetPhi;
   vector<int>     *AK8_JetGenPartonID;
   vector<float>   *AK8_JetGenEn;
   vector<float>   *AK8_JetGenPt;
   vector<float>   *AK8_JetGenEta;
   vector<float>   *AK8_JetGenPhi;
   vector<int>     *AK8_JetGenPartonMomID;
   vector<int>     *n_AK8SDSJ;
   vector<vector<float> > *AK8_SDSJPt;
   vector<vector<float> > *AK8_SDSJEta;
   vector<vector<float> > *AK8_SDSJPhi;
   vector<vector<float> > *AK8_SDSJMass;
   vector<vector<float> > *AK8_SDSJE;
   vector<vector<int> > *AK8_SDSJCharge;
   vector<vector<int> > *AK8_SDSJFlavour;
   vector<vector<float> > *AK8_SDSJCSV;
   vector<int>     *n_AK8puppiSDSJ;
   vector<vector<float> > *AK8_puppiSDSJPt;
   vector<vector<float> > *AK8_puppiSDSJEta;
   vector<vector<float> > *AK8_puppiSDSJPhi;
   vector<vector<float> > *AK8_puppiSDSJMass;
   vector<vector<float> > *AK8_puppiSDSJE;
   vector<vector<int> > *AK8_puppiSDSJCharge;
   vector<vector<int> > *AK8_puppiSDSJFlavour;
   vector<vector<float> > *AK8_puppiSDSJCSV;

   // List of branches
   TBranch        *b_v_event;   //!
   TBranch        *b_n_Vtx;   //!
   TBranch        *b_n_GoodVtx;   //!
   TBranch        *b_n_TrksPV;   //!
   TBranch        *b_is_PVGood;   //!
   TBranch        *b_v_vtx;   //!
   TBranch        *b_v_vty;   //!
   TBranch        *b_v_vtz;   //!
   TBranch        *b_rho_Central;   //!
   TBranch        *b_n_MC;   //!
   TBranch        *b_mc_PID;   //!
   TBranch        *b_mc_Vtx;   //!
   TBranch        *b_mc_Vty;   //!
   TBranch        *b_mc_Vtz;   //!
   TBranch        *b_mc_Pt;   //!
   TBranch        *b_mc_Mass;   //!
   TBranch        *b_mc_Eta;   //!
   TBranch        *b_mc_Phi;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_Et;   //!
   TBranch        *b_mc_GMomPID;   //!
   TBranch        *b_mc_MomPID;   //!
   TBranch        *b_mc_MomPt;   //!
   TBranch        *b_mc_MomMass;   //!
   TBranch        *b_mc_MomEta;   //!
   TBranch        *b_mc_MomPhi;   //!
   TBranch        *b_mc_StatusFlag;   //!
   TBranch        *b_mc_Parentage;   //!
   TBranch        *b_mc_Status;   //!
   TBranch        *b_mc_CalIsoDR03;   //!
   TBranch        *b_mc_TrkIsoDR03;   //!
   TBranch        *b_mc_CalIsoDR04;   //!
   TBranch        *b_mc_TrkIsoDR04;   //!
   TBranch        *b_gen_MET;   //!
   TBranch        *b_gen_METPhi;   //!
   TBranch        *b_pf_MET;   //!
   TBranch        *b_pf_METPhi;   //!
   TBranch        *b_pf_METsumEt;   //!
   TBranch        *b_pf_METmEtSig;   //!
   TBranch        *b_pf_METSig;   //!
   TBranch        *b_n_Ele;   //!
   TBranch        *b_ele_Charge;   //!
   TBranch        *b_ele_En;   //!
   TBranch        *b_ele_D0;   //!
   TBranch        *b_ele_Dz;   //!
   TBranch        *b_ele_Pt;   //!
   TBranch        *b_ele_Eta;   //!
   TBranch        *b_ele_Phi;   //!
   TBranch        *b_ele_R9;   //!
   TBranch        *b_ele_SCEta;   //!
   TBranch        *b_ele_SCPhi;   //!
   TBranch        *b_ele_HoverE;   //!
   TBranch        *b_ele_EoverP;   //!
   TBranch        *b_ele_EoverPout;   //!
   TBranch        *b_ele_EoverPInv;   //!
   TBranch        *b_ele_dEtaAtVtx;   //!
   TBranch        *b_ele_dPhiAtVtx;   //!
   TBranch        *b_ele_SigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_ele_SigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_ele_ConvVeto;   //!
   TBranch        *b_ele_MissHits;   //!
   TBranch        *b_ele_PFChIso;   //!
   TBranch        *b_ele_PFPhoIso;   //!
   TBranch        *b_ele_PFNeuIso;   //!
   TBranch        *b_ele_PFMiniIso;   //!
   TBranch        *b_ele_dEtaseedAtVtx;   //!
   TBranch        *b_N_Mu;   //!
   TBranch        *b_mu_Pt;   //!
   TBranch        *b_mu_En;   //!
   TBranch        *b_mu_Eta;   //!
   TBranch        *b_mu_Phi;   //!
   TBranch        *b_mu_Charge;   //!
   TBranch        *b_mu_IDbit;   //!
   TBranch        *b_mu_D0;   //!
   TBranch        *b_mu_Dz;   //!
   TBranch        *b_mu_Chi2NDF;   //!
   TBranch        *b_mu_InnerD0;   //!
   TBranch        *b_mu_InnerDz;   //!
   TBranch        *b_mu_InnervalidFraction;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_chi2LocalPosition;   //!
   TBranch        *b_mu_trkKink;   //!
   TBranch        *b_mu_PFChIso;   //!
   TBranch        *b_mu_PFPhoIso;   //!
   TBranch        *b_mu_PFNeuIso;   //!
   TBranch        *b_mu_PFMiniIso;   //!
   TBranch        *b_mu_TrkLayers;   //!
   TBranch        *b_mu_PixelLayers;   //!
   TBranch        *b_mu_PixelHits;   //!
   TBranch        *b_mu_MuonHits;   //!
   TBranch        *b_mu_Stations;   //!
   TBranch        *b_mu_Matches;   //!
   TBranch        *b_mu_TrkQuality;   //!
   TBranch        *b_mu_IsoTrk;   //!
   TBranch        *b_n_Jet;   //!
   TBranch        *b_jet_Pt;   //!
   TBranch        *b_jet_En;   //!
   TBranch        *b_jet_Eta;   //!
   TBranch        *b_jet_Phi;   //!
   TBranch        *b_jet_Area;   //!
   TBranch        *b_jet_Mt;   //!
   TBranch        *b_jet_CSV2BJetTags;   //!
   TBranch        *b_jet_JetProbabilityBJetTags;   //!
   TBranch        *b_jet_pfCombinedMVAV2BJetTags;   //!
   TBranch        *b_jet_DeepCSVTags_b;   //!
   TBranch        *b_jet_DeepCSVTags_bb;   //!
   TBranch        *b_jet_DeepCSVTags_c;   //!
   TBranch        *b_jet_DeepCSVTags_cc;   //!
   TBranch        *b_jet_DeepCSVTags_udsg;   //!
   TBranch        *b_jet_PartonID;   //!
   TBranch        *b_jet_HadFlvr;   //!
   TBranch        *b_jet_GenJetEn;   //!
   TBranch        *b_jet_GenJetPt;   //!
   TBranch        *b_jet_GenJetEta;   //!
   TBranch        *b_jet_GenJetPhi;   //!
   TBranch        *b_jet_GenPartonID;   //!
   TBranch        *b_jet_GenEn;   //!
   TBranch        *b_jet_GenPt;   //!
   TBranch        *b_jet_GenEta;   //!
   TBranch        *b_jet_GenPhi;   //!
   TBranch        *b_jet_GenPartonMomID;   //!
   TBranch        *b_jet_PFLooseId;   //!
   TBranch        *b_jet_ID;   //!
   TBranch        *b_jet_PUID;   //!
   TBranch        *b_jet_PUFullID;   //!
   TBranch        *b_jet_CHF;   //!
   TBranch        *b_jet_NHF;   //!
   TBranch        *b_jet_CEF;   //!
   TBranch        *b_jet_NEF;   //!
   TBranch        *b_jet_NCH;   //!
   TBranch        *b_jet_NNP;   //!
   TBranch        *b_jet_MUF;   //!
   TBranch        *b_N_AK8Jet;   //!
   TBranch        *b_n_AK8Jetpuppi;   //!
   TBranch        *b_AK8_JetPt;   //!
   TBranch        *b_AK8_JetEn;   //!
   TBranch        *b_AK8_JetEta;   //!
   TBranch        *b_AK8_JetPhi;   //!
   TBranch        *b_AK8_JetMass;   //!
   TBranch        *b_AK8_Jet_tau1;   //!
   TBranch        *b_AK8_Jet_tau2;   //!
   TBranch        *b_AK8_Jet_tau3;   //!
   TBranch        *b_AK8_Jet_CHStau1;   //!
   TBranch        *b_AK8_Jet_CHStau2;   //!
   TBranch        *b_AK8_Jet_CHStau3;   //!
   TBranch        *b_AK8_Jet_CHStau4;   //!
   TBranch        *b_AK8_JetCHF;   //!
   TBranch        *b_AK8_JetNHF;   //!
   TBranch        *b_AK8_JetCEF;   //!
   TBranch        *b_AK8_JetNEF;   //!
   TBranch        *b_AK8_JetNCH;   //!
   TBranch        *b_AK8_JetNNP;   //!
   TBranch        *b_AK8_JetMUF;   //!
   TBranch        *b_AK8_Jetnconstituents;   //!
   TBranch        *b_AK8_JetPFLooseId;   //!
   TBranch        *b_AK8_JetPFTightLepVetoId;   //!
   TBranch        *b_AK8_JetSoftDropMass;   //!
   TBranch        *b_AK8_JetSoftDropMassCorr;   //!
   TBranch        *b_AK8_JetPrunedMassCorr;   //!
   TBranch        *b_AK8_JetL2L3corr;   //!
   TBranch        *b_AK8_puppiSDMassL2L3Corr;   //!
   TBranch        *b_AK8_puppiSDL2L3Corr;   //!
   TBranch        *b_AK8_JetPrunedMass;   //!
   TBranch        *b_AK8_JetpfBoostedDSVBTag;   //!
   TBranch        *b_AK8_JetCSV;   //!
   TBranch        *b_AK8_JetDSVnewV4;   //!
   TBranch        *b_AK8_puppiPt;   //!
   TBranch        *b_AK8_puppiMass;   //!
   TBranch        *b_AK8_puppiEta;   //!
   TBranch        *b_AK8_puppiPhi;   //!
   TBranch        *b_AK8_puppiTau1;   //!
   TBranch        *b_AK8_puppiTau2;   //!
   TBranch        *b_AK8_puppiTau3;   //!
   TBranch        *b_AK8_puppiTau4;   //!
   TBranch        *b_AK8_puppiSDMass;   //!
   TBranch        *b_AK8_JetPartonID;   //!
   TBranch        *b_AK8_JetHadFlvr;   //!
   TBranch        *b_AK8_JetGenJetIndex;   //!
   TBranch        *b_AK8_JetGenJetEn;   //!
   TBranch        *b_AK8_JetGenJetPt;   //!
   TBranch        *b_AK8_JetGenJetEta;   //!
   TBranch        *b_AK8_JetGenJetPhi;   //!
   TBranch        *b_AK8_JetGenPartonID;   //!
   TBranch        *b_AK8_JetGenEn;   //!
   TBranch        *b_AK8_JetGenPt;   //!
   TBranch        *b_AK8_JetGenEta;   //!
   TBranch        *b_AK8_JetGenPhi;   //!
   TBranch        *b_AK8_JetGenPartonMomID;   //!
   TBranch        *b_n_AK8SDSJ;   //!
   TBranch        *b_AK8_SDSJPt;   //!
   TBranch        *b_AK8_SDSJEta;   //!
   TBranch        *b_AK8_SDSJPhi;   //!
   TBranch        *b_AK8_SDSJMass;   //!
   TBranch        *b_AK8_SDSJE;   //!
   TBranch        *b_AK8_SDSJCharge;   //!
   TBranch        *b_AK8_SDSJFlavour;   //!
   TBranch        *b_AK8_SDSJCSV;   //!
   TBranch        *b_n_AK8puppiSDSJ;   //!
   TBranch        *b_AK8_puppiSDSJPt;   //!
   TBranch        *b_AK8_puppiSDSJEta;   //!
   TBranch        *b_AK8_puppiSDSJPhi;   //!
   TBranch        *b_AK8_puppiSDSJMass;   //!
   TBranch        *b_AK8_puppiSDSJE;   //!
   TBranch        *b_AK8_puppiSDSJCharge;   //!
   TBranch        *b_AK8_puppiSDSJFlavour;   //!
   TBranch        *b_AK8_puppiSDSJCSV;   //!

   TTbar_PUPPItau_loc_190204(TString inputFile);
   virtual ~TTbar_PUPPItau_loc_190204();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(TString OutputFileName, int y);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

//===================User Defined Function List=========================
	virtual  void    Clear_Vector() ;
	virtual  float  delta_phi(float phi1,float phi2);

	// Histogram Defining & Filling Functions
	//1. For Gen Level Objects
	virtual void      DefineMC_NPtEta_Histo();
	virtual void      dRHisto_MCObject() ;
	virtual void      Fill_MC_PtEta_Histo(int id, int en);
	// virtual  float    Fill_mcMu_recoMu_plots(int Y, int Z, int lvl);  
	virtual void     dR_Plots_Genlvl() ;
	virtual void     GenvsRecoMass(int obj, int ID ) ;
	virtual void     GenvsRecoObj(int obj, int ID ) ;

	// 2. For Reco Objects  at Preselection
	virtual void  	  Define_2DMass_Histo() ;
	virtual void      Define_NPtEta_Histo();  // for reco objects definitions at preselection level
	virtual void      dRHisto_RecoObject();  // for reco objects dR Plots  at preselection level
	virtual void      dRHisto_MCRecoObject();
	virtual void      Fill_NPtEta_Histo(int id, int en, int idx);
	virtual void      RecoPlots_dRHisto();
	virtual void      Fill_RecoObject();
	virtual float    Fill_MET_var(int lep, int lvl);  // X dlt
	//3. For Tag jets before Category Selection
	virtual void      Define_Tag_Jet_Histo() ;       // for  tag jets before Category Selection
	virtual void      dR_tagjetHisto();
	virtual void      Define_Reco_tagjetHisto() ; 
	virtual void       Fill_Puppi_jet(int jet, int idx);

	//---- Object Selection Functions------
	virtual  bool   Cut_Muon(int c_muon);
	virtual  void   Cut_bjet(int c_jet, int wp);
	virtual  void   Cut_AK8jet(int c_jet, int var);

	// -----dR & dPt Calculating Functions--------
	//1. For MC Objects
	virtual  float   dR_mcbjet(int Y, int Z, int X, int idx, int lvl);
	virtual  float   dR_mc_mub(int Y, int Z, int A);
	virtual  void   dR_mc_bqi(int Y, int Z, int A);
	virtual  void   dR_mc_qiqj(int Y, int Z, int A);
	virtual  float  dR_mcAK8(int Y, int Z, int X, int idx);
	virtual  float   dPt_Gen_jetlep(int jet, int lep, int idx);

	//2. For Reco Objects
	virtual  float   dR_mu(int Y, int Z_jet, int idx, int lvl);
	virtual  float   dR_mu_AK8(int Y, int Z_jet, int idx, int lvl );
	virtual  float   dRPlots_bjet(int Y, int Z, int lvl) ;
	virtual  float   dR_AK8jet(int Y, int Z, int id8, int lvl);
	virtual  float   dR_AK8bjet(int Y, int Z, int id8, int idb, int lvl);
	virtual  void   dR_qjet_objects(int qjet, int obj, int idq, int idob);  
	virtual  float   dPt_lep(int jet, int lep, int idx, int lvl);

	//---- Jet Tagging Selections & Plots Functions------
	// 0. Optimisation
	virtual  void    Higgs_Optimisation(int var) ;
	virtual  void    Top_Optimisation(int var) ;

	//1. Selections
	virtual  void    Higgs_selection(int var) ;
	virtual  void    Top_selection(int var) ;
	virtual  void    Wjet_selection(int var) ;
	virtual  void    Fatjet_selection() ;
	virtual  int     WtagMatch( int Wjet);
	// 2. Plot Function
	virtual  void     Wjet_Plots( int bjet );
	virtual  void     Topjet_Plots(int var ) ;
	virtual  void      Higgsjet_Plots( int bjet) ;
	virtual  void     Fatjet_Plots(int topbjet) ; 
	virtual void      TagJets_dRPlots() ;

	// ---------Signal Categorization , Histogram & Plot Function----------
	//1. Categorization
	virtual  void     Wtag0_Category(); 
	virtual  void     Wtag1_Category(); 
	// 2. Histogram definition
	virtual void      Category_Object_Histo() ;
	virtual void      Category_Object_MtHisto() ;  
	virtual void      Category_Object_dRHisto() ;  
	//3. Plot Functions  
	virtual  void     top_fatjet_Plots();
	virtual  void     top_Wjet_Plots() ; 
	virtual  void     W_fatjet_Plots() ;
	virtual  void     WW_lvbjet_Plots() ;
	virtual void      Higgs_lbjet_Plots() ;
	virtual void      CatI_Objects_Plots() ;  
	virtual void      CatII_Objects_Plots() ;    
	virtual void      CatIII_Objects_Plots() ;
	virtual void      CatIV_Objects_Plots() ;  
	virtual void      CatV_Objects_Plots() ;    
	virtual void      CatVI_Objects_Plots() ;
	virtual void      Category_Wjet_Plot( int cat, int W, int idx) ;
	virtual void      Category_Top_Plot( int cat,int top, int idx) ;
	virtual void      Category_Fat_Plot( int cat,int fat, int idx) ;
	virtual void      Category_Higgs_Plot( int cat, int higg, int idx) ;
	// 3. Transverse Mass Calculation   
	virtual  void    toptag_MTCalculation(int Cat ) ;
	virtual  void    Higgstag_MTCalculation()  ;
	virtual  void     genlvl_MTCalculation(int type ) ;
	virtual  void     WWtag_MTCalculation() ;
	virtual  void     WfatV_MTCalculation( ) ;
	virtual  void     WfatIII_MTCalculation( ) ;          
	virtual  void  	  Group5_WbH() ;
	virtual	 void  	  Group3_TopbH() ;
	//////////////
	virtual void      Define_Mt_Histo();  // for transverse mass histo only...delete it as soon for new histo

	// ===================create an array of Histograms======================
	// Category wise tag object Histogram

	std::array< std::array< TH1F*, 5> , 6> h_Histo_Pt;
	std::array< std::array< TH1F*, 5> , 6> h_Histo_Mt ;  
	std::array< std::array< TH1F*, 5> , 6> h_Histo_Eta ;                          
	std::array< std::array< TH1F*, 2> , 6> h_Histo_Mass ;  
	std::array< std::array< TH1F*, 2> , 6> h_Histo_SD ;  
	std::array< std::array< TH1F*, 2> , 6> h_Histo_tau ; 

	std::array< std::array< TH1F*, 6> , 6> h_Histo_dR ;   
	std::array< std::array< TH1F*, 6> , 6> h_Histo_dPt ;   

	// tag object histogram before Category selection
	std::array< TH1F*, 4> h_tag_N ;  
	std::array< TH1F*, 6> h_tag_Pt ;  
	std::array< TH1F*, 6> h_tag_Eta ;                          
	std::array< TH1F*, 6> h_tag_Mass ;  
	std::array< TH1F*, 6> h_tag_PUPPImass ;
	std::array< TH1F*, 6> h_tag_SD ;  
	std::array< TH1F*, 6> h_tag_Pruned;
	std::array< TH1F*, 6> h_tag_tau21 ;  
	std::array< TH1F*, 6> h_tag_tau32 ;  
	std::array< TH1F*, 6> h_tag_Puppitau21 ;  
	std::array< TH1F*, 6> h_tag_Puppitau32 ;  
	std::array< TH1F*, 6> h_tag_Puppitau42 ;  
	std::array< TH1F*, 6> h_tag_CHStau42 ;  
	std::array< std::array< TH1F*,  6> ,  6> h_dR_tagjet ;

	// dR Histogram for tagjet - reco object
	std::array< TH1F*, 6>  h_dR_Recomu_tagjet ;
	std::array< TH1F*, 6>  h_dPt_lep_tagjet ;
	std::array< TH1F*, 6>  h_dR_Recob1_tagjet ;
	std::array< TH1F*, 6>  h_dR_Recob2_tagjet ;

	// mc level object histogram
	std::array< TH1F*, 7>  h_mcobject_pt;
	std::array< TH1F*, 7>  h_mcobject_eta;
	std::array< TH1F*, 7>  h_mcobject_mass;		
	std::array< std::array< TH1F*, 7> , 7>  h_dR_MC ;
	std::array< std::array< TH1F*, 7> , 7>  h_dPt_MC ;
	std::array< std::array< TH1F*, 5> , 5>  h_dR_MCReco ;

	// reco object variables histogram
	std::array< TH1F*, 9>  h_object_pt ;
	std::array< TH1F*, 9>  h_object_eta ;
	std::array< TH1F*, 5>  h_object_no ;
	std::array< TH1F*, 3>  h_AK8_Jetmass ;
	std::array< TH1F*, 3>  h_AK8_PUPPImass ;
	std::array< TH1F*, 3>  h_AK8_PUPPISDmass;
	std::array< TH2F*, 3>  h_AK8_PUPPIvsmass ;
	std::array< TH2F*, 3>  h_AK8_CHSvsmass ;
	std::array< TH1F*, 3>  h_AK8_CHSmass ;
	std::array< TH1F*, 3>  h_AK8_PUPPItau21 ; 
	std::array< TH1F*, 3>  h_AK8_PUPPItau31 ; 
	std::array< TH1F*, 3>  h_AK8_PUPPItau32 ; 
	std::array< TH1F*, 3>  h_AK8_PUPPItau41 ; 
	std::array< TH1F*, 3>  h_AK8_PUPPItau42 ; 
	std::array< TH1F*, 3>  h_AK8_PUPPItau43 ; 
	std::array< TH1F*, 3>  h_AK8_tau21 ; 
	std::array< TH1F*, 3>  h_AK8_tau32 ; 
	std::array< TH1F*, 3>  h_AK8_CHStau42 ; 
	std::array< std::array< TH1F*, 8> , 8>  h_dR_Reco ;
	std::array< std::array< TH1F*, 8> , 8>  h_dPt_Reco ;
	std::array< TH1F*, 3>  h_MET_var ;
	std::array< TH1F*, 9>  h_object_MT ;
	TH2F* Ptbjet_dRW;
	TH2F* Ptbjet_dRmu ;
	std::array< TH2F*, 4>  h_AK81vsGenObj ;
	std::array< TH2F*, 4>  h_AK82vsGenObj ;
	std::array< TH2F*, 4>  h_AK83vsGenObj ;

	//======Global Variables =================

	float tauval[6] = {0.65 , 0.65, 0.65, 0.65, 0.54, 0.54} ;
	float tau_21[6] = {0.55 , 0.55, 0.40, 0.40, 0.35, 0.35} ;
	float PUPPItau[4] = { 0.70 , 0.60 ,0.50, 0.40 } ;
	int vr = -1 ;
	int  T_top = -1;
	int  T_higgs = -1;
	int  Higgs_W = -1 ;
	int  Top_W = -1 ;
	float METCut_cat12[5] = { 100.0,150.0,200.0,250.0,300.0} ;
	float METCut_cat3[5] = { 100.0,150.0,200.0,0.0,0.0} ;

	int event_Higgs 		= -1 ;
	int event_top 			= -1 ;	
	int event_bjet 			= -1 ;	
	int event_dR_Hmu  		= -1 ;
	int event_LepIso_bmu 		= -1 ;
	int event_dR_Wt_tmu 		= -1 ;
	int event_LepIso_Wmu_cat2 	= -1 ;
	int event_bjet_cat3 		= -1 ;
	int event_2W 			= -1 ;
        int event_0W 			= -1 ;
        int event_1W 			= -1 ;
	int event_dR_Wb_bmu 		= -1 ;
	int event_LepIso_Wmu_cat3 		= -1 ;
	int event_topPt			= 0 ;	
	int event_HiggsPt 		= 0 ;
	int event_W_Pt 			= 0 ;

	int  b_top = -1 ;
	int  b_asso = -1 ;
	int  q_forw = -1 ;
	int  Higgs_indx = - 1 ;
	int  top_indx = -1 ;


	vector<int>   puppi_jet ;
	vector <int> n_ele;
	vector <int> b_jet_loose;
	vector <int> b_jet;
	vector <int> n_Mu;
	vector <int> top_Mu;
	vector <int> Higgs_Mu;                
	vector <int> b_jet_tight;
	vector <int> n_jet;
	vector <int> n_forwjet;
	vector <int> n_AK8Jet;
	vector <int> b_jet_medium;
	vector <int> fat_jet;
	vector <int> Higgs_Jet;
	vector <int> Higgsjets;                
	vector <int> Higgs_GenJet;
	vector <int> Top_GenJet;
	vector <int> topW_q;
	vector <int> W_boson;
	vector <int> bjet_match ; 
	vector <int> topjet ;

	vector <int> CatI_Objects;              // for top, W, muon & MET
	vector <int> CatII_Objects;             // for top, fat, muon & MET
	vector <int> CatIII_Objects;            // for W, b, fat, muon & MET
	vector <int> CatIV_Objects;            // for  W, b, W, muon & MET
	vector <int> CatV_Objects;             // for  muon, MET, b, W & fat
	vector <int> CatVI_Objects;            // for Higgs, b, l & MET 

};

#endif

#ifdef TTbar_PUPPItau_loc_190204_cxx
TTbar_PUPPItau_loc_190204::TTbar_PUPPItau_loc_190204(TString inputFile)
{
	TChain *tree = new TChain("v_tree");
	//        add input files here


	ifstream datafile;
	datafile.open(inputFile.Data(), ifstream::in );

	TString datafilename;

	//            for Qjet sample
//	string location ="/home/logan/VLQAnalysis/New_Method/WjetsQQ_files/";
	 string location ="root://cmseos.fnal.gov//store/user/achhetri/";

	while (true) {
		datafile >> datafilename;
		string fname = datafilename.Data();
		string nameFile = fname;
		string string_search (".root");
		size_t found =nameFile.find(string_search);

		if(found != string::npos){
			tree->Add((location+fname).c_str());
			cout<<" FileName =" <<fname<<" , Events= "<<tree->GetEntries()<<"  file= "<<((location+fname).c_str())<<endl;
		}
		if( datafile.eof() ) break;
	}
	Init(tree);
}

TTbar_PUPPItau_loc_190204::~TTbar_PUPPItau_loc_190204()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//====================User Defined Functions' Definitions====================

void  TTbar_PUPPItau_loc_190204::Clear_Vector()
{
	//   for clearing objects vector
	puppi_jet.clear();
	n_ele.clear();
	n_Mu.clear();
	top_Mu.clear();
	topjet.clear();        
	Higgs_Mu.clear();        
	n_jet.clear();
	n_forwjet.clear();
	b_jet_loose.clear();
	b_jet.clear();
	b_jet_medium.clear();
	b_jet_tight.clear();
	fat_jet.clear();
	Higgs_Jet.clear();
	Higgsjets.clear();        
	Higgs_GenJet.clear();
	Top_GenJet.clear();
	n_AK8Jet.clear();
	topW_q.clear();
	W_boson.clear();
	bjet_match.clear();
	CatI_Objects.clear() ;    
	CatII_Objects.clear() ;
	CatIII_Objects.clear() ;
	CatIV_Objects.clear() ;
	CatV_Objects.clear() ;
	CatVI_Objects.clear() ;        
	//     CategoryII.clear();             
}


//========== Histogram Definations======================================


void TTbar_PUPPItau_loc_190204::Define_Mt_Histo()
{
	//   New defination of histogram
	char  dr_MTname[100] , dr_MTtitle[100] ;

	int s = -1;
	int rk = -1 ;
	string MT_comp[9] = {"Mt_muMeT", "Mt_WmuMeT", "Mt_Tprime", "Mt_topjet", "Mt_Wjet", "Mt_Wmc", "MT_HiggsMC", "MT_topmc", "MT_Tmc"} ;
	//cout<<"\nFunc 1" ;
	// for Transverse mass plots, we have 6 combinations for objects
	for( int Sk = 0 ; Sk < 9 ; Sk ++ )
	{
		sprintf( dr_MTname, "%s",  MT_comp[Sk].c_str() ) ;
		sprintf( dr_MTtitle, "%s Distribution",  MT_comp[Sk].c_str() ) ;		
		h_object_MT.at(Sk) = new TH1F(dr_MTname,dr_MTtitle,  600, 0.0 , 3000.0 );
		h_object_MT.at(Sk)->GetYaxis()->SetTitle("Events");
		//cout<< " , h_object_Mt["<<Sk<<"].name() : "<< h_object_MT.at(Sk) ->GetName()<<endl; 
	}      
}

// ---------------------------------------------------mc level objects ----------------------------------------------

void TTbar_PUPPItau_loc_190204::DefineMC_NPtEta_Histo()
{
	//   New defination of histogram
	//        TString   mc_pT , mc_pT_title, gh_name, gh_title;
	char mc_pT[100], mc_pT_title[100], gh_name[100], gh_title[100] ;

	std::string mc_objectI[7] = {"topb", "topWl","higgWq1", "higgWq2", "higgWl", "assob", "forwq" };
	//std::string mc_objectI[3] = {"topb", "topWl","higgWq1"};
	// for object pT histograms
	//cout << "\nFunc 2" ;
	for(Int_t k = 0 ; k < 7; ++k) {
		sprintf (mc_pT, "Pt(%s)", mc_objectI[k].c_str() )  ;
		sprintf (mc_pT_title, "Pt(%s) Distribution", mc_objectI[k].c_str() )  ;               
		h_mcobject_pt.at(k) = new TH1F(mc_pT,mc_pT_title, 200, 0, 1000.0);
		h_mcobject_pt.at(k) ->GetYaxis()->SetTitle("Events");        
		//	cout<< " , h_mcobject_pt["<<k<<"].name() : "<< h_mcobject_pt.at(k) ->GetName()<<endl; 	
	}


	for(Int_t h = 0 ; h < 7; ++h) {
		sprintf (gh_name, "Eta(%s)", mc_objectI[h].c_str() )  ;
		sprintf (gh_title, "Eta(%s) Distribution", mc_objectI[h].c_str() )  ;                            
		h_mcobject_eta.at(h) = new  TH1F(gh_name,gh_title, 200, -5.0, 5.0);
		h_mcobject_eta.at(h) ->GetYaxis()->SetTitle("Events");
		//cout<<" , h_mcobject_eta["<<h<<"].name() : "<< h_mcobject_eta.at(h) ->GetName()<<endl; 
	}

	for(Int_t h = 0 ; h < 7; ++h) {
		sprintf (gh_name, "Mass(%s)", mc_objectI[h].c_str() )  ;
		sprintf (gh_title, "Mass(%s) Distribution", mc_objectI[h].c_str() )  ;                            
		h_mcobject_mass.at(h) = new  TH1F(gh_name,gh_title, 600, 0.0, 3000.0);
		h_mcobject_mass.at(h) ->GetYaxis()->SetTitle("Events");
		//cout<<" , h_mcobject_eta["<<h<<"].name() : "<< h_mcobject_eta.at(h) ->GetName()<<endl; 
	}
}



void TTbar_PUPPItau_loc_190204::dRHisto_MCObject()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string mc_object[7] = {"topb", "topWl","higgWq1", "higgWq2", "higgWl", "assob", "forwq" };
	// Here, i = 0-1 top decayed, i = 2-4 Higgs decayed, i = 5 assob, i = 6 forwq
	Int_t s = -1;
	//cout << "\nFunc 3" ;	
	for (Int_t i = 0; i < 7; i++) {     
		for (Int_t k = i+1; k < 7; k++) {
			// for mc  objects dR histograms
			sprintf(dR_name, "DeltaR(%s_%s)", mc_object[i].c_str(), mc_object[k].c_str()) ;
			sprintf(dR_title, "DeltaR(%s_%s) Distribution", mc_object[i].c_str(), mc_object[k].c_str()) ;
			h_dR_MC.at(i).at(k) = new TH1F(dR_name,dR_title, 10000, 0, 10.0);
			h_dR_MC.at(i).at(k) ->GetYaxis()->SetTitle("Events");			
			//cout<<" , h_dR_MC["<<i<<"]["<<k<<"].name() : "<< h_dR_MC.at(i).at(k) ->GetName()<<endl; 
			// for mc Object dPt histograms     
			sprintf(dPt_name, "DeltaPt(%s_%s)", mc_object[i].c_str(), mc_object[k].c_str()) ;
			sprintf(dPt_title, "DeltaPt(%s_%s) Distribution", mc_object[i].c_str(), mc_object[k].c_str()) ;
			h_dPt_MC.at(i).at(k) = new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
			h_dPt_MC.at(i).at(k) ->GetYaxis()->SetTitle("Events");			
			//cout<<" h_dPt_MC["<<i<<"]["<<k<<"].name() : "<< h_dPt_MC.at(i).at(k) ->GetName()<<endl;			
		}
	}
}     


//------------------------------------------reco level objects - Preselection-----------------------------------------------

void TTbar_PUPPItau_loc_190204::Define_2DMass_Histo()
{
	//   New defination of histogram
	char   mass[100] , massT[100] ; 
	string object_no[4] = {"Higgs", "WHiggs","Top","WTop"};
	string reco_object[3] = {"AK8jet1", "AK8jet2","AK8jet3"};

	for(Int_t k = 0 ; k <  3; k++) {
		for(Int_t j = 0; j < 4; j++) {

			sprintf( mass , "%svs%s_Mass" , reco_object[k].c_str(),  object_no[j].c_str()  ) ;
			sprintf( massT , "%svs%s_Mass Distribution" , reco_object[k].c_str(),  object_no[j].c_str()) ;			

			if ( k == 0) {
				h_AK81vsGenObj.at(j) = new TH2F(mass,massT, 200, 0.0, 1000.0, 200, 0.0, 1000.0);
				h_AK81vsGenObj.at(j) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_pt["<<k<<"].name() : "<< h_object_pt.at(k) ->GetName()<<endl; 
			}

			if ( k == 1) {
				h_AK82vsGenObj.at(j) = new TH2F(mass,massT, 200, 0.0, 1000.0, 200, 0.0, 1000.0);
				h_AK82vsGenObj.at(j) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_pt["<<k<<"].name() : "<< h_object_pt.at(k) ->GetName()<<endl; 
			}

			if ( k == 2) {
				h_AK83vsGenObj.at(j) = new TH2F(mass,massT, 200, 0.0, 1000.0, 200, 0.0, 1000.0);
				h_AK83vsGenObj.at(j) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_pt["<<k<<"].name() : "<< h_object_pt.at(k) ->GetName()<<endl; 
			}

		}
	}
}


void TTbar_PUPPItau_loc_190204::Define_NPtEta_Histo()
{
	//   New defination of histogram
	char   pT_name[100] , pT_title[100] ; 
	string object_no[5] = {"mu", "bjet", "AK8jet","qjet","forwjet"};
	string reco_object[9] = {"MET", "mu1", "bjet1", "bjet2","fjet1", "fjet2","AK8jet1", "AK8jet2","AK8jet3"};
	string variable[17] = {"pt", "eta", "puppimass", "chsprunedmass", "puppitau21", "puppitau32","tau21", "tau32", "chstau42","puppisdvsmass", "chsprunedvsmass","jetmass","puppisdmass","puppitau42","puppitau41", "puppitau43","puppitau31"} ;
	string ar;
	//cout << "\nFunc 4" ;
	int s = -1;

	// for object pT histograms
	for(Int_t k = 0 ; k <  9; k++) {
		for(Int_t j = 0; j < 17; j++) {

			sprintf( pT_name , "%s(%s)" , variable[j].c_str(), reco_object[k].c_str() ) ;
			sprintf( pT_title , "%s(%s) Distribution" , variable[j].c_str(), reco_object[k].c_str() ) ;			

			if ( j == 0) {
				h_object_pt.at(k) = new TH1F(pT_name,pT_title, 200, 0, 1000.0);
				h_object_pt.at(k) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_pt["<<k<<"].name() : "<< h_object_pt.at(k) ->GetName()<<endl; 
			}
			if ( j == 1) {
				h_object_eta.at(k) = new TH1F(pT_name,pT_title, 200, -5.0, 5.0);
				h_object_eta.at(k) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_eta["<<k<<"].name() : "<< h_object_eta.at(k) ->GetName()<<endl; 	
			}
			if(!( k > 5) ) continue ;
			if( j == 2) {
				s = k +j - 8 ;
				h_AK8_PUPPImass.at(s) =  new TH1F(pT_name,pT_title, 200, 0.0,  1000.0);
				h_AK8_PUPPImass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPImass["<<s<<"].name() : "<< h_AK8_PUPPImass.at(s)->GetName()<<endl; 	
			}
			if( j == 3) {
				s = k +j - 9 ;
				h_AK8_CHSmass.at(s) =  new TH1F(pT_name,pT_title, 102, -10.0, 500.0);
				h_AK8_CHSmass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_CHSmass["<<s<<"].name() : "<< h_AK8_CHSmass.at(s)->GetName()<<endl; 	
			}

			if( j == 4) {
				s = k +j - 10;
				h_AK8_PUPPItau21.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau21.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau21["<<s<<"].name() : "<< h_AK8_PUPPItau21.at(s)->GetName()<<endl; 					
			}
			if( j == 5) {
				s = k +j - 11 ;
				h_AK8_PUPPItau32.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau32.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau32["<<s<<"].name() : "<< h_AK8_PUPPItau32.at(s)->GetName()<<endl; 									
			}
			if( j == 6) {
				s = k +j - 12 ;
				h_AK8_tau21.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_tau21.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_tau21["<<s<<"].name() : "<< h_AK8_tau21.at(s)->GetName()<<endl; 					
			}
			if( j == 7) {
				s = k +j - 13 ;
				h_AK8_tau32.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_tau32.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_tau32["<<s<<"].name() : "<< h_AK8_tau32.at(s)->GetName()<<endl; 									
			}

			if( j == 8) {
				s = k +j - 14 ;
				h_AK8_CHStau42.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_CHStau42.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_CHStau42["<<s<<"].name() : "<< h_AK8_CHStau42.at(s)->GetName()<<endl; 									
			}
			if( j == 9) {
				s = k +j - 15 ;
				h_AK8_PUPPIvsmass.at(s) =  new TH2F(pT_name,pT_title, 202, -10.0, 1000.0, 202, -10.0, 1000.0) ;
				h_AK8_PUPPIvsmass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_Puppivsmass["<<s<<"].name() : "<< h_AK8_PUPPIvsmass.at(s)->GetName()<<endl; 									
			}
			if( j == 10) {
				s = k +j - 16 ;
				h_AK8_CHSvsmass.at(s) =  new TH2F(pT_name,pT_title, 202, -10.0, 1000.0, 202, -10.0, 1000.0) ;
				h_AK8_CHSvsmass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_CHSvsmass["<<s<<"].name() : "<< h_AK8_CHSvsmass.at(s)->GetName()<<endl; 									
			}

			if( j == 11) {
				s = k +j - 17 ;
				h_AK8_Jetmass.at(s) =  new TH1F(pT_name,pT_title, 200, 0.0, 1000.0);
				h_AK8_Jetmass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPImass["<<s<<"].name() : "<< h_AK8_Jetmass.at(s)->GetName()<<endl; 	
			}

			if( j == 12) {
				s = k +j - 18 ;
				h_AK8_PUPPISDmass.at(s) =  new TH1F(pT_name,pT_title, 102,  -10.0, 500.0) ;
				h_AK8_PUPPISDmass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPImass["<<s<<"].name() : "<< h_AK8_PUPPImass.at(s)->GetName()<<endl; 	
			}
			if( j == 13) {
				s = k +j - 19;
				h_AK8_PUPPItau42.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau42.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau42["<<s<<"].name() : "<< h_AK8_PUPPItau42.at(s)->GetName()<<endl; 					
			}
	if( j == 14) {
				s = k +j - 20;
				h_AK8_PUPPItau41.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau41.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau42["<<s<<"].name() : "<< h_AK8_PUPPItau42.at(s)->GetName()<<endl; 					
			}
	if( j == 15) {
				s = k +j - 21;
				h_AK8_PUPPItau43.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau43.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau42["<<s<<"].name() : "<< h_AK8_PUPPItau42.at(s)->GetName()<<endl; 					
			}
	if( j == 16) {
				s = k +j - 22;
				h_AK8_PUPPItau31.at(s) =  new TH1F(pT_name,pT_title, 200, -1.0, 1.0) ;
				h_AK8_PUPPItau31.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau42["<<s<<"].name() : "<< h_AK8_PUPPItau42.at(s)->GetName()<<endl; 					
			}

		}
	}

	// for No. object  histograms
	for(Int_t j = 0; j < 5; j++) {	
		sprintf( pT_name , "N(%s)" , object_no[j].c_str() ) ;
		sprintf( pT_title , "N(%s) Distribution" , object_no[j].c_str() ) ;	
		h_object_no.at(j) = new TH1F(pT_name,pT_title, 10, 0, 10);  
		h_object_no.at(j) ->GetYaxis()->SetTitle("Events");
		//cout<<"h_object_no["<<j<<"].name() : "<< h_object_no.at(j)->GetName()<<endl; 			
	}

}


void TTbar_PUPPItau_loc_190204::dRHisto_RecoObject() 
{
	// New defination of histogram for reco level at Preselection 
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string reco_object[7] = {"mu1", "bjet1", "bjet2","fjet1", "fjet2","AK8jet1", "AK8jet2"};
	// Here i -> 0 muon, i -> 1-2 bjet & i -> 3-4 forwjet, i -> 5-6 AK8jet
	//cout<<"\nFunc 5" ;	
	for (Int_t i = 0; i < 7; i++) {          
		for (Int_t k = i+1; k < 7; k++) {      

			// for reco preselected objects dR histograms
			sprintf( dR_name , "DeltaR(%s_%s)" , reco_object[i].c_str(), reco_object[k].c_str() ) ;
			sprintf( dR_title ,"DeltaR(%s_%s) Distribution",reco_object[i].c_str(),reco_object[k].c_str()) ;
			h_dR_Reco.at(i).at(k) = new TH1F(dR_name,dR_title, 10000, 0, 10.0);
			h_dR_Reco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dR_Reco[" <<i<<"]["<<k<<"].name() : "<< h_dR_Reco.at(i).at(k)->GetName()<<endl; 									


			// for lep-bjet dPt histograms  
			sprintf( dPt_name , "DeltaPt(%s_%s)" , reco_object[i].c_str(), reco_object[k].c_str() ) ;
			sprintf( dPt_title ,"DeltaPt(%s_%s) Distribution",reco_object[i].c_str(),reco_object[k].c_str());
			h_dPt_Reco.at(i).at(k) = new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
			h_dPt_Reco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dPt_Reco[" <<i<<"]["<<k<<"].name() : "<< h_dPt_Reco.at(i).at(k)->GetName()<<endl; 			
		}
	}
}     

//------------------------------------------mc-Reco dR Histogram--------------------------------------

void TTbar_PUPPItau_loc_190204::dRHisto_MCRecoObject()
{
	// New defination of histogram
	char dR_name[100], dR_title[100];
	string mc_object[4] = {"Higgs","W_Higgs", "Top", "W_Top"};
	string reco_object[3] = {"AK8jet1", "AK8jet2","AK8jet3"};
	//cout<<"\nFunc 6" ;
	for (Int_t i = 0; i < 3; i++) {          
		for (Int_t k = 0; k < 4; k++) {
			// for mc  objects dR histograms
			sprintf( dR_name, "DeltaR(%s_%s)" ,reco_object[i].c_str(), mc_object[k].c_str()) ;
			sprintf( dR_title,"DeltaR(%s_%s) Distribution" , reco_object[i].c_str(), mc_object[k].c_str()) ;
			h_dR_MCReco.at(i).at(k) = new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_MCReco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dR_MCReco[" <<i<<"]["<<k<<"].name() : "<< h_dR_MCReco.at(i).at(k)->GetName()<<endl; 			
		}
	}     
}
//-----------------------------Tagged jets Histogram ---------------------------------------------------

void  TTbar_PUPPItau_loc_190204::Define_Tag_Jet_Histo()
{
	char tag_name[100] , tag_title[100] ;
	int x1 = -1 ;
	int rk = -1 ;

	string jet_tag[4] = { "Wjet", "topjet", "Higgsjet", "fatjet" } ; 
	string  Wjet_var[12]  = { "Pt", "Eta", "Mass", "PuppiSDMass","CHSPrunedMass", "tau21","tau32","Puppitau21", "PUPPItau32", "CHStau42","PUPPIMass", "PUPPItau42"} ;
	string   top_var[12]  = { "Pt", "Eta", "Mass", "PuppiSDMass","CHSPrunedMass", "tau21","tau32","Puppitau21", "PUPPItau32", "CHStau42","PUPPIMass", "PUPPItau42"} ;
	string Higgs_var[12]  = {"Pt", "Eta", "Mass", "PuppiSDMass","CHSPrunedMass", "tau21","tau32","Puppitau21", "PUPPItau32", "CHStau42","PUPPIMass", "PUPPItau42"} ;
	//cout<<"\nFunc 7" ;
	for( int h =0 ; h < 4 ; h ++) {       
		sprintf( tag_name, "N(%s)" , jet_tag[h].c_str() ) ;
		sprintf( tag_title, "N(%s) Distribution" , jet_tag[h].c_str() ) ;
		h_tag_N.at(h)  =  new TH1F(tag_name,tag_title, 10, 0, 10);
		h_tag_N.at(h)  -> GetYaxis()->SetTitle("Events");
		//cout<<"h_tag_N[" <<h<<"].name() : "<< h_tag_N.at(h)->GetName()<<endl; 		
	}

	// Here i = 0-1 Wjets, i -> 2 topjet, i = 3 Higgsjet, i = 4-5  fatjet,Use this index cycle for filling plots  
	for(     int Sj = 0 ; Sj <  6; Sj ++) {
		for( int Tk = 0 ; Tk < 12; Tk ++ ) {

			if ( Sj < 2  ) {
				x1 = Sj + 1 ;                // for Wjets index   0-1     
				sprintf( tag_name, "%s(%s%d)",  Wjet_var[Tk].c_str(), jet_tag[0].c_str(), x1 ) ;  
				sprintf( tag_title,"%s(%s%d) Distribution", Wjet_var[Tk].c_str(), jet_tag[0].c_str(), x1);  			
			}
			if ( Sj == 2 ) {
				x1 = Sj - 1 ;           // for topjet index  2
				sprintf( tag_name, "%s(%s%d)",  top_var[Tk].c_str(), jet_tag[1].c_str(), x1 ) ;  
				sprintf( tag_title,"%s(%s%d) Distribution", top_var[Tk].c_str(), jet_tag[1].c_str(), x1); 
			}
			if ( Sj == 3 ) {
				x1 = Sj - 2 ;           // for Higgsjet index   3 
				sprintf( tag_name, "%s(%s%d)",  Higgs_var[Tk].c_str(), jet_tag[2].c_str(), x1 ) ;     	
				sprintf( tag_title,"%s(%s%d) Distribution", Higgs_var[Tk].c_str(), jet_tag[2].c_str(), x1); 	
			}          
			if ( Sj > 3)      {
				x1 = Sj - 3 ;            // for fatjet index 4-5
				sprintf( tag_name, "%s(%s%d)",  Wjet_var[Tk].c_str(), jet_tag[3].c_str(), x1 ) ;  		
				sprintf( tag_title,"%s(%s%d) Distribution", Wjet_var[Tk].c_str(), jet_tag[3].c_str(), x1); 	
			}         

			if(Tk == 0 ) {
				h_tag_Pt.at(Sj)             = new TH1F(tag_name,tag_title, 200, 0, 1000.0);
				h_tag_Pt.at(Sj)             ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Pt[" <<Sj<<"].name() : "<< h_tag_Pt.at(Sj)->GetName()<<endl; 						
			}      
			if(Tk == 1 )  {
				h_tag_Eta.at(Sj)           = new TH1F(tag_name,tag_title, 200, -5.0, 5.0);                   
				h_tag_Eta.at(Sj)            ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Eta[" <<Sj<<"].name() : "<< h_tag_Eta.at(Sj)->GetName()<<endl; 			
			}      
			if(Tk == 2 ) {
				h_tag_Mass.at(Sj)         = new TH1F(tag_name,tag_title, 200, 0, 1000.0);         
				h_tag_Mass.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Mass[" <<Sj<<"].name() : "<< h_tag_Mass.at(Sj)->GetName()<<endl; 						
			}                                
			if(Tk == 3  )  {
				h_tag_SD.at(Sj)      = new TH1F(tag_name,tag_title, 102, -10.0, 500.0);       
				h_tag_SD.at(Sj)       ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_SD[" <<Sj<<"].name() : "<< h_tag_SD.at(Sj)->GetName()<<endl;				
			}            
			if( Tk == 4 )  {
				h_tag_Pruned.at(Sj)      = new TH1F(tag_name,tag_title, 102, -10.0, 500.0);       
				h_tag_Pruned.at(Sj)       ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_SD[" <<Sj<<"].name() : "<< h_tag_SD.at(Sj)->GetName()<<endl;				
			}            

			if(Tk == 5 )  {
				h_tag_tau21.at(Sj)         = new TH1F(tag_name,tag_title, 200, -1.0, 1.0);               
				h_tag_tau21.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_tau21.at(Sj)->GetName()<<endl;				
			}                        
			if(Tk == 6 )  {
				h_tag_tau32.at(Sj)         = new TH1F(tag_name,tag_title,200, -1.0, 1.0);               
				h_tag_tau32.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_tau32.at(Sj)->GetName()<<endl;				
			}   
			if(Tk == 7 )  {
				h_tag_Puppitau21.at(Sj)         = new TH1F(tag_name,tag_title, 200, -1.0, 1.0);               
				h_tag_Puppitau21.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_Puppitau21.at(Sj)->GetName()<<endl;				
			}   
			if(Tk == 8 )  {
				h_tag_Puppitau32.at(Sj)         = new TH1F(tag_name,tag_title, 200, -1.0, 1.0);               
				h_tag_Puppitau32.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_Puppitau32.at(Sj)->GetName()<<endl;				
			}   
			if(Tk == 9 )  {
				h_tag_CHStau42.at(Sj)         = new TH1F(tag_name,tag_title, 200, -1.0, 1.0);               
				h_tag_CHStau42.at(Sj)         ->GetYaxis()->SetTitle("Events");
			} 
			if(Tk == 10){ 
				h_tag_PUPPImass.at(Sj)         = new TH1F(tag_name,tag_title, 200, 0, 1000.0);       
				h_tag_PUPPImass.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_PUPPImass[" <<Sj<<"].name() : "<< h_tag_PUPPImass.at(Sj)->GetName()<<endl;				
			} 
			if(Tk == 11 )  {
				h_tag_Puppitau42.at(Sj)         = new TH1F(tag_name,tag_title, 200, -1.0, 1.0);               
				h_tag_Puppitau42.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_Puppitau32.at(Sj)->GetName()<<endl;				
			}   


		}      
	}
}


void TTbar_PUPPItau_loc_190204::dR_tagjetHisto()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string tag_jet[6] = {"Wjet1", "Wjet2", "topjet1",  "Higgsjet1", "fatjet1", "fatjet2"} ;                   
	//for tag jet   dR histograms
	// Here i = 0-1 Wjets, i -> 2 topjet, i = 3 Higgsjet, i = 4-5 fatjet, Use this index cycle for illing plots
	//cout<<"\nFunc 8" ;
	for (int i = 0; i < 6; i++) {
		for(int j = i+1; j < 6; j++) {

			sprintf(dR_name,"DeltaR(%s_%s)" , tag_jet[i].c_str(), tag_jet[j].c_str() ) ;
			sprintf(dR_title,"DeltaR(%s_%s) Distribution", tag_jet[i].c_str(), tag_jet[j].c_str() ) ;	
			h_dR_tagjet.at(i).at(j)   =    new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_tagjet.at(i).at(j)   ->   GetYaxis()->SetTitle("Events");         
			//cout<<"h_dR_tagjet["<<i<<"]["<<j<<"].name() : "<< h_dR_tagjet.at(i).at(j)->GetName()<<endl;			
		}
	}
}


void TTbar_PUPPItau_loc_190204::Define_Reco_tagjetHisto()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string reco_object[3]  =  {"muon","bjet1", "bjet2"} ;
	string tag_jet[6] = {"Wjet1", "Wjet2", "topjet1", "Higgsjet1", "fatjet1", "fatjet2" } ;  

	//for tag jet - muon  dR histograms
	// Here i = 0-1 Wjets, i -> 2 topjet, i = 3 Higgsjet, i = 4-5 fatjet, Use this index cycle for illing plots
	//cout<<"\nFunc 9" ;	
	int s = -1 ;
	for (Int_t i = 0; i < 3; i++) {
		for(Int_t j = 0; j < 6; j++) {
			//cout << "\n value i = " << i << " && value j << " << j ;        

			if( i == 0){			
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str()) ;
				h_dR_Recomu_tagjet.at(j)= new TH1F(dR_name,dR_title, 5000, 0, 5.0);
				h_dR_Recomu_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recomu_tagjet]"<<j<<"].name() : "<< h_dR_Recomu_tagjet.at(j)->GetName()<<endl;					

				sprintf( dPt_name,"DeltaPt(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dPt_title,"DeltaPt(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str());
				h_dPt_lep_tagjet.at(j)= new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
				h_dPt_lep_tagjet.at(j)->GetYaxis()->SetTitle("Events");				
				//cout<<"h_dPt_lep_tagjet]"<<j<<"].name() : "<< h_dPt_lep_tagjet.at(j)->GetName()<<endl;						
			}
			if( i == 1){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str());	
				h_dR_Recob1_tagjet.at(j)= new TH1F(dR_name,dR_title, 5000, 0, 5.0);
				h_dR_Recob1_tagjet.at(j)->GetYaxis()->SetTitle("Events");      
				//cout<<"h_dR_Recob1_tagjet]"<<j<<"].name() : "<< h_dR_Recob1_tagjet.at(j)->GetName()<<endl;										
			}
			if( i == 2){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str() );
				h_dR_Recob2_tagjet.at(j)= new TH1F(dR_name,dR_title, 5000, 0, 5.0);
				h_dR_Recob2_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recob2_tagjet]"<<j<<"].name() : "<< h_dR_Recob2_tagjet.at(j)->GetName()<<endl;				
			}
		}//END OF J
	}// END OF I
	//	cout << "Size in loop = " <<  sizeof(h_dR_Recomu_tagjet) <<endl;
} // END OF FUNC


//================Signal Category Objects Histogram==================


void TTbar_PUPPItau_loc_190204::Category_Object_Histo()
{
	char Histo_name[100], Histo_title[100] ;
	string  obj[5], var  ;
	string  reco_var[2] = { "Pt", "Eta" } ;
	string Wjet_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau21"} ;
	string top_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau32"} ;
	string Higgs_var[5] = { "Pt", "Eta", "TransMass", "CHSPrunedMass",  "CHStau42"} ;
	int x ;
	int f = 0 ;
	//cout<<"\nFunc 10" ;	
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Topjet" ; obj[1] = "HWjet"; obj[2] = "muon" ; obj[3] = "MET";
			f = 4;
		}
		if ( u == 1 )  {
			obj[0] = "Topjet" ; obj[1] = "HFatjet"; obj[2] = "muon" ; obj[3] = "MET";
			f =4 ;
		}
		if ( u == 2 )  {
			obj[0] = "TWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 3 )  {
			obj[0] = "TWjet" ; obj[1] = "HWjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 4 )  {
			obj[0] = "HWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgsjet" ; obj[2] = "bjet"; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}

		for ( int j = 0; j < f; j ++ )
		{
			if ( obj[j] == "no" ) continue;
			for ( int k = 0; k < 5 ; k++) 
			{
				if ( j > 1 && k > 1) continue ;

				if( obj[j] == "Topjet" )                                var =  top_var[k]  ;
				if(obj[j] == "Higgsjet")                                var =  Higgs_var[k]  ;
				if( obj[j] == "muon" )                                 var =  reco_var[k]  ;
				if( obj[j] == "MET" )                                 var =  reco_var[k]  ;
				if( obj[j] == "bjet" )                                   var =  reco_var[k] ;
				if( obj[j] == "TWjet" || obj[j] == "HWjet" )		var  =  Wjet_var[k]  ;				
				if( obj[j] == "TFatjet" || obj[j] == "HFatjet")     var  =  Wjet_var[k]  ;				

				sprintf( Histo_name, "Category%d_%s(%s)", x, obj[j].c_str(), var.c_str() ) ; 
				sprintf( Histo_title, "Category%d_%s(%s) Distribution", x, obj[j].c_str(), var.c_str() ) ; 				
				if(k == 0 ) {
					h_Histo_Pt.at(u).at(j)       = new TH1F(Histo_name,Histo_title, 200, 0, 1000.0);
					h_Histo_Pt.at(u).at(j)             ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Pt["<<u<<"]["<<j<<"].name() : "<< h_Histo_Pt.at(u).at(j)->GetName()<<endl;     				
				}      
				if(k == 1 )  {
					h_Histo_Eta.at(u).at(j)      = new TH1F(Histo_name,Histo_title, 200, -5.0, 5.0);    
					h_Histo_Eta.at(u).at(j)      ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Eta["<<u<<"]["<<j<<"].name() : "<< h_Histo_Eta.at(u).at(j)->GetName()<<endl;					
				}      
				if(k == 2 ) {
					h_Histo_Mass.at(u).at(j)     = new TH1F(Histo_name,Histo_title, 100, 0, 1000.0);   
					h_Histo_Mass.at(u).at(j)      ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Mass["<<u<<"]["<<j<<"].name() : "<< h_Histo_Mass.at(u).at(j)->GetName()<<endl;						
				}                                
				if(k == 3 )  {
					h_Histo_SD.at(u).at(j)        = new TH1F(Histo_name,Histo_title, 100, 0, 1000.0);
					h_Histo_SD.at(u).at(j)        ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_SD["<<u<<"]["<<j<<"].name() : "<< h_Histo_SD.at(u).at(j)->GetName()<<endl;						
				}            
				if(k == 4 )  {
					h_Histo_tau.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 200, 0, 1.0);  
					h_Histo_tau.at(u).at(j)         ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_tau["<<u<<"]["<<j<<"].name() : "<< h_Histo_tau.at(u).at(j)->GetName()<<endl;
				}  // end of k ==4                       

			}  //end of kth loop
		}    //end of jth loop

	}    // end of uth loop
}  //  end of function


void TTbar_PUPPItau_loc_190204::Category_Object_MtHisto()
{
	char Histo_name[100], Histo_title[100];
	string TransM, transM1 , transM2 ;
	string  obj[5]  ;
	int x ;
	int f = 0 ;
	//cout<<"\nFunc 11" ;		
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Top" ; obj[1] = "HW"; obj[2] = "mu" ; obj[3] = "MET";
			f = 4;
		}
		if ( u == 1 )  {
			obj[0] = "Top" ; obj[1] = "HFat"; obj[2] = "mu" ; obj[3] = "MET";
			f =4 ;
		}
		if ( u == 2 )  {
			obj[0] = "TW" ; obj[1] = "b"; obj[2] = "HFat" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 3 )  {
			obj[0] = "TW" ; obj[1] = "b"; obj[2] = "HW" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 4 )  {
			obj[0] = "HW" ; obj[1] = "HFat"; obj[2] = "b" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgs" ; obj[2] = "b"; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}

		TransM = TransM + "MET" ;
		for ( int j = f-2; j > -1; j -- )
		{
			transM1 = "" ;		
			transM2 = "" ;		               
			if ( obj[j] == "no" )                                continue;    
			if( obj[j] == "Top" )                               TransM = TransM + obj[j] ;
			if(obj[j] == "Higgs")                               TransM = TransM + obj[j] ;
			if( obj[j] == "mu" )                                TransM = TransM + obj[j] ;   
			if(  obj[j] == "b" )                                 TransM = TransM + obj[j] ;
			if(  obj[j] == "HFat" )                             TransM = TransM + obj[j] ;     
			if( obj[j] == "TW" || obj[j] == "HW" ) {
				if ( u == 2 || u == 3 || u == 4 )   TransM = TransM + obj[j+1] + obj[j];
				else {
					TransM = TransM + obj[j] ;
				}
			}
			transM2 = TransM ;
			if(  obj[j] == "b" && ( u == 2 || u == 3)  )  transM2 = obj[j] + obj[j-1]; 
			if( obj[j] == "HFat" && u == 4 )               transM2 = obj[j] + obj[j-1];    

			sprintf( Histo_name, "Cat%d_Mt(%s)", x , transM2.c_str() ) ;
			sprintf( Histo_title, "Cat%d_Mt(%s) Distribution", x , transM2.c_str() ) ;			
			if(u <= 4){
			if ( j == 0 )	h_Histo_Mt.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 60, 0, 3000.0);
			if ( j != 0 )	h_Histo_Mt.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 300, 0, 3000.0);
			}
			if(u == 5 ) { 
			if ( j == 1 )	h_Histo_Mt.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 60, 0, 3000.0);
			if ( j != 1 )	h_Histo_Mt.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 300, 0, 3000.0);
			}
			h_Histo_Mt		  .at(u).at(j)         ->GetYaxis()->SetTitle("Events");   
			//	cout<<"h_Histo_Mt["<<u<<"]["<<j<<"].name() : "<< h_Histo_Mt.at(u).at(j)->GetName()<<endl;			
		} 
		TransM = "" ;
	}
}


void TTbar_PUPPItau_loc_190204::Category_Object_dRHisto()
{
	char Histo_name[100], Histo_title[100] , dp_name[100], dp_title[100];
	string  obj[5]  ;
	int x ;
	int f = 0 , g  ;
	//cout<<"\nFunc 12" ;		
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		g = -1 ; 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Topjet" ; obj[1] = "HWjet"; obj[2] = "muon" ; obj[3] = "MET";
			f = 3;
		}
		if ( u == 1 )  {
			obj[0] = "Topjet" ; obj[1] = "HFatjet"; obj[2] = "muon" ; obj[3] = "MET";
			f =3 ;
		}
		if ( u == 2 )  {
			obj[0] = "TWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 3 )  {
			obj[0] = "TWjet" ; obj[1] = "HWjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 4 )  {
			obj[0] = "HWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgsjet" ; obj[2] = "bjet"; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}

		for ( int j = 0; j < f; j ++ )
		{
			if ( obj[j] == "no" ) continue;
			for ( int k = j+1; k < f ; k++) 
			{
				g ++ ;
				if( obj[j] == "Topjet" ){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if(obj[j] == "Higgsjet") {
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "bjet" ) {
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "TWjet" || obj[j] == "HWjet" ){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "TFatjet" || obj[j] == "HFatjet"){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				} 

				sprintf(Histo_title,"Category%d_dR(%s_%s) Distribution",x ,obj[j].c_str(), obj[k].c_str() );
				sprintf(dp_title, "Category%d_dPt(%s_%s) Distribution",x ,obj[j].c_str(), obj[k].c_str() ) ;

				h_Histo_dR.at(u).at(g)             = new TH1F(Histo_name,Histo_title, 5000, 0.0, 5.0);
				h_Histo_dR.at(u).at(g)             ->GetYaxis()->SetTitle("Events");
				//cout<<"h_Histo_dR["<<u<<"]["<<g<<"].name() : "<< h_Histo_dR.at(u).at(g)->GetName()<<endl;				

				if( obj[k] == "muon"){
					h_Histo_dPt.at(u).at(g)        = new TH1F(dp_name,dp_title, 500, 0.0, 500.0);
					h_Histo_dPt.at(u).at(g)        ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_dPt["<<u<<"]["<<g<<"].name() : "<< h_Histo_dPt.at(u).at(g)->GetName()<<endl;					
				}

			} // end of kth loop 
		} // end of jth loop
	}  //  end of uth loop
}

//===============================================================


void  TTbar_PUPPItau_loc_190204::dR_Plots_Genlvl()
{
	//  for pt & eta topb= 0, topWl = 1, higgWq1 = 2, higgWq2= 3, higgWl=4, assob = 5, forwq = 6
	//  for dR , Here, i = 0-1 top decayed, i = 2-4 Higgs decayed, i = 5 assob, i = 6 forwq

	// int it = topW_q[0] ;
	//  int it = T_top ;
	//  int jet1 = topW_q[0] ;
	//   int j =  Higgs_Mu[0] ;    

	//int it = Higgs_Mu[0] ;  

	//int jet1 =  Higgs_W ;
	//int j = top_Mu[0] ;
	int jet1 =  Higgs_Jet[0] ;
	//int qj = Higgs_Mu[0] ;
	int j = topW_q[0] ;
	//Fill_MC_PtEta_Histo(4, qj) ;
	int A = -1 ; // for mom ID
	float dR = 0.0 ;

	h_mcobject_pt.at(0) ->Fill((*mc_MomPt)[b_top]);       
	h_mcobject_eta.at(0) ->Fill((*mc_MomEta)[b_top]); 
	h_mcobject_mass.at(0) ->Fill((*mc_MomMass)[b_top]);	       
	h_mcobject_pt.at(1) ->Fill((*mc_MomPt)[j]);       
	h_mcobject_eta.at(1) ->Fill((*mc_MomEta)[j]);
	h_mcobject_mass.at(1) ->Fill((*mc_MomMass)[j]);	     

	h_mcobject_pt.at(5) ->Fill((*mc_Pt)[b_asso]);       
	h_mcobject_eta.at(5) ->Fill((*mc_Eta)[b_asso]);   
	h_mcobject_pt.at(6) ->Fill((*mc_Pt)[q_forw]);       
	h_mcobject_eta.at(6) ->Fill((*mc_Eta)[q_forw]);   

	h_mcobject_pt.at(2) ->Fill((*mc_Pt)[T_higgs]);       
	h_mcobject_eta.at(2) ->Fill((*mc_Eta)[T_higgs]);
	h_mcobject_mass.at(2) ->Fill((*mc_MomMass)[T_higgs]);

	/*h_mcobject_pt.at(3) ->Fill((*mc_MomPt)[it]);       
	  h_mcobject_eta.at(3) ->Fill((*mc_MomEta)[it]);
	  h_mcobject_mass.at(3) ->Fill((*mc_MomMass)[it]);

	  h_mcobject_pt.at(4) ->Fill((*mc_Pt)[qj]);       
	  h_mcobject_eta.at(4) ->Fill((*mc_Eta)[qj]);
	  h_mcobject_mass.at(4) ->Fill((*mc_MomMass)[qj]);
	 */

	//for Higgs -> hadronically & top ->  leptonically
	/*
	   dR = dR_mc_mub( j, b_top , 1);   
	   h_dR_MC.at(0).at(1)  -> Fill(dR) ;       // dR (top muon, btop)	  	  	     
	   dR = dPt_Gen_jetlep(j, b_top , 1);
	   h_dPt_MC.at(0).at(1) -> Fill(dR) ; 

	   dR = dR_mc_mub(jet1, b_top, A) ;     
	   h_dR_MC.at(0).at(2)  -> Fill(dR) ;	   // (btop & higgs )

	   dR = dR_mc_mub(it, b_top, A) ;     
	   h_dR_MC.at(0).at(3)  -> Fill(dR) ;	   // (btop & higgs W)

	   dR = dR_mc_mub(b_top , b_asso, 1); 
	   h_dR_MC.at(0).at(5)  -> Fill(dR) ; 		

	   dR = dR_mc_mub(jet1, j, A) ;     
	   h_dR_MC.at(1).at(2)  -> Fill(dR) ;	   // (top muon & higgs )

	   dR = dR_mc_mub(it, j, A) ;     
	   h_dR_MC.at(1).at(3)  -> Fill(dR) ;	     // (W higgs & top muon)

	   dR = dR_mc_mub( j, b_asso , 1);   // dR (Wtop, basso)	 
	   h_dR_MC.at(1).at(5)  -> Fill(dR) ;  

	   dR = dR_mc_mub(it, jet1, A) ;     
	   h_dR_MC.at(2).at(3)  -> Fill(dR) ; 

	   dR = dR_mc_mub( jet1, b_asso , A);   // dR (WHiggs, basso)	 
	   h_dR_MC.at(2).at(5)  -> Fill(dR) ;  



	//for top -> hadronically & Higgs -> semi leptonically
	dR = dR_mc_mub( j, b_top , A);   
	h_dR_MC.at(0).at(1)  -> Fill(dR) ;       // dR (Wtop, btop)	  	  
	dR = dR_mc_mub(b_top, jet1, A) ;     
	h_dR_MC.at(0).at(2)  -> Fill(dR) ;	   // (top & higgs W)

	dR = dR_mc_mub(jet1, b_top, A) ;     
	h_dR_MC.at(0).at(3)  -> Fill(dR) ;	   // (btop & higgs W)

	dR = dR_mc_mub(b_top , b_asso, 1); 
	h_dR_MC.at(0).at(5)  -> Fill(dR) ; 	

	dR = dR_mc_mub(b_top , q_forw, 1); 
	h_dR_MC.at(0).at(6)  -> Fill(dR) ; 	

	dR = dR_mc_mub(b_top, qj, A) ;     
	h_dR_MC.at(0).at(4)  -> Fill(dR) ;	   // (btop & higgs muon)

	dR = dR_mc_mub(j, jet1, A) ;     
	h_dR_MC.at(1).at(2)  -> Fill(dR) ;	   // (btop & higgs W)

	dR = dR_mc_mub(j, qj, A) ;     
	h_dR_MC.at(1).at(4)  -> Fill(dR) ;	     // (Wtop & higgs muon)

	dR = dR_mc_mub( j, b_asso , A);   // dR (Wtop, basso)	 
	h_dR_MC.at(1).at(5)  -> Fill(dR) ;  
	dR = dR_mc_mub(j , q_forw, 1); 
	h_dR_MC.at(1).at(6)  -> Fill(dR) ; 	

	dR = dR_mc_mub(jet1, qj, A) ;     
	h_dR_MC.at(2).at(4)  -> Fill(dR) ;

	dR = dPt_Gen_jetlep(jet1 , qj, A);
	h_dPt_MC.at(2).at(4) -> Fill(dR) ;  

	dR = dR_mc_mub( jet1, b_asso , A);   // dR (WHiggs, basso) 
	h_dR_MC.at(2).at(5)  -> Fill(dR) ;  

	dR = dR_mc_mub(jet1 , q_forw, A); 
	h_dR_MC.at(2).at(6)  -> Fill(dR) ; 	
	*/

}



void   TTbar_PUPPItau_loc_190204::Fill_NPtEta_Histo(int id, int en, int idx)
{
	float ptau1 = 0.0 ;
	float puppi_mass = 0.0 , CHS_mass = 0.0 ; 
	//   id for defining particle type, ele =1, mu= 3, bjet = 05, AK8jet= 24(W) , qjet & forward jets = 13
	if (id == 0) {
		h_object_no.at(0) ->Fill(n_Mu.size());
		h_object_no.at(1) ->Fill(b_jet.size());
		h_object_no.at(2) ->Fill(n_AK8Jet.size()) ;
		h_object_no.at(3) ->Fill(n_jet.size()) ;
		h_object_no.at(4) ->Fill(n_forwjet.size() ) ;
		h_object_pt.at(0) ->Fill(pf_MET); 
	}

	if ( id == 3){
		h_object_pt.at(idx+1)    ->Fill((*mu_Pt)[en]); 
		h_object_eta.at(idx+1)   ->Fill((*mu_Eta)[en]);
	}

	if( id == 2){

		//h_object_eta.at(id) ->Fill(pf_METEta);
	}
	if( id == 5) {
		h_object_pt.at(idx+2) ->Fill((*jet_Pt)[en]); 
		h_object_eta.at(idx+2) ->Fill((*jet_Eta)[en]);
	}
	if (id == 24 ) {
		if (vr%2 == 1 ){
			puppi_mass = (*AK8_puppiSDMassL2L3Corr)[en] ;
			CHS_mass = (*AK8_JetPrunedMassCorr)[en] ; 
		}
		if (vr%2 == 0){
			puppi_mass = (*AK8_puppiSDMass)[en] ;
			CHS_mass = (*AK8_JetPrunedMass)[en] ;
		}
		if ( puppi_mass > 00.0) {
			h_object_pt.at(idx+6) ->Fill((*AK8_JetPt)[en]); 
			h_object_eta.at(idx+6) ->Fill((*AK8_JetEta)[en]);

			h_AK8_Jetmass.at(idx) -> Fill((*AK8_JetMass)[en]) ;
			h_AK8_PUPPImass.at(idx) -> Fill((*AK8_puppiMass)[en]);

			h_AK8_PUPPISDmass.at(idx) -> Fill(puppi_mass);
			h_AK8_CHSmass.at(idx) -> Fill(CHS_mass) ;		

			h_AK8_PUPPIvsmass.at(idx) -> Fill(puppi_mass,(*AK8_puppiMass)[en]) ;
			h_AK8_CHSvsmass.at(idx) -> Fill(CHS_mass,(*AK8_JetMass)[en]) ;

			//		h_AK8_PUPPIvsmass.at(idx) -> Fill((*mc_MomMass)[jet],(*AK8_JetMass)[en]) ;
			//	h_AK8_CHSvsmass.at(idx) -> Fill((*AK8_JetMass)[en],(*AK8_JetPt)[en]) ;
			//h_AK8_CHSvsmass.at(idx) -> Fill((*mc_MomMass)[jet],(*AK8_JetPt)[en]) ;

			ptau1 = ((*AK8_puppiTau2)[en]/(*AK8_puppiTau1)[en]) ;
			h_AK8_PUPPItau21.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_puppiTau3)[en]/(*AK8_puppiTau2)[en]) ;
			h_AK8_PUPPItau32.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_puppiTau3)[en]/(*AK8_puppiTau1)[en]) ;
			h_AK8_PUPPItau31.at(idx) -> Fill(ptau1);

			ptau1 = ((*AK8_puppiTau4)[en]/(*AK8_puppiTau1)[en]) ;
			h_AK8_PUPPItau41.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_puppiTau4)[en]/(*AK8_puppiTau2)[en]) ;
			h_AK8_PUPPItau42.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_puppiTau4)[en]/(*AK8_puppiTau3)[en]) ;
			h_AK8_PUPPItau43.at(idx) -> Fill(ptau1);

			ptau1 = ((*AK8_Jet_tau2)[en]/(*AK8_Jet_tau1)[en]) ;
			h_AK8_tau21.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_Jet_tau3)[en]/(*AK8_Jet_tau2)[en]) ;
			h_AK8_tau32.at(idx) -> Fill(ptau1);
			ptau1 = ((*AK8_Jet_CHStau4)[en]/(*AK8_Jet_CHStau2)[en]) ;
			h_AK8_CHStau42.at(idx) -> Fill(ptau1);
		}


	}  
	if(id == 14){
		h_object_pt.at(idx+4) ->Fill((*jet_Pt)[en]);
		h_object_eta.at(idx+4) ->Fill((*jet_Eta)[en]);
	}
}


void   TTbar_PUPPItau_loc_190204::Fill_MC_PtEta_Histo(int id, int en)
{
	//   id for defining particle type, 6 = asso bquark, 7= top bquark, 8-9 Wdeacay(Higgs) qi, 10-11 Wdeacay(top) qi, 12 = muon, 13= forw quark
	h_mcobject_pt.at(en) ->Fill((*mc_Pt)[id]); 
	//cout<<", ID = "<< id;
	h_mcobject_eta.at(en) ->Fill((*mc_Eta)[id]);
}


float TTbar_PUPPItau_loc_190204::Fill_MET_var(int lep, int lvl)
{
	float pass = 1.0 ;
	//bool pass = true ;
	float dPhi = 0.0 ;
	for(int h =0; h< n_MC; h++){
		if( abs((*mc_PID)[h]) == 14  && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 6) {
			dPhi = delta_phi((*mc_Phi)[h] , pf_METPhi);
		}
	}
	if (lvl == 1){
		h_MET_var.at(0) ->Fill(pf_MET);
		float sum_HT = pf_MET + (*mu_Pt)[lep] ;
		h_MET_var.at(1) ->Fill(sum_HT);
		h_MET_var.at(2) ->Fill(dPhi);
	}
	if ( dPhi < 0.6 ) pass = dPhi;
	return pass ;
}

void   TTbar_PUPPItau_loc_190204::Fill_Puppi_jet(int jet, int idx)
{// here idx = 0 for toptag, idx= 1,2 for Wtag, idx =3 for fatjet, idx =4,5 for Higgtag

	float ptau3 =   ((*AK8_puppiTau3)[jet]/(*AK8_puppiTau2)[jet]);
	float ptau2 =   ((*AK8_puppiTau2)[jet]/(*AK8_puppiTau1)[jet]) ;
	float ptau4 =   ((*AK8_Jet_CHStau4)[jet]/(*AK8_Jet_CHStau2)[jet]) ;
	float mass_puppi = 0.0;
	TLorentzVector v1 ;

	v1.SetPtEtaPhiE((*AK8_JetPt)[jet],(*AK8_JetEta)[jet],(*AK8_JetPhi)[jet],(*AK8_JetEn)[jet]);
	mass_puppi = v1.M();

}


bool TTbar_PUPPItau_loc_190204::Cut_Muon(int c_muon){

	bool pass_muon = true;
	if ( (*mu_Pt)[c_muon] <= 40.0 ) pass_muon = false;
	if ( fabs((*mu_Eta)[c_muon]) >= 2.1 ) pass_muon = false;
	bool tight_muon = (*mu_IDbit)[c_muon]>>2 & 1;    //for tight muon ID
	if ( tight_muon == 0 ) pass_muon = false;

	//cout << "Tight mu = " << tight_muon << endl ;
	// for medium ID        
	/*   //      if ((*mu_Chi2NDF)[c_muon] >= 3) pass_muon = false;
	     if ((*mu_chi2LocalPosition)[c_muon] >= 12 ) pass_muon = false;
	     if ((*mu_trkKink)[c_muon] >= 20 ) pass_muon = false;
	     if ((*mu_InnervalidFraction)[c_muon] <= 0.8) pass_muon = false;
	     if ((*mu_segmentCompatibility)[c_muon] <= 0.303 ) pass_muon = false;
	// for tight id
	//     if ((*mu_Chi2NDF)[c_muon] >= 10) pass_muon = false   ;
	if ((*mu_TrkLayers)[c_muon] <= 5  ) pass_muon = false ;
	if ((*mu_PixelHits)[c_muon] <= 0) pass_muon = false ;
	if ((*mu_MuonHits)[c_muon] <= 0 ) pass_muon = false ;
	if ((*mu_Stations)[c_muon] <= 1) pass_muon = false ;
	if((*mu_D0)[c_muon] >= 0.2) pass_muon = false  ;
	if((*mu_Dz)[c_muon] >= 0.5) pass_muon = false   ;  
	 */

	return pass_muon;

}


void TTbar_PUPPItau_loc_190204::Cut_bjet(int c_jet, int wp)
{     
	if( (*jet_Pt)[c_jet] >= 30.0){
		if (fabs((*jet_Eta)[c_jet]) < 2.4 ){ 
			if (wp == 1 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.2219) {
					b_jet_loose.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

			if (wp == 2 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.6324) {
					b_jet_medium.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

			if (wp == 3 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.8958) {
					b_jet_tight.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

		}
		else{
			n_forwjet.push_back(c_jet) ;

		}}}


void   TTbar_PUPPItau_loc_190204::Cut_AK8jet(int c_jet, int var)
{
	float cut1 = 0.0, cut2 = 1000.0 ;
	if (var <= 6) cut1 = 40.0 ;
	if (var == 0) {
	cut1 = 105.0;
	cut2 = 140.0;
	}
	if((*AK8_puppiPt)[c_jet] >= 200.0 && fabs((*AK8_puppiEta)[c_jet]) <= 2.4 && (*AK8_puppiSDMassL2L3Corr)[c_jet] >= cut1 && (*AK8_puppiSDMassL2L3Corr)[c_jet] <= cut2) n_AK8Jet.push_back(c_jet) ;
}


//===================dR distribution functions===================

float  TTbar_PUPPItau_loc_190204::dR_AK8jet(int Y, int Z, int id8, int lvl)
{ // here lvl value -1 for Puppi & CHS jets dR calculation & lvl = 1 for  dR with tagged jets
	bool pass_dr = false;
	float dPhi ; 
	float dEta;
	float dR   ;        
	if(lvl == 1){
		dPhi = delta_phi((*AK8_JetPhi)[Y], (*AK8_JetPhi)[Z]) ;
		dEta = fabs((*AK8_JetEta)[Y] - (*AK8_JetEta)[Z]);
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	}       
	else {
		dPhi = delta_phi((*AK8_JetPhi)[Y], (*AK8_puppiPhi)[Z]) ;
		dEta = fabs((*AK8_JetEta)[Y] - (*AK8_puppiEta)[Z]);
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	}
	return dR;
}


float  TTbar_PUPPItau_loc_190204::dRPlots_bjet(int Y, int Z, int lvl)
{
	// for bjets only  
	float dPhi = delta_phi((*jet_Phi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*jet_Eta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;     
	return dR;
}


float  TTbar_PUPPItau_loc_190204::dR_AK8bjet(int Y, int Z, int id8, int idb, int lvl)
{
	// for topjet id8 = 0, fatjet id8 = 1, Wtag id8 = 2,3
	bool pass_dr = false;
	int d =  id8 +4*idb;

	float dPhi = delta_phi((*AK8_JetPhi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*AK8_JetEta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;     
	return dR;
}

float  TTbar_PUPPItau_loc_190204::dR_mc_mub(int Y, int Z, int A)
{
	// here Z for muon/bquark & Y for top/Wtop & Higgs/Whiggs
	float phi, eta ;
	if ( A == -1 ) {
		phi = (*mc_MomPhi)[Y] ;
		eta = (*mc_MomEta)[Y] ;
	}  
	else {
		phi = (*mc_Phi)[Y] ;
		eta = (*mc_Eta)[Y] ;  
	}
	float dPhi = delta_phi(phi , (*mc_Phi)[Z]) ;
	float dEta = eta - (*mc_Eta)[Z];
	float dR   = Sqrt(dPhi*dPhi + dEta*dEta) ;               
	return dR ; 
}


void  TTbar_PUPPItau_loc_190204::dR_mc_qiqj(int Y, int Z,int A)
{
	bool pass_dr = true;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*mc_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*mc_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
	pass_dr = false;
}

void  TTbar_PUPPItau_loc_190204::dR_mc_bqi(int Y, int Z,int A)
{
	bool pass_dr = true;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*mc_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*mc_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
	pass_dr = false;
}


float  TTbar_PUPPItau_loc_190204::dR_mcbjet(int Y, int Z, int X, int idx, int lvl)
{
	float pass_dr = 1.0;
	int d = 2*X + Z;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR; 
}


float TTbar_PUPPItau_loc_190204::dR_mcAK8(int Y, int Z, int X, int idx)
{
	bool pass_dr = true;
	int d = 2*X + idx;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*AK8_JetPhi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*AK8_JetEta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
}

float  TTbar_PUPPItau_loc_190204::dR_mu(int Y, int Z_jet, int idx, int lvl)
{
	float pass_dr = 2.0;
	double dPhi = delta_phi((*mu_Phi)[Y] , (*jet_Phi)[Z_jet]) ;
	double dEta = (*mu_Eta)[Y] - (*jet_Eta)[Z_jet];
	double dR   = sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR;
}

float TTbar_PUPPItau_loc_190204::dR_mu_AK8(int Y, int Z_jet, int idx, int lvl)
{
	bool pass_dr = true;
	double dPhi = delta_phi((*mu_Phi)[Y] , (*AK8_JetPhi)[Z_jet]) ;
	double dEta = (*mu_Eta)[Y] - (*AK8_JetEta)[Z_jet];
	double dR   = sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR;
}


float TTbar_PUPPItau_loc_190204::dPt_lep(int jet, int lep, int idx, int lvl)
{
	float  pass_dPt =  0.0;
	float  vec_mag = 0.0 ;
	TVector3 vec_Xprod, vec_lep, vec_jet;
	float dPt = 0.0 ;

	vec_lep.SetPtEtaPhi( (*mu_Pt)[lep], (*mu_Eta)[lep] , (*mu_Phi)[lep] ) ;
	if(idx == 0) vec_jet.SetPtEtaPhi( (*jet_Pt)[jet], (*jet_Eta)[jet], (*jet_Phi)[jet] ) ;
	else {
		vec_jet.SetPtEtaPhi( (*AK8_JetPt)[jet], (*AK8_JetEta)[jet], (*AK8_JetPhi)[jet] ) ; 
	}    

	vec_Xprod = vec_lep. Cross(vec_jet);
	vec_mag = vec_jet.Mag();
	dPt =  (vec_Xprod.Mag() ) / vec_mag ;
	return dPt ;
}


float TTbar_PUPPItau_loc_190204::dPt_Gen_jetlep(int jet, int lep, int idx)
{
	// idx = -1 for Mom & idx = 1 for daughter
	float  pass_dPt =  0.0;
	float  vec_mag = 0.0 ;
	TVector3 vec_Xprod, vec_lep, vec_jet;
	float dPt = 0.0 ;

	vec_lep.SetPtEtaPhi( (*mc_Pt)[lep], (*mc_Eta)[lep] , (*mc_Phi)[lep] ) ;
	if(idx == -1) vec_jet.SetPtEtaPhi( (*mc_MomPt)[jet], (*mc_MomEta)[jet] , (*mc_MomPhi)[jet] ) ;
	if(idx == 1 ) vec_jet.SetPtEtaPhi( (*mc_Pt)[jet], (*mc_Eta)[jet] , (*mc_Phi)[jet] ) ;

	vec_Xprod = vec_lep. Cross(vec_jet);
	vec_mag = vec_jet.Mag();
	dPt =  (vec_Xprod.Mag() ) / vec_mag ;                      
	return dPt ;
}


void  TTbar_PUPPItau_loc_190204::dR_qjet_objects(int qjet, int obj, int idq, int idob)   
{
	// for reco, mu = 0, bjets = 1&2, ak8jet = 3&4, qjet = 13. for gen, assob =5, topb = 6, q from W(higgs) = 7 & 8, q from W(top) = 9 & 10, mu = 11, forw q = 12,
	bool pass_dr = false;
	float dR = 0.0;
	float dPhi = 0.0 , dEta  = 0.0 ;        
	if( idob == 0){         
		dPhi = delta_phi((*mu_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*mu_Eta)[obj] - (*jet_Eta)[qjet];
	}   

	if( idob >= 5 && idob <=12 ){
		dPhi = delta_phi((*mc_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*mc_Eta)[obj] - (*jet_Eta)[qjet];       
	}

	if( idob ==1 || idob == 2){
		dPhi = delta_phi((*jet_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*jet_Eta)[obj] - (*jet_Eta)[qjet];
	}

	if( idob ==3|| idob == 4){
		dPhi = delta_phi((*AK8_JetPhi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*AK8_JetEta)[obj] - (*jet_Eta)[qjet];
	}

	if( idob <= 12 ) {
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;

	}

	//        return pass_dr;
}


float TTbar_PUPPItau_loc_190204::delta_phi(float phi1,float phi2)
{
	const float PI=2.0*acos(0.);
	const float TWOPI=2.0*PI;
	float PHI=fabs( phi1 - phi2 ) ;
	return (PHI<=PI)? PHI : TWOPI-PHI;

}
// ---------------------------------------Reco Object Plots----------------------------------------

void TTbar_PUPPItau_loc_190204::GenvsRecoMass(int obj, int ID )
{
	//   New defination of histogram
	float dR = 0.0 ;	
	int jet = -1 ;
	float obj_mass = (*mc_Mass)[obj];
	float jet_mass = 0.0;
	// h_AK81vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
	for(Int_t k = 0 ; k < n_AK8Jet.size(); k++) {
		jet = 	n_AK8Jet[k] ;
		//continue ;
		jet_mass = (*AK8_puppiMass)[jet];
		dR = dR_mcAK8( obj, jet, 0 , 0); 
		if ( dR > 0.1 ) continue;		

		if( k < 3) h_dR_MCReco[k][ID] -> Fill(dR);
		if ( k == 0){
		 h_AK81vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
		 if( ID == 0) Fill_NPtEta_Histo(24,jet,0) ;
		 if( ID == 1 || ID == 3) Fill_NPtEta_Histo(24,jet,1) ;
		 if( ID == 2) Fill_NPtEta_Histo(24,jet,2) ;
		}
		if ( k == 1) {
		 h_AK82vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
		 if( ID == 0) Fill_NPtEta_Histo(24,jet,0) ;
		 if( ID == 1 || ID == 3) Fill_NPtEta_Histo(24,jet,1) ;
		 if( ID == 2) Fill_NPtEta_Histo(24,jet,2) ;
		}

		if ( k == 2) h_AK83vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
	}
}


void TTbar_PUPPItau_loc_190204::GenvsRecoObj(int obj, int ID )
{
	//   New defination of histogram
	float dR = 0.0 ;	
	int jet = -1, quark = -1  ;
	int pass = 0 ;
	//float obj_mass = (*mc_Mass)[obj];
	float jet_mass = 0.0;
	int Msize = -1 ;
	if (ID == 0 ) Msize = Higgs_Jet.size() ;
	if (ID == 2 ) Msize = topW_q.size() ; 
	// h_AK81vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
	int AK8size  = ( n_AK8Jet.size() >= 3) ? 3 : n_AK8Jet.size() ; 
	for(Int_t k = 0 ; k < AK8size; k++) {
		jet = 	n_AK8Jet[k] ;
		//continue ;
		//jet_mass = (*AK8_puppiMass)[jet];
		pass = 0 ;
		for (Int_t m = 0 ; m < Msize ; m ++ ) {

	 	if ( ID == 0 ) quark = Higgs_Jet[m] ;
		if ( ID == 2 ) quark = topW_q[m] ;

		dR = dR_mcAK8( quark, jet, 0 , 0); 
		if ( dR < 0.1 ) pass ++ ;		
		}
		
		if ( pass == 4){
		 //if( ID == 0) Fill_NPtEta_Histo(24,jet,0) ;
		 if( ID == 0) Higgs_GenJet.push_back(jet);
		 }
		if ( pass == 2 && ID == 2) {	
	//	 if( ID == 1 || ID == 3) Fill_NPtEta_Histo(24,jet,1) ;
		 if (dR_mcAK8( b_top, jet, 0 , 0) < 0.1 ){
		 Top_GenJet.push_back(jet);	
		 //Fill_NPtEta_Histo(24,jet,2) ;
		}
		}
	//	if ( k == 2) h_AK83vsGenObj.at(ID) ->Fill(jet_mass, obj_mass);
	}
}



void TTbar_PUPPItau_loc_190204::RecoPlots_dRHisto()
{
	float dR = 0.0;
	int b1 = -1, b2 =-1 ;
	int Msize  = ( n_Mu.size() >= 1) ? 1 : n_Mu.size() ;
	int Asize  = ( n_AK8Jet.size() >= 2) ? 2 : n_AK8Jet.size() ;
	int Bsize  = ( b_jet.size() >= 2) ? 2 : b_jet.size() ;
	int Fsize  = ( n_forwjet.size() >= 2) ? 2 : n_forwjet.size() ;                
	// for muon w.r.t other jets dR plots
	for ( int h = 0 ; h < Msize ; h ++ ){
		b1 = n_Mu[h];   
		for( int g = 0 ; g < Bsize ; g ++ ) {
			b2 = b_jet[g];
			dR = dR_mu(b1, b2, b1, b1);
			h_dR_Reco.at(h).at(g+1) ->Fill(dR) ;
			dR = dPt_lep(b2, b1, 0, b1);
			h_dPt_Reco.at(h).at(g+1) ->Fill(dR) ;
		}

		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = n_forwjet[g];
			dR = dR_mu(b1, b2, b1, b1);
			h_dR_Reco.at(h).at(g+3) ->Fill(dR) ;
			dR = dPt_lep(b2, b1, 0, b1);
			h_dPt_Reco.at(h).at(g+3) ->Fill(dR) ;
		}

		for( int g = 0 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_mu_AK8(b1, b2, b1, b2);
			h_dR_Reco.at(h).at(g+5) ->Fill(dR) ;
			dR = dPt_lep(b2, b1, 1, b1);
			h_dPt_Reco.at(h).at(g+5) ->Fill(dR) ;
		}  
	}

	// for bjet w.r.t other jets dR plots
	for( int h = 0 ; h < Bsize ; h ++ ) {
		b1 = b_jet[h];
		for( int g = h+1 ; g < Bsize ; g ++ ) {
			b2 = b_jet[g];  
			dR = dRPlots_bjet(b1, b2, -1);
			h_dR_Reco.at(h+1).at(g+2) ->Fill(dR) ;  
		}

		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = n_forwjet[g];
			dR = dRPlots_bjet(b1, b2, -1);
			h_dR_Reco.at(h+1).at(g+3) ->Fill(dR) ;
		}

		for( int g = 0 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_AK8bjet(b2, b1, b1, 0, 0);
			h_dR_Reco.at(h+1).at(g+5) ->Fill(dR) ; 
		}  
	}

	// for AK8jet dR Plots
	for( int h = 0 ; h < Asize ; h ++ ) {
		b1 = n_AK8Jet[h];
		for( int g = h+1 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_AK8jet(b1, b2, 0, 1);
			h_dR_Reco.at(h+5).at(g+5) ->Fill(dR) ; 
		}  
		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = n_forwjet[g];
			dR = dR_AK8bjet(b1, b2, b1, 0, 0);
			h_dR_Reco.at(g+3).at(h+5) ->Fill(dR) ;
		}

	}
}



void TTbar_PUPPItau_loc_190204::Fill_RecoObject()
{
	int Msize  = ( n_Mu.size() >= 1) ? 1 : n_Mu.size() ;             // id = 3 for muon
	int Asize  = ( n_AK8Jet.size() >= 3) ? 3 : n_AK8Jet.size() ;    // id = 24 for AK8jet
	int Bsize  = ( b_jet.size() >= 2) ? 2 : b_jet.size() ;                // id = 5 for bje
	int Fsize  = ( n_forwjet.size() >= 2) ? 2 : n_forwjet.size() ;     
	int b2 = -1 ;

	Fill_NPtEta_Histo(0, 0, 0) ;                                     // for MET & object population

	for ( int i = 0; i < Msize; i ++)
	{    b2 = n_Mu[i] ;
		Fill_NPtEta_Histo(3, b2, i) ;
	}

	for ( int i = 0; i < Bsize; i ++)
	{    b2 = b_jet[i] ;
		Fill_NPtEta_Histo(5, b2, i) ;
	}

	for ( int i = 0; i < Asize; i ++)
	{    b2 = n_AK8Jet[i] ;
		if ( WtagMatch(b2) == -1) ;
		if ( (*AK8_JetPt)[b2] > 20.0)Fill_NPtEta_Histo(24,b2,i) ;
	}

	for ( int i = 0; i < Fsize; i ++)
	{    b2 = n_forwjet[i] ;
		Fill_NPtEta_Histo(14, b2, i) ;
	}

	RecoPlots_dRHisto() ;  
}
// ================= Jet Selection optimisation========================

void  TTbar_PUPPItau_loc_190204::Higgs_Optimisation(int var)  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ; 
	int sel = -1 ;
	int pm  = -1 ;
	float Pruned_mass;
	float min = 0.0 ;
	float tau42 = 0.0, tau32 = 0.0 ;
	float t42 = PUPPItau[var-1] ;
	float max = 0.0 ;
	float pt = 0.0 ;
	int Wsize  = -1 ;
	if(var <= 3 ) Wsize  = (Higgs_GenJet.size() >= 2) ? 2 : Higgs_GenJet.size() ;
	if(var >= 4 ) Wsize  = (Top_GenJet.size() >= 2)   ? 2 : Top_GenJet.size() ;
	//t42 = 0.75 ;
	//if (var%2 == 0) t42 = 0.65 ;

	for (int g = 0 ; g < Wsize ; g ++ ){
		//if( g==0) continue ;
		//if ( Higgs_GenJet[g] == -1 )   continue;
		//cout << "\n Higgs=" << g;
		pm = -1 ;
		sel = 1 ;

		if (var <= 3){
		//	Pruned_mass 	= 	(*AK8_puppiSDMass)[jt] ;
			jt = Higgs_GenJet[g] ;

			Pruned_mass 	= 	(*AK8_puppiSDMassL2L3Corr)[jt] ;
		 	tau42 		= 	((*AK8_puppiTau4)[jt]/(*AK8_puppiTau2)[jt]) ;
			pt     	        = 	(*AK8_puppiPt)[jt] ; 
			if( Pruned_mass >=  90.0 && Pruned_mass <= 140.0 ) pm = 1 ;
			//if( g == 0 )min = 	fabs(Pruned_mass - 125.0) ;
	        	//max 	     	=  	fabs(Pruned_mass - 125.0) ;
			sel 		= 	-1 ;
		/*for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.8484) sel = 1;
			} */

		}
	
		if (var >= 4 ) {
			jt 		= 	Top_GenJet[g]   ;
	        	Pruned_mass  	= 	(*AK8_puppiSDMass)[jt] ;		
	        	tau42 	     	= 	((*AK8_puppiTau4)[jt]/(*AK8_puppiTau2)[jt]);
	        	pt 	     	=	(*AK8_puppiPt)[jt] ;
			if( Pruned_mass > 105.0 && Pruned_mass < 140.0 ) pm = 1;
			tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;
			if( tau32 > 0.65 ) sel = -1 ; 
	        }


		//for tag jet selection     

		if 	( tau42 > t42   ) 	                       continue ;
		if      ( pt    < 300.0 )   			       continue ;
		if	( pm  == -1 ) 				       continue ;
		//if 	( sel ==  1 ) 				       continue ;
		if( var <= 3) Higgsjets.push_back(jt);			
		if( var >= 4) {
		if 	( sel ==  1 ) topjet.push_back(jt) ;
		}	  
	}

	 
	
//v_tree ->Draw("AK8_puppiSDMassL2L3Corr[0]: AK8_puppiTau4[0]/AK8_puppiTau2[0]", "AK8_puppiPt[0] > 300.0 && abs(AK8_puppiEta[0]) < 2.4","colz");
//v_tree ->Draw("AK8_JetPrunedMassCorr[0]: AK8_Jet_CHStau4[0]/AK8_Jet_CHStau2[0]", "AK8_JetPt[0] > 300.0 && abs(AK8_JetEta[0]) < 2.4","colz"); 
}


void  TTbar_PUPPItau_loc_190204::Top_Optimisation( int var )  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ;
	int sel = -1 ;
	float tau32 = 0.0 ;
	float SDmass = 0.0  ;
	int t32 = -1 ;
	int SD = -1 ;
	float min = 0.0 ;
	float max = 0.0 ;
	int Itop = -1 ;
	int pt = -1;
//	float tau1 = taIuval[var-1] ;
	float tau1 = 0.65 ;
	int Wsize  = (Top_GenJet.size() >= 2) ? 2 : Top_GenJet.size() ;

	for (int g = 0 ; g < Wsize ; g ++ ){
		//if( g ==0) continue ;
		//cout << "\n Tsize = " << Wsize ;
		if ( Top_GenJet[g] == -1 )   continue;
		//cout << "\n Topp=" << g;
		sel = -1 ;
		jt = Top_GenJet[g] ;
		t32 = -1;
		SD = -1 ;
		pt = -1 ;
		if( var <= 6) {
			tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;

			if(tau32 < tau1 ) t32 = 1 ;
		//	SDmass = (*AK8_puppiSDMassL2L3Corr)[jt];
			SDmass = (*AK8_puppiSDMass)[jt];
			if( g == 0 ) min = fabs(SDmass - 172.5) ;
			max = fabs(SDmass - 172.5) ;
			if ((*AK8_puppiPt)[jt] > 400.0) pt = 1;
			if ( SDmass > 140.0 && SDmass <= 210.0 ) SD = 1;
			for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.5426) sel = 1;
			}
		} 

		if( var%2 == 100 ){
			tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;
			if(tau32 < tau1 ) t32 = 1 ;
			SDmass = (*AK8_puppiSDMassL2L3Corr)[jt];
			if( g == 0 ) min = fabs(SDmass - 172.5) ;
			max = fabs(SDmass - 172.5) ;
			if ((*AK8_puppiPt)[jt] > 400.0) pt = 1;
			if ( SDmass > 140.0 && SDmass < 210.0 ) SD = 1;
			for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.5426) sel = 1;
			}
		}

		/*if( var == 44 ){
			tau32 = ((*AK8_Jet_tau3)[jt]/(*AK8_Jet_tau2)[jt]) ;
			if(tau32 < 0.67 ) t32 = 1 ;
			//SDmass = (*AK8_JetPrunedMass)[jt];
			SDmass = (*AK8_JetPrunedMassCorr)[jt];
			if ((*AK8_JetPt)[jt] > 400.0) pt = 1 ;
			if( g == 0 ) min = fabs(SDmass -172.5) ;
			max = fabs(SDmass - 172.5) ;
			if ( SDmass >105.0 && SDmass < 220.0 ) SD = 1;
			for(int jg = 0 ; jg < (*n_AK8SDSJ)[jt] ; jg ++ )
			{
			if ((*AK8_SDSJCSV)[jt][jg] > 0.5426) sel = 1;
			}

		} */


		// for top jet selection     
		if (t32 == -1 ) continue ;
		if (pt ==  -1) continue ;
		if (SD == -1 ) continue ;
		if (sel == -1) continue ;
		//if( min < max) continue ;
		//min = max ;
		Itop = g ;
		topjet.push_back(jt) ;
		//n_AK8Jet[g] = -1 ;


	}
	if ( Itop == -10) {
		jt = n_AK8Jet[Itop] ;
		topjet.push_back(jt) ;
		n_AK8Jet[Itop] = -1 ;

	}

}

//===============Jet Category Selection=================================
void  TTbar_PUPPItau_loc_190204::Higgs_selection(int var)  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ; 
	int sel = -1 ;
	int pm  = -1 ;
	float Pruned_mass;
	float min = 0.0 ;
	float tau42 = 0.0 ;
	//float t42 = PUPPItau[var-1] ;
	float t42 = 0.70 ;
	float max = 0.0 ;
	float pt = 0.0 ;
	int Wsize  = ( n_AK8Jet.size() >= 3) ? 3 :  n_AK8Jet.size() ;
	
	//t42 = 0.75 ;
	//if (var%2 == 0) t42 = 0.65 ;

	for (int g = 0 ; g < Wsize ; g ++ ){
		//if( g==0) continue ;
		if ( n_AK8Jet[g] == -1 )   continue;
		//cout << "\n Higgs=" << g;
		jt = n_AK8Jet[g] ;
		pm = -1 ;

		if (var <= 6){
		//	Pruned_mass 	= 	(*AK8_puppiSDMass)[jt] ;
			Pruned_mass 	= 	(*AK8_puppiSDMassL2L3Corr)[jt] ;
		 	tau42 		= 	((*AK8_puppiTau4)[jt]/(*AK8_puppiTau2)[jt]) ;
			pt     	        = 	(*AK8_puppiPt)[jt] ; 
			if( Pruned_mass >=  90.0 && Pruned_mass <= 140.0 ) pm = 1 ;
			if( g == 0 )min = 	fabs(Pruned_mass - 125.0) ;
	        	max 	     	=  	fabs(Pruned_mass - 125.0) ;
			sel 		= 	-1 ;
		/*for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.9535) sel = 1;
			}
		*/
		}
		
	
		if (var%2 == 100) {
	        	Pruned_mass  	= 	(*AK8_JetPrunedMassCorr)[jt] ;		
	        	tau42 	     	= 	((*AK8_Jet_CHStau4)[jt]/(*AK8_Jet_CHStau2)[jt]);
	        	pt 	     	=	(*AK8_JetPt)[jt] ;
			if( Pruned_mass > 110.0 && Pruned_mass < 135.0 ) pm = 1;
	        	if( g == 0 )min = 	fabs(Pruned_mass - 125.0) ;
	        	max 	     	=  	fabs(Pruned_mass - 125.0) ;
		}


		//for higgs jet selection     

		if 	( tau42 > t42   ) 	                       continue ;
		if      ( pt    < 300.0 )   			       continue ;
		if	( pm  == -1 ) 				       continue ;
		if 	( sel ==  1 ) 				       continue ;
		//if	( min < max) 				       continue ;
		//min = max ;
			//sel = g ;
			Higgsjets.push_back(jt);
			n_AK8Jet[g] = -1 ;    
	}

		if ( sel == -10){
		jt = n_AK8Jet[sel] ;
		Higgsjets.push_back(jt);
		n_AK8Jet[sel] = -1 ;    
		}
//v_tree ->Draw("AK8_puppiSDMassL2L3Corr[0]: AK8_puppiTau4[0]/AK8_puppiTau2[0]", "AK8_puppiPt[0] > 300.0 && abs(AK8_puppiEta[0]) < 2.4","colz");
//v_tree ->Draw("AK8_JetPrunedMassCorr[0]: AK8_Jet_CHStau4[0]/AK8_Jet_CHStau2[0]", "AK8_JetPt[0] > 300.0 && abs(AK8_JetEta[0]) < 2.4","colz"); 
}


void  TTbar_PUPPItau_loc_190204::Top_selection( int var )  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ;
	int sel = -1 ;
	float tau32 = 0.0 ;
	float SDmass = 0.0  ;
	int t32 = -1 ;
	int SD = -1 ;
	float min = 0.0 ;
	float max = 0.0 ;
	int Itop = -1 ;
	int pt = -1;
//	float tau1 = taIuval[var-1] ;
	float tau1 = 0.80 ;
	int Wsize  = ( n_AK8Jet.size() >= 3) ? 3 :  n_AK8Jet.size() ;

	for (int g = 0 ; g < Wsize ; g ++ ){
		//if( g ==0) continue ;
		//cout << "\n Tsize = " << Wsize ;
		if ( n_AK8Jet[g] == -1 )   continue;
		//cout << "\n Topp=" << g;
		sel = -1 ;
		jt = n_AK8Jet[g] ;
		t32 = -1;
		SD = -1 ;
		pt = -1 ;
		if( var <= 6) {
			tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;

			if(tau32 < tau1 ) t32 = 1 ;
		//	SDmass = (*AK8_puppiSDMassL2L3Corr)[jt];
			SDmass = (*AK8_puppiSDMass)[jt];
			if( g == 0 ) min = fabs(SDmass - 172.5) ;
			max = fabs(SDmass - 172.5) ;
			if ((*AK8_puppiPt)[jt] > 400.0) pt = 1;
			if ( SDmass > 105.0 && SDmass <= 210.0 ) SD = 1;
			for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.5426) sel = 1;
			}
		}

		if( var%2 == 100 ){
			tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;
			if(tau32 < tau1 ) t32 = 1 ;
			SDmass = (*AK8_puppiSDMassL2L3Corr)[jt];
			if( g == 0 ) min = fabs(SDmass - 172.5) ;
			max = fabs(SDmass - 172.5) ;
			if ((*AK8_puppiPt)[jt] > 400.0) pt = 1;
			if ( SDmass > 140.0 && SDmass < 210.0 ) SD = 1;
			for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
			{
				if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.5426) sel = 1;
			}
		}

	
		// for top jet selection     
		if (t32 == -1 ) continue ;
		if (pt ==  -1) continue ;
		if (SD == -1 ) continue ;
		if (sel == -1) continue ;
		//if( min < max) continue ;
		//min = max ;
		Itop = g ;
		topjet.push_back(jt) ;
		n_AK8Jet[g] = -1 ;


	}
	if ( Itop == -10) {
		jt = n_AK8Jet[Itop] ;
		topjet.push_back(jt) ;
		n_AK8Jet[Itop] = -1 ;

	}

}


void  TTbar_PUPPItau_loc_190204::Wjet_selection(int var)   // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1  ; 
	int sel = -1 ;
	float tau21 ;
	float SDmass, SDmass1 ;
	int t21 = -1;
	int SD  = -1 ;
	int SD1 = -1 ;
	float min = 0.0 ;
	float max = 0.0 ;
	int pt = -1; 
	int Wsize  = ( n_AK8Jet.size() >= 3) ? 3 :  n_AK8Jet.size() ;
	//int Wsize  =  n_AK8Jet.size()  ;
	float tau1 = tau_21[var-1] ;
	tau1	= 0.55;	

	for (int g = 0 ; g < Wsize ; g ++ ){
		//if( g ==0) continue ;
		//cout << "\n Wiggs=" << Wsize;
		if ( n_AK8Jet[g] == -1 )   continue;
		//cout << "\n Wjet=" << g;
		jt = n_AK8Jet[g] ;
		t21 = -1 ;
		SD = -1 ;
		pt = -1 ;
		if( var <= 6) {
			tau21 = ((*AK8_puppiTau2)[jt]/(*AK8_puppiTau1)[jt]) ;
			if(tau21 < tau1) t21 = 1 ;
			SDmass  = (*AK8_puppiSDMassL2L3Corr)[jt];
			SDmass1 = (*AK8_JetPrunedMassCorr)[jt];
			if( g == 0 ) min = fabs(SDmass - 80.5) ;
			max = fabs(SDmass -  80.5) ;
			if ( SDmass >=  65.0 && SDmass <   105.0 ) SD = 1;
		//	if ( SDmass1 <   105.0 ) SD1 = 1;
			if ((*AK8_puppiPt)[jt] > 200.0) pt = 1;
		}


		if( var%2 == 100) {
			tau21 = ((*AK8_puppiTau2)[jt]/(*AK8_puppiTau1)[jt]) ;
			if(tau21 < tau1) t21 = 1 ;
			SDmass = (*AK8_puppiSDMassL2L3Corr)[jt];
			if( g == 0 ) min = fabs(SDmass - 80.5) ;
			max = fabs(SDmass -  80.5) ;
			if ( SDmass >  65.0 && SDmass <  105.0 ) SD = 1;
			if ((*AK8_puppiPt)[jt] > 200.0) pt = 1;
		}

		// for Wjet selection     
		if ( t21 == -1 ) continue ;
		if ( pt  == -1 ) continue ;
		if ( SD  == -1 ) continue ;
		//if ( SD1 == -1 ) continue ;
		//if( min < max) continue ;
		//min = max ;
		sel = g ;
		W_boson.push_back(jt) ;
		n_AK8Jet[g] = -1 ;

	}

	if (sel == -10){
		jt = n_AK8Jet[sel] ;
		W_boson.push_back(jt) ;
		n_AK8Jet[sel] = -1 ;

	}

}


void  TTbar_PUPPItau_loc_190204::Fatjet_selection()
{
	int g3 = -1 ;
	float dR2 = 0.0 ;  

	for (int g = 0 ; g < n_AK8Jet.size() ; g ++ )
	{
		g3 = n_AK8Jet[g]  ;
		if( g3 == -1 ) continue ;
		if ((*AK8_JetPt)[g3] < 170.0) continue ;
		dR2 =((*AK8_puppiTau2)[g3]/(*AK8_puppiTau1)[g3]) ;
		if ( ! ( dR2 > 0.00 && dR2 < 0.75) ) continue ;
		//if ( ! (  dR2 < 0.80) ) continue ;
		if ( (*AK8_puppiSDMass)[g3] > 65.0 && (*AK8_puppiSDMass)[g3] < 210.0 ) continue ;
		//if ( (*AK8_puppiSDMass)[g3] <=  1.0 ) continue ;  
		//if ( (*AK8_JetPrunedMass)[g3] > 65.0 && (*AK8_puppiSDMass)[g3] < 210.0 ) continue ;
		fat_jet.push_back(g3) ;
	}

}

int  TTbar_PUPPItau_loc_190204::WtagMatch(int Wjet)
{
	//int Wjet = W_boson[0] ;
	//int q = Higgs_Jet[0] ;
	int q = Higgs_W ;
	float dPhi = delta_phi((*mc_MomPhi)[q] , (*AK8_JetPhi)[Wjet]) ;
	float dEta = (*mc_MomEta)[q] - (*AK8_JetEta)[Wjet];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       

//	float dR = dR_mcAK8(q , Wjet, 0 , 0) ;
	int match = - 1;
	if (dR < 0.1 ) match = 1 ;	
	return match ;	

}

// ============================Category Study =================
void     TTbar_PUPPItau_loc_190204::Wtag0_Category()
{		
	if (Higgsjets.size() != 0 ) {
			event_Higgs ++ ;
			if (b_jet.size() != 0 ) {
				Higgs_lbjet_Plots() ;	
				event_bjet ++ ;		                 
			}
		}

/*	else{
	if (topjet.size() != 0) {
		if ( fat_jet.size() != 0 ) {
			//top_fatjet_Plots() ;	
			eventII ++ ;		                 
		}
	} 
	}*/
}

void     TTbar_PUPPItau_loc_190204::Wtag1_Category()
{
	if (topjet.size() != 0 ) {
		top_Wjet_Plots() ; 
		event_top ++ ;		                 		}
/*	else{
		if (b_jet.size() != 0 ) {

			if ( fat_jet.size() != 0 ) W_fatjet_Plots() ;   
		}
	} 
*/
}
// ======================Signal Profile Plots ==================
void  TTbar_PUPPItau_loc_190204::Wjet_Plots(int bjet ) 
{
	int j = -1 ;
	int jet = -1;
	//	int jet = Higgs_Jet[0];
	int fill = -1 ;
	int Pt_idx = 0; // for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	float mass_chs = 0.0 ;
	float max = 0.0 ;
	float min = 0.0 ;
	TLorentzVector  Ws ;
	int Wsize        =   W_boson.size()  ;
	if ( Wsize != 0 )  h_tag_N.at(0)    ->Fill(Wsize) ;

	for (int k = 0 ; k < Wsize ; k ++ ){
		t2		=  	W_boson[k] ;
		if(k == 1 )t3	=	t2 ;
		mass_puppi 	=	(*AK8_puppiSDMassL2L3Corr)[t2] ;	
		if( k == 0 ){
		min 		= 	fabs(mass_puppi - 80.0);
		j   		=  	t2 ;
		}
	        max 	     	=  	fabs(mass_puppi - 80.0) ;		
		if( min > max) {
		min 		= 	max ;
		t3 		=	j   ;
		j   		= 	t2  ;
		}		
		}
	

	if (t3 != -1 && Wsize != 0 ){
		W_boson[0] = j  ;
		W_boson[1] = t3 ;
	}

	Wsize  = ( W_boson.size() >= 2) ? 2 : W_boson.size() ;

	for (int k = 0 ; k < Wsize ; k ++ ) {
		

		t2		=  	W_boson[k] ;

		fill =   -1 ;	
		if (bjet == 1)   {		
		tau =  ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		if (tau >= 0.00) fill = 1 ;
		}

		if (bjet == 2 ) {		
		tau =  ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		if (tau < 0.00) fill = 1 ;
		}
		
		//if ( fill == -1 ) continue ;

		if( bjet <= 6 ) {
			//tau  = (*AK8_puppiPt)[t2]; 
			//tau1 = (*AK8_puppiEta)[t2]; 
			mass_chs = (*AK8_JetPrunedMassCorr)[t2];
			mass_puppi = (*AK8_puppiSDMassL2L3Corr)[t2]; 
		}

		if( bjet%2 ==10){
			mass_chs = (*AK8_JetPrunedMass)[t2];
			mass_puppi = (*AK8_puppiSDMass)[t2]; 		
			//tau  = (*AK8_JetPt)[t2]; 
			//tau1 = (*AK8_JetEta)[t2]; 
		}

		tau  = (*AK8_puppiPt)[t2]; 
		tau1 = (*AK8_puppiEta)[t2]; 
		h_tag_Pt.	at(k)    -> Fill(tau) ;
		h_tag_Eta.	at(k)  	 -> Fill(tau1) ;

		h_tag_Pruned.	at(k)  	 ->Fill(mass_chs);
		h_tag_SD.	at(k) 	 ->Fill(mass_puppi);

		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_Puppitau21.at(k) 	 -> Fill(tau) ;
		tau =  ((*AK8_puppiTau3)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau32.at(k)   -> Fill(tau) ;
		tau =  1.0 - ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau42.at(k)   -> Fill(tau) ;

		tau =  ((*AK8_Jet_tau2)[t2]/(*AK8_Jet_tau1)[t2]) ;
		h_tag_tau21.	at(k)  	 -> Fill(tau) ;
		tau =  ((*AK8_Jet_tau3)[t2]/(*AK8_Jet_tau2)[t2]) ;
		h_tag_tau32.	at(k)    -> Fill(tau) ;

		tau =  ((*AK8_Jet_CHStau4)[t2]/(*AK8_Jet_CHStau2)[t2]) ;
		h_tag_CHStau42.	at(k)    -> Fill(tau) ;

		mass_puppi = (*AK8_puppiMass)[t2];
		h_tag_PUPPImass.at(k)   ->Fill(mass_puppi);  

		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.     at(k)    ->Fill(mass_puppi);      

		//if(k == 0 ) h_AK8_PUPPIvsmass.at(1) -> Fill((*mc_MomMass)[jet], mass_puppi) ;
		//		if(k == 1) h_AK8_CHSvsmass.at(1) -> Fill((*mc_MomMass)[jet], mass_puppi) ;


		// for muon dR plots
		if ( n_Mu.size() != 0) {
			j = n_Mu[0] ;
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dR_Recomu_tagjet.at(k) ->Fill(tau1);
			h_dPt_lep_tagjet.at(k) ->Fill(tau);
		}

		// for bjet dR plots
		if ( bjet != -1 ) {
			tau = dR_AK8bjet( t2, bjet, k, 0, fill);
			h_dR_Recob1_tagjet.at(k) -> Fill(tau) ;
		}
		else{
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k) -> Fill(tau) ;
			}
		}
	}

	// for dR plots in Wjets
	for ( int h = 0 ; h < W_boson.size() ; h ++ ){
		t2 = W_boson[h];
		for ( int g = h+1; g < W_boson.size(); g++) {
			t3 = W_boson[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h).at(g) ->Fill(tau) ;
		}  
	}
}

void  TTbar_PUPPItau_loc_190204::Topjet_Plots(int var ) 
{
	int j = -1 ;
	int fill = -1 ;
	int Pt_idx = 2; //for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float max = 0.0 ;
	float min = 0.0 ;
	TLorentzVector  Ws ;
	int Tsize        =   topjet.size()  ;
	float mass_puppi = 0.0;
	float mass_chs = 0.0 ;

	if ( Tsize != 0 )  h_tag_N.at(1)     ->Fill(Tsize) ;
	for (int k = 0 ; k < Tsize ; k ++ ){
		t2		=  	topjet[k] ;
		mass_puppi 	=	(*AK8_puppiSDMassL2L3Corr)[t2] ;	
		if( k == 0 ) {
		min 		= 	fabs(mass_puppi - 172.5) ;
		j 		=  	t2 ;
		}
	        max 	     	=  	fabs(mass_puppi - 172.5) ;		
		if( min > max) {
		min = max ;
		j   = t2 ;
		}
	//cout << "\nValue["<< k << "] = " << j ;
	}

	if ( j != -1 && Tsize != 0 ) topjet[0] = j ;

	Tsize  = ( topjet.size() >= 1) ? 1 : topjet.size() ;
	for (int k = 0 ; k < Tsize ; k ++ ) {

		t2   =   topjet[k] ;
		fill =   -1 ;	
		if (var == 1)   {		
		tau =  ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		if (tau >= 0.00) fill = 1 ;
		}

		if (var == 2 ) {		
		tau =  ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		if (tau < 0.00) fill = 1 ;
		}
		
		//if ( fill == -1 ) continue ;
		if( vr   <= 6){
			//tau  = (*AK8_puppiPt)[t2]; 
			//tau1 = (*AK8_puppiEta)[t2]; 
			mass_chs = (*AK8_JetPrunedMassCorr)[t2];
			mass_puppi = (*AK8_puppiSDMassL2L3Corr)[t2];
		}

		if( vr %2 ==10){
			mass_chs = (*AK8_JetPrunedMass)[t2];
			mass_puppi = (*AK8_puppiSDMass)[t2];
			//tau  = (*AK8_JetPt)[t2]; 
			//tau1 = (*AK8_JetEta)[t2]; 
		}

		tau  = (*AK8_puppiPt)[t2];
		tau1 = (*AK8_puppiEta)[t2];
		h_tag_Pt.       at(k+2)    -> Fill(tau) ;
		h_tag_Eta.      at(k+2)          -> Fill(tau1) ;

		h_tag_Pruned.   at(k+2)          ->Fill(mass_chs);
		h_tag_SD.       at(k+2)          ->Fill(mass_puppi);

		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_Puppitau21.at(k+2)         -> Fill(tau) ;
		tau =  ((*AK8_puppiTau3)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau32.at(k+2)   -> Fill(tau) ;
		tau =  ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau42.at(k+2)   -> Fill(tau) ;

		tau =  ((*AK8_Jet_tau2)[t2]/(*AK8_Jet_tau1)[t2]) ;
		h_tag_tau21.    at(k+2)          -> Fill(tau) ;
		tau =  ((*AK8_Jet_tau3)[t2]/(*AK8_Jet_tau2)[t2]) ;
		h_tag_tau32.    at(k+2)    -> Fill(tau) ;

		tau =  ((*AK8_Jet_CHStau4)[t2]/(*AK8_Jet_CHStau2)[t2]) ;
		h_tag_CHStau42. at(k+2)    -> Fill(tau) ;

		mass_puppi = (*AK8_puppiMass)[t2];
		h_tag_PUPPImass.at(k+2) ->Fill(mass_puppi);

		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.     at(k+2)    ->Fill(mass_puppi);


		if(n_Mu.size() != 0) {
			j = n_Mu[0] ;
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+2) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+2) ->Fill(tau1);
		}     

		// for bjet dR plots
		for(int l = 0 ; l < b_jet.size() ; l ++){
			j = b_jet[l] ;
			tau = dR_AK8bjet( t2, j, k, l, fill);
			if( l == 0) h_dR_Recob1_tagjet.at(k+2) -> Fill(tau) ;
			if( l == 1) h_dR_Recob2_tagjet.at(k+2) -> Fill(tau) ;
		}    
	}  

	// for dR plots in Wjets
	for ( int h = 0 ; h < Tsize ; h ++ ){
		t2 = topjet[h];
		for ( int g = h+1; g < Tsize; g++) {
			t3 = topjet[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
				h_dR_tagjet.at(h+2).at(g+2) ->Fill(tau) ;
		}  
	}
}


void  TTbar_PUPPItau_loc_190204::Higgsjet_Plots( int bjet) 
{
	int j = -1 ;
	int fill = -1 ;
	int Pt_idx = 3; //for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;     
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	float max = 0.0 ;
	float min = 0.0 ;
	TLorentzVector  Ws ;
	int Hsize        =   Higgsjets.size()  ;
	if ( Hsize != 0 )  h_tag_N.at(2)     ->Fill(Hsize) ;
	for (int k = 0 ; k < Hsize ; k ++ ){
		t2		=  	Higgsjets[k] ;
		mass_puppi 	=	(*AK8_puppiSDMassL2L3Corr)[t2] ;	
		if( k == 0 ) {
		min 		= 	fabs(mass_puppi - 125.0) ;
		j 		= 	t2 ;
		}

	        max 	     	=  	fabs(mass_puppi - 125.0) ;		
		if( min > max) {
		min = max ;
		j   = t2 ;
		}
	}

	if ( j != -1 && Hsize != 0 ) Higgsjets[0] = j ;
	Hsize  = ( Higgsjets.size() >= 1) ? 1 : Higgsjets.size() ;
	for (int k = 0 ; k < Hsize ; k ++ ) {

		t2   =   Higgsjets[k] ;
		h_tag_Pt.at(k+3)  -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k+3) -> Fill((*AK8_JetEta)[t2]) ; 

		h_tag_Pruned.at(k+3)  	 ->Fill((*AK8_JetPrunedMassCorr)[t2]);
		h_tag_SD.at(k+3) 	 ->Fill((*AK8_puppiSDMassL2L3Corr)[t2]);
		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_Puppitau21.at(k+3) -> Fill(tau) ;
		tau =  ((*AK8_puppiTau3)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau32.at(k+3) -> Fill(tau) ;
		tau = ((*AK8_puppiTau4)[t2]/(*AK8_puppiTau2)[t2]) ;
		h_tag_Puppitau42.at(k+3) -> Fill(tau) ;

		tau =  ((*AK8_Jet_tau2)[t2]/(*AK8_Jet_tau1)[t2]) ;
		h_tag_tau21.at(k+3)  -> Fill(tau) ;
		tau =  ((*AK8_Jet_tau3)[t2]/(*AK8_Jet_tau2)[t2]) ;
		h_tag_tau32.at(k+3)  -> Fill(tau) ;

		tau =  ((*AK8_Jet_CHStau4)[t2]/(*AK8_Jet_CHStau2)[t2]) ;
		//	if (tau < 0.005 ) cout << "\n Tau42 = " << tau ;
		h_tag_CHStau42.at(k+3)  -> Fill(tau) ;

		mass_puppi = (*AK8_puppiMass)[t2];
		h_tag_PUPPImass.at(k+3) ->Fill(mass_puppi);

		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k+3) ->Fill(mass_puppi);

		if ( n_Mu.size() != 0) {
			j = n_Mu[0] ;
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+3) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+3) ->Fill(tau1);

		}   

		if ( bjet != -1 ) {
			// for bjet dR plots
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k+3) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k+3) -> Fill(tau) ;
			}
		} 
	}

	// for dR plots in Wjets
	for ( int h = 0 ; h < Hsize ; h ++ ){
		t2 = Higgsjets[h];
		for ( int g = h+1; g < Hsize; g++) {
			t3 = Higgsjets[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h+3).at(g+3) ->Fill(tau) ;
		}  
	}
}


void  TTbar_PUPPItau_loc_190204::Fatjet_Plots(int topbjet ) 
{
	int j =  -1 ;
	int fill = -1 ;
	int Pt_idx = 4; //for dPt w.r.t muon plots     
	int t2 = -1 ;
	int t3 = -1 ;     
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	TLorentzVector  Ws ;
	int Fsize        =   fat_jet.size()  ;
	if ( Fsize != 0 )  h_tag_N.at(3)     ->Fill(Fsize) ;
	Fsize  = ( fat_jet.size() >= 2) ? 2 : fat_jet.size() ;
	for (int k = 0 ; k < Fsize ; k ++ ) {

		t2   =   fat_jet[k] ;
		h_tag_Pt.at(k+4)  -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k+4) -> Fill((*AK8_JetEta)[t2]) ; 
		h_tag_SD.at(k+4)  ->Fill((*AK8_puppiSDMass)[t2]);
		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_Puppitau21.at(k+4)  -> Fill(tau) ;

		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k+4) ->Fill(mass_puppi);

		if ( n_Mu.size() != 0) {
			j = n_Mu[0] ;
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+4) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+4) ->Fill(tau1);
		}     

		// for bjet dR plots
		if ( topbjet != -1 ) {
			tau = dR_AK8bjet( t2, topbjet, k, 0, fill);
			h_dR_Recob1_tagjet.at(k+4) -> Fill(tau) ;
		}
		else {
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k+4) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k+4) -> Fill(tau) ;
			}
		} 
	}

	// for dR plots in fatjet
	for ( int h = 0 ; h < Fsize ; h ++ ){
		t2 = fat_jet[h];
		for ( int g = h+1; g < Fsize; g++) {
			t3 = fat_jet[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h+4).at(g+4) ->Fill(tau) ;
		}  
	}
}


void  TTbar_PUPPItau_loc_190204::Category_Wjet_Plot( int cat, int W, int idx) 
{    
	// II is for whether, W is subleading Wboson belonging to category IV
	int fill = 0 ;
	float  tau = 0.0, pt = 0.0, eta = 0.0 ;
	TLorentzVector  Ws ;
	int Wsize        =   W_boson.size()  ;     
	float mass_puppi = 0.0;
		if( vr <= 6 ) {
			mass_puppi =  (*AK8_puppiMass)[W]; 
			tau 	   =  ((*AK8_puppiTau2)[W]/(*AK8_puppiTau1)[W]) ;
			pt  	   =  (*AK8_puppiPt)[W]; 
			eta 	   =  (*AK8_puppiEta)[W]; 
			Ws.SetPtEtaPhiM((*AK8_puppiPt)[W],(*AK8_puppiEta)[W],(*AK8_puppiPhi)[W],(*AK8_puppiMass)[W]);

		}

		if( vr%2 == 10){
			mass_puppi =  (*AK8_JetMass)[W]; 		
			tau 	   =  ((*AK8_puppiTau2)[W]/(*AK8_puppiTau1)[W]) ;
			pt         =  (*AK8_JetPt)[W] ; 	
			eta	   =  (*AK8_JetEta)[W] ;  
			Ws.SetPtEtaPhiE((*AK8_JetPt)[W],(*AK8_JetEta)[W],(*AK8_JetPhi)[W],(*AK8_JetEn)[W]);
		}

	h_Histo_Pt.at(cat).at(idx)  -> Fill(pt) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill(eta) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill(mass_puppi);
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;
	float   Et_W        	  =  Ws.Et() ;
	mass_puppi  		  =  Sqrt (fabs ((Et_W * Et_W) - (pt * pt) ) );
        //mass_puppi  		  =  Ws.M() ;
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);
}


void  TTbar_PUPPItau_loc_190204::Category_Top_Plot( int cat,int top, int idx) 
{    
	int fill = -1 ;
	float  tau = 0.0, pt = 0.0, eta = 0.0 ;     
	TLorentzVector  Ws ;
	int Tsize        =   topjet.size()  ;
	float mass_puppi = 0.0;
		if( vr <= 6 ) {
			mass_puppi =  (*AK8_puppiMass)[top]; 
			tau 	   =  ((*AK8_puppiTau3)[top]/(*AK8_puppiTau2)[top]) ;
			pt  	   =  (*AK8_puppiPt)[top]; 
			eta 	   =  (*AK8_puppiEta)[top]; 
			Ws.SetPtEtaPhiM((*AK8_puppiPt)[top],(*AK8_puppiEta)[top],(*AK8_puppiPhi)[top],(*AK8_puppiMass)[top]);

		}

		if( vr%2 == 10){
			mass_puppi =  (*AK8_JetMass)[top]; 		
			tau 	   =  ((*AK8_puppiTau3)[top]/(*AK8_puppiTau2)[top]) ;
			pt         =  (*AK8_JetPt)[top] ; 	
			eta	   =  (*AK8_JetEta)[top] ;  
			Ws.SetPtEtaPhiE((*AK8_JetPt)[top],(*AK8_JetEta)[top],(*AK8_JetPhi)[top],(*AK8_JetEn)[top]);
		}

	h_Histo_Pt.at(cat).at(idx)  -> Fill(pt) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill(eta) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill(mass_puppi);
	  
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;

	float   Et_top        =  Ws.Et() ;
	mass_puppi  = Sqrt (fabs ((Et_top * Et_top) - (pt * pt) ) );    
        //mass_puppi  		  =  Ws.M() ;
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);

}


void  TTbar_PUPPItau_loc_190204::Category_Higgs_Plot( int cat, int higg, int idx) 
{      
	float  tau = 0.0, pt = 0.0, eta = 0.0 ;
	TLorentzVector  Ws ;
	int Hsize        =   Higgsjets.size()  ;
	float mass_puppi = 0.0;
		if( vr <= 6 ) {
			mass_puppi =  (*AK8_puppiMass)[higg]; 
			tau 	   =  ((*AK8_puppiTau4)[higg]/(*AK8_puppiTau2)[higg]) ;
			pt  	   =  (*AK8_puppiPt)[higg]; 
			eta 	   =  (*AK8_puppiEta)[higg]; 
			Ws.SetPtEtaPhiM((*AK8_puppiPt)[higg],(*AK8_puppiEta)[higg],(*AK8_puppiPhi)[higg],(*AK8_puppiMass)[higg]);

		}

		if( vr%2 == 10){
			mass_puppi =  (*AK8_JetMass)[higg]; 		
			tau 	   =  ((*AK8_Jet_CHStau4)[higg]/(*AK8_Jet_CHStau2)[higg]) ;
			pt         =  (*AK8_JetPt)[higg] ; 	
			eta	   =  (*AK8_JetEta)[higg] ;  
			Ws.SetPtEtaPhiE((*AK8_JetPt)[higg],(*AK8_JetEta)[higg],(*AK8_JetPhi)[higg],(*AK8_JetEn)[higg]);
		}

	h_Histo_Pt.at(cat).at(idx)  -> 	Fill(pt) ;
	h_Histo_Eta.at(cat).at(idx) -> 	Fill(eta) ; 
	h_Histo_SD.at(cat).at(idx)  -> 	Fill(mass_puppi) ;               
	h_Histo_tau.at(cat).at(idx) ->  Fill(tau) ;

	float   Et_top        	  =   Ws.Et() ;	
	mass_puppi  		  =   Sqrt (fabs ((Et_top * Et_top) - (pt * pt) ) );    
        //mass_puppi  		  =   Ws.M() ;
	h_Histo_Mass.at(cat).at(idx)->  Fill(mass_puppi);      
}        


void  TTbar_PUPPItau_loc_190204::Category_Fat_Plot( int cat,int fat, int idx) 
{    
	int fill = -1 ;
	float  tau = 0.0 ;
	TLorentzVector  Ws ;
	int Tsize        =   fat_jet.size()  ;

	h_Histo_Pt.at(cat).at(idx)  -> Fill((*AK8_JetPt)[fat]) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill((*AK8_JetEta)[fat]) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill((*AK8_puppiSDMass)[fat]);

	tau =  ((*AK8_puppiTau2)[fat]/(*AK8_puppiTau1)[fat]) ;
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;
	Ws.SetPtEtaPhiE((*AK8_JetPt)[fat],(*AK8_JetEta)[fat],(*AK8_JetPhi)[fat],(*AK8_JetEn)[fat]);
	float   Et_top        =  Ws.Et() ;
	float   pT_top  =  (*AK8_JetPt)[fat] ;
	double mass_puppi  = Sqrt (fabs ((Et_top * Et_top) - (pT_top * pT_top) ) );    
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);

}


void TTbar_PUPPItau_loc_190204::TagJets_dRPlots()
{
	float dR = 0.0;
	int b1 = -1, b2 =-1 ;
//	int Fsize  = ( fat_jet.size() >= 4) ? 4 : fat_jet.size() ;
	int Wsize  = ( W_boson.size() >= 2) ? 2 : W_boson.size() ;
	int Tsize  = ( topjet.size() >= 1) ? 1 : topjet.size() ;
	int Hsize  = ( Higgsjets.size() >= 1 ) ? 1 : Higgsjets.size() ;
	// for Wjet w.r.t tag jets dR plots
	for ( int h = 0 ; h < Wsize ; h ++ ){
		b1 = W_boson[h];  
		for( int g = 0 ; g < Tsize ; g ++ ) {
			b2 = topjet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+2) ->Fill(dR) ;
		}
		for( int g = 0 ; g < Hsize ; g ++ ) {
			b2 = Higgsjets[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+3) ->Fill(dR) ;
		}
	/*	for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+8) ->Fill(dR) ;
		}  */
	} 
	// for topjet w.r.t tag jets dR plots
	for( int h = 0 ; h < Tsize ; h ++ ) {
		b1 = topjet[h]; 
		for( int g = 0 ; g < Hsize ; g ++ ) {
			b2 = Higgsjets[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+2).at(g+3) ->Fill(dR) ;
		}
	/*	for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+4).at(g+8) ->Fill(dR) ;
		}  */
	}
	// for Higgsjet w.r.t tag jets dR plots  
	for( int h = 0 ; h < Hsize ; h ++ ) {
		b1 = Higgsjets[h];  
		/*for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+6).at(g+8) ->Fill(dR) ;
		}  */  
	}

}

// ============Mt Calculation ========================

void  TTbar_PUPPItau_loc_190204::WfatIII_MTCalculation( ) 
{
	int top          =   CatIII_Objects[0] ;  	
	int higgs2      =   CatIII_Objects[1] ;  
	int topbjet     =   CatIII_Objects[2] ;  	
	int mu         =   CatIII_Objects[3] ;  	

	TLorentzVector v_higgs2, v_mu, v_topbjet, v_top ;
	float MET_phi[4] ;
	float Obj_phi[6] ;

	v_top.SetPtEtaPhiE((*AK8_JetPt)[top], (*AK8_JetEta)[top], (*AK8_JetPhi)[top], (*AK8_JetEn)[top] ) ;
	v_higgs2.SetPtEtaPhiE((*AK8_JetPt)[higgs2], (*AK8_JetEta)[higgs2], (*AK8_JetPhi)[higgs2], (*AK8_JetEn)[higgs2] ) ;
	v_topbjet.SetPtEtaPhiE((*jet_Pt)[topbjet], (*jet_Eta)[topbjet], (*jet_Phi)[topbjet], (*jet_En)[topbjet] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_top        =  v_top.Et() ;	
	float Et_higgs2        =  v_higgs2.Et() ;
	float Et_bjet          =   v_topbjet.Et() ;
	float Et_mu          =    v_mu.Et() ;
	double Et_sum       =   0.0 ;

	float pT_top       =  (*AK8_JetPt)[top] ;
	float pT_higgs2       =  (*AK8_JetPt)[higgs2] ;
	float pT_bjet         =  (*jet_Pt)[topbjet] ;
	float pT_mu         =  (*mu_Pt)[mu] ;
	double pT_sum      =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  top_Trans            = 0.0 ;  
	double  higgs2_Trans            = 0.0 ;
	double  TprimEt_Trans         = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  

	//============ w.r.t MET ======
	float  dPhi            =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	MET_phi[0]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet] , pf_METPhi)  ;   //  for topbjet & MET dphi
	MET_phi[1]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs2] , pf_METPhi)  ;   //  for higgs2 & MET dphi
	MET_phi[2]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[top] , pf_METPhi)  ;   //  for top & MET dphi
	MET_phi[3]          =  pf_MET * Cos(dPhi) ;


	// ======= w.r.t each other ========
	dPhi                   =  delta_phi((*AK8_JetPhi)[top], (*mu_Phi)[mu]) ;  ;   //  for top & muon dphi
	Obj_phi[0]            =   pT_mu * Cos(dPhi) ;  

	dPhi                   =  delta_phi((*AK8_JetPhi)[top], (*AK8_JetPhi)[higgs2])  ; //  for top & higgs dphi
	Obj_phi[1]            =   pT_higgs2 * Cos(dPhi) ;

	dPhi                   =  delta_phi( (*AK8_JetPhi)[top], (*jet_Phi)[topbjet])  ; //  for W & higgs dphi
	Obj_phi[2]            =   pT_bjet * Cos(dPhi) ;	

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*mu_Phi)[mu])  ;   //  for topbjet & muon dphi
	Obj_phi[3]            =  pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*AK8_JetPhi)[higgs2])  ; //  for topbjet & higgs2 dphi
	Obj_phi[4]            =  pT_higgs2 * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs2], (*mu_Phi)[mu]) ; //  for higgs & muon dphi
	Obj_phi[5]            =  pT_mu * Cos(dPhi) ;


	float  Dot_prod = 0.0   ;

	// ---------------Transverse mass---------------
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi[0] ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_Histo_Mt.at(2).at(3) -> Fill( trans_mass) ;

	//  for higgs2 transverse mass 
	Et_sum               =      Et_sum +   Et_higgs2 ;
	pT_sum              =       pT_sum + ( pT_higgs2* pT_higgs2 ) ;
	Dot_prod             =      Dot_prod + 2.0 * pT_higgs2 *( MET_phi[2] + Obj_phi[5] ) ;
	higgs2_Trans         =      fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass           =      Sqrt ( higgs2_Trans) ;
	h_Histo_Mt.at(2).at(2) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum              =    Et_sum + Et_top + Et_bjet ;
	pT_sum             =     pT_sum + (pT_top * pT_top) + (pT_bjet * pT_bjet);    
	Dot_prod            =    Dot_prod +2.0 *pT_top* (MET_phi[3]+ Obj_phi[0]+ Obj_phi[1]+ Obj_phi[2]) +2.0 *pT_bjet* (MET_phi[1]+ Obj_phi[3]+ Obj_phi[4]) ;
	TprimEt_Trans     =    fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         =     Sqrt ( TprimEt_Trans ) ;
	//	cout << "\nTop mass = " << trans_mass;
	h_Histo_Mt.at(2).at(0) -> Fill( trans_mass) ;

	// for Top transverse mass
	Et_sum                        =     0.0 + Et_top + Et_bjet ;
	pT_sum                       =     0.0 + (pT_top * pT_top) + (pT_bjet * pT_bjet);            
	Dot_prod                     =     0.0  +  2.0 *pT_top* Obj_phi[2] ;
	top_Trans                   =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass                 =    Sqrt ( top_Trans) ;
	h_Histo_Mt.at(2).at(1)   -> Fill( trans_mass) ;

}


void  TTbar_PUPPItau_loc_190204::WfatV_MTCalculation( ) 
{
	int higgs1      =  CatV_Objects[0] ;  	
	int higgs2      =  CatV_Objects[1] ;  
	int topbjet     =  CatV_Objects[2] ;  	
	int mu         =     CatV_Objects[3] ;  	
	TLorentzVector v_higgs2, v_mu, v_topbjet, v_higgs1 ;
	float MET_phi[4] ;
	float Obj_phi[6] ;

	v_higgs1.SetPtEtaPhiE((*AK8_JetPt)[higgs1], (*AK8_JetEta)[higgs1], (*AK8_JetPhi)[higgs1], (*AK8_JetEn)[higgs1] ) ;
	v_higgs2.SetPtEtaPhiE((*AK8_JetPt)[higgs2], (*AK8_JetEta)[higgs2], (*AK8_JetPhi)[higgs2], (*AK8_JetEn)[higgs2] ) ;
	v_topbjet.SetPtEtaPhiE((*jet_Pt)[topbjet], (*jet_Eta)[topbjet], (*jet_Phi)[topbjet], (*jet_En)[topbjet] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_higgs1        =  v_higgs1.Et() ;	
	float Et_higgs2        =  v_higgs2.Et() ;
	float Et_bjet          =   v_topbjet.Et() ;
	float Et_mu          =    v_mu.Et() ;
	double Et_sum       =   0.0 ;

	float pT_higgs1       =  (*AK8_JetPt)[higgs1] ;
	float pT_higgs2       =  (*AK8_JetPt)[higgs2] ;
	float pT_bjet         =  (*jet_Pt)[topbjet] ;
	float pT_mu         =  (*mu_Pt)[mu] ;
	double pT_sum      =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  higgs1_Trans            = 0.0 ;  
	double  higgs2_Trans            = 0.0 ;
	double  TprimEt_Trans         = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  

	//============ w.r.t MET ======
	float  dPhi            =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	MET_phi[0]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet] , pf_METPhi)  ;   //  for topbjet & MET dphi
	MET_phi[1]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs2] , pf_METPhi)  ;   //  for higgs2 & MET dphi
	MET_phi[2]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs1] , pf_METPhi)  ;   //  for top & MET dphi
	MET_phi[3]          =  pf_MET * Cos(dPhi) ;


	// ======= w.r.t each other ========
	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs1], (*mu_Phi)[mu]) ;  ;   //  for top & muon dphi
	Obj_phi[0]            =   pT_mu * Cos(dPhi) ;  

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs1], (*AK8_JetPhi)[higgs2])  ; //  for top & higgs dphi
	Obj_phi[1]            =   pT_higgs2 * Cos(dPhi) ;

	dPhi                   =  delta_phi( (*AK8_JetPhi)[higgs1], (*jet_Phi)[topbjet])  ; //  for W & higgs dphi
	Obj_phi[2]            =   pT_bjet * Cos(dPhi) ;	

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*mu_Phi)[mu])  ;   //  for topbjet & muon dphi
	Obj_phi[3]            =  pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs2], (*jet_Phi)[topbjet])  ; //  for topbjet & higgs2 dphi
	Obj_phi[4]            =  pT_bjet * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs2], (*mu_Phi)[mu]) ; //  for higgs & muon dphi
	Obj_phi[5]            =  pT_mu * Cos(dPhi) ;


	float  Dot_prod = 0.0   ;

	// ---------------Transverse mass---------------
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi[0] ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_Histo_Mt.at(4).at(3) -> Fill( trans_mass) ;

	//  for Top transverse mass 
	Et_sum               =      Et_sum +   Et_bjet ;
	pT_sum              =       pT_sum + ( pT_bjet* pT_bjet ) ;
	Dot_prod             =      Dot_prod + 2.0 * pT_bjet *( MET_phi[1] + Obj_phi[3] ) ;
	higgs1_Trans         =      fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass           =      Sqrt ( higgs1_Trans) ;
	h_Histo_Mt.at(4).at(1) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum              =    Et_sum + Et_higgs1 + Et_higgs2 ;
	pT_sum             =     pT_sum + (pT_higgs1 * pT_higgs1) + (pT_higgs2 * pT_higgs2);    
	Dot_prod            =  Dot_prod +2.0*pT_higgs1*(MET_phi[3]+ Obj_phi[0]+ Obj_phi[1]+ Obj_phi[2]) +2.0*pT_higgs2*(MET_phi[2]+ Obj_phi[4]+ Obj_phi[5]) ;
	TprimEt_Trans     =    fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         =     Sqrt ( TprimEt_Trans ) ;
	//	cout << "\nTop mass = " << trans_mass;
	h_Histo_Mt.at(4).at(0) -> Fill( trans_mass) ;

	// for Higgs transverse mass
	Et_sum                        =     0.0 + Et_higgs1 + Et_higgs2 ;
	pT_sum                       =     0.0 + (pT_higgs1 * pT_higgs1) + (pT_higgs2 * pT_higgs2);            
	Dot_prod                     =     0.0  +  2.0 *pT_higgs1* Obj_phi[1] ;
	higgs2_Trans                   =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass                 =    Sqrt ( higgs2_Trans) ;
	h_Histo_Mt.at(4).at(2)   -> Fill( trans_mass) ;

}


void  TTbar_PUPPItau_loc_190204::WWtag_MTCalculation() 
{
	int top = CatIV_Objects[0] ;  	
	int higgs = CatIV_Objects[1] ;  
	int topbjet = CatIV_Objects[2] ;  	
	int mu = CatIV_Objects[3] ;  	
	TLorentzVector v_higgs, v_mu, v_topbjet, v_top ;
	float MET_phi[4] ;
	float Obj_phi[6] ;
	float pT_top = 0.0, pT_higgs = 0.0 ;

	if( vr < 6 ) {
 	pT_top 	       =  (*AK8_puppiPt)[top] ;
	pT_higgs      =  (*AK8_puppiPt)[higgs] ;

	v_top.SetPtEtaPhiM((*AK8_puppiPt)[top],(*AK8_puppiEta)[top],(*AK8_puppiPhi)[top],(*AK8_puppiMass)[top]);
	v_higgs.SetPtEtaPhiM((*AK8_puppiPt)[higgs],(*AK8_puppiEta)[higgs],(*AK8_puppiPhi)[higgs],(*AK8_puppiMass)[higgs]);
	}


	v_topbjet.SetPtEtaPhiE((*jet_Pt)[topbjet], (*jet_Eta)[topbjet], (*jet_Phi)[topbjet], (*jet_En)[topbjet] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_top          =  v_top.Et() ;	
	float Et_higgs        =  v_higgs.Et() ;
	float Et_bjet          = v_topbjet.Et() ;
	float Et_mu          = v_mu.Et() ;
	double Et_sum       = 0.0 ;

	float pT_bjet         =  (*jet_Pt)[topbjet] ;
	float pT_mu         =  (*mu_Pt)[mu] ;
	double pT_sum      =  0.0 ;
	float  Dot_prod     =  0.0   ;

	double  W_Trans                = 0.0 ;
	double  top_Trans                = 0.0 ;  
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans        = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  

	//============ w.r.t MET ======
	float  dPhi            =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	MET_phi[0]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet] , pf_METPhi)  ;   //  for topbjet & MET dphi
	MET_phi[1]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[higgs] , pf_METPhi)  ;   //  for higgs & MET dphi
	MET_phi[2]          =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[top] , pf_METPhi)  ;   //  for top & MET dphi
	MET_phi[3]          =  pf_MET * Cos(dPhi) ;


	// ======= w.r.t each other ========
	dPhi                   =  delta_phi((*AK8_puppiPhi)[top], (*mu_Phi)[mu]) ;  ;   //  for top & muon dphi
	Obj_phi[0]            =   pT_mu * Cos(dPhi) ;  

	dPhi                   =  delta_phi((*AK8_puppiPhi)[top], (*AK8_puppiPhi)[higgs])  ; //  for top & higgs dphi
	Obj_phi[1]            =   pT_higgs * Cos(dPhi) ;

	dPhi                   =  delta_phi( (*AK8_puppiPhi)[top], (*jet_Phi)[topbjet])  ; //  for W & higgs dphi
	Obj_phi[2]            =   pT_bjet * Cos(dPhi) ;	

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*mu_Phi)[mu])  ;   //  for topbjet & muon dphi
	Obj_phi[3]            =  pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*AK8_puppiPhi)[higgs])  ; //  for topbjet & higgs dphi
	Obj_phi[4]            =  pT_higgs * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[higgs], (*mu_Phi)[mu]) ; //  for higgs & muon dphi
	Obj_phi[5]            =  pT_mu * Cos(dPhi) ;


	// ---------------Transverse mass---------------
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi[0] ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_Histo_Mt.at(3).at(3) -> Fill( trans_mass) ;

	//  for Higgs transverse mass 
	Et_sum               =     Et_sum +   Et_higgs ;
	pT_sum              =     pT_sum + ( pT_higgs* pT_higgs ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_higgs *( MET_phi[2] + Obj_phi[5] ) ;
	Higgs_Trans           =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( Higgs_Trans) ;
	h_Histo_Mt.at(3).at(2) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum              =    Et_sum + Et_top + Et_bjet ;
	pT_sum             =     pT_sum + (pT_top * pT_top) + (pT_bjet * pT_bjet);    
	Dot_prod            =    Dot_prod +2.0 *pT_top* (MET_phi[3]+ Obj_phi[0]+ Obj_phi[1]+ Obj_phi[2]) +2.0 *pT_bjet* (MET_phi[1]+ Obj_phi[3]+ Obj_phi[4]) ;
	TprimEt_Trans     =    fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         =     Sqrt ( TprimEt_Trans ) ;
	h_Histo_Mt.at(3).at(0) -> Fill( trans_mass) ;

	// for Top transverse mass
	Et_sum                        =     0.0 + Et_top + Et_bjet ;
	pT_sum                       =     0.0 + (pT_top * pT_top) + (pT_bjet * pT_bjet);            
	Dot_prod                     =     0.0  +  2.0 *pT_top* Obj_phi[2] ;
	top_Trans                   =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass                 =    Sqrt ( top_Trans) ;
	h_Histo_Mt.at(3).at(1)   -> Fill( trans_mass) ;
}


void  TTbar_PUPPItau_loc_190204::toptag_MTCalculation( int Cat ) 
{
	int top = -1, W = -1, mu  = -1 ;
	if ( Cat == 1) {
		top = CatI_Objects[0];
		W  = CatI_Objects[1];          
		mu = CatI_Objects[2];          
	}
	else{
		top = CatII_Objects[0];
		W  = CatII_Objects[1];          
		mu = CatII_Objects[2];                
	}      


	TLorentzVector v_top, v_mu, v_W ;
	float pT_top = 0.0, pT_W = 0.0 ;

	if( vr < 6 ) {
 	pT_top 	       =  (*AK8_puppiPt)[top] ;
	pT_W  	       =  (*AK8_puppiPt)[W] ;

	v_top.SetPtEtaPhiM((*AK8_puppiPt)[top],(*AK8_puppiEta)[top],(*AK8_puppiPhi)[top],(*AK8_puppiMass)[top]);
	v_W.SetPtEtaPhiM((*AK8_puppiPt)[W],(*AK8_puppiEta)[W],(*AK8_puppiPhi)[W],(*AK8_puppiMass)[W]);
	}

/*	if( vr%2 == 0 ) {
	pT_top 	       =  (*AK8_JetPt)[top] ;
	pT_W  	       =  (*AK8_JetPt)[W] ;
	v_top.SetPtEtaPhiE((*AK8_JetPt)[top], (*AK8_JetEta)[top], (*AK8_JetPhi)[top], (*AK8_JetEn)[top] ) ;
	v_W.SetPtEtaPhiE((*AK8_JetPt)[W], (*AK8_JetEta)[W], (*AK8_JetPhi)[W], (*AK8_JetEn)[W] ) ;

	}
*/
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_top        =  v_top.Et() ;
	float Et_W         = v_W.Et() ;
	float Et_mu        = v_mu.Et() ;
	double Et_sum    = 0.0 ;

	float pT_mu   =  (*mu_Pt)[mu] ;
	double pT_sum  = 0.0 ;

	double  W_Trans                = 0.0 ;
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans          = 0.0 ;
	double  trans_mass            = 0.0 ;

	// for cos(deltaPhi) calculation  
	float  dPhi           =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	float  MET_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[W] , pf_METPhi)  ;   //  for W & MET dphi
	float  MET2_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[top] , pf_METPhi)  ;   //  for top & MET dphi
	float  MET3_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[W] , (*mu_Phi)[mu])  ;   //  for W & muon dphi
	float  mu_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[top] , (*mu_Phi)[mu])  ;   //  for top & muon dphi
	float  mu2_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[W] , (*AK8_puppiPhi)[top])  ;   //  for W & top dphi
	float  W_phi       =   pT_W * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	if (Cat == 1 ) h_Histo_Mt.at(0).at(2) -> Fill( trans_mass) ;
	if (Cat == 2 ) h_Histo_Mt.at(1).at(2) -> Fill( trans_mass) ;

	// for Higgs Transverse mass
	Et_sum               =     Et_sum +   Et_W ;
	pT_sum              =     pT_sum + ( pT_W* pT_W ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_W *( MET2_phi + mu_phi ) ;
	Higgs_Trans        =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( Higgs_Trans) ;
	if (Cat == 1 ) h_Histo_Mt.at(0).at(1) -> Fill( trans_mass) ;
	if (Cat == 2 ) h_Histo_Mt.at(1).at(1) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum          = Et_sum + Et_top ;
	pT_sum        =    pT_sum + (pT_top * pT_top) ;    
	Dot_prod       =   Dot_prod + 2.0 * pT_top * ( MET3_phi + mu2_phi + W_phi ) ;
	TprimEt_Trans =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( TprimEt_Trans ) ;
	if (Cat == 1 ) h_Histo_Mt.at(0).at(0) -> Fill( trans_mass) ;
	if (Cat == 2 ) h_Histo_Mt.at(1).at(0) -> Fill( trans_mass) ;

}


// =========for Higgs tag Mt calculation =======================

void  TTbar_PUPPItau_loc_190204::Higgstag_MTCalculation() 
{
	int higgs    = -1 ;  
	int topbjet  = -1 ;  	
	int mu       = -1 ; 
	TLorentzVector v_higgs, v_mu, v_topbjet ;
	float pT_higgs = 0.0 ;

	
	
	higgs 		= CatVI_Objects[1] ;  
	topbjet		= CatVI_Objects[2] ;  	
	mu		= CatVI_Objects[3] ;  
	v_topbjet.SetPtEtaPhiE((*jet_Pt)[topbjet], (*jet_Eta)[topbjet], (*jet_Phi)[topbjet], (*jet_En)[topbjet] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);



	if( vr< 6 ) {
	pT_higgs  	   =  (*AK8_puppiPt)[higgs]; 
	v_higgs.SetPtEtaPhiM((*AK8_puppiPt)[higgs],(*AK8_puppiEta)[higgs],(*AK8_puppiPhi)[higgs],(*AK8_puppiMass)[higgs]);
	}

/*	if( vr == 0 ) {
	v_higgs.SetPtEtaPhiE((*AK8_JetPt)[higgs], (*AK8_JetEta)[higgs], (*AK8_JetPhi)[higgs], (*AK8_JetEn)[higgs] ) ;
	pT_higgs       =  (*AK8_JetPt)[higgs] ;
	} */


	float Et_higgs        =  v_higgs.Et() ;
	float Et_bjet          = v_topbjet.Et() ;
	float Et_mu          = v_mu.Et() ;
	double Et_sum       = 0.0 ;


	float pT_bjet         =  (*jet_Pt)[topbjet] ;
	float pT_mu         =  (*mu_Pt)[mu] ;
	double pT_sum      =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  top_Trans                = 0.0 ;  
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans        = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  
	float  dPhi           =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	float  MET_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet] , pf_METPhi)  ;   //  for topbjet & MET dphi
	float  MET2_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[higgs] , pf_METPhi)  ;   //  for higgs & MET dphi
	float  MET3_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*mu_Phi)[mu])  ;   //  for topbjet & muon dphi
	float  mu_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[higgs], (*mu_Phi)[mu]) ; //  for higgs & muon dphi
	float  mu2_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*AK8_puppiPhi)[higgs])  ; //  for W & higgs dphi
	float  bjet_phi       =   pT_bjet * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	 h_Histo_Mt.at(5).at(3) -> Fill( trans_mass) ;

	//  for Top transverse mass 
	Et_sum               =     Et_sum +   Et_bjet ;
	pT_sum              =     pT_sum + ( pT_bjet* pT_bjet ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_bjet *( MET2_phi + mu_phi ) ;
	top_Trans           =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( top_Trans) ;
	 h_Histo_Mt.at(5).at(2) -> Fill( trans_mass) ;
	
	// for Tprime Transverse mass
	Et_sum              = Et_sum + Et_higgs ;
	pT_sum             =    pT_sum + (pT_higgs * pT_higgs) ;    
	Dot_prod            =   Dot_prod + 2.0 * pT_higgs * ( MET3_phi + mu2_phi + bjet_phi ) ;
	TprimEt_Trans     =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         = Sqrt ( TprimEt_Trans ) ;

	 h_Histo_Mt.at(5).at(1) -> Fill( trans_mass) ;
}


void  TTbar_PUPPItau_loc_190204::genlvl_MTCalculation(int type ) 
{
	TLorentzVector mc_Tprime, mc_mu ;
	int ID = -1 ;
	float Et = 0.0 ;
	float Pt2 = 0.0 ;
	float MT = 0.0 ;

	float  dPhi      =  delta_phi((*mc_Phi)[T_higgs] , (*mc_Phi)[T_top]) ;   // for muon & MET dphi

	mc_Tprime.SetPtEtaPhiM((*mc_MomPt)[T_higgs],(*mc_MomEta)[T_higgs],(*mc_MomPhi)[T_higgs], (*mc_MomMass)[T_higgs] ) ;

	Et  =  (*mc_Et)[T_higgs] + (*mc_Et)[T_top] ;
	Pt2  = Power( (*mc_Pt)[T_higgs], 2)  + Power ((*mc_Pt)[T_top], 2 ) + 2.0 * (*mc_Pt)[T_higgs] * (*mc_Pt)[T_top] * Cos(dPhi) ;
	MT =  Sqrt (fabs( Power( Et , 2 ) -  Pt2  ) ) ;
	h_object_MT.at(8) -> Fill (MT );

	MT = Sqrt (fabs( Power( (*mc_Et)[T_higgs] , 2 ) - Power ( (*mc_Pt)[T_higgs] , 2 ) ) );
	h_object_MT.at(6) -> Fill ( MT );

	MT = Sqrt( fabs( Power( (*mc_Et)[T_top] , 2 ) - Power ( (*mc_Pt)[T_top] , 2 ) ) );
	h_object_MT.at(7) -> Fill ( MT );	

	if ( type  == -1 ) {
		if ( Higgs_Mu.size() != 0 ) ID = Higgs_Mu[0] ;
	}
	if ( type  == 1 ) {
		if ( top_Mu.size() != 0 )  ID = top_Mu[0] ;
	}
	if ( ID != -1 ) {  
		mc_mu.SetPtEtaPhiM((*mc_MomPt)[ID],(*mc_MomEta)[ID],(*mc_MomPhi)[ID], (*mc_MomMass)[ID] ) ;
		Et  =  mc_mu.Et();
		Pt2  = (*mc_MomPt)[ID] ;
		MT =  Sqrt (fabs( Power( Et , 2 ) - Power( Pt2 , 2 ) ) ) ;
		h_object_MT.at(5) -> Fill (MT );  
	}   
}

// =====Auxillary group functions =============
void  TTbar_PUPPItau_loc_190204::Group3_TopbH()
{
	int Hjet 	= Higgsjets[0];
	int Tjet 	= topjet[0];
	int bjet5	= b_jet[0] ; 
	TLorentzVector v_higgs, v_mu, v_topbjet ;
	float pT_higgs = 0.0 ;
	
	float dR_Tb 	= dR_AK8bjet( Tjet, bjet5, 0, 0, 0) ;
	float dR_Hb 	= dR_AK8bjet( Hjet, bjet5, 0, 0, 0) ;
	float dR_TH	= dR_AK8jet( Hjet, Tjet, 0, 0) ;

	if (  dR_TH > 2.0 ) {

	v_topbjet.SetPtEtaPhiM((*AK8_puppiPt)[Tjet], (*AK8_puppiEta)[Tjet], (*AK8_puppiPhi)[Tjet], (*AK8_puppiMass)[Tjet] ) ;
	//v_mu.SetPtEtaPhiE( (*jet_Pt)[bjet5], (*jet_Eta)[bjet5], (*jet_Phi)[bjet5], (*jet_En)[bjet5]);
	v_mu = v_topbjet ;
	
	pT_higgs  	   =  (*AK8_puppiPt)[Hjet]; 
	v_higgs.SetPtEtaPhiM((*AK8_puppiPt)[Hjet],(*AK8_puppiEta)[Hjet],(*AK8_puppiPhi)[Hjet],(*AK8_puppiMass)[Hjet]);

	float Et_Hjet        	=  v_higgs.Et() ;
	float Et_Tjet           = v_topbjet.Et() ;
	float Et_mu             = v_mu.Et() ;
	double Et_sum           = 0.0 ;


	float pT_Tjet           =  (*AK8_puppiPt)[Tjet] ;
	float pT_mu         	=  (*jet_Pt)[bjet5] ;
	double pT_sum           =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  top_Trans                = 0.0 ;  
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans        = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  

	float  dPhi                   =  delta_phi((*jet_Phi)[bjet5], (*AK8_puppiPhi)[Tjet])  ;   //  for topbjet & muon dphi
	float  mu_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[Hjet], (*jet_Phi)[bjet5]) ; //  for higgs & muon dphi
	float  mu2_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[Tjet],(*AK8_puppiPhi)[Hjet])  ; //  for W & higgs dphi
	float  bjet_phi       =   pT_Tjet * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for W transverse mass    
	W_Trans = fabs (( Et_Tjet * Et_Tjet) - (pT_Tjet * pT_Tjet) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_Histo_Mt.at(2).at(3) -> Fill( trans_mass) ;

	//  for Top transverse mass 
	top_Trans           =     fabs (( Et_Hjet * Et_Hjet)  - ( pT_higgs * pT_higgs) )  ;
	trans_mass          =    Sqrt ( top_Trans) ;
	 h_Histo_Mt.at(2).at(2) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum              = Et_sum + Et_Hjet + Et_Tjet ;
	pT_sum             =    pT_sum + (pT_higgs * pT_higgs) + (pT_Tjet * pT_Tjet) ;    
	Dot_prod            =   Dot_prod + 2.0 * pT_higgs * ( bjet_phi ) ;
	TprimEt_Trans     =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         = Sqrt ( TprimEt_Trans ) ;
	
	h_Histo_Mt.at(2).at(1) -> Fill( trans_mass) ;

}
}



void  TTbar_PUPPItau_loc_190204::Group5_WbH()
{
	int Hjet 	= Higgsjets[0];
	int Wjet 	= W_boson[0];
	int bjet5	= b_jet[0] ; 
	TLorentzVector v_higgs, v_bjet5, v_Wjet ;
	float pT_higgs = 0.0 ;
	
	float dR_Wb 	= dR_AK8bjet( Wjet, bjet5, 0, 0, 0) ;
	float dR_Hb 	= dR_AK8bjet( Hjet, bjet5, 0, 0, 0) ;
	float dR_WH	= dR_AK8jet( Hjet, Wjet, 0, 0) ;

	if ( dR_Wb < 2.0 &&  dR_Hb > 2.0 &&  dR_WH > 2.0 ) {

	v_Wjet.SetPtEtaPhiM((*AK8_puppiPt)[Wjet], (*AK8_puppiEta)[Wjet], (*AK8_puppiPhi)[Wjet], (*AK8_puppiMass)[Wjet] ) ;
	v_bjet5.SetPtEtaPhiE( (*jet_Pt)[bjet5], (*jet_Eta)[bjet5], (*jet_Phi)[bjet5], (*jet_En)[bjet5]);

	pT_higgs  	   =  (*AK8_puppiPt)[Hjet]; 
	v_higgs.SetPtEtaPhiM((*AK8_puppiPt)[Hjet],(*AK8_puppiEta)[Hjet],(*AK8_puppiPhi)[Hjet],(*AK8_puppiMass)[Hjet]);

	float  Et_Hjet        	  =  v_higgs.Et() ;
	float  Et_Wjet            =  v_Wjet.Et() ;
	float  Et_bjet            =  v_bjet5.Et() ;
	double Et_sum             =  0.0 ;
	double Mt_Wjet	 	  =  v_Wjet.Mt() ; 


	float pT_Wjet           =  v_Wjet.Pt() ;
	float pT_bjet         	=  (*jet_Pt)[bjet5] ;
	double pT_sum           =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  top_Trans                = 0.0 ;  
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans        = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  

	float  dPhi                   =  delta_phi((*jet_Phi)[bjet5], (*AK8_puppiPhi)[Wjet])  ;   //  for topbjet & muon dphi
	float  mu_phi       =   pT_bjet * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[Hjet], (*jet_Phi)[bjet5]) ; //  for higgs & muon dphi
	float  mu2_phi       =   pT_bjet * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_puppiPhi)[Wjet],(*AK8_puppiPhi)[Hjet])  ; //  for W & higgs dphi
	float  bjet_phi       =   pT_Wjet * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for Higgs transverse mass    
	W_Trans = fabs ((Et_Hjet * Et_Hjet) - ( pT_higgs * pT_higgs ) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
//	cout << " Mt(W) = " << trans_mass << " --  = " << Mt_Wjet << endl ;
//	cout << " Eta = " << (*AK8_puppiEta)[Wjet] << "  , Phi = " << (*AK8_puppiPhi)[Wjet] << "  , Mass = " << (*AK8_puppiMass)[Wjet] << "  , pt = " << pT_Wjet << " , Et = " << Et_Wjet <<endl ;
	h_Histo_Mt.at(4).at(3) -> Fill( trans_mass) ;

	//  for Top transverse mass 
	Et_sum               =     Et_sum +   Et_bjet + Et_Wjet ;
	pT_sum              =     pT_sum + ( pT_bjet* pT_bjet ) + ( pT_Wjet * pT_Wjet ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_bjet *(  mu_phi ) ;
	top_Trans           =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( top_Trans) ;
	
	 h_Histo_Mt.at(4).at(2) -> Fill( trans_mass) ;
	// for Tprime Transverse mass
	Et_sum              = Et_sum + Et_Hjet ;
	pT_sum             =    pT_sum + (pT_higgs * pT_higgs) ;    
	Dot_prod            =   Dot_prod + 2.0 * pT_higgs * (  mu2_phi + bjet_phi ) ;
	TprimEt_Trans     =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         = Sqrt ( TprimEt_Trans ) ;
	
	h_Histo_Mt.at(4).at(1) -> Fill( trans_mass) ;

}
}
//===========================Category Plots ===============

void     TTbar_PUPPItau_loc_190204::top_fatjet_Plots()
{
	TString Cat_pass  = "no" ;
	int jet1 = topjet[0] ;
	int jet2 = fat_jet[0] ;
	int j2 = n_Mu[0] ;
	int lvl = -1 ;
	int mu = -1 ;
	int top = -1 ;
	int W = -1;  
	int muW = -1 ;               
	float dR  =   dR_mu_AK8(j2, jet2, j2, lvl) ;
	float dPt =    dPt_lep( jet2, j2, 1, lvl);
	float 	ptcut = METCut_cat12[vr - 1] ;


	Ptbjet_dRW -> Fill( (*AK8_JetPt)[jet2] , dR) ;
	if ( dR_mu_AK8(j2, jet1, j2, lvl) >  2.5)  mu = 1 ;
	if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.5) top = 1 ;
	if ( (( dR < 1.5 && dR > 0.05) || dPt > 10.0 ) ) muW = 1 ;     
	if ( dR < 1.5 && dR > 0.00) W = 1 ;
	if ( mu != -1 && top != -1 && W != -1 && muW != -1 ) {
		Cat_pass  = "yes" ;
		CatII_Objects.push_back(jet1);               // topjet at 0th place
		CatII_Objects.push_back(jet2);			  //  fatjet at 1st place	
		CatII_Objects.push_back(j2);				  //  muon at 2nd place & MET is global
		CatII_Objects_Plots();
	}				

	//	if ( Cat_pass  == "yes" ) { }

	// genlvl_MTCalculation(1) ;    			// for muon from Higgs decayed W type
}			


void  TTbar_PUPPItau_loc_190204::top_Wjet_Plots() {
	TString Cat_pass  = "no" ;
	int jet1 = topjet[0] ;
	int jet2 = W_boson[0] ;
	int j = n_Mu[0] ;
	int lvl = -1 ;
	int mu = -1 ;
	int top = -1 ;
	int W = -1;                 
	int muW = -1 ;
	float dR  =   dR_mu_AK8(j, jet2, j, lvl) ;
	float dPt =    dPt_lep( jet2, j, 1, lvl);
	float jetPt	= (*AK8_puppiPt)[jet1] ;

	float ptcut = METCut_cat12[vr - 1] ;
	if ( dR_mu_AK8(j, jet1, j, lvl) >  2.5)  mu = 1 ;
	if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.5) top = 1 ;
	if ( dR < 1.5  ) W =1 ;
	if ( ( (dR > 0.1  ) || dPt > 20.0 ) ) muW = 1 ;     

	if ( mu != -1 && top != -1 ) {
		event_dR_Wt_tmu ++;
	if ( W != -1 && muW != -1 )  {
		event_LepIso_Wmu_cat2 ++ ;
	if ( pf_MET > ptcut ) {
		event_topPt ++ ;

		CatI_Objects.push_back(jet1);               // topjet at 0th place
		CatI_Objects.push_back(jet2);			  //  Wjet at 1st place
		CatI_Objects.push_back(j);				  //  muon at 2nd place & MET is global
		CatI_Objects_Plots();   
	}
	}
	}
	}


void  TTbar_PUPPItau_loc_190204::Higgs_lbjet_Plots() {

	TString Cat_pass  = "no" ;
	int jet1 	= Higgsjets[0] ;
	float jetPt     = (*AK8_puppiPt)[jet1];
	int jet2 	= -1 ; 
	int b2   	= -1 ;
	int j    	= n_Mu[0] ;
	int lvl  	= -1 ;
	bool pass 	= false ;
	int mu = -1 , top = -1, higg = -1 , bjet = -1 ;
	float dR = 0.0, dPt = 0.0 ;
	float ptcut = METCut_cat12[vr - 1] ;

	if ( dR_mu_AK8(j, jet1, j, lvl) < 2.0 ) {
	 pass  = true ;
	event_dR_Hmu  ++ ;
	}
	int size4  = ( b_jet.size() >= 3) ? 3 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){
		if ( jet2 != -1 ) continue ;
		b2  	= 	b_jet[h] ;                
		dR 	= 	dR_mu(j, b2, 0, lvl);
		dPt 	= 	dPt_lep(j, b2, 0, lvl);   
		if ( dR < 1.5 && dR > 0.1) top = 1  ;    
		if(((dR > 0.4 && dR < 1.5)  || dPt > 40.0))  mu = 1 ;                                            
		if ( (*jet_Pt)[b2] > 50.0 ) bjet = 1 ;
		if ( dR_AK8bjet( jet1, b2, 6, 0, lvl) > 2.0 )  higg = 1 ;           
		if ( top == 1&& mu == 1 && higg == 1 && bjet == 1)  jet2  = b2 ;
		mu = -1;
		top = -1;
		higg = -1;
	}    
	if ( pass == true && jet2 != -1 && (*mu_Pt)[j] >= 40.0 ) {
		event_LepIso_bmu ++ ;
	if ( pf_MET > ptcut ) {
		Cat_pass  = "yes"  ; 
		event_HiggsPt ++ ;
		CatVI_Objects.push_back(-1);    //vague value as all histogram defined from index =1                 
		CatVI_Objects.push_back(jet1);               // Higgsjet at 0th place
		CatVI_Objects.push_back(jet2);			  //  bjet at 1st place	
		CatVI_Objects.push_back(j);				  //  muon at 2nd place & MET is global
		CatVI_Objects_Plots();    
	}
	}

} 



void  TTbar_PUPPItau_loc_190204::W_fatjet_Plots() {
	TString Cat_pass  = "no" ;
	int jet1 = W_boson[0] ;
	int jet2mu = -1 ;  // jet2 for lvfatjet
	int jet3mu = -1 ;  // jet3 for Wfatjet
	int jet2W = -1 ;  // jet2 for lvfatjet
	int jet3W = -1 ;  // jet3 for Wfatjet	
	int j    = n_Mu[0] ;
	float dR, dPt ;
	int lvl = -1 ;
	int mu = -1 , top = -1 , higg = -1 , topbjet = -1 , hW = -1, hfat = -1, tb = -1, tmu = -1 , bj = -1;                   
	int b2  = -1;  
	/*if ( dR_mu_AK8(j, jet1, j, lvl) > 2.0) {   
	  TagJets_dRPlots() ;

	  }*/
	// only N(fatjet) =2 relevent for studying
	int size4  = ( fat_jet.size() >= 2) ? 2 : fat_jet.size() ;
	for ( int g = 0; g < fat_jet.size() ; g ++ ) {
		b2  = fat_jet[0] ;		
		dPt = dR_AK8jet(b2, jet1, jet1, 1) ;
		dR =  dR_mu_AK8(j, b2, j, lvl) ;
		jet2W = -1 ;
		jet3W = -1 ;
		jet2mu = -1 ;		
		jet3mu = -1 ;	
		//if ( dPt > 2.5 && dR < 1.5 ) Ptbjet_dRW ->Fill (dPt, dR) ;
		//if ( dPt < 2.5 && dR > 2.0 ) 	Ptbjet_dRmu ->Fill (dPt, dR) ;		

		if (  dPt > 2.5 ) jet2W =  b2 ;
		if (  dPt < 2.5 ) jet3W =  b2 ;
		if (  dR < 1.5  && dR > 0.00) jet2mu =  b2 ;
		if (  dR > 2.0 ) jet3mu =  b2 ;
		if ( jet2W != -1 && jet2mu != -1 ) jet2mu = b2 ;
		if ( jet3W != -1 && jet3mu != -1 ) jet3mu = b2 ;

	}		


	size4  = ( b_jet.size() >= 3) ? 3 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){
		b2  = b_jet[h] ;
		//         if ( dR_AK8bjet( jet2, b2, h, 1, lvl) > 2.0 ) higg = 1 ; 
		//         if ( dR_AK8bjet( jet3, b2, h, 1, lvl) > 2.0 ) hfat = 1 ;          // for Cat V
		dR = dR_AK8bjet( jet1, b2, h, 0, lvl) ;
		dPt = dR_mu(j, b2, h, lvl) ;
		if ( (*jet_Pt)[b2] > 50.0 ) bj = 1;
		if ( dR  < 1.5 && dR > 0.02) top = 1 ;
		if ( dR  > 2.0 ) hW = 1 ;
		if ( dPt > 2.5 ) mu = 1 ;
		if ( dPt < 1.5 ) tmu = 1 ;
		if ( tmu == 1 && hW == 1 && bj!= -1 ) tb = b2 ;
		if ( mu == 1 && top == 1 && bj != -1 ) {
			if ( dR > 0.00 ) topbjet = b2 ;
		}
		mu = -1;  tmu = -1;
		higg = -1 ;  top = - 1 ;
		hW = -1;  hfat = -1 ;
		bj = -1 ;
	}


	if ( dR_mu_AK8(j, jet1, j, lvl) > 2.0) mu = 1;
	//    if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.5 )higg = 1;
	if ( topbjet != -1 && mu == 1 && jet2mu != -1 )  {
		Cat_pass  = "yes"  ;
		CatIII_Objects.push_back(jet1) ;
		CatIII_Objects.push_back(jet2mu) ;
		CatIII_Objects.push_back(topbjet) ;
		CatIII_Objects.push_back(j) ;        
	}               
	//   if ( dR_AK8jet(jet3, jet1, jet1, 1) < 2.5 ) hW = 1;
	if ( tb != -1 && mu == 1 && jet3mu != -1 )  {
		CatV_Objects.push_back(jet1) ;
		CatV_Objects.push_back(jet3mu) ;
		CatV_Objects.push_back(tb) ;
		CatV_Objects.push_back(j) ;        
	}

	if (  CatV_Objects.size() != 0 )  	CatV_Objects_Plots() ;       
	if( CatIII_Objects.size() != 0 && CatV_Objects.size() == 0)     CatIII_Objects_Plots() ;

}


void  TTbar_PUPPItau_loc_190204::WW_lvbjet_Plots() 
{
	TString Cat_pass  = "no" ;
	int jet1 = -1;   // jet1 = top decayed jet
	int jet2 = -1 ;   // jet2 = higgs decayed jet
	int Wtop = W_boson[0] ;
	int Whiggs = W_boson[1] ;
	int j = n_Mu[0] ;
	float dR, dPt ;
	int lvl = -1 ;
	int b2 = -1 ;   
	int topbjet  = -1 ;    
	int mu = -1, top = -1 , higg = -1, bj = -1 ;
	float jetPt = - 1.0;
	float ptcut = METCut_cat12[vr - 1] ;

	vector <int> Higg_WJet;
	vector <int> Top_WJet;
	//  for differentiating top W from higgs W  
	dR  =   dR_mu_AK8(j, Wtop, j, lvl) ;
	dPt  =  dR_mu_AK8(j, Whiggs, j, lvl)  ;    
	if ( dR > 2.5 ) {
	Top_WJet.push_back(Wtop);
	}
	else {
	if ( dR < 1.5 ) {	
	Higg_WJet.push_back(Wtop);
	}
	}

	if ( dPt > 2.5)  {
	Top_WJet.push_back(Whiggs);
	}
	else {
	if ( dPt < 1.5 ) {
	Higg_WJet.push_back(Whiggs);
	}
	}

	if ( ( Top_WJet.size() == 1 ) && (Higg_WJet.size() == 1 ) ) {	
	jet1 	= Top_WJet[0]  ;
	jet2 	= Higg_WJet[0] ;
	jetPt     = (*AK8_puppiPt)[jet1];
	}
 
	//  for selecting top bjet                   
	int size4  = ( b_jet.size() >= 3) ? 3 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){ 
		if ( jet1 == -1 || jet2 == -1) continue ;
		//if (size4 >= 2 ) b2  = b_jet[1] ;  		             
		b2  = b_jet[h] ;  
		dR = dR_AK8bjet( jet1, b2, h, 0, lvl) ;
		dPt = dR_AK8bjet( jet2, b2, h, 1, lvl) ;
		if ( topbjet != -1) continue ;
		if ( dR_mu(j, b2, h, lvl) > 2.0 ) mu = 1 ;
		if ( dR > 0.05 &&  dR < 1.5 ) top = 1 ;  				
		if ( dPt > 2.0 ) higg = 1 ; 		
		if ( (*jet_Pt)[b2] > 50.0 ) bj = 1;   				        
		if ( mu == 1 && higg == 1 && top == 1 && bj == 1) topbjet = b2 ;

		mu = -1; 
		higg = -1;
		top = -1 ; 
		bj = -1 ;
	}

	if ( topbjet != -1 && jet1 != -1 ) {
		event_dR_Wb_bmu ++ ;
		if( jet2 != -1 ) {
		event_LepIso_Wmu_cat3 ++ ;
		if (pf_MET > ptcut ) {
		event_W_Pt ++ ;
		CatIV_Objects.push_back(jet1);               // Wjet at 0th place
		CatIV_Objects.push_back(jet2);			  //  Wjet at 1st place
		CatIV_Objects.push_back(topbjet);				//	bjet at 2nd place	
		CatIV_Objects.push_back(j);				  //  muon at 3rd place & MET is global
		CatIV_Objects_Plots();              
	}
	}
	}
	Top_WJet.clear() ;
	Higg_WJet.clear() ;
}

//=======================Plotting Functions=====================

void     TTbar_PUPPItau_loc_190204::CatI_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(0).at(3)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 0 ; g < CatI_Objects.size(); g ++ )
	{   
		prtcl = CatI_Objects[g] ;
		if (g == 0)       	  Category_Top_Plot( 0,prtcl, g) ; 
		if (g == 1)    		  Category_Wjet_Plot( 0,prtcl, g) ; 
		if (g == 2)  {
			h_Histo_Pt.at(0).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(0).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}        
		for( int h = g+1; h < CatI_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatI_Objects[h] ;          
			if ( j == 0)   {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(0).at(j) ->Fill(dR) ;
			}
			if ( j >= 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(0).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(0).at(j) ->Fill(dR) ;         
			}             				  
		}
	}
	//         genlvl_MTCalculation(1) ;
	toptag_MTCalculation(1) ; // for category I plots
}     	


void     TTbar_PUPPItau_loc_190204::CatII_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(1).at(3)  -> Fill(pf_MET) ;    // for MET being  global & [cat)(histono.]				     

	for(int g = 0 ; g < CatII_Objects.size(); g ++ )
	{   
		prtcl = CatII_Objects[g] ;
		if (g == 0)       	  Category_Top_Plot( 1,prtcl, g) ; 
		if (g == 1)          Category_Fat_Plot( 1,prtcl, g) ; 
		if (g == 2)  {
			h_Histo_Pt.at(1).at(g) -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(1).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatII_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatII_Objects[h] ;          
			if ( j == 0)   {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(1).at(j) ->Fill(dR) ;
			}
			if ( j >= 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(1).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(1).at(j) ->Fill(dR) ;         
			}             

		}
	}
	toptag_MTCalculation(2) ; // for category II plots
}     	


void     TTbar_PUPPItau_loc_190204::CatIII_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(2).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     
	WfatIII_MTCalculation( ) ;
	for(int g = 0 ; g < CatIII_Objects.size(); g ++ )
	{   
		prtcl = CatIII_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 2,prtcl, g) ; 
		if (g == 1)       	  Category_Fat_Plot( 2,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(2).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(2).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(2).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(2).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatIII_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatIII_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(2).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(2).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}  


void     TTbar_PUPPItau_loc_190204::CatIV_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;            
	WWtag_MTCalculation() ;        
	h_Histo_Pt.at(3).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     


	// for MET & object population
	for(int g = 0 ; g < CatIV_Objects.size(); g ++ )
	{   
		prtcl = CatIV_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 3,prtcl, g) ; 
		if (g == 1)       	  Category_Wjet_Plot( 3,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(3).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(3).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(3).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(3).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatIV_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatIV_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(3).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(3).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}     	


void     TTbar_PUPPItau_loc_190204::CatV_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(4).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     
	WfatV_MTCalculation( ) ;

	// for MET & object population

	for(int g = 0 ; g < CatV_Objects.size(); g ++ )
	{   
		prtcl = CatV_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 4,prtcl, g) ; 
		if (g == 1)       	  Category_Fat_Plot( 4,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(4).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(4).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(4).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(4).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatV_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatV_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(4).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(4).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}  


void     TTbar_PUPPItau_loc_190204::CatVI_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(5).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 1 ; g < CatVI_Objects.size(); g ++ )
	{   
		prtcl = CatVI_Objects[g] ;
		if (g == 1)       	  Category_Higgs_Plot( 5,prtcl, g) ; 
		if (g == 2)    {
			h_Histo_Pt.at(5).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(5).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(5).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(5).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}
		for( int h = g+1; h < CatVI_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatVI_Objects[h] ;          
			if ( j == 0 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
			}
			if ( j == 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(5).at(j) ->Fill(dR) ;         
			}
			if ( j == 2 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(5).at(j) ->Fill(dR) ;   
			}             
		}
	}

	Higgstag_MTCalculation() ;
	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}   

//=======================Default Funtions=========================



Int_t TTbar_PUPPItau_loc_190204::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TTbar_PUPPItau_loc_190204::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TTbar_PUPPItau_loc_190204::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_PID = 0;
   mc_Vtx = 0;
   mc_Vty = 0;
   mc_Vtz = 0;
   mc_Pt = 0;
   mc_Mass = 0;
   mc_Eta = 0;
   mc_Phi = 0;
   mc_E = 0;
   mc_Et = 0;
   mc_GMomPID = 0;
   mc_MomPID = 0;
   mc_MomPt = 0;
   mc_MomMass = 0;
   mc_MomEta = 0;
   mc_MomPhi = 0;
   mc_StatusFlag = 0;
   mc_Parentage = 0;
   mc_Status = 0;
   mc_CalIsoDR03 = 0;
   mc_TrkIsoDR03 = 0;
   mc_CalIsoDR04 = 0;
   mc_TrkIsoDR04 = 0;
   ele_Charge = 0;
   ele_En = 0;
   ele_D0 = 0;
   ele_Dz = 0;
   ele_Pt = 0;
   ele_Eta = 0;
   ele_Phi = 0;
   ele_R9 = 0;
   ele_SCEta = 0;
   ele_SCPhi = 0;
   ele_HoverE = 0;
   ele_EoverP = 0;
   ele_EoverPout = 0;
   ele_EoverPInv = 0;
   ele_dEtaAtVtx = 0;
   ele_dPhiAtVtx = 0;
   ele_SigmaIEtaIEtaFull5x5 = 0;
   ele_SigmaIPhiIPhiFull5x5 = 0;
   ele_ConvVeto = 0;
   ele_MissHits = 0;
   ele_PFChIso = 0;
   ele_PFPhoIso = 0;
   ele_PFNeuIso = 0;
   ele_PFMiniIso = 0;
   ele_dEtaseedAtVtx = 0;
   mu_Pt = 0;
   mu_En = 0;
   mu_Eta = 0;
   mu_Phi = 0;
   mu_Charge = 0;
   mu_IDbit = 0;
   mu_D0 = 0;
   mu_Dz = 0;
   mu_Chi2NDF = 0;
   mu_InnerD0 = 0;
   mu_InnerDz = 0;
   mu_InnervalidFraction = 0;
   mu_segmentCompatibility = 0;
   mu_chi2LocalPosition = 0;
   mu_trkKink = 0;
   mu_PFChIso = 0;
   mu_PFPhoIso = 0;
   mu_PFNeuIso = 0;
   mu_PFMiniIso = 0;
   mu_TrkLayers = 0;
   mu_PixelLayers = 0;
   mu_PixelHits = 0;
   mu_MuonHits = 0;
   mu_Stations = 0;
   mu_Matches = 0;
   mu_TrkQuality = 0;
   mu_IsoTrk = 0;
   jet_Pt = 0;
   jet_En = 0;
   jet_Eta = 0;
   jet_Phi = 0;
   jet_Area = 0;
   jet_Mt = 0;
   jet_CSV2BJetTags = 0;
   jet_JetProbabilityBJetTags = 0;
   jet_pfCombinedMVAV2BJetTags = 0;
   jet_DeepCSVTags_b = 0;
   jet_DeepCSVTags_bb = 0;
   jet_DeepCSVTags_c = 0;
   jet_DeepCSVTags_cc = 0;
   jet_DeepCSVTags_udsg = 0;
   jet_PartonID = 0;
   jet_HadFlvr = 0;
   jet_GenJetEn = 0;
   jet_GenJetPt = 0;
   jet_GenJetEta = 0;
   jet_GenJetPhi = 0;
   jet_GenPartonID = 0;
   jet_GenEn = 0;
   jet_GenPt = 0;
   jet_GenEta = 0;
   jet_GenPhi = 0;
   jet_GenPartonMomID = 0;
   jet_PFLooseId = 0;
   jet_ID = 0;
   jet_PUID = 0;
   jet_PUFullID = 0;
   jet_CHF = 0;
   jet_NHF = 0;
   jet_CEF = 0;
   jet_NEF = 0;
   jet_NCH = 0;
   jet_NNP = 0;
   jet_MUF = 0;
   AK8_JetPt = 0;
   AK8_JetEn = 0;
   AK8_JetEta = 0;
   AK8_JetPhi = 0;
   AK8_JetMass = 0;
   AK8_Jet_tau1 = 0;
   AK8_Jet_tau2 = 0;
   AK8_Jet_tau3 = 0;
   AK8_Jet_CHStau1 = 0;
   AK8_Jet_CHStau2 = 0;
   AK8_Jet_CHStau3 = 0;
   AK8_Jet_CHStau4 = 0;
   AK8_JetCHF = 0;
   AK8_JetNHF = 0;
   AK8_JetCEF = 0;
   AK8_JetNEF = 0;
   AK8_JetNCH = 0;
   AK8_JetNNP = 0;
   AK8_JetMUF = 0;
   AK8_Jetnconstituents = 0;
   AK8_JetPFLooseId = 0;
   AK8_JetPFTightLepVetoId = 0;
   AK8_JetSoftDropMass = 0;
   AK8_JetSoftDropMassCorr = 0;
   AK8_JetPrunedMassCorr = 0;
   AK8_JetL2L3corr = 0;
   AK8_puppiSDMassL2L3Corr = 0;
   AK8_puppiSDL2L3Corr = 0;
   AK8_JetPrunedMass = 0;
   AK8_JetpfBoostedDSVBTag = 0;
   AK8_JetCSV = 0;
   AK8_JetDSVnewV4 = 0;
   AK8_puppiPt = 0;
   AK8_puppiMass = 0;
   AK8_puppiEta = 0;
   AK8_puppiPhi = 0;
   AK8_puppiTau1 = 0;
   AK8_puppiTau2 = 0;
   AK8_puppiTau3 = 0;
   AK8_puppiTau4 = 0;
   AK8_puppiSDMass = 0;
   AK8_JetPartonID = 0;
   AK8_JetHadFlvr = 0;
   AK8_JetGenJetIndex = 0;
   AK8_JetGenJetEn = 0;
   AK8_JetGenJetPt = 0;
   AK8_JetGenJetEta = 0;
   AK8_JetGenJetPhi = 0;
   AK8_JetGenPartonID = 0;
   AK8_JetGenEn = 0;
   AK8_JetGenPt = 0;
   AK8_JetGenEta = 0;
   AK8_JetGenPhi = 0;
   AK8_JetGenPartonMomID = 0;
   n_AK8SDSJ = 0;
   AK8_SDSJPt = 0;
   AK8_SDSJEta = 0;
   AK8_SDSJPhi = 0;
   AK8_SDSJMass = 0;
   AK8_SDSJE = 0;
   AK8_SDSJCharge = 0;
   AK8_SDSJFlavour = 0;
   AK8_SDSJCSV = 0;
   n_AK8puppiSDSJ = 0;
   AK8_puppiSDSJPt = 0;
   AK8_puppiSDSJEta = 0;
   AK8_puppiSDSJPhi = 0;
   AK8_puppiSDSJMass = 0;
   AK8_puppiSDSJE = 0;
   AK8_puppiSDSJCharge = 0;
   AK8_puppiSDSJFlavour = 0;
   AK8_puppiSDSJCSV = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("v_event", &v_event, &b_v_event);
   fChain->SetBranchAddress("n_Vtx", &n_Vtx, &b_n_Vtx);
   fChain->SetBranchAddress("n_GoodVtx", &n_GoodVtx, &b_n_GoodVtx);
   fChain->SetBranchAddress("n_TrksPV", &n_TrksPV, &b_n_TrksPV);
   fChain->SetBranchAddress("is_PVGood", &is_PVGood, &b_is_PVGood);
   fChain->SetBranchAddress("v_vtx", &v_vtx, &b_v_vtx);
   fChain->SetBranchAddress("v_vty", &v_vty, &b_v_vty);
   fChain->SetBranchAddress("v_vtz", &v_vtz, &b_v_vtz);
   fChain->SetBranchAddress("rho_Central", &rho_Central, &b_rho_Central);
   fChain->SetBranchAddress("n_MC", &n_MC, &b_n_MC);
   fChain->SetBranchAddress("mc_PID", &mc_PID, &b_mc_PID);
   fChain->SetBranchAddress("mc_Vtx", &mc_Vtx, &b_mc_Vtx);
   fChain->SetBranchAddress("mc_Vty", &mc_Vty, &b_mc_Vty);
   fChain->SetBranchAddress("mc_Vtz", &mc_Vtz, &b_mc_Vtz);
   fChain->SetBranchAddress("mc_Pt", &mc_Pt, &b_mc_Pt);
   fChain->SetBranchAddress("mc_Mass", &mc_Mass, &b_mc_Mass);
   fChain->SetBranchAddress("mc_Eta", &mc_Eta, &b_mc_Eta);
   fChain->SetBranchAddress("mc_Phi", &mc_Phi, &b_mc_Phi);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_Et", &mc_Et, &b_mc_Et);
   fChain->SetBranchAddress("mc_GMomPID", &mc_GMomPID, &b_mc_GMomPID);
   fChain->SetBranchAddress("mc_MomPID", &mc_MomPID, &b_mc_MomPID);
   fChain->SetBranchAddress("mc_MomPt", &mc_MomPt, &b_mc_MomPt);
   fChain->SetBranchAddress("mc_MomMass", &mc_MomMass, &b_mc_MomMass);
   fChain->SetBranchAddress("mc_MomEta", &mc_MomEta, &b_mc_MomEta);
   fChain->SetBranchAddress("mc_MomPhi", &mc_MomPhi, &b_mc_MomPhi);
   fChain->SetBranchAddress("mc_StatusFlag", &mc_StatusFlag, &b_mc_StatusFlag);
   fChain->SetBranchAddress("mc_Parentage", &mc_Parentage, &b_mc_Parentage);
   fChain->SetBranchAddress("mc_Status", &mc_Status, &b_mc_Status);
   fChain->SetBranchAddress("mc_CalIsoDR03", &mc_CalIsoDR03, &b_mc_CalIsoDR03);
   fChain->SetBranchAddress("mc_TrkIsoDR03", &mc_TrkIsoDR03, &b_mc_TrkIsoDR03);
   fChain->SetBranchAddress("mc_CalIsoDR04", &mc_CalIsoDR04, &b_mc_CalIsoDR04);
   fChain->SetBranchAddress("mc_TrkIsoDR04", &mc_TrkIsoDR04, &b_mc_TrkIsoDR04);
   fChain->SetBranchAddress("gen_MET", &gen_MET, &b_gen_MET);
   fChain->SetBranchAddress("gen_METPhi", &gen_METPhi, &b_gen_METPhi);
   fChain->SetBranchAddress("pf_MET", &pf_MET, &b_pf_MET);
   fChain->SetBranchAddress("pf_METPhi", &pf_METPhi, &b_pf_METPhi);
   fChain->SetBranchAddress("pf_METsumEt", &pf_METsumEt, &b_pf_METsumEt);
   fChain->SetBranchAddress("pf_METmEtSig", &pf_METmEtSig, &b_pf_METmEtSig);
   fChain->SetBranchAddress("pf_METSig", &pf_METSig, &b_pf_METSig);
   fChain->SetBranchAddress("n_Ele", &n_Ele, &b_n_Ele);
   fChain->SetBranchAddress("ele_Charge", &ele_Charge, &b_ele_Charge);
   fChain->SetBranchAddress("ele_En", &ele_En, &b_ele_En);
   fChain->SetBranchAddress("ele_D0", &ele_D0, &b_ele_D0);
   fChain->SetBranchAddress("ele_Dz", &ele_Dz, &b_ele_Dz);
   fChain->SetBranchAddress("ele_Pt", &ele_Pt, &b_ele_Pt);
   fChain->SetBranchAddress("ele_Eta", &ele_Eta, &b_ele_Eta);
   fChain->SetBranchAddress("ele_Phi", &ele_Phi, &b_ele_Phi);
   fChain->SetBranchAddress("ele_R9", &ele_R9, &b_ele_R9);
   fChain->SetBranchAddress("ele_SCEta", &ele_SCEta, &b_ele_SCEta);
   fChain->SetBranchAddress("ele_SCPhi", &ele_SCPhi, &b_ele_SCPhi);
   fChain->SetBranchAddress("ele_HoverE", &ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_EoverP", &ele_EoverP, &b_ele_EoverP);
   fChain->SetBranchAddress("ele_EoverPout", &ele_EoverPout, &b_ele_EoverPout);
   fChain->SetBranchAddress("ele_EoverPInv", &ele_EoverPInv, &b_ele_EoverPInv);
   fChain->SetBranchAddress("ele_dEtaAtVtx", &ele_dEtaAtVtx, &b_ele_dEtaAtVtx);
   fChain->SetBranchAddress("ele_dPhiAtVtx", &ele_dPhiAtVtx, &b_ele_dPhiAtVtx);
   fChain->SetBranchAddress("ele_SigmaIEtaIEtaFull5x5", &ele_SigmaIEtaIEtaFull5x5, &b_ele_SigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("ele_SigmaIPhiIPhiFull5x5", &ele_SigmaIPhiIPhiFull5x5, &b_ele_SigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("ele_ConvVeto", &ele_ConvVeto, &b_ele_ConvVeto);
   fChain->SetBranchAddress("ele_MissHits", &ele_MissHits, &b_ele_MissHits);
   fChain->SetBranchAddress("ele_PFChIso", &ele_PFChIso, &b_ele_PFChIso);
   fChain->SetBranchAddress("ele_PFPhoIso", &ele_PFPhoIso, &b_ele_PFPhoIso);
   fChain->SetBranchAddress("ele_PFNeuIso", &ele_PFNeuIso, &b_ele_PFNeuIso);
   fChain->SetBranchAddress("ele_PFMiniIso", &ele_PFMiniIso, &b_ele_PFMiniIso);
   fChain->SetBranchAddress("ele_dEtaseedAtVtx", &ele_dEtaseedAtVtx, &b_ele_dEtaseedAtVtx);
   fChain->SetBranchAddress("N_Mu", &N_Mu, &b_N_Mu);
   fChain->SetBranchAddress("mu_Pt", &mu_Pt, &b_mu_Pt);
   fChain->SetBranchAddress("mu_En", &mu_En, &b_mu_En);
   fChain->SetBranchAddress("mu_Eta", &mu_Eta, &b_mu_Eta);
   fChain->SetBranchAddress("mu_Phi", &mu_Phi, &b_mu_Phi);
   fChain->SetBranchAddress("mu_Charge", &mu_Charge, &b_mu_Charge);
   fChain->SetBranchAddress("mu_IDbit", &mu_IDbit, &b_mu_IDbit);
   fChain->SetBranchAddress("mu_D0", &mu_D0, &b_mu_D0);
   fChain->SetBranchAddress("mu_Dz", &mu_Dz, &b_mu_Dz);
   fChain->SetBranchAddress("mu_Chi2NDF", &mu_Chi2NDF, &b_mu_Chi2NDF);
   fChain->SetBranchAddress("mu_InnerD0", &mu_InnerD0, &b_mu_InnerD0);
   fChain->SetBranchAddress("mu_InnerDz", &mu_InnerDz, &b_mu_InnerDz);
   fChain->SetBranchAddress("mu_InnervalidFraction", &mu_InnervalidFraction, &b_mu_InnervalidFraction);
   fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
   fChain->SetBranchAddress("mu_chi2LocalPosition", &mu_chi2LocalPosition, &b_mu_chi2LocalPosition);
   fChain->SetBranchAddress("mu_trkKink", &mu_trkKink, &b_mu_trkKink);
   fChain->SetBranchAddress("mu_PFChIso", &mu_PFChIso, &b_mu_PFChIso);
   fChain->SetBranchAddress("mu_PFPhoIso", &mu_PFPhoIso, &b_mu_PFPhoIso);
   fChain->SetBranchAddress("mu_PFNeuIso", &mu_PFNeuIso, &b_mu_PFNeuIso);
   fChain->SetBranchAddress("mu_PFMiniIso", &mu_PFMiniIso, &b_mu_PFMiniIso);
   fChain->SetBranchAddress("mu_TrkLayers", &mu_TrkLayers, &b_mu_TrkLayers);
   fChain->SetBranchAddress("mu_PixelLayers", &mu_PixelLayers, &b_mu_PixelLayers);
   fChain->SetBranchAddress("mu_PixelHits", &mu_PixelHits, &b_mu_PixelHits);
   fChain->SetBranchAddress("mu_MuonHits", &mu_MuonHits, &b_mu_MuonHits);
   fChain->SetBranchAddress("mu_Stations", &mu_Stations, &b_mu_Stations);
   fChain->SetBranchAddress("mu_Matches", &mu_Matches, &b_mu_Matches);
   fChain->SetBranchAddress("mu_TrkQuality", &mu_TrkQuality, &b_mu_TrkQuality);
   fChain->SetBranchAddress("mu_IsoTrk", &mu_IsoTrk, &b_mu_IsoTrk);
   fChain->SetBranchAddress("n_Jet", &n_Jet, &b_n_Jet);
   fChain->SetBranchAddress("jet_Pt", &jet_Pt, &b_jet_Pt);
   fChain->SetBranchAddress("jet_En", &jet_En, &b_jet_En);
   fChain->SetBranchAddress("jet_Eta", &jet_Eta, &b_jet_Eta);
   fChain->SetBranchAddress("jet_Phi", &jet_Phi, &b_jet_Phi);
   fChain->SetBranchAddress("jet_Area", &jet_Area, &b_jet_Area);
   fChain->SetBranchAddress("jet_Mt", &jet_Mt, &b_jet_Mt);
   fChain->SetBranchAddress("jet_CSV2BJetTags", &jet_CSV2BJetTags, &b_jet_CSV2BJetTags);
   fChain->SetBranchAddress("jet_JetProbabilityBJetTags", &jet_JetProbabilityBJetTags, &b_jet_JetProbabilityBJetTags);
   fChain->SetBranchAddress("jet_pfCombinedMVAV2BJetTags", &jet_pfCombinedMVAV2BJetTags, &b_jet_pfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("jet_DeepCSVTags_b", &jet_DeepCSVTags_b, &b_jet_DeepCSVTags_b);
   fChain->SetBranchAddress("jet_DeepCSVTags_bb", &jet_DeepCSVTags_bb, &b_jet_DeepCSVTags_bb);
   fChain->SetBranchAddress("jet_DeepCSVTags_c", &jet_DeepCSVTags_c, &b_jet_DeepCSVTags_c);
   fChain->SetBranchAddress("jet_DeepCSVTags_cc", &jet_DeepCSVTags_cc, &b_jet_DeepCSVTags_cc);
   fChain->SetBranchAddress("jet_DeepCSVTags_udsg", &jet_DeepCSVTags_udsg, &b_jet_DeepCSVTags_udsg);
   fChain->SetBranchAddress("jet_PartonID", &jet_PartonID, &b_jet_PartonID);
   fChain->SetBranchAddress("jet_HadFlvr", &jet_HadFlvr, &b_jet_HadFlvr);
   fChain->SetBranchAddress("jet_GenJetEn", &jet_GenJetEn, &b_jet_GenJetEn);
   fChain->SetBranchAddress("jet_GenJetPt", &jet_GenJetPt, &b_jet_GenJetPt);
   fChain->SetBranchAddress("jet_GenJetEta", &jet_GenJetEta, &b_jet_GenJetEta);
   fChain->SetBranchAddress("jet_GenJetPhi", &jet_GenJetPhi, &b_jet_GenJetPhi);
   fChain->SetBranchAddress("jet_GenPartonID", &jet_GenPartonID, &b_jet_GenPartonID);
   fChain->SetBranchAddress("jet_GenEn", &jet_GenEn, &b_jet_GenEn);
   fChain->SetBranchAddress("jet_GenPt", &jet_GenPt, &b_jet_GenPt);
   fChain->SetBranchAddress("jet_GenEta", &jet_GenEta, &b_jet_GenEta);
   fChain->SetBranchAddress("jet_GenPhi", &jet_GenPhi, &b_jet_GenPhi);
   fChain->SetBranchAddress("jet_GenPartonMomID", &jet_GenPartonMomID, &b_jet_GenPartonMomID);
   fChain->SetBranchAddress("jet_PFLooseId", &jet_PFLooseId, &b_jet_PFLooseId);
   fChain->SetBranchAddress("jet_ID", &jet_ID, &b_jet_ID);
   fChain->SetBranchAddress("jet_PUID", &jet_PUID, &b_jet_PUID);
   fChain->SetBranchAddress("jet_PUFullID", &jet_PUFullID, &b_jet_PUFullID);
   fChain->SetBranchAddress("jet_CHF", &jet_CHF, &b_jet_CHF);
   fChain->SetBranchAddress("jet_NHF", &jet_NHF, &b_jet_NHF);
   fChain->SetBranchAddress("jet_CEF", &jet_CEF, &b_jet_CEF);
   fChain->SetBranchAddress("jet_NEF", &jet_NEF, &b_jet_NEF);
   fChain->SetBranchAddress("jet_NCH", &jet_NCH, &b_jet_NCH);
   fChain->SetBranchAddress("jet_NNP", &jet_NNP, &b_jet_NNP);
   fChain->SetBranchAddress("jet_MUF", &jet_MUF, &b_jet_MUF);
   fChain->SetBranchAddress("N_AK8Jet", &N_AK8Jet, &b_N_AK8Jet);
   fChain->SetBranchAddress("n_AK8Jetpuppi", &n_AK8Jetpuppi, &b_n_AK8Jetpuppi);
   fChain->SetBranchAddress("AK8_JetPt", &AK8_JetPt, &b_AK8_JetPt);
   fChain->SetBranchAddress("AK8_JetEn", &AK8_JetEn, &b_AK8_JetEn);
   fChain->SetBranchAddress("AK8_JetEta", &AK8_JetEta, &b_AK8_JetEta);
   fChain->SetBranchAddress("AK8_JetPhi", &AK8_JetPhi, &b_AK8_JetPhi);
   fChain->SetBranchAddress("AK8_JetMass", &AK8_JetMass, &b_AK8_JetMass);
   fChain->SetBranchAddress("AK8_Jet_tau1", &AK8_Jet_tau1, &b_AK8_Jet_tau1);
   fChain->SetBranchAddress("AK8_Jet_tau2", &AK8_Jet_tau2, &b_AK8_Jet_tau2);
   fChain->SetBranchAddress("AK8_Jet_tau3", &AK8_Jet_tau3, &b_AK8_Jet_tau3);
   fChain->SetBranchAddress("AK8_Jet_CHStau1", &AK8_Jet_CHStau1, &b_AK8_Jet_CHStau1);
   fChain->SetBranchAddress("AK8_Jet_CHStau2", &AK8_Jet_CHStau2, &b_AK8_Jet_CHStau2);
   fChain->SetBranchAddress("AK8_Jet_CHStau3", &AK8_Jet_CHStau3, &b_AK8_Jet_CHStau3);
   fChain->SetBranchAddress("AK8_Jet_CHStau4", &AK8_Jet_CHStau4, &b_AK8_Jet_CHStau4);
   fChain->SetBranchAddress("AK8_JetCHF", &AK8_JetCHF, &b_AK8_JetCHF);
   fChain->SetBranchAddress("AK8_JetNHF", &AK8_JetNHF, &b_AK8_JetNHF);
   fChain->SetBranchAddress("AK8_JetCEF", &AK8_JetCEF, &b_AK8_JetCEF);
   fChain->SetBranchAddress("AK8_JetNEF", &AK8_JetNEF, &b_AK8_JetNEF);
   fChain->SetBranchAddress("AK8_JetNCH", &AK8_JetNCH, &b_AK8_JetNCH);
   fChain->SetBranchAddress("AK8_JetNNP", &AK8_JetNNP, &b_AK8_JetNNP);
   fChain->SetBranchAddress("AK8_JetMUF", &AK8_JetMUF, &b_AK8_JetMUF);
   fChain->SetBranchAddress("AK8_Jetnconstituents", &AK8_Jetnconstituents, &b_AK8_Jetnconstituents);
   fChain->SetBranchAddress("AK8_JetPFLooseId", &AK8_JetPFLooseId, &b_AK8_JetPFLooseId);
   fChain->SetBranchAddress("AK8_JetPFTightLepVetoId", &AK8_JetPFTightLepVetoId, &b_AK8_JetPFTightLepVetoId);
   fChain->SetBranchAddress("AK8_JetSoftDropMass", &AK8_JetSoftDropMass, &b_AK8_JetSoftDropMass);
   fChain->SetBranchAddress("AK8_JetSoftDropMassCorr", &AK8_JetSoftDropMassCorr, &b_AK8_JetSoftDropMassCorr);
   fChain->SetBranchAddress("AK8_JetPrunedMassCorr", &AK8_JetPrunedMassCorr, &b_AK8_JetPrunedMassCorr);
   fChain->SetBranchAddress("AK8_JetL2L3corr", &AK8_JetL2L3corr, &b_AK8_JetL2L3corr);
   fChain->SetBranchAddress("AK8_puppiSDMassL2L3Corr", &AK8_puppiSDMassL2L3Corr, &b_AK8_puppiSDMassL2L3Corr);
   fChain->SetBranchAddress("AK8_puppiSDL2L3Corr", &AK8_puppiSDL2L3Corr, &b_AK8_puppiSDL2L3Corr);
   fChain->SetBranchAddress("AK8_JetPrunedMass", &AK8_JetPrunedMass, &b_AK8_JetPrunedMass);
   fChain->SetBranchAddress("AK8_JetpfBoostedDSVBTag", &AK8_JetpfBoostedDSVBTag, &b_AK8_JetpfBoostedDSVBTag);
   fChain->SetBranchAddress("AK8_JetCSV", &AK8_JetCSV, &b_AK8_JetCSV);
   fChain->SetBranchAddress("AK8_JetDSVnewV4", &AK8_JetDSVnewV4, &b_AK8_JetDSVnewV4);
   fChain->SetBranchAddress("AK8_puppiPt", &AK8_puppiPt, &b_AK8_puppiPt);
   fChain->SetBranchAddress("AK8_puppiMass", &AK8_puppiMass, &b_AK8_puppiMass);
   fChain->SetBranchAddress("AK8_puppiEta", &AK8_puppiEta, &b_AK8_puppiEta);
   fChain->SetBranchAddress("AK8_puppiPhi", &AK8_puppiPhi, &b_AK8_puppiPhi);
   fChain->SetBranchAddress("AK8_puppiTau1", &AK8_puppiTau1, &b_AK8_puppiTau1);
   fChain->SetBranchAddress("AK8_puppiTau2", &AK8_puppiTau2, &b_AK8_puppiTau2);
   fChain->SetBranchAddress("AK8_puppiTau3", &AK8_puppiTau3, &b_AK8_puppiTau3);
   fChain->SetBranchAddress("AK8_puppiTau4", &AK8_puppiTau4, &b_AK8_puppiTau4);
   fChain->SetBranchAddress("AK8_puppiSDMass", &AK8_puppiSDMass, &b_AK8_puppiSDMass);
   fChain->SetBranchAddress("AK8_JetPartonID", &AK8_JetPartonID, &b_AK8_JetPartonID);
   fChain->SetBranchAddress("AK8_JetHadFlvr", &AK8_JetHadFlvr, &b_AK8_JetHadFlvr);
   fChain->SetBranchAddress("AK8_JetGenJetIndex", &AK8_JetGenJetIndex, &b_AK8_JetGenJetIndex);
   fChain->SetBranchAddress("AK8_JetGenJetEn", &AK8_JetGenJetEn, &b_AK8_JetGenJetEn);
   fChain->SetBranchAddress("AK8_JetGenJetPt", &AK8_JetGenJetPt, &b_AK8_JetGenJetPt);
   fChain->SetBranchAddress("AK8_JetGenJetEta", &AK8_JetGenJetEta, &b_AK8_JetGenJetEta);
   fChain->SetBranchAddress("AK8_JetGenJetPhi", &AK8_JetGenJetPhi, &b_AK8_JetGenJetPhi);
   fChain->SetBranchAddress("AK8_JetGenPartonID", &AK8_JetGenPartonID, &b_AK8_JetGenPartonID);
   fChain->SetBranchAddress("AK8_JetGenEn", &AK8_JetGenEn, &b_AK8_JetGenEn);
   fChain->SetBranchAddress("AK8_JetGenPt", &AK8_JetGenPt, &b_AK8_JetGenPt);
   fChain->SetBranchAddress("AK8_JetGenEta", &AK8_JetGenEta, &b_AK8_JetGenEta);
   fChain->SetBranchAddress("AK8_JetGenPhi", &AK8_JetGenPhi, &b_AK8_JetGenPhi);
   fChain->SetBranchAddress("AK8_JetGenPartonMomID", &AK8_JetGenPartonMomID, &b_AK8_JetGenPartonMomID);
   fChain->SetBranchAddress("n_AK8SDSJ", &n_AK8SDSJ, &b_n_AK8SDSJ);
   fChain->SetBranchAddress("AK8_SDSJPt", &AK8_SDSJPt, &b_AK8_SDSJPt);
   fChain->SetBranchAddress("AK8_SDSJEta", &AK8_SDSJEta, &b_AK8_SDSJEta);
   fChain->SetBranchAddress("AK8_SDSJPhi", &AK8_SDSJPhi, &b_AK8_SDSJPhi);
   fChain->SetBranchAddress("AK8_SDSJMass", &AK8_SDSJMass, &b_AK8_SDSJMass);
   fChain->SetBranchAddress("AK8_SDSJE", &AK8_SDSJE, &b_AK8_SDSJE);
   fChain->SetBranchAddress("AK8_SDSJCharge", &AK8_SDSJCharge, &b_AK8_SDSJCharge);
   fChain->SetBranchAddress("AK8_SDSJFlavour", &AK8_SDSJFlavour, &b_AK8_SDSJFlavour);
   fChain->SetBranchAddress("AK8_SDSJCSV", &AK8_SDSJCSV, &b_AK8_SDSJCSV);
   fChain->SetBranchAddress("n_AK8puppiSDSJ", &n_AK8puppiSDSJ, &b_n_AK8puppiSDSJ);
   fChain->SetBranchAddress("AK8_puppiSDSJPt", &AK8_puppiSDSJPt, &b_AK8_puppiSDSJPt);
   fChain->SetBranchAddress("AK8_puppiSDSJEta", &AK8_puppiSDSJEta, &b_AK8_puppiSDSJEta);
   fChain->SetBranchAddress("AK8_puppiSDSJPhi", &AK8_puppiSDSJPhi, &b_AK8_puppiSDSJPhi);
   fChain->SetBranchAddress("AK8_puppiSDSJMass", &AK8_puppiSDSJMass, &b_AK8_puppiSDSJMass);
   fChain->SetBranchAddress("AK8_puppiSDSJE", &AK8_puppiSDSJE, &b_AK8_puppiSDSJE);
   fChain->SetBranchAddress("AK8_puppiSDSJCharge", &AK8_puppiSDSJCharge, &b_AK8_puppiSDSJCharge);
   fChain->SetBranchAddress("AK8_puppiSDSJFlavour", &AK8_puppiSDSJFlavour, &b_AK8_puppiSDSJFlavour);
   fChain->SetBranchAddress("AK8_puppiSDSJCSV", &AK8_puppiSDSJCSV, &b_AK8_puppiSDSJCSV);
   Notify();
}

Bool_t TTbar_PUPPItau_loc_190204::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TTbar_PUPPItau_loc_190204::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TTbar_PUPPItau_loc_190204::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TTbar_PUPPItau_loc_190204_cxx
