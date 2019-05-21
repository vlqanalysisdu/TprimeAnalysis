#ifndef TTbar_mc_bkg_PUPPI_030219_h
#define TTbar_mc_bkg_PUPPI_030219_h

// Header file for the classes stored in the TTree if any.
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
//#include "TDCacheFile.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include <map>
//#include <TRFIOFile.h>
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <TLatex.h>
#include <Riostream.h>
#include <set>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TKDE.h"


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
using namespace ROOT;

class TTbar_mc_bkg_PUPPI_030219 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   //Int_t           nTrksPV;
   Bool_t          isPVGood;
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         genPho1;
   Float_t         genPho2;
   TString         *EventTag;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR04;
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           metFilters;
   Float_t         pfMET;
   //Float_t         pfMETPhi;
   //Float_t         pfMETsumEt;
   //Float_t         pfMETmEtSig;
   //Float_t         pfMETSig;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   //vector<float>   *phoE1x3;
   //vector<float>   *phoE1x5;
   //vector<float>   *phoE2x2;
   //vector<float>   *phoE2x5Max;
   //vector<float>   *phoE5x5;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   //vector<float>   *phoE1x3Full5x5;
   //vector<float>   *phoE1x5Full5x5;
   vector<float>   *phoE2x2Full5x5;
   //vector<float>   *phoE2x5MaxFull5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   //vector<float>   *phoPFRandConeChIso;
   //vector<float>   *phoCITKChIso;
   //vector<float>   *phoCITKPhoIso;
   //vector<float>   *phoCITKNeuIso;
   vector<float>   *phoIDMVA;
   vector<ULong64_t> *phoFiredSingleTrgs;
   vector<ULong64_t> *phoFiredDoubleTrgs;
   vector<ULong64_t> *phoFiredL1Trgs;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<unsigned short> *phoxtalBits;
   vector<unsigned short> *phoIDbit;
   vector<float>   *phoScale_stat_up;
   vector<float>   *phoScale_stat_dn;
   vector<float>   *phoScale_syst_up;
   vector<float>   *phoScale_syst_dn;
   vector<float>   *phoScale_gain_up;
   vector<float>   *phoScale_gain_dn;
   vector<float>   *phoResol_rho_up;
   vector<float>   *phoResol_rho_dn;
   vector<float>   *phoResol_phi_up;
   vector<float>   *phoResol_phi_dn;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleEcalEn;
   vector<float>   *eleESEnP1;
   vector<float>   *eleESEnP2;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   //vector<float>   *eledEtaAtCalo;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   //vector<float>   *elePFMiniIso;
   //vector<float>   *eleIDMVA;
   //vector<float>   *eleIDMVAHZZ;
   //vector<float>   *eledEtaseedAtVtx;
   //vector<float>   *eleE1x5;
   //vector<float>   *eleE2x5;
   //vector<float>   *eleE5x5;
   //vector<float>   *eleE1x5Full5x5;
   //vector<float>   *eleE2x5Full5x5;
   //vector<float>   *eleE5x5Full5x5;
   vector<float>   *eleR9Full5x5;
   vector<int>     *eleEcalDrivenSeed;
   //vector<float>   *eleDr03EcalRecHitSumEt;
   //vector<float>   *eleDr03HcalDepth1TowerSumEt;
   //vector<float>   *eleDr03HcalDepth2TowerSumEt;
   //vector<float>   *eleDr03HcalTowerSumEt;
   //vector<float>   *eleDr03TkSumPt;
   //vector<float>   *elecaloEnergy;
   vector<float>   *eleTrkdxy;
   vector<float>   *eleKFHits;
   vector<float>   *eleKFChi2;
   vector<float>   *eleGSFChi2;
   vector<vector<float> > *eleGSFPt;
   vector<vector<float> > *eleGSFEta;
   vector<vector<float> > *eleGSFPhi;
   vector<vector<float> > *eleGSFCharge;
   vector<vector<int> > *eleGSFHits;
   vector<vector<int> > *eleGSFMissHits;
   vector<vector<int> > *eleGSFNHitsMax;
   vector<vector<float> > *eleGSFVtxProb;
   vector<vector<float> > *eleGSFlxyPV;
   vector<vector<float> > *eleGSFlxyBS;
   //vector<vector<float> > *eleBCEn;
   //vector<vector<float> > *eleBCEta;
   //vector<vector<float> > *eleBCPhi;
   //vector<vector<float> > *eleBCS25;
   //vector<vector<float> > *eleBCS15;
   //vector<vector<float> > *eleBCSieie;
   //vector<vector<float> > *eleBCSieip;
   //vector<vector<float> > *eleBCSipip;
   vector<ULong64_t> *eleFiredSingleTrgs;
   vector<ULong64_t> *eleFiredDoubleTrgs;
   vector<ULong64_t> *eleFiredL1Trgs;
   vector<unsigned short> *eleIDbit;
   vector<float>   *eleScale_stat_up;
   vector<float>   *eleScale_stat_dn;
   vector<float>   *eleScale_syst_up;
   vector<float>   *eleScale_syst_dn;
   vector<float>   *eleScale_gain_up;
   vector<float>   *eleScale_gain_dn;
   vector<float>   *eleResol_rho_up;
   vector<float>   *eleResol_rho_dn;
   vector<float>   *eleResol_phi_up;
   vector<float>   *eleResol_phi_dn;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<unsigned short> *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<float>   *muPFChIso03;
   vector<float>   *muPFPhoIso03;
   vector<float>   *muPFNeuIso03;
   vector<float>   *muPFPUIso03;
   //vector<float>   *muPFMiniIso;
   vector<ULong64_t> *muFiredTrgs;
   vector<ULong64_t> *muFiredL1Trgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   vector<int>     *muBestTrkType;
   Int_t           npfPho;
   //vector<float>   *pfphoEt;
   //vector<float>   *pfphoEta;
   //vector<float>   *pfphoPhi;
   Int_t           npfHF;
   vector<float>   *pfHFEn;
   vector<float>   *pfHFECALEn;
   vector<float>   *pfHFHCALEn;
   vector<float>   *pfHFPt;
   vector<float>   *pfHFEta;
   vector<float>   *pfHFPhi;
   vector<float>   *pfHFIso;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetLeadTrackEta;
   vector<float>   *jetLeadTrackPhi;
   vector<int>     *jetLepTrackPID;
   vector<float>   *jetLepTrackPt;
   vector<float>   *jetLepTrackEta;
   vector<float>   *jetLepTrackPhi;
   vector<float>   *jetCSV2BJetTags;
   //vector<float>   *jetJetProbabilityBJetTags;
   //vector<float>   *jetpfCombinedMVAV2BJetTags;
   vector<float>   *jetDeepCSVTags_b;
   vector<float>   *jetDeepCSVTags_bb;
   vector<float>   *jetDeepCSVTags_c;
   //vector<float>   *jetDeepCSVTags_cc;
   vector<float>   *jetDeepCSVTags_udsg;
   vector<int>     *jetPartonID;
   vector<int>     *jetHadFlvr;
   vector<float>   *jetGenJetEn;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenEn;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   vector<float>   *jetP4Smear;
   vector<float>   *jetP4SmearUp;
   vector<float>   *jetP4SmearDo;
   vector<bool>    *jetPFLooseId;
   vector<int>     *jetID;
   vector<float>   *jetPUID;
   vector<int>     *jetPUFullID;
   vector<float>   *jetJECUnc;
   vector<ULong64_t> *jetFiredTrgs;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<int>     *jetNCH;
   vector<int>     *jetNNP;
   vector<float>   *jetMUF;
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtxNtrks;
   vector<float>   *jetVtx3DVal;
   vector<float>   *jetVtx3DSig;
   //Int_t           nAK8Jetpuppi;
   Int_t           nAK8Jet;
   vector<float>   *AK8JetPt;
   vector<float>   *AK8JetEn;
   vector<float>   *AK8JetRawPt;
   vector<float>   *AK8JetRawEn;
   vector<float>   *AK8JetEta;
   vector<float>   *AK8JetPhi;
   vector<float>   *AK8JetMass;
   vector<float>   *AK8Jet_tau1;
   vector<float>   *AK8Jet_tau2;
   vector<float>   *AK8Jet_tau3;
   //vector<float>   *AK8Jet_CHStau1;
   //vector<float>   *AK8Jet_CHStau2;
   //vector<float>   *AK8Jet_CHStau3;
   //vector<float>   *AK8Jet_CHStau4;
   vector<float>   *AK8JetCHF;
   vector<float>   *AK8JetNHF;
   vector<float>   *AK8JetCEF;
   vector<float>   *AK8JetNEF;
   vector<int>     *AK8JetNCH;
   vector<int>     *AK8JetNNP;
   vector<float>   *AK8JetMUF;
   vector<int>     *AK8Jetnconstituents;
   vector<bool>    *AK8JetPFLooseId;
   vector<bool>    *AK8JetPFTightLepVetoId;
   vector<float>   *AK8JetSoftDropMass;
   vector<float>   *AK8JetSoftDropMassCorr;
   vector<float>   *AK8JetPrunedMass;
   vector<float>   *AK8JetPrunedMassCorr;
   vector<float>   *AK8JetpfBoostedDSVBTag;
   vector<float>   *AK8JetDSVnewV4;
   vector<float>   *AK8JetCSV;
   vector<float>   *AK8JetJECUnc;
   vector<float>   *AK8JetL2L3corr;
   vector<float>   *AK8puppiPt;
   vector<float>   *AK8puppiMass;
   vector<float>   *AK8puppiEta;
   vector<float>   *AK8puppiPhi;
   vector<float>   *AK8puppiTau1;
   vector<float>   *AK8puppiTau2;
   vector<float>   *AK8puppiTau3;
   //vector<float>   *AK8puppiTau4;
   vector<float>   *AK8puppiSDL2L3corr;
   vector<float>   *AK8puppiSDMass;
   vector<float>   *AK8puppiSDMassL2L3Corr;
   vector<int>     *AK8JetPartonID;
   vector<int>     *AK8JetHadFlvr;
   vector<int>     *AK8JetGenJetIndex;
   vector<float>   *AK8JetGenJetEn;
   vector<float>   *AK8JetGenJetPt;
   vector<float>   *AK8JetGenJetEta;
   vector<float>   *AK8JetGenJetPhi;
   vector<int>     *AK8JetGenPartonID;
   vector<float>   *AK8JetGenEn;
   vector<float>   *AK8JetGenPt;
   vector<float>   *AK8JetGenEta;
   vector<float>   *AK8JetGenPhi;
   vector<int>     *AK8JetGenPartonMomID;
   vector<float>   *AK8JetP4Smear;
   vector<float>   *AK8JetP4SmearUp;
   vector<float>   *AK8JetP4SmearDo;
   vector<int>     *nAK8SDSJ;
   vector<vector<float> > *AK8SDSJPt;
   vector<vector<float> > *AK8SDSJEta;
   vector<vector<float> > *AK8SDSJPhi;
   vector<vector<float> > *AK8SDSJMass;
   vector<vector<float> > *AK8SDSJE;
   vector<vector<int> > *AK8SDSJCharge;
   vector<vector<int> > *AK8SDSJFlavour;
   vector<vector<float> > *AK8SDSJCSV;
   vector<int>     *nAK8puppiSDSJ;
   vector<vector<float> > *AK8puppiSDSJPt;
   vector<vector<float> > *AK8puppiSDSJEta;
   vector<vector<float> > *AK8puppiSDSJPhi;
   vector<vector<float> > *AK8puppiSDSJMass;
   vector<vector<float> > *AK8puppiSDSJE;
   vector<vector<int> > *AK8puppiSDSJCharge;
   vector<vector<int> > *AK8puppiSDSJFlavour;
   vector<vector<float> > *AK8puppiSDSJCSV;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   //TBranch        *b_nTrksPV;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_genPho1;   //!
   TBranch        *b_genPho2;   //!
   TBranch        *b_EventTag;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   //TBranch        *b_pfMETPhi;   //!
   //TBranch        *b_pfMETsumEt;   //!
   //TBranch        *b_pfMETmEtSig;   //!
   //TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   //TBranch        *b_phoE1x3;   //!
   //TBranch        *b_phoE1x5;   //!
   //TBranch        *b_phoE2x2;   //!
   //TBranch        *b_phoE2x5Max;   //!
   //TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   //TBranch        *b_phoE1x3Full5x5;   //!
   //TBranch        *b_phoE1x5Full5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   //TBranch        *b_phoE2x5MaxFull5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   //TBranch        *b_phoPFRandConeChIso;   //!
   //TBranch        *b_phoCITKChIso;   //!
   //TBranch        *b_phoCITKPhoIso;   //!
   //TBranch        *b_phoCITKNeuIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoxtalBits;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_phoScale_stat_up;   //!
   TBranch        *b_phoScale_stat_dn;   //!
   TBranch        *b_phoScale_syst_up;   //!
   TBranch        *b_phoScale_syst_dn;   //!
   TBranch        *b_phoScale_gain_up;   //!
   TBranch        *b_phoScale_gain_dn;   //!
   TBranch        *b_phoResol_rho_up;   //!
   TBranch        *b_phoResol_rho_dn;   //!
   TBranch        *b_phoResol_phi_up;   //!
   TBranch        *b_phoResol_phi_dn;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleESEnP1;   //!
   TBranch        *b_eleESEnP2;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   //TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   //TBranch        *b_elePFMiniIso;   //!
   //TBranch        *b_eleIDMVA;   //!
   //TBranch        *b_eleIDMVAHZZ;   //!
   //TBranch        *b_eledEtaseedAtVtx;   //!
   //TBranch        *b_eleE1x5;   //!
   //TBranch        *b_eleE2x5;   //!
   //TBranch        *b_eleE5x5;   //!
   //TBranch        *b_eleE1x5Full5x5;   //!
   //TBranch        *b_eleE2x5Full5x5;   //!
   //TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   //TBranch        *b_eleDr03EcalRecHitSumEt;   //!
   //TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
   //TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
   //TBranch        *b_eleDr03HcalTowerSumEt;   //!
   //TBranch        *b_eleDr03TkSumPt;   //!
   //TBranch        *b_elecaloEnergy;   //!
   TBranch        *b_eleTrkdxy;   //!
   TBranch        *b_eleKFHits;   //!
   TBranch        *b_eleKFChi2;   //!
   TBranch        *b_eleGSFChi2;   //!
   TBranch        *b_eleGSFPt;   //!
   TBranch        *b_eleGSFEta;   //!
   TBranch        *b_eleGSFPhi;   //!
   TBranch        *b_eleGSFCharge;   //!
   TBranch        *b_eleGSFHits;   //!
   TBranch        *b_eleGSFMissHits;   //!
   TBranch        *b_eleGSFNHitsMax;   //!
   TBranch        *b_eleGSFVtxProb;   //!
   TBranch        *b_eleGSFlxyPV;   //!
   TBranch        *b_eleGSFlxyBS;   //!
   //TBranch        *b_eleBCEn;   //!
   //TBranch        *b_eleBCEta;   //!
   //TBranch        *b_eleBCPhi;   //!
   //TBranch        *b_eleBCS25;   //!
   //TBranch        *b_eleBCS15;   //!
   //TBranch        *b_eleBCSieie;   //!
   //TBranch        *b_eleBCSieip;   //!
   //TBranch        *b_eleBCSipip;   //!
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleScale_stat_up;   //!
   TBranch        *b_eleScale_stat_dn;   //!
   TBranch        *b_eleScale_syst_up;   //!
   TBranch        *b_eleScale_syst_dn;   //!
   TBranch        *b_eleScale_gain_up;   //!
   TBranch        *b_eleScale_gain_dn;   //!
   TBranch        *b_eleResol_rho_up;   //!
   TBranch        *b_eleResol_rho_dn;   //!
   TBranch        *b_eleResol_phi_up;   //!
   TBranch        *b_eleResol_phi_dn;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muPFChIso03;   //!
   TBranch        *b_muPFPhoIso03;   //!
   TBranch        *b_muPFNeuIso03;   //!
   TBranch        *b_muPFPUIso03;   //!
   //TBranch        *b_muPFMiniIso;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_muBestTrkType;   //!
   TBranch        *b_npfPho;   //!
   //TBranch        *b_pfphoEt;   //!
   //TBranch        *b_pfphoEta;   //!
   //TBranch        *b_pfphoPhi;   //!
   TBranch        *b_npfHF;   //!
   TBranch        *b_pfHFEn;   //!
   TBranch        *b_pfHFECALEn;   //!
   TBranch        *b_pfHFHCALEn;   //!
   TBranch        *b_pfHFPt;   //!
   TBranch        *b_pfHFEta;   //!
   TBranch        *b_pfHFPhi;   //!
   TBranch        *b_pfHFIso;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetLeadTrackEta;   //!
   TBranch        *b_jetLeadTrackPhi;   //!
   TBranch        *b_jetLepTrackPID;   //!
   TBranch        *b_jetLepTrackPt;   //!
   TBranch        *b_jetLepTrackEta;   //!
   TBranch        *b_jetLepTrackPhi;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   //TBranch        *b_jetJetProbabilityBJetTags;   //!
   //TBranch        *b_jetpfCombinedMVAV2BJetTags;   //!
   TBranch        *b_jetDeepCSVTags_b;   //!
   TBranch        *b_jetDeepCSVTags_bb;   //!
   TBranch        *b_jetDeepCSVTags_c;   //!
   //TBranch        *b_jetDeepCSVTags_cc;   //!
   TBranch        *b_jetDeepCSVTags_udsg;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetHadFlvr;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_jetP4Smear;   //!
   TBranch        *b_jetP4SmearUp;   //!
   TBranch        *b_jetP4SmearDo;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetPUFullID;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetFiredTrgs;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetNNP;   //!
   TBranch        *b_jetMUF;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtxNtrks;   //!
   TBranch        *b_jetVtx3DVal;   //!
   TBranch        *b_jetVtx3DSig;   //!
   //TBranch        *b_nAK8Jetpuppi;   //!
   TBranch        *b_nAK8Jet;   //!
   TBranch        *b_AK8JetPt;   //!
   TBranch        *b_AK8JetEn;   //!
   TBranch        *b_AK8JetRawPt;   //!
   TBranch        *b_AK8JetRawEn;   //!
   TBranch        *b_AK8JetEta;   //!
   TBranch        *b_AK8JetPhi;   //!
   TBranch        *b_AK8JetMass;   //!
   TBranch        *b_AK8Jet_tau1;   //!
   TBranch        *b_AK8Jet_tau2;   //!
   TBranch        *b_AK8Jet_tau3;   //!
   //TBranch        *b_AK8Jet_CHStau1;   //!
   //TBranch        *b_AK8Jet_CHStau2;   //!
   //TBranch        *b_AK8Jet_CHStau3;   //!
   //TBranch        *b_AK8Jet_CHStau4;   //!
   TBranch        *b_AK8JetCHF;   //!
   TBranch        *b_AK8JetNHF;   //!
   TBranch        *b_AK8JetCEF;   //!
   TBranch        *b_AK8JetNEF;   //!
   TBranch        *b_AK8JetNCH;   //!
   TBranch        *b_AK8JetNNP;   //!
   TBranch        *b_AK8JetMUF;   //!
   TBranch        *b_AK8Jetnconstituents;   //!
   TBranch        *b_AK8JetPFLooseId;   //!
   TBranch        *b_AK8JetPFTightLepVetoId;   //!
   TBranch        *b_AK8JetSoftDropMass;   //!
   TBranch        *b_AK8JetSoftDropMassCorr;   //!
   TBranch        *b_AK8JetPrunedMass;   //!
   TBranch        *b_AK8JetPrunedMassCorr;   //!
   TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
   TBranch        *b_AK8JetDSVnewV4;   //!
   TBranch        *b_AK8JetCSV;   //!
   TBranch        *b_AK8JetJECUnc;   //!
   TBranch        *b_AK8JetL2L3corr;   //!
   TBranch        *b_AK8puppiPt;   //!
   TBranch        *b_AK8puppiMass;   //!
   TBranch        *b_AK8puppiEta;   //!
   TBranch        *b_AK8puppiPhi;   //!
   TBranch        *b_AK8puppiTau1;   //!
   TBranch        *b_AK8puppiTau2;   //!
   TBranch        *b_AK8puppiTau3;   //!
   //TBranch        *b_AK8puppiTau4;   //!
   TBranch        *b_AK8puppiSDL2L3corr;   //!
   TBranch        *b_AK8puppiSDMass;   //!
   TBranch        *b_AK8puppiSDMassL2L3Corr;   //!
   TBranch        *b_AK8JetPartonID;   //!
   TBranch        *b_AK8JetHadFlvr;   //!
   TBranch        *b_AK8JetGenJetIndex;   //!
   TBranch        *b_AK8JetGenJetEn;   //!
   TBranch        *b_AK8JetGenJetPt;   //!
   TBranch        *b_AK8JetGenJetEta;   //!
   TBranch        *b_AK8JetGenJetPhi;   //!
   TBranch        *b_AK8JetGenPartonID;   //!
   TBranch        *b_AK8JetGenEn;   //!
   TBranch        *b_AK8JetGenPt;   //!
   TBranch        *b_AK8JetGenEta;   //!
   TBranch        *b_AK8JetGenPhi;   //!
   TBranch        *b_AK8JetGenPartonMomID;   //!
   TBranch        *b_AK8JetP4Smear;   //!
   TBranch        *b_AK8JetP4SmearUp;   //!
   TBranch        *b_AK8JetP4SmearDo;   //!
   TBranch        *b_nAK8SDSJ;   //!
   TBranch        *b_AK8SDSJPt;   //!
   TBranch        *b_AK8SDSJEta;   //!
   TBranch        *b_AK8SDSJPhi;   //!
   TBranch        *b_AK8SDSJMass;   //!
   TBranch        *b_AK8SDSJE;   //!
   TBranch        *b_AK8SDSJCharge;   //!
   TBranch        *b_AK8SDSJFlavour;   //!
   TBranch        *b_AK8SDSJCSV;   //!
   TBranch        *b_nAK8puppiSDSJ;   //!
   TBranch        *b_AK8puppiSDSJPt;   //!
   TBranch        *b_AK8puppiSDSJEta;   //!
   TBranch        *b_AK8puppiSDSJPhi;   //!
   TBranch        *b_AK8puppiSDSJMass;   //!
   TBranch        *b_AK8puppiSDSJE;   //!
   TBranch        *b_AK8puppiSDSJCharge;   //!
   TBranch        *b_AK8puppiSDSJFlavour;   //!
   TBranch        *b_AK8puppiSDSJCSV;   //!

	//---Vector elements --------------------------------------------------------
		//1. Extra variables
		Int_t           v_event;
		Int_t           n_Vtx;
		Int_t           n_GoodVtx;
		Int_t           n_TrksPV;
		Bool_t          is_PVGood;
		Float_t         v_vtx;
		Float_t         v_vty;
		Float_t         v_vtz;
		Float_t         rho_Central;

		// mc variables

		Int_t           n_MC;
		vector<int>     mc_PID;
		vector<float>   mc_Vtx;
		vector<float>   mc_Vty;
		vector<float>   mc_Vtz;
		vector<float>   mc_Pt;
		vector<float>   mc_Mass;
		vector<float>   mc_Eta;
		vector<float>   mc_Phi;
		vector<float>   mc_E;
		vector<float>   mc_Et;
		vector<int>     mc_GMomPID;
		vector<int>     mc_MomPID;
		vector<float>   mc_MomPt;
		vector<float>   mc_MomMass;
		vector<float>   mc_MomEta;
		vector<float>   mc_MomPhi;
		vector<unsigned short> mc_StatusFlag;
		vector<int>     mc_Parentage;
		vector<int>     mc_Status;
		vector<float>   mc_CalIsoDR03;
		vector<float>   mc_TrkIsoDR03;
		vector<float>   mc_CalIsoDR04;
		vector<float>   mc_TrkIsoDR04;
		Float_t         gen_MET;
		Float_t         gen_METPhi;


		//3. Electron 
		Int_t           n_Ele;
		vector<int>     ele_Charge;
		vector<float>   ele_En;
		vector<float>   ele_D0;
		vector<float>   ele_Dz;
		vector<float>   ele_Pt;
		vector<float>   ele_Eta;
		vector<float>   ele_Phi;
		vector<float>   ele_R9;
		vector<float>   ele_SCEta;
		vector<float>   ele_SCPhi;		
		vector<float>   ele_HoverE;
		vector<float>   ele_EoverP;
		vector<float>   ele_EoverPout;
		vector<float>   ele_EoverPInv;
		vector<float>   ele_dEtaAtVtx;
		vector<float>   ele_dPhiAtVtx;

		vector<float>   ele_SigmaIEtaIEtaFull5x5;
		vector<float>   ele_SigmaIPhiIPhiFull5x5;
		vector<int>     ele_ConvVeto;
		vector<int>     ele_MissHits;

		vector<float>   ele_PFChIso;
		vector<float>   ele_PFPhoIso;
		vector<float>   ele_PFNeuIso;
		vector<float>   ele_PFMiniIso;
		vector<float>   ele_dEtaseedAtVtx; 


		//4. Muon 

		Int_t           N_Mu;
		vector<float>   mu_Pt;
		vector<float>   mu_En;
		vector<float>   mu_Eta;
		vector<float>   mu_Phi;
		vector<int>     mu_Charge;
		vector<unsigned short> mu_IDbit;
		vector<float>   mu_D0;
		vector<float>   mu_Dz;		
		vector<float>   mu_Chi2NDF;
		vector<float>   mu_InnerD0;
		vector<float>   mu_InnerDz;

		vector<float>   mu_PFChIso;
		vector<float>   mu_PFPhoIso;
		vector<float>   mu_PFNeuIso;
		vector<float>   mu_PFMiniIso;

		vector<float>   mu_InnervalidFraction;
		vector<float>   mu_segmentCompatibility;
		vector<float>   mu_chi2LocalPosition;
		vector<float>   mu_trkKink;

		vector<int>     mu_TrkLayers;
		vector<int>     mu_PixelLayers;
		vector<int>     mu_PixelHits;
		vector<int>     mu_MuonHits;
		vector<int>     mu_Stations;
		vector<int>     mu_Matches;
		vector<int>     mu_TrkQuality;
		vector<float>   mu_IsoTrk;
		//5. MET Variables

		Float_t         pf_MET;
		Float_t         pf_METPhi;
		Float_t         pf_METsumEt;
		Float_t         pf_METmEtSig;
		//Float_t         pf_METSig;


		//6. Bjet 

		Int_t           n_Jet;
		vector<float>   jet_Pt;
		vector<float>   jet_En;
		vector<float>   jet_Eta;
		vector<float>   jet_Phi;
		vector<float>   jet_Area;
		vector<float>   jet_Mt;

		vector<float>   jet_CSV2BJetTags;
		vector<float>   jet_JetProbabilityBJetTags;
		vector<float>   jet_pfCombinedMVAV2BJetTags;
		vector<float>   jet_DeepCSVTags_b;
		vector<float>   jet_DeepCSVTags_bb;
		vector<float>   jet_DeepCSVTags_c;
		vector<float>   jet_DeepCSVTags_cc;
		vector<float>   jet_DeepCSVTags_udsg;

		vector<int>     jet_PartonID;
		vector<int>     jet_HadFlvr;
		vector<float>   jet_GenJetEn;
		vector<float>   jet_GenJetPt;
		vector<float>   jet_GenJetEta;
		vector<float>   jet_GenJetPhi;
		vector<int>     jet_GenPartonID;
		vector<float>   jet_GenEn;
		vector<float>   jet_GenPt;
		vector<float>   jet_GenEta;
		vector<float>   jet_GenPhi;
		vector<int>     jet_GenPartonMomID;

		vector<bool>    jet_PFLooseId;
		vector<int>     jet_ID;
		vector<float>   jet_PUID;
		vector<int>     jet_PUFullID;
		vector<float>   jet_CHF;
		vector<float>   jet_NHF;
		vector<float>   jet_CEF;
		vector<float>   jet_NEF;
		vector<int>     jet_NCH;
		vector<int>     jet_NNP;
		vector<float>   jet_MUF;

		//7. Fat jets

		Int_t           N_AK8Jet;
	 	Int_t 		n_AK8Jetpuppi;
		vector<float>   AK8_JetPt;
		vector<float>   AK8_JetEn;

		vector<float>   AK8_JetEta;
		vector<float>   AK8_JetPhi;
		vector<float>   AK8_JetMass;
		vector<float>   AK8_Jet_tau1;
		vector<float>   AK8_Jet_tau2;
		vector<float>   AK8_Jet_tau3;
		vector<float>   AK8_Jet_CHStau1;
		vector<float>   AK8_Jet_CHStau2;
		vector<float>   AK8_Jet_CHStau3;
		vector<float>   AK8_Jet_CHStau4;
		vector<float>   AK8_JetCHF;
		vector<float>   AK8_JetNHF;
		vector<float>   AK8_JetCEF;
		vector<float>   AK8_JetNEF;
		vector<int>     AK8_JetNCH;
		vector<int>     AK8_JetNNP;
		vector<float>   AK8_JetMUF;
		vector<int>     AK8_Jetnconstituents;
		vector<bool>    AK8_JetPFLooseId;
		vector<bool>    AK8_JetPFTightLepVetoId;
		vector<float>   AK8_JetSoftDropMass;
		vector<float>   AK8_JetSoftDropMassCorr;

		vector<float>   AK8_JetPrunedMass;
		vector<float>   AK8_JetPrunedMassCorr;
		vector<float>   AK8_JetpfBoostedDSVBTag;
		vector<float>   AK8_JetDSVnewV4;
		vector<float>   AK8_JetCSV;
		vector<float>   AK8_JetL2L3corr;

		vector<float>   AK8_puppiPt;
		vector<float>   AK8_puppiMass;
		vector<float>   AK8_puppiEta;
		vector<float>   AK8_puppiPhi;
		vector<float>   AK8_puppiTau1;
		vector<float>   AK8_puppiTau2;
		vector<float>   AK8_puppiTau3;
		vector<float>   AK8_puppiTau4;
		vector<float>   AK8_puppiSDMass;
		vector<float>   AK8_puppiSDMassL2L3Corr;
		vector<float>   AK8_puppiSDL2L3Corr;

		vector<int>     AK8_JetPartonID;
		vector<int>     AK8_JetHadFlvr;
		vector<int>     AK8_JetGenJetIndex;
		vector<float>   AK8_JetGenJetEn;
		vector<float>   AK8_JetGenJetPt;
		vector<float>   AK8_JetGenJetEta;
		vector<float>   AK8_JetGenJetPhi;
		vector<int>     AK8_JetGenPartonID;
		vector<float>   AK8_JetGenEn;
		vector<float>   AK8_JetGenPt;
		vector<float>   AK8_JetGenEta;
		vector<float>   AK8_JetGenPhi;
		vector<int>     AK8_JetGenPartonMomID;

		vector<int>     n_AK8SDSJ;
		vector<vector<float> > AK8_SDSJPt;
		vector<vector<float> > AK8_SDSJEta;
		vector<vector<float> > AK8_SDSJPhi;
		vector<vector<float> > AK8_SDSJMass;
		vector<vector<float> > AK8_SDSJE;
		vector<vector<int> > AK8_SDSJCharge;
		vector<vector<int> > AK8_SDSJFlavour;
		vector<vector<float> > AK8_SDSJCSV;
		vector<int>     n_AK8puppiSDSJ;
		vector<vector<float> > AK8_puppiSDSJPt;
		vector<vector<float> > AK8_puppiSDSJEta;
		vector<vector<float> > AK8_puppiSDSJPhi;
		vector<vector<float> > AK8_puppiSDSJMass;
		vector<vector<float> > AK8_puppiSDSJE;
		vector<vector<int> > AK8_puppiSDSJCharge;
		vector<vector<int> > AK8_puppiSDSJFlavour;
		vector<vector<float> > AK8_puppiSDSJCSV;

		//--vector elements for objects-----------

		vector <int> n_ele;
		vector <int> b_jet_loose;
		vector <int> n_Mu;
		vector <int> b_jet_tight;
		vector <int> n_jet;
		vector <int> n_AK8Jet;
		vector <int> b_jet_medium;
		vector <int> forward_jet;
		vector <int> Higgs_Jet;
		vector <int> Higgs_W;
		vector <int> W_Off_Shell;
		float W_mass = 80.385 ;
		float H_mass = 125.09 ;
		float tau21 = 0.0 ;	


   TTbar_mc_bkg_PUPPI_030219(TString inputFile);
   virtual ~TTbar_mc_bkg_PUPPI_030219();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(TString OutputFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

 //..........Added Functions..............................................
		//TTree *v_tree = new TTree("v_tree","Tree with vectors");
		virtual  void        Define_tree(TTree *v_tree);
                virtual  void        Clear_Branches();
		virtual  bool       Cut_Muon(int c_muon) ;
		virtual   bool     Cut_Jet(int c_jet);
		virtual   bool     Cut_AK8Jet(int c_jet8);
		virtual   float    DeltaR_CHS_PUPPI( int chs, int puppi) ;	
		virtual   float    delta_phi(float phi1,float phi2) ;
		virtual   bool     Cut_AK8PUPPIJet(int c_jet8);
		virtual   bool     Cut_Electron(int c_ele);
		virtual   void     Fill_Histo_Ele(int k_ele );
		virtual   void     Fill_Histo_Muon(int k_mu );
		virtual   void     Fill_Histo_Jets(int l_jet);
		virtual   void     Fill_Histo_AK8Jets(int m_jet);
		virtual   void     Fill_Histo_AK8PUPPIJets(int m_jet);
                virtual   void     Fill_MCevent(int MC);
		virtual  float 	   getPUPPIweight(float puppipt, float puppieta) ;
	
		TFile  *fi ;
		TF1* puppisd_corrGEN	  ; 
		TF1* puppisd_corrRECO_cen ; 
		TF1* puppisd_corrRECO_for ; 

// ...................CHS Subjet Vectors....................
		std::vector<float> vecSDSJcsv ;
	        std::vector<float> vecSDSJpt ;
	        std::vector<float> vecSDSJeta ;
	        std::vector<float> vecSDSJmass ;
	        std::vector<float> vecSDSJphi ;
	        std::vector<float> vecSDSJe ;
	        std::vector<int > vecSDSJcharge ;
	        std::vector<int > vecSDSJflavour;
	    	    
// .................PUPPI Subjet Vectors....................
	        std::vector<float> vecPuppiSDSJcsv ;
	        std::vector<float> vecPuppiSDSJpt ;
	        std::vector<float> vecPuppiSDSJeta ;
	        std::vector<float> vecPuppiSDSJmass ;
	        std::vector<float> vecPuppiSDSJphi ;
	        std::vector<float> vecPuppiSDSJe ;
	        std::vector<int > vecPuppiSDSJcharge ;
	        std::vector<int > vecPuppiSDSJflavour;


};

#endif

#ifdef TTbar_mc_bkg_PUPPI_030219_cxx
TTbar_mc_bkg_PUPPI_030219::TTbar_mc_bkg_PUPPI_030219(TString inputFile)
{

	TChain *tree = new TChain("ggNtuplizer/EventTree");
	//        add input files here

	ifstream datafile;
	datafile.open(inputFile.Data(), ifstream::in );

	TString datafilename;

	// for signal
//	string location ="root://cmseos.fnal.gov//store/user/achhetri/Tb_tH_single_analysis/TT_bkg/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_Job_Tb_tH_MC_bkg_TTbar/";
//string location ="/uscms/homes/a/ajay/work/Tprime/CMSSW_9_4_13/src/ggAnalysis/ggNtuplizer/test/";
string location ="";


	//string location ="/uscms_data/d3/achhetri/Tb_tH_DeepCSV/CMSSW_8_0_26_patch1/src/Tb_tH_lvqqbqq/Tb_tH_LH_1500/";
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

TTbar_mc_bkg_PUPPI_030219::~TTbar_mc_bkg_PUPPI_030219()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// =========User Added Functions ==============================

void TTbar_mc_bkg_PUPPI_030219::Define_tree(TTree *v_tree)
{
	v_tree->Branch("v_event" ,                        &v_event);
	v_tree->Branch("n_Vtx",                           &n_Vtx);
	v_tree->Branch("n_GoodVtx" ,                      &n_GoodVtx);
	v_tree->Branch("n_TrksPV" ,                       &n_TrksPV);
	v_tree->Branch("is_PVGood" ,                      &is_PVGood);
	v_tree->Branch("v_vtx" ,                          &v_vtx);
	v_tree->Branch("v_vty" ,                          &v_vty);
	v_tree->Branch("v_vtz" ,                          &v_vtz);
	v_tree->Branch("rho_Central" ,                    &rho_Central);

	// 1. MC variables 

	v_tree->Branch("n_MC" ,                          &n_MC);
	v_tree->Branch("mc_PID" ,                        &mc_PID);
	v_tree->Branch("mc_Vtx" ,                        &mc_Vtx);
	v_tree->Branch("mc_Vty" ,                        &mc_Vty);
	v_tree->Branch("mc_Vtz" ,                        &mc_Vtz);
	v_tree->Branch("mc_Pt" ,                         &mc_Pt);
	v_tree->Branch("mc_Mass" ,                       &mc_Mass);
	v_tree->Branch("mc_Eta" ,                        &mc_Eta);
	v_tree->Branch("mc_Phi" ,                        &mc_Phi);
	v_tree->Branch("mc_E" ,                          &mc_E);
	v_tree->Branch("mc_Et" ,                         &mc_Et);
	v_tree->Branch("mc_GMomPID" ,                    &mc_GMomPID);
	v_tree->Branch("mc_MomPID" ,                     &mc_MomPID);
	v_tree->Branch("mc_MomPt" ,                      &mc_MomPt);
	v_tree->Branch("mc_MomMass" ,                    &mc_MomMass);
	v_tree->Branch("mc_MomEta" ,                     &mc_MomEta);
	v_tree->Branch("mc_MomPhi" ,                     &mc_MomPhi);
	v_tree->Branch("mc_StatusFlag" ,                 &mc_StatusFlag);
	v_tree->Branch("mc_Parentage" ,                  &mc_Parentage);
	v_tree->Branch("mc_Status" ,                     &mc_Status);
	v_tree->Branch("mc_CalIsoDR03" ,                 &mc_CalIsoDR03);
	v_tree->Branch("mc_TrkIsoDR03" ,                 &mc_TrkIsoDR03);
	v_tree->Branch("mc_CalIsoDR04" ,                 &mc_CalIsoDR04);
	v_tree->Branch("mc_TrkIsoDR04" ,                 &mc_TrkIsoDR04);
	v_tree->Branch("gen_MET" ,                       &gen_MET);
	v_tree->Branch("gen_METPhi" ,                    &gen_METPhi);


	//2. MET Variables
	v_tree->Branch("pf_MET" ,                         &pf_MET);         
	v_tree->Branch("pf_METPhi",                       &pf_METPhi);
	v_tree->Branch("pf_METsumEt",                     &pf_METsumEt);
	v_tree->Branch("pf_METmEtSig",                    &pf_METmEtSig);
	//v_tree->Branch("pf_METSig",                       &pf_METSig);



	//3. Electron 
	v_tree->Branch("n_Ele",                           &n_Ele);
	v_tree->Branch("ele_Charge",                      &ele_Charge);
	v_tree->Branch("ele_En",                          &ele_En);
	v_tree->Branch("ele_D0",                          &ele_D0);
	v_tree->Branch("ele_Dz",                          &ele_Dz);
	v_tree->Branch("ele_Pt",                          &ele_Pt);
	v_tree->Branch("ele_Eta",                         &ele_Eta);
	v_tree->Branch("ele_Phi",                         &ele_Phi);
	v_tree->Branch("ele_R9",                          &ele_R9);
	v_tree->Branch("ele_SCEta",                       &ele_SCEta);
	v_tree->Branch("ele_SCPhi",                       &ele_SCPhi);		
	v_tree->Branch("ele_HoverE",                      &ele_HoverE);
	v_tree->Branch("ele_EoverP",                      &ele_EoverP);
	v_tree->Branch("ele_EoverPout",                   &ele_EoverPout);
	v_tree->Branch("ele_EoverPInv",                   &ele_EoverPInv);
	v_tree->Branch("ele_dEtaAtVtx",                   &ele_dEtaAtVtx);
	v_tree->Branch("ele_dPhiAtVtx",                   &ele_dPhiAtVtx);

	v_tree->Branch("ele_SigmaIEtaIEtaFull5x5",        &ele_SigmaIEtaIEtaFull5x5);
	v_tree->Branch("ele_SigmaIPhiIPhiFull5x5",        &ele_SigmaIPhiIPhiFull5x5);
	v_tree->Branch("ele_ConvVeto",                    &ele_ConvVeto);
	v_tree->Branch("ele_MissHits",                    &ele_MissHits);

	v_tree->Branch("ele_PFChIso",                     &ele_PFChIso);
	v_tree->Branch("ele_PFPhoIso",                    &ele_PFPhoIso);
	v_tree->Branch("ele_PFNeuIso",                    &ele_PFNeuIso);
	v_tree->Branch("ele_PFMiniIso" ,                  &ele_PFMiniIso);
	v_tree->Branch("ele_dEtaseedAtVtx",               &ele_dEtaseedAtVtx);


	//4. Muon 
	v_tree->Branch("N_Mu",           &N_Mu);
	v_tree->Branch("mu_Pt",          &mu_Pt);
	v_tree->Branch("mu_En",          &mu_En);
	v_tree->Branch("mu_Eta",         &mu_Eta);
	v_tree->Branch("mu_Phi",         &mu_Phi);
	v_tree->Branch("mu_Charge",      &mu_Charge);

	v_tree->Branch("mu_IDbit",       &mu_IDbit);
	v_tree->Branch("mu_D0",          &mu_D0);
	v_tree->Branch("mu_Dz",          &mu_Dz);

	v_tree->Branch("mu_Chi2NDF",     &mu_Chi2NDF);
	v_tree->Branch("mu_InnerD0",     &mu_InnerD0);
	v_tree->Branch("mu_InnerDz",     &mu_InnerDz);
	v_tree->Branch("mu_InnervalidFraction",   &mu_InnervalidFraction);
	v_tree->Branch("mu_segmentCompatibility", &mu_segmentCompatibility);
	v_tree->Branch("mu_chi2LocalPosition",    &mu_chi2LocalPosition);
	v_tree->Branch("mu_trkKink",              &mu_trkKink);

	v_tree->Branch("mu_PFChIso",     &mu_PFChIso);
	v_tree->Branch("mu_PFPhoIso",    &mu_PFPhoIso);
	v_tree->Branch("mu_PFNeuIso",    &mu_PFNeuIso);
	v_tree->Branch("mu_PFMiniIso",   &mu_PFMiniIso);

	v_tree->Branch("mu_TrkLayers" ,                          &mu_TrkLayers);
	v_tree->Branch("mu_PixelLayers" ,                         &mu_PixelLayers);
	v_tree->Branch("mu_PixelHits" ,                            &mu_PixelHits);
	v_tree->Branch("mu_MuonHits" ,                          &mu_MuonHits);
	v_tree->Branch("mu_Stations" ,                              &mu_Stations);
	v_tree->Branch("mu_Matches" ,                              &mu_Matches);
	v_tree->Branch("mu_TrkQuality" ,                         &mu_TrkQuality);
	v_tree->Branch("mu_IsoTrk" ,                              &mu_IsoTrk);

	//4. Bjet 


	v_tree->Branch("n_Jet",            &n_Jet);
	v_tree->Branch("jet_Pt",           &jet_Pt);
	v_tree->Branch("jet_En",           &jet_En);
	v_tree->Branch("jet_Eta",          &jet_Eta);
	v_tree->Branch("jet_Phi",          &jet_Phi);
	v_tree->Branch("jet_Area",         &jet_Area);
	v_tree->Branch("jet_Mt",           &jet_Mt);

	v_tree->Branch("jet_CSV2BJetTags",              &jet_CSV2BJetTags);
	v_tree->Branch("jet_JetProbabilityBJetTags",    &jet_JetProbabilityBJetTags);
	v_tree->Branch("jet_pfCombinedMVAV2BJetTags",   &jet_pfCombinedMVAV2BJetTags);
	v_tree->Branch("jet_DeepCSVTags_b",             &jet_DeepCSVTags_b);
	v_tree->Branch("jet_DeepCSVTags_bb",            &jet_DeepCSVTags_bb);
	v_tree->Branch("jet_DeepCSVTags_c",             &jet_DeepCSVTags_c);
	v_tree->Branch("jet_DeepCSVTags_cc",            &jet_DeepCSVTags_cc);
	v_tree->Branch("jet_DeepCSVTags_udsg",          &jet_DeepCSVTags_udsg);

	v_tree->Branch("jet_PartonID",                  &jet_PartonID);
	v_tree->Branch("jet_HadFlvr",     	        &jet_HadFlvr);
	v_tree->Branch("jet_GenJetEn",                	&jet_GenJetEn);
	v_tree->Branch("jet_GenJetPt",                	&jet_GenJetPt);
	v_tree->Branch("jet_GenJetEta",                	&jet_GenJetEta);
	v_tree->Branch("jet_GenJetPhi",                	&jet_GenJetPhi);
	v_tree->Branch("jet_GenPartonID",               &jet_GenPartonID);
	v_tree->Branch("jet_GenEn",                	&jet_GenEn);
	v_tree->Branch("jet_GenPt",                	&jet_GenPt);
	v_tree->Branch("jet_GenEta",                	&jet_GenEta);
	v_tree->Branch("jet_GenPhi",                	&jet_GenPhi);
	v_tree->Branch("jet_GenPartonMomID",            &jet_GenPartonMomID);

	v_tree->Branch("jet_PFLooseId", &jet_PFLooseId);
	v_tree->Branch("jet_ID",        &jet_ID);
	v_tree->Branch("jet_PUID",      &jet_PUID);
	v_tree->Branch("jet_PUFullID",  &jet_PUFullID);
	v_tree->Branch("jet_CHF",       &jet_CHF);
	v_tree->Branch("jet_NHF",       &jet_NHF);
	v_tree->Branch("jet_CEF",       &jet_CEF);
	v_tree->Branch("jet_NEF",       &jet_NEF);
	v_tree->Branch("jet_NCH",       &jet_NCH);
	v_tree->Branch("jet_NNP",       &jet_NNP);
	v_tree->Branch("jet_MUF",       &jet_MUF);


	//5. Fat jets

	v_tree->Branch("N_AK8Jet",                  &N_AK8Jet);
	v_tree->Branch("n_AK8Jetpuppi",             &n_AK8Jetpuppi);
	v_tree->Branch("AK8_JetPt",                 &AK8_JetPt);
	v_tree->Branch("AK8_JetEn",                 &AK8_JetEn);
	v_tree->Branch("AK8_JetEta",                &AK8_JetEta);
	v_tree->Branch("AK8_JetPhi",                &AK8_JetPhi);
	v_tree->Branch("AK8_JetMass",               &AK8_JetMass);

	v_tree->Branch("AK8_Jet_tau1",              &AK8_Jet_tau1);
	v_tree->Branch("AK8_Jet_tau2",              &AK8_Jet_tau2);
	v_tree->Branch("AK8_Jet_tau3",              &AK8_Jet_tau3);
	v_tree->Branch("AK8_Jet_CHStau1",              &AK8_Jet_CHStau1);
	v_tree->Branch("AK8_Jet_CHStau2",              &AK8_Jet_CHStau2);
	v_tree->Branch("AK8_Jet_CHStau3",              &AK8_Jet_CHStau3);
	v_tree->Branch("AK8_Jet_CHStau4",              &AK8_Jet_CHStau4);
	v_tree->Branch("AK8_JetCHF",                &AK8_JetCHF);
	v_tree->Branch("AK8_JetNHF",                &AK8_JetNHF);
	v_tree->Branch("AK8_JetCEF",                &AK8_JetCEF);
	v_tree->Branch("AK8_JetNEF",                &AK8_JetNEF);
	v_tree->Branch("AK8_JetNCH",                &AK8_JetNCH);
	v_tree->Branch("AK8_JetNNP",                &AK8_JetNNP);
	v_tree->Branch("AK8_JetMUF",                &AK8_JetMUF);

	v_tree->Branch("AK8_Jetnconstituents",      &AK8_Jetnconstituents);
	v_tree->Branch("AK8_JetPFLooseId",          &AK8_JetPFLooseId);
	v_tree->Branch("AK8_JetPFTightLepVetoId",   &AK8_JetPFTightLepVetoId);
	v_tree->Branch("AK8_JetSoftDropMass",       &AK8_JetSoftDropMass);

	v_tree->Branch("AK8_JetSoftDropMassCorr", 		&AK8_JetSoftDropMassCorr ) ;
	v_tree->Branch("AK8_JetPrunedMassCorr",			&AK8_JetPrunedMassCorr ) ; 
	v_tree->Branch("AK8_JetL2L3corr",   			&AK8_JetL2L3corr ) ; 
	v_tree->Branch("AK8_puppiSDMassL2L3Corr",   		&AK8_puppiSDMassL2L3Corr ) ; 
	v_tree->Branch("AK8_puppiSDL2L3Corr",       		&AK8_puppiSDL2L3Corr ) ;
 
	v_tree->Branch("AK8_JetPrunedMass",         &AK8_JetPrunedMass);
	v_tree->Branch("AK8_JetpfBoostedDSVBTag",   &AK8_JetpfBoostedDSVBTag);
	v_tree->Branch("AK8_JetCSV",                &AK8_JetCSV);
	v_tree->Branch("AK8_JetDSVnewV4",           &AK8_JetDSVnewV4);

	v_tree->Branch("AK8_puppiPt",               &AK8_puppiPt);
	v_tree->Branch("AK8_puppiMass",             &AK8_puppiMass);
	v_tree->Branch("AK8_puppiEta",              &AK8_puppiEta);
	v_tree->Branch("AK8_puppiPhi",              &AK8_puppiPhi);
	v_tree->Branch("AK8_puppiTau1",             &AK8_puppiTau1);
	v_tree->Branch("AK8_puppiTau2",             &AK8_puppiTau2);
	v_tree->Branch("AK8_puppiTau3",             &AK8_puppiTau3);
	v_tree->Branch("AK8_puppiTau4",             &AK8_puppiTau4);
	v_tree->Branch("AK8_puppiSDMass",           &AK8_puppiSDMass);

	v_tree->Branch("AK8_JetPartonID",                	&AK8_JetPartonID);
	v_tree->Branch("AK8_JetHadFlvr",                	&AK8_JetHadFlvr);
	v_tree->Branch("AK8_JetGenJetIndex",                	&AK8_JetGenJetIndex);
	v_tree->Branch("AK8_JetGenJetEn",                	&AK8_JetGenJetEn);
	v_tree->Branch("AK8_JetGenJetPt",                	&AK8_JetGenJetPt);
	v_tree->Branch("AK8_JetGenJetEta",                	&AK8_JetGenJetEta);
	v_tree->Branch("AK8_JetGenJetPhi",                	&AK8_JetGenJetPhi);
	v_tree->Branch("AK8_JetGenPartonID",              	&AK8_JetGenPartonID);
	v_tree->Branch("AK8_JetGenEn",                		&AK8_JetGenEn);
	v_tree->Branch("AK8_JetGenPt",                		&AK8_JetGenPt);
	v_tree->Branch("AK8_JetGenEta",                		&AK8_JetGenEta);
	v_tree->Branch("AK8_JetGenPhi",                		&AK8_JetGenPhi);
	v_tree->Branch("AK8_JetGenPartonMomID",          	&AK8_JetGenPartonMomID);

	v_tree->Branch("n_AK8SDSJ",                		&n_AK8SDSJ);

	v_tree->Branch("AK8_SDSJPt",                		&AK8_SDSJPt);
	v_tree->Branch("AK8_SDSJEta",                		&AK8_SDSJEta);
	v_tree->Branch("AK8_SDSJPhi",                		&AK8_SDSJPhi);
	v_tree->Branch("AK8_SDSJMass",                		&AK8_SDSJMass);
	v_tree->Branch("AK8_SDSJE",                		&AK8_SDSJE);
	v_tree->Branch("AK8_SDSJCharge",                	&AK8_SDSJCharge);
	v_tree->Branch("AK8_SDSJFlavour",                	&AK8_SDSJFlavour);
	v_tree->Branch("AK8_SDSJCSV",                		&AK8_SDSJCSV);

	v_tree->Branch("n_AK8puppiSDSJ",                	&n_AK8puppiSDSJ);

	v_tree->Branch("AK8_puppiSDSJPt",                	&AK8_puppiSDSJPt);
	v_tree->Branch("AK8_puppiSDSJEta",                	&AK8_puppiSDSJEta);
	v_tree->Branch("AK8_puppiSDSJPhi",                	&AK8_puppiSDSJPhi);
	v_tree->Branch("AK8_puppiSDSJMass",                	&AK8_puppiSDSJMass);
	v_tree->Branch("AK8_puppiSDSJE",                	&AK8_puppiSDSJE);
	v_tree->Branch("AK8_puppiSDSJCharge",                	&AK8_puppiSDSJCharge);
	v_tree->Branch("AK8_puppiSDSJFlavour",                	&AK8_puppiSDSJFlavour);
	v_tree->Branch("AK8_puppiSDSJCSV",                	&AK8_puppiSDSJCSV);


}



void TTbar_mc_bkg_PUPPI_030219::Clear_Branches()
{


	mc_PID                          .clear();
	mc_Vtx                          .clear();
	mc_Vty                          .clear();
	mc_Vtz                          .clear();
	mc_Pt                         	.clear();
	mc_Mass                         .clear();
	mc_Eta                          .clear();
	mc_Phi                          .clear();
	mc_E                          	.clear();
	mc_Et                          	.clear();
	mc_GMomPID                      .clear();
	mc_MomPID                       .clear();
	mc_MomPt                        .clear();
	mc_MomMass                      .clear();
	mc_MomEta                       .clear();
	mc_MomPhi                       .clear();
	mc_StatusFlag                   .clear();
	mc_Parentage                    .clear();
	mc_Status                       .clear();
	mc_CalIsoDR03                   .clear();
	mc_TrkIsoDR03                   .clear();
	mc_CalIsoDR04                   .clear();
	mc_TrkIsoDR04                   .clear();


	//3. Electron 

	ele_Charge                          	.clear();
	ele_En                          	.clear();
	ele_D0                          	.clear();
	ele_Dz                          	.clear();
	ele_Pt                          	.clear();
	ele_Eta                          	.clear();
	ele_Phi                          	.clear();
	ele_R9                          	.clear();
	ele_SCEta                          	.clear();
	ele_SCPhi                          	.clear();		
	ele_HoverE                          	.clear();
	ele_EoverP                          	.clear();
	ele_EoverPout                          	.clear();
	ele_EoverPInv                           .clear();
	ele_dEtaAtVtx                           .clear();
	ele_dPhiAtVtx                           .clear();

	ele_SigmaIEtaIEtaFull5x5          	.clear();
	ele_SigmaIPhiIPhiFull5x5          	.clear();
	ele_ConvVeto                          	.clear();
	ele_MissHits                          	.clear();

	ele_PFChIso                          	.clear();
	ele_PFPhoIso                          	.clear();
	ele_PFNeuIso                          	.clear();
	ele_PFMiniIso                          	.clear();
	ele_dEtaseedAtVtx                       .clear(); 


	//4. Muon 

	mu_Pt                          		.clear();
	mu_En                          		.clear();
	mu_Eta                          	.clear();
	mu_Phi                          	.clear();
	mu_Charge                          	.clear();
	mu_IDbit                          	.clear();
	mu_D0                          		.clear();
	mu_Dz                          		.clear();		
	mu_Chi2NDF                          	.clear();
	mu_InnerD0                          	.clear();
	mu_InnerDz                          	.clear();

	mu_PFChIso                          	.clear();
	mu_PFPhoIso                          	.clear();
	mu_PFNeuIso                          	.clear();
	mu_PFMiniIso                          	.clear();

	mu_InnervalidFraction             	.clear();
	mu_segmentCompatibility        		.clear();
	mu_chi2LocalPosition           		.clear();
	mu_trkKink                          	.clear();

	mu_TrkLayers                            .clear();
	mu_PixelLayers                          .clear();
	mu_PixelHits                             .clear();
	mu_MuonHits                            .clear();
	mu_Stations                              .clear();
	mu_Matches                              .clear();
	mu_TrkQuality                           .clear();
	mu_IsoTrk                                .clear();

	//6. Bjet 

	jet_Pt                          	.clear();
	jet_En                          	.clear();
	jet_Eta                          	.clear();
	jet_Phi                          	.clear();
	jet_Area                          	.clear();
	jet_Mt                          	.clear();

	jet_CSV2BJetTags                        .clear();
	jet_JetProbabilityBJetTags      	.clear();
	jet_pfCombinedMVAV2BJetTags      	.clear();
	jet_DeepCSVTags_b      	 		.clear();
	jet_DeepCSVTags_bb       	 	.clear();
	jet_DeepCSVTags_c           		.clear();
	jet_DeepCSVTags_cc        	 	.clear();
	jet_DeepCSVTags_udsg   	 		.clear();

	jet_PartonID                          	.clear();
	jet_HadFlvr                          	.clear();
	jet_GenJetEn                          	.clear();
	jet_GenJetPt                          	.clear();
	jet_GenJetEta                          	.clear();
	jet_GenJetPhi                          	.clear();
	jet_GenPartonID                         .clear();
	jet_GenEn                          	.clear();
	jet_GenPt                          	.clear();
	jet_GenEta                          	.clear();
	jet_GenPhi                          	.clear();
	jet_GenPartonMomID                      .clear();

	jet_PFLooseId                          	.clear();
	jet_ID                          	.clear();
	jet_PUID                          	.clear();
	jet_PUFullID                          	.clear();
	jet_CHF                          	.clear();
	jet_NHF                          	.clear();
	jet_CEF                          	.clear();
	jet_NEF                          	.clear();
	jet_NCH                          	.clear();
	jet_NNP                          	.clear();
	jet_MUF                          	.clear();

	//7. Fat jets

	AK8_JetPt                          	.clear();
	AK8_JetEn                          	.clear();

	AK8_JetEta                          	.clear();
	AK8_JetPhi                          	.clear();
	AK8_JetMass                          	.clear();
	AK8_Jet_tau1                          	.clear();
	AK8_Jet_tau2                          	.clear();
	AK8_Jet_tau3                          	.clear();
	AK8_Jet_CHStau1                          	.clear();
	AK8_Jet_CHStau2                          	.clear();
	AK8_Jet_CHStau3                          	.clear();
	AK8_Jet_CHStau4                          	.clear();
	AK8_JetCHF                          	.clear();
	AK8_JetNHF                          	.clear();
	AK8_JetCEF                          	.clear();
	AK8_JetNEF                          	.clear();
	AK8_JetNCH                          	.clear();
	AK8_JetNNP                          	.clear();
	AK8_JetMUF                          	.clear();
	AK8_Jetnconstituents      	 	.clear();
	AK8_JetPFLooseId         		.clear();
	AK8_JetPFTightLepVetoId             	.clear();
	AK8_JetSoftDropMass         		.clear();
	AK8_JetSoftDropMassCorr 		.clear();
	AK8_JetPrunedMassCorr			 .clear(); 
	AK8_JetL2L3corr   			 .clear(); 
	AK8_puppiSDMassL2L3Corr   		 .clear(); 
	AK8_puppiSDL2L3Corr       		 .clear(); 

	AK8_JetPrunedMass          		.clear();
	AK8_JetpfBoostedDSVBTag        		.clear();
	AK8_JetDSVnewV4                         .clear();
	AK8_JetCSV                          	.clear();

	AK8_puppiPt                          	.clear();
	AK8_puppiMass                          	.clear();
	AK8_puppiEta                          	.clear();
	AK8_puppiPhi                          	.clear();
	AK8_puppiTau1                          	.clear();
	AK8_puppiTau2                          	.clear();
	AK8_puppiTau3                          	.clear();
	AK8_puppiTau4                          	.clear();
	AK8_puppiSDMass                         .clear();

	AK8_JetPartonID                         .clear();
	AK8_JetHadFlvr                          .clear();
	AK8_JetGenJetIndex                      .clear();
	AK8_JetGenJetEn                         .clear();
	AK8_JetGenJetPt                         .clear();
	AK8_JetGenJetEta                        .clear();
	AK8_JetGenJetPhi                        .clear();
	AK8_JetGenPartonID                      .clear();
	AK8_JetGenEn                          	.clear();
	AK8_JetGenPt                          	.clear();
	AK8_JetGenEta                          	.clear();
	AK8_JetGenPhi                          	.clear();
	AK8_JetGenPartonMomID                   .clear();

	n_AK8SDSJ                          	.clear();
	AK8_SDSJPt                          	.clear();
	AK8_SDSJEta                          	.clear();
	AK8_SDSJPhi                          	.clear();
	AK8_SDSJMass                          	.clear();
	AK8_SDSJE                          	.clear();
	AK8_SDSJCharge                          .clear();
	AK8_SDSJFlavour                         .clear();
	AK8_SDSJCSV                          	.clear();

	n_AK8puppiSDSJ                          .clear();
	AK8_puppiSDSJPt                         .clear();
	AK8_puppiSDSJEta                        .clear();
	AK8_puppiSDSJPhi                        .clear();
	AK8_puppiSDSJMass                       .clear();
	AK8_puppiSDSJE                          .clear();
	AK8_puppiSDSJCharge                     .clear();
	AK8_puppiSDSJFlavour                    .clear();
	AK8_puppiSDSJCSV                        .clear();


}

//=================Object Selection Functions ==========================

bool TTbar_mc_bkg_PUPPI_030219::Cut_Muon(int c_muon){

	bool pass_muon = true;

	if ((*muPt)[c_muon  ] <= 25.0) pass_muon = false;
	if (fabs((*muEta)[c_muon  ]) >= 3.0) pass_muon = false;

	return pass_muon;
}

bool TTbar_mc_bkg_PUPPI_030219::Cut_Jet(int c_jet){

	bool pass_jet = true;
	if((*jetPt)[c_jet] <= 30.0) pass_jet = false;
	if( fabs((*jetEta)[c_jet] ) >= 5.0) pass_jet = false;
	if (!((*jetPFLooseId)[c_jet ])) pass_jet = false;
	return pass_jet;
}

bool TTbar_mc_bkg_PUPPI_030219::Cut_AK8Jet(int c_jet8){

	bool pass_AK8jet = true;
	if((*AK8JetPt)[c_jet8] <= 190.0) pass_AK8jet = false;
	if(! (*AK8JetPFLooseId)[c_jet8]  ) pass_AK8jet = false;
	if( fabs((*AK8JetEta)[c_jet8] ) >= 5.0) pass_AK8jet = false;

	return pass_AK8jet ;
}

bool TTbar_mc_bkg_PUPPI_030219::Cut_AK8PUPPIJet(int c_jet8){

	bool pass_AK8jet = true;
	if((*AK8JetPt)[c_jet8] <= 150.0) pass_AK8jet = false;
	if(! (*AK8JetPFLooseId)[c_jet8]  ) pass_AK8jet = false;
	if( fabs((*AK8JetEta)[c_jet8] ) >= 5.0) pass_AK8jet = false;

	return pass_AK8jet ;
}


float TTbar_mc_bkg_PUPPI_030219::DeltaR_CHS_PUPPI( int chs, int puppi) 
{
	float dEta = (*AK8JetEta)[chs] - (*AK8puppiEta)[puppi] ; 
	float dPhi = delta_phi((*AK8JetPhi)[chs] , (*AK8puppiPhi)[puppi]) ;
	float dR   = TMath::Sqrt(dPhi*dPhi + dEta*dEta) ;

	return dR ;
}
float TTbar_mc_bkg_PUPPI_030219::delta_phi(float phi1,float phi2)
{
	const float PI=2.0*acos(0.);
	const float TWOPI=2.0*PI;
	float PHI=fabs( phi1 - phi2 ) ;
	return (PHI<=PI)? PHI : TWOPI-PHI;

}



bool TTbar_mc_bkg_PUPPI_030219::Cut_Electron(int c_ele){

	bool pass_ele1 = true ;

	if((*elePt)[c_ele ] <= 30) pass_ele1 = false ;
	if(fabs((*eleSCEta)[c_ele ])>= 2.5) pass_ele1 = false ;


	if (fabs((*eleSCEta)[c_ele ]) <= 1.479){
		//if(fabs((*eledEtaseedAtVtx)[c_ele]) >=  0.00477 ) pass_ele1 = false ;
		if(fabs((*eledPhiAtVtx)[c_ele]) >=  0.222  ) pass_ele1 = false;
		if( (*eleSigmaIEtaIEtaFull5x5)[c_ele] >=  0.011) pass_ele1 = false;
		if( (*eleHoverE)[c_ele] >=  0.298  ) pass_ele1 = false;
		if (fabs((*eleEoverPInv)[c_ele]) >=   0.241  ) pass_ele1 = false;
		if ((*eleConvVeto)[c_ele] == 0 ) pass_ele1 = false;
		if ((*eleMissHits)[c_ele] > 1 ) pass_ele1 = false;
	}

	else {

		//if(fabs((*eledEtaseedAtVtx)[c_ele]) >= 0.00868 ) pass_ele1 = false ;
		if(fabs((*eledPhiAtVtx)[c_ele]) >=  0.213 ) pass_ele1 = false;
		if( (*eleSigmaIEtaIEtaFull5x5)[c_ele] >=  0.0314  ) pass_ele1 = false;
		if( (*eleHoverE)[c_ele] >=  0.101 ) pass_ele1 = false;
		if (fabs((*eleEoverPInv)[c_ele]) >=   0.14 ) pass_ele1 = false;
		if ((*eleConvVeto)[c_ele] == 0 ) pass_ele1 = false;
		if ((*eleMissHits)[c_ele] > 1 ) pass_ele1 = false;

	}


	return pass_ele1 ;

}

float  TTbar_mc_bkg_PUPPI_030219::getPUPPIweight(float puppipt, float puppieta ){
//TFile *file = TFile::Open( "puppiCorr.root","READ");
//TFile* file = TFile::Open( puppiFile_.c_str(),"READ" );
float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  
  totalWeight = genCorr * recoCorr;

  return totalWeight;

}


// ============Filling Elements =========================================

void  TTbar_mc_bkg_PUPPI_030219::Fill_Histo_Ele(int k_ele) {

	n_Ele ++ ; // not forget to assign 0 value at start of each event


	ele_Charge                 .push_back((*eleCharge)[k_ele]);
	ele_En                     .push_back((*eleEn)[k_ele]);

	ele_D0                     .push_back((*eleD0)[k_ele]);
	ele_Dz                     .push_back((*eleDz)[k_ele]);

	ele_Pt                     .push_back((*elePt)[k_ele]);
	ele_Eta                    .push_back((*eleEta)[k_ele]);
	ele_Phi                    .push_back((*elePhi)[k_ele]);
	ele_R9                     .push_back((*eleR9)[k_ele]);
	ele_SCEta                  .push_back((*eleSCEta)[k_ele]);
	ele_SCPhi                  .push_back((*eleSCPhi)[k_ele]);
	ele_HoverE                 .push_back((*eleHoverE)[k_ele]);
	ele_EoverP                 .push_back((*eleEoverP)[k_ele]);
	ele_EoverPout              .push_back((*eleEoverPout)[k_ele]);
	ele_EoverPInv              .push_back((*eleEoverPInv)[k_ele]);
	ele_dEtaAtVtx              .push_back((*eledEtaAtVtx)[k_ele]);
	ele_dPhiAtVtx              .push_back((*eledPhiAtVtx)[k_ele]);

	ele_SigmaIEtaIEtaFull5x5   .push_back((*eleSigmaIEtaIEtaFull5x5)[k_ele]);
	ele_SigmaIPhiIPhiFull5x5   .push_back((*eleSigmaIPhiIPhiFull5x5)[k_ele]);
	ele_ConvVeto               .push_back((*eleConvVeto)[k_ele]);
	ele_MissHits               .push_back((*eleMissHits)[k_ele]);
	ele_PFChIso                .push_back((*elePFChIso)[k_ele]);
	ele_PFPhoIso               .push_back((*elePFPhoIso)[k_ele]);
	ele_PFNeuIso               .push_back((*elePFNeuIso)[k_ele]);
	//ele_PFMiniIso              .push_back((*elePFMiniIso)[k_ele]);
	//ele_dEtaseedAtVtx          .push_back((*eledEtaseedAtVtx)[k_ele]);

	return ;
}


void TTbar_mc_bkg_PUPPI_030219::Fill_Histo_Muon(int k_mu)
{

	N_Mu ++ ;

	mu_Pt        .push_back((*muPt)[k_mu]);
	mu_En        .push_back((*muEn)[k_mu]);
	mu_Eta       .push_back((*muEta)[k_mu]);
	mu_Phi       .push_back((*muPhi)[k_mu]);
	mu_Charge    .push_back((*muCharge)[k_mu]);

	mu_IDbit     .push_back((*muIDbit)[k_mu]);
	mu_D0        .push_back((*muD0)[k_mu]);
	mu_Dz        .push_back((*muDz)[k_mu]);

	mu_Chi2NDF   .push_back((*muChi2NDF)[k_mu]);
	mu_InnerD0   .push_back((*muInnerD0)[k_mu]);
	mu_InnerDz   .push_back((*muInnerDz)[k_mu]);
	mu_InnervalidFraction    .push_back((*muInnervalidFraction)[k_mu]);
	mu_segmentCompatibility .push_back((*musegmentCompatibility)[k_mu]);
	mu_chi2LocalPosition    .push_back((*muchi2LocalPosition)[k_mu]);
	mu_trkKink              .push_back((*mutrkKink)[k_mu]);
	mu_PFChIso   .push_back((*muPFChIso)[k_mu]);
	mu_PFPhoIso  .push_back((*muPFPhoIso)[k_mu]);
	mu_PFNeuIso  .push_back((*muPFNeuIso)[k_mu]);
	//mu_PFMiniIso .push_back((*muPFMiniIso)[k_mu]);

	mu_TrkLayers                            .push_back((*muTrkLayers)[k_mu]);
	mu_PixelLayers                          .push_back((*muPixelLayers)[k_mu]);
	mu_PixelHits                             .push_back((*muPixelHits)[k_mu]);
	mu_MuonHits                            .push_back((*muMuonHits)[k_mu]);
	mu_Stations                              .push_back((*muStations)[k_mu]);
	mu_Matches                              .push_back((*muMatches)[k_mu]);
	mu_TrkQuality                           .push_back((*muTrkQuality)[k_mu]);
	mu_IsoTrk                                .push_back((*muIsoTrk)[k_mu]);


}


void TTbar_mc_bkg_PUPPI_030219::Fill_Histo_AK8Jets(int m_jet)
{

	//===================for CHS jets=======================================

	N_AK8Jet ++;
	n_AK8Jetpuppi ++;

	AK8_JetPt          .push_back((*AK8JetPt)[m_jet]);
	AK8_JetEn          .push_back((*AK8JetEn)[m_jet]);
	AK8_JetEta         .push_back((*AK8JetEta)[m_jet]);
	AK8_JetPhi         .push_back((*AK8JetPhi)[m_jet]);
	AK8_JetMass        .push_back((*AK8JetMass)[m_jet]);
	AK8_Jet_tau1       .push_back((*AK8Jet_tau1)[m_jet]);
	AK8_Jet_tau2       .push_back((*AK8Jet_tau2)[m_jet]);
	AK8_Jet_tau3       .push_back((*AK8Jet_tau3)[m_jet]);
	//AK8_Jet_CHStau1       .push_back((*AK8Jet_CHStau1)[m_jet]);
	//AK8_Jet_CHStau2       .push_back((*AK8Jet_CHStau2)[m_jet]);
	//AK8_Jet_CHStau3       .push_back((*AK8Jet_CHStau3)[m_jet]);
	//AK8_Jet_CHStau4       .push_back((*AK8Jet_CHStau4)[m_jet]);
	AK8_JetCHF         .push_back((*AK8JetCHF)[m_jet]);
	AK8_JetNHF         .push_back((*AK8JetNHF)[m_jet]);
	AK8_JetCEF         .push_back((*AK8JetCEF)[m_jet]);
	AK8_JetNEF         .push_back((*AK8JetNEF)[m_jet]);
	AK8_JetNCH         .push_back((*AK8JetNCH)[m_jet]);
	AK8_JetNNP         .push_back((*AK8JetNNP)[m_jet]);
	AK8_Jetnconstituents .push_back((*AK8Jetnconstituents)[m_jet]);
	AK8_JetMUF         .push_back((*AK8JetMUF)[m_jet]);
	AK8_JetPFLooseId        .push_back((*AK8JetPFLooseId)[m_jet]);
	AK8_JetPFTightLepVetoId .push_back((*AK8JetPFTightLepVetoId)[m_jet]);

	AK8_JetSoftDropMass  	 .push_back((*AK8JetSoftDropMass)[m_jet]);
	AK8_JetSoftDropMassCorr  .push_back((*AK8JetSoftDropMassCorr)[m_jet]);
	AK8_JetPrunedMass 	 .push_back((*AK8JetPrunedMass)[m_jet]);
	AK8_JetPrunedMassCorr	 .push_back((*AK8JetPrunedMassCorr)[m_jet]);
	AK8_JetpfBoostedDSVBTag  .push_back((*AK8JetpfBoostedDSVBTag)[m_jet]);
	AK8_JetCSV         .push_back((*AK8JetCSV)[m_jet]);
	AK8_JetL2L3corr    .push_back((*AK8JetL2L3corr)[m_jet]); ;

	AK8_JetDSVnewV4                         .push_back((*AK8JetDSVnewV4)[m_jet]);

	AK8_JetPartonID                         .push_back((*AK8JetPartonID )[m_jet]);
	AK8_JetHadFlvr                          .push_back((*AK8JetHadFlvr )[m_jet]);
	AK8_JetGenJetIndex                      .push_back((*AK8JetGenJetIndex )[m_jet]);
	AK8_JetGenJetEn                         .push_back((*AK8JetGenJetEn )[m_jet]);
	AK8_JetGenJetPt                         .push_back((*AK8JetGenJetPt )[m_jet]);
	AK8_JetGenJetEta                        .push_back((*AK8JetGenJetEta )[m_jet]);
	AK8_JetGenJetPhi                        .push_back((*AK8JetGenJetPhi )[m_jet]);
	AK8_JetGenPartonID                      .push_back((*AK8JetGenPartonID )[m_jet]);
	AK8_JetGenEn                          	.push_back((*AK8JetGenEn )[m_jet]);
	AK8_JetGenPt                          	.push_back((*AK8JetGenPt )[m_jet]);
	AK8_JetGenEta                          	.push_back((*AK8JetGenEta )[m_jet]);
	AK8_JetGenPhi                          	.push_back((*AK8JetGenPhi )[m_jet]);
	AK8_JetGenPartonMomID                   .push_back((*AK8JetGenPartonMomID )[m_jet]);


        vecSDSJcsv.clear();
        vecSDSJpt.clear();
        vecSDSJeta.clear();
        vecSDSJmass.clear();
        vecSDSJphi.clear();
        vecSDSJe.clear();
        vecSDSJcharge.clear();
        vecSDSJflavour.clear();


	n_AK8SDSJ                          	.push_back((*nAK8SDSJ )[m_jet]);

        for (int sb = 0 ; sb < (*nAK8SDSJ)[m_jet]; sb ++) {

	vecSDSJpt                     	 	.push_back((*AK8SDSJPt)		[m_jet][sb]);
	vecSDSJeta                      	.push_back((*AK8SDSJEta)	[m_jet][sb]);
	vecSDSJphi                          	.push_back((*AK8SDSJPhi)	[m_jet][sb]);
	vecSDSJmass                          	.push_back((*AK8SDSJMass)	[m_jet][sb]);
	vecSDSJe                          	.push_back((*AK8SDSJE)		[m_jet][sb]);
	vecSDSJcharge                           .push_back((*AK8SDSJCharge)	[m_jet][sb]);
	vecSDSJflavour                          .push_back((*AK8SDSJFlavour)	[m_jet][sb]);
	vecSDSJcsv                        	.push_back((*AK8SDSJCSV)	[m_jet][sb]);
      }

     	AK8_SDSJPt                          	.push_back(vecSDSJpt);
	AK8_SDSJEta                          	.push_back(vecSDSJeta); 
	AK8_SDSJPhi                          	.push_back(vecSDSJphi); 
	AK8_SDSJMass                          	.push_back(vecSDSJmass); 
	AK8_SDSJE                          	.push_back(vecSDSJe); 		
	AK8_SDSJCharge                          .push_back(vecSDSJcharge);	
	AK8_SDSJFlavour                         .push_back(vecSDSJflavour);  	
	AK8_SDSJCSV                          	.push_back(vecSDSJcsv); 	

	Fill_Histo_AK8PUPPIJets(m_jet);
}


void TTbar_mc_bkg_PUPPI_030219::Fill_Histo_AK8PUPPIJets(int m_jet)
{

	//============================For Puppi variables ==============================



	AK8_puppiMass      .push_back((*AK8puppiMass)[m_jet]);
	AK8_puppiPhi       .push_back((*AK8puppiPhi)[m_jet]);
	AK8_puppiTau1      .push_back((*AK8puppiTau1)[m_jet]);
	AK8_puppiTau2      .push_back((*AK8puppiTau2)[m_jet]);
	AK8_puppiTau3      .push_back((*AK8puppiTau3)[m_jet]);
	//AK8_puppiTau4      .push_back((*AK8puppiTau4)[m_jet]);
	AK8_puppiSDMass    .push_back((*AK8puppiSDMass)[m_jet]);


	float puppiPt 				= (*AK8puppiPt)[m_jet] ;
	float puppiEta				= (*AK8puppiEta)[m_jet];
	AK8_puppiPt     		   	.push_back(puppiPt);
	AK8_puppiEta  			        .push_back(puppiEta);
	float corr_SD			 	=  getPUPPIweight(puppiPt, puppiEta) ;
	float corr_SDMass			= corr_SD * (*AK8puppiSDMass)[m_jet] ;
	AK8_puppiSDMassL2L3Corr   		.push_back(corr_SDMass);
	AK8_puppiSDL2L3Corr       		.push_back(corr_SD) ; 


	vecPuppiSDSJcsv.clear();
        vecPuppiSDSJpt.clear();
        vecPuppiSDSJeta.clear();
        vecPuppiSDSJmass.clear();
        vecPuppiSDSJphi.clear();
        vecPuppiSDSJe.clear();
        vecPuppiSDSJcharge.clear();
        vecPuppiSDSJflavour.clear();


	n_AK8puppiSDSJ                          .push_back((*nAK8puppiSDSJ )[m_jet]);

	for (int sbn = 0 ; sbn < (*nAK8puppiSDSJ)[m_jet] ; sbn ++) {
	vecPuppiSDSJpt                         .push_back((*AK8puppiSDSJPt) 		  [m_jet][sbn]); 
	vecPuppiSDSJeta                        .push_back((*AK8puppiSDSJEta)		  [m_jet][sbn]);      
	vecPuppiSDSJphi                        .push_back((*AK8puppiSDSJPhi)		  [m_jet][sbn]);      
	vecPuppiSDSJmass                       .push_back((*AK8puppiSDSJMass)		  [m_jet][sbn]);      
	vecPuppiSDSJe                          .push_back((*AK8puppiSDSJE)		  [m_jet][sbn]);      
	vecPuppiSDSJcharge                     .push_back((*AK8puppiSDSJCharge)	 	  [m_jet][sbn]);      
	vecPuppiSDSJflavour                    .push_back((*AK8puppiSDSJFlavour) 	  [m_jet][sbn]);       
	vecPuppiSDSJcsv                        .push_back((*AK8puppiSDSJCSV)		  [m_jet][sbn]);      
}
	//int sbn = 0 ;
	AK8_puppiSDSJPt                         .push_back(vecPuppiSDSJpt)		;//  [m_jet][sbn]); 
	AK8_puppiSDSJEta                        .push_back(vecPuppiSDSJeta)	 	;//[m_jet][sbn]);      
	AK8_puppiSDSJPhi                        .push_back(vecPuppiSDSJphi) 		;//[m_jet][sbn]);      
	AK8_puppiSDSJMass                       .push_back(vecPuppiSDSJmass)		;//[m_jet][sbn]);      
	AK8_puppiSDSJE                          .push_back(vecPuppiSDSJe)		;//[m_jet][sbn]);      
	AK8_puppiSDSJCharge                     .push_back(vecPuppiSDSJcharge)		;//[m_jet][sbn]);      
	AK8_puppiSDSJFlavour                    .push_back(vecPuppiSDSJflavour) 	;//[m_jet][sbn]);       
	AK8_puppiSDSJCSV                        .push_back(vecPuppiSDSJcsv) 		;//[m_jet][sbn]);     


}


void TTbar_mc_bkg_PUPPI_030219::Fill_Histo_Jets(int l_jet )
{


	n_Jet ++ ;

	jet_Pt                                 .push_back((*jetPt)[l_jet]);
	jet_En                                 .push_back((*jetEn)[l_jet]);
	jet_Eta                                .push_back((*jetEta)[l_jet]);
	jet_Phi                                .push_back((*jetPhi)[l_jet]);
	jet_CSV2BJetTags                       .push_back((*jetCSV2BJetTags)[l_jet]);
	//jet_JetProbabilityBJetTags             .push_back((*jetJetProbabilityBJetTags)[l_jet]);
	//jet_pfCombinedMVAV2BJetTags            .push_back((*jetpfCombinedMVAV2BJetTags)[l_jet]);

	jet_PFLooseId                          .push_back((*jetPFLooseId)[l_jet]);
	jet_ID                                 .push_back((*jetID)[l_jet]);
	jet_PUID                               .push_back((*jetPUID)[l_jet]);
	jet_PUFullID                           .push_back((*jetPUFullID)[l_jet]);
	jet_CHF                                .push_back((*jetCHF)[l_jet]);
	jet_NHF                                .push_back((*jetNHF)[l_jet]);
	jet_CEF                                .push_back((*jetCEF)[l_jet]);
	jet_NEF                                .push_back((*jetNEF)[l_jet]);
	jet_NCH                                .push_back((*jetNCH)[l_jet]);
	jet_NNP                                .push_back((*jetNNP)[l_jet]);
	jet_MUF                                .push_back((*jetMUF)[l_jet]);


	jet_Area                          	.push_back((*jetArea)[l_jet]);
	jet_Mt                          	.push_back((*jetMt)[l_jet]);

	jet_DeepCSVTags_b      	 		.push_back((*jetDeepCSVTags_b)[l_jet]);
	jet_DeepCSVTags_bb       	 	.push_back((*jetDeepCSVTags_bb)[l_jet]);
	jet_DeepCSVTags_c           		.push_back((*jetDeepCSVTags_c)[l_jet]);
	//jet_DeepCSVTags_cc        	 	.push_back((*jetDeepCSVTags_cc)[l_jet]);
	jet_DeepCSVTags_udsg   	 		.push_back((*jetDeepCSVTags_udsg)[l_jet]);

	jet_PartonID                          	.push_back((*jetPartonID)[l_jet]);
	jet_HadFlvr                          	.push_back((*jetHadFlvr)[l_jet]);
	jet_GenJetEn                          	.push_back((*jetGenJetEn)[l_jet]);
	jet_GenJetPt                          	.push_back((*jetGenJetPt)[l_jet]);
	jet_GenJetEta                          	.push_back((*jetGenJetEta)[l_jet]);
	jet_GenJetPhi                          	.push_back((*jetGenJetPhi)[l_jet]);
	jet_GenPartonID                         .push_back((*jetGenPartonID)[l_jet]);
	jet_GenEn                          	.push_back((*jetGenEn)[l_jet]);
	jet_GenPt                          	.push_back((*jetGenPt)[l_jet]);
	jet_GenEta                          	.push_back((*jetGenEta)[l_jet]);
	jet_GenPhi                          	.push_back((*jetGenPhi)[l_jet]);
	jet_GenPartonMomID                      .push_back((*jetGenPartonMomID)[l_jet]);


}


void     TTbar_mc_bkg_PUPPI_030219::Fill_MCevent(int MC)
{
	n_MC ++ ;

	mc_PID                      .push_back((*mcPID)[MC]);
	mc_Vtx                      .push_back((*mcVtx)[MC]);
	mc_Vty                      .push_back((*mcVty)[MC]);
	mc_Vtz                      .push_back((*mcVtz)[MC]);
	mc_Pt                       .push_back((*mcPt)[MC]);
	mc_Mass                     .push_back((*mcMass)[MC]);
	mc_Eta                      .push_back((*mcEta)[MC]);
	mc_Phi                      .push_back((*mcPhi)[MC]);
	mc_E                        .push_back((*mcE)[MC]);
	mc_Et                       .push_back((*mcEt)[MC]);
	mc_GMomPID                  .push_back((*mcGMomPID)[MC]);
	mc_MomPID                   .push_back((*mcMomPID)[MC]);
	mc_MomPt                    .push_back((*mcMomPt)[MC]);
	mc_MomMass                  .push_back((*mcMomMass)[MC]);
	mc_MomEta                   .push_back((*mcMomEta)[MC]);
	mc_MomPhi                   .push_back((*mcMomPhi)[MC]);
	mc_StatusFlag               .push_back((*mcStatusFlag)[MC]);
	mc_Parentage                .push_back((*mcParentage)[MC]);
	mc_Status                   .push_back((*mcStatus)[MC]);
	mc_CalIsoDR03               .push_back((*mcCalIsoDR03)[MC]);
	mc_TrkIsoDR03               .push_back((*mcTrkIsoDR03)[MC]);
	mc_CalIsoDR04               .push_back((*mcCalIsoDR04)[MC]);
	mc_TrkIsoDR04               .push_back((*mcTrkIsoDR04)[MC]);


}


//////Default Functions ============


Int_t TTbar_mc_bkg_PUPPI_030219::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TTbar_mc_bkg_PUPPI_030219::LoadTree(Long64_t entry)
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

void TTbar_mc_bkg_PUPPI_030219::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdf = 0;
   EventTag = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   //phoE1x3 = 0;
   //phoE1x5 = 0;
   //phoE2x2 = 0;
   //phoE2x5Max = 0;
   //phoE5x5 = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   //phoE1x3Full5x5 = 0;
   //phoE1x5Full5x5 = 0;
   phoE2x2Full5x5 = 0;
   //phoE2x5MaxFull5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   //phoPFRandConeChIso = 0;
   //phoCITKChIso = 0;
   //phoCITKPhoIso = 0;
   //phoCITKNeuIso = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoxtalBits = 0;
   phoIDbit = 0;
   phoScale_stat_up = 0;
   phoScale_stat_dn = 0;
   phoScale_syst_up = 0;
   phoScale_syst_dn = 0;
   phoScale_gain_up = 0;
   phoScale_gain_dn = 0;
   phoResol_rho_up = 0;
   phoResol_rho_dn = 0;
   phoResol_phi_up = 0;
   phoResol_phi_dn = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleEcalEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   //eledEtaAtCalo = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   //elePFMiniIso = 0;
   //eleIDMVA = 0;
   //eleIDMVAHZZ = 0;
   //eledEtaseedAtVtx = 0;
   //eleE1x5 = 0;
   //eleE2x5 = 0;
   //eleE5x5 = 0;
   //eleE1x5Full5x5 = 0;
   //eleE2x5Full5x5 = 0;
   //eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   //eleDr03EcalRecHitSumEt = 0;
   //eleDr03HcalDepth1TowerSumEt = 0;
   //eleDr03HcalDepth2TowerSumEt = 0;
   //eleDr03HcalTowerSumEt = 0;
   //eleDr03TkSumPt = 0;
   //elecaloEnergy = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFChi2 = 0;
   eleGSFPt = 0;
   eleGSFEta = 0;
   eleGSFPhi = 0;
   eleGSFCharge = 0;
   eleGSFHits = 0;
   eleGSFMissHits = 0;
   eleGSFNHitsMax = 0;
   eleGSFVtxProb = 0;
   eleGSFlxyPV = 0;
   eleGSFlxyBS = 0;
   //eleBCEn = 0;
   //eleBCEta = 0;
   //eleBCPhi = 0;
   //eleBCS25 = 0;
   //eleBCS15 = 0;
   //eleBCSieie = 0;
   //eleBCSieip = 0;
   //eleBCSipip = 0;
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleIDbit = 0;
   eleScale_stat_up = 0;
   eleScale_stat_dn = 0;
   eleScale_syst_up = 0;
   eleScale_syst_dn = 0;
   eleScale_gain_up = 0;
   eleScale_gain_dn = 0;
   eleResol_rho_up = 0;
   eleResol_rho_dn = 0;
   eleResol_phi_up = 0;
   eleResol_phi_dn = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muPFChIso03 = 0;
   muPFPhoIso03 = 0;
   muPFNeuIso03 = 0;
   muPFPUIso03 = 0;
   //muPFMiniIso = 0;
   muFiredTrgs = 0;
   muFiredL1Trgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   muBestTrkType = 0;
   //pfphoEt = 0;
   //pfphoEta = 0;
   //pfphoPhi = 0;
   pfHFEn = 0;
   pfHFECALEn = 0;
   pfHFHCALEn = 0;
   pfHFPt = 0;
   pfHFEta = 0;
   pfHFPhi = 0;
   pfHFIso = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetCSV2BJetTags = 0;
   //jetJetProbabilityBJetTags = 0;
   //jetpfCombinedMVAV2BJetTags = 0;
   jetDeepCSVTags_b = 0;
   jetDeepCSVTags_bb = 0;
   jetDeepCSVTags_c = 0;
   //jetDeepCSVTags_cc = 0;
   jetDeepCSVTags_udsg = 0;
   jetPartonID = 0;
   jetHadFlvr = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetP4Smear = 0;
   jetP4SmearUp = 0;
   jetP4SmearDo = 0;
   jetPFLooseId = 0;
   jetID = 0;
   jetPUID = 0;
   jetPUFullID = 0;
   jetJECUnc = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetNNP = 0;
   jetMUF = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtxNtrks = 0;
   jetVtx3DVal = 0;
   jetVtx3DSig = 0;
   AK8JetPt = 0;
   AK8JetEn = 0;
   AK8JetRawPt = 0;
   AK8JetRawEn = 0;
   AK8JetEta = 0;
   AK8JetPhi = 0;
   AK8JetMass = 0;
   AK8Jet_tau1 = 0;
   AK8Jet_tau2 = 0;
   AK8Jet_tau3 = 0;
   //AK8Jet_CHStau1 = 0;
   //AK8Jet_CHStau2 = 0;
   //AK8Jet_CHStau3 = 0;
   //AK8Jet_CHStau4 = 0;
   AK8JetCHF = 0;
   AK8JetNHF = 0;
   AK8JetCEF = 0;
   AK8JetNEF = 0;
   AK8JetNCH = 0;
   AK8JetNNP = 0;
   AK8JetMUF = 0;
   AK8Jetnconstituents = 0;
   AK8JetPFLooseId = 0;
   AK8JetPFTightLepVetoId = 0;
   AK8JetSoftDropMass = 0;
   AK8JetSoftDropMassCorr = 0;
   AK8JetPrunedMass = 0;
   AK8JetPrunedMassCorr = 0;
   AK8JetpfBoostedDSVBTag = 0;
   AK8JetDSVnewV4 = 0;
   AK8JetCSV = 0;
   AK8JetJECUnc = 0;
   AK8JetL2L3corr = 0;
   AK8puppiPt = 0;
   AK8puppiMass = 0;
   AK8puppiEta = 0;
   AK8puppiPhi = 0;
   AK8puppiTau1 = 0;
   AK8puppiTau2 = 0;
   AK8puppiTau3 = 0;
   //AK8puppiTau4 = 0;
   AK8puppiSDL2L3corr = 0;
   AK8puppiSDMass = 0;
   AK8puppiSDMassL2L3Corr = 0;
   AK8JetPartonID = 0;
   AK8JetHadFlvr = 0;
   AK8JetGenJetIndex = 0;
   AK8JetGenJetEn = 0;
   AK8JetGenJetPt = 0;
   AK8JetGenJetEta = 0;
   AK8JetGenJetPhi = 0;
   AK8JetGenPartonID = 0;
   AK8JetGenEn = 0;
   AK8JetGenPt = 0;
   AK8JetGenEta = 0;
   AK8JetGenPhi = 0;
   AK8JetGenPartonMomID = 0;
   AK8JetP4Smear = 0;
   AK8JetP4SmearUp = 0;
   AK8JetP4SmearDo = 0;
   nAK8SDSJ = 0;
   AK8SDSJPt = 0;
   AK8SDSJEta = 0;
   AK8SDSJPhi = 0;
   AK8SDSJMass = 0;
   AK8SDSJE = 0;
   AK8SDSJCharge = 0;
   AK8SDSJFlavour = 0;
   AK8SDSJCSV = 0;
   nAK8puppiSDSJ = 0;
   AK8puppiSDSJPt = 0;
   AK8puppiSDSJEta = 0;
   AK8puppiSDSJPhi = 0;
   AK8puppiSDSJMass = 0;
   AK8puppiSDSJE = 0;
   AK8puppiSDSJCharge = 0;
   AK8puppiSDSJFlavour = 0;
   AK8puppiSDSJCSV = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   //fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("genPho1", &genPho1, &b_genPho1);
   fChain->SetBranchAddress("genPho2", &genPho2, &b_genPho2);
   fChain->SetBranchAddress("EventTag", &EventTag, &b_EventTag);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   //fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   //fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   //fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   //fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   //fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   //fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   //fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   //fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   //fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   //fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   //fChain->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   //fChain->SetBranchAddress("phoE1x5Full5x5", &phoE1x5Full5x5, &b_phoE1x5Full5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   //fChain->SetBranchAddress("phoE2x5MaxFull5x5", &phoE2x5MaxFull5x5, &b_phoE2x5MaxFull5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   //fChain->SetBranchAddress("phoPFRandConeChIso", &phoPFRandConeChIso, &b_phoPFRandConeChIso);
   //fChain->SetBranchAddress("phoCITKChIso", &phoCITKChIso, &b_phoCITKChIso);
   //fChain->SetBranchAddress("phoCITKPhoIso", &phoCITKPhoIso, &b_phoCITKPhoIso);
   //fChain->SetBranchAddress("phoCITKNeuIso", &phoCITKNeuIso, &b_phoCITKNeuIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoxtalBits", &phoxtalBits, &b_phoxtalBits);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("phoScale_stat_up", &phoScale_stat_up, &b_phoScale_stat_up);
   fChain->SetBranchAddress("phoScale_stat_dn", &phoScale_stat_dn, &b_phoScale_stat_dn);
   fChain->SetBranchAddress("phoScale_syst_up", &phoScale_syst_up, &b_phoScale_syst_up);
   fChain->SetBranchAddress("phoScale_syst_dn", &phoScale_syst_dn, &b_phoScale_syst_dn);
   fChain->SetBranchAddress("phoScale_gain_up", &phoScale_gain_up, &b_phoScale_gain_up);
   fChain->SetBranchAddress("phoScale_gain_dn", &phoScale_gain_dn, &b_phoScale_gain_dn);
   fChain->SetBranchAddress("phoResol_rho_up", &phoResol_rho_up, &b_phoResol_rho_up);
   fChain->SetBranchAddress("phoResol_rho_dn", &phoResol_rho_dn, &b_phoResol_rho_dn);
   fChain->SetBranchAddress("phoResol_phi_up", &phoResol_phi_up, &b_phoResol_phi_up);
   fChain->SetBranchAddress("phoResol_phi_dn", &phoResol_phi_dn, &b_phoResol_phi_dn);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   //fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   //fChain->SetBranchAddress("elePFMiniIso", &elePFMiniIso, &b_elePFMiniIso);
   //fChain->SetBranchAddress("eleIDMVA", &eleIDMVA, &b_eleIDMVA);
   //fChain->SetBranchAddress("eleIDMVAHZZ", &eleIDMVAHZZ, &b_eleIDMVAHZZ);
   //fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   //fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   //fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   //fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   //fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   //fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   //fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   //fChain->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   //fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   //fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   //fChain->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   //fChain->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   //fChain->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFChi2", &eleGSFChi2, &b_eleGSFChi2);
   fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   fChain->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   fChain->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   fChain->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   fChain->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   fChain->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   fChain->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   //fChain->SetBranchAddress("eleBCEn", &eleBCEn, &b_eleBCEn);
   //fChain->SetBranchAddress("eleBCEta", &eleBCEta, &b_eleBCEta);
   //fChain->SetBranchAddress("eleBCPhi", &eleBCPhi, &b_eleBCPhi);
   //fChain->SetBranchAddress("eleBCS25", &eleBCS25, &b_eleBCS25);
   //fChain->SetBranchAddress("eleBCS15", &eleBCS15, &b_eleBCS15);
   //fChain->SetBranchAddress("eleBCSieie", &eleBCSieie, &b_eleBCSieie);
   //fChain->SetBranchAddress("eleBCSieip", &eleBCSieip, &b_eleBCSieip);
   //fChain->SetBranchAddress("eleBCSipip", &eleBCSipip, &b_eleBCSipip);
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleScale_stat_up", &eleScale_stat_up, &b_eleScale_stat_up);
   fChain->SetBranchAddress("eleScale_stat_dn", &eleScale_stat_dn, &b_eleScale_stat_dn);
   fChain->SetBranchAddress("eleScale_syst_up", &eleScale_syst_up, &b_eleScale_syst_up);
   fChain->SetBranchAddress("eleScale_syst_dn", &eleScale_syst_dn, &b_eleScale_syst_dn);
   fChain->SetBranchAddress("eleScale_gain_up", &eleScale_gain_up, &b_eleScale_gain_up);
   fChain->SetBranchAddress("eleScale_gain_dn", &eleScale_gain_dn, &b_eleScale_gain_dn);
   fChain->SetBranchAddress("eleResol_rho_up", &eleResol_rho_up, &b_eleResol_rho_up);
   fChain->SetBranchAddress("eleResol_rho_dn", &eleResol_rho_dn, &b_eleResol_rho_dn);
   fChain->SetBranchAddress("eleResol_phi_up", &eleResol_phi_up, &b_eleResol_phi_up);
   fChain->SetBranchAddress("eleResol_phi_dn", &eleResol_phi_dn, &b_eleResol_phi_dn);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muPFChIso03", &muPFChIso03, &b_muPFChIso03);
   fChain->SetBranchAddress("muPFPhoIso03", &muPFPhoIso03, &b_muPFPhoIso03);
   fChain->SetBranchAddress("muPFNeuIso03", &muPFNeuIso03, &b_muPFNeuIso03);
   fChain->SetBranchAddress("muPFPUIso03", &muPFPUIso03, &b_muPFPUIso03);
   //fChain->SetBranchAddress("muPFMiniIso", &muPFMiniIso, &b_muPFMiniIso);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("muBestTrkType", &muBestTrkType, &b_muBestTrkType);
   fChain->SetBranchAddress("npfPho", &npfPho, &b_npfPho);
   //fChain->SetBranchAddress("pfphoEt", &pfphoEt, &b_pfphoEt);
   //fChain->SetBranchAddress("pfphoEta", &pfphoEta, &b_pfphoEta);
   //fChain->SetBranchAddress("pfphoPhi", &pfphoPhi, &b_pfphoPhi);
   fChain->SetBranchAddress("npfHF", &npfHF, &b_npfHF);
   fChain->SetBranchAddress("pfHFEn", &pfHFEn, &b_pfHFEn);
   fChain->SetBranchAddress("pfHFECALEn", &pfHFECALEn, &b_pfHFECALEn);
   fChain->SetBranchAddress("pfHFHCALEn", &pfHFHCALEn, &b_pfHFHCALEn);
   fChain->SetBranchAddress("pfHFPt", &pfHFPt, &b_pfHFPt);
   fChain->SetBranchAddress("pfHFEta", &pfHFEta, &b_pfHFEta);
   fChain->SetBranchAddress("pfHFPhi", &pfHFPhi, &b_pfHFPhi);
   fChain->SetBranchAddress("pfHFIso", &pfHFIso, &b_pfHFIso);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   //fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   //fChain->SetBranchAddress("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags, &b_jetpfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("jetDeepCSVTags_b", &jetDeepCSVTags_b, &b_jetDeepCSVTags_b);
   fChain->SetBranchAddress("jetDeepCSVTags_bb", &jetDeepCSVTags_bb, &b_jetDeepCSVTags_bb);
   fChain->SetBranchAddress("jetDeepCSVTags_c", &jetDeepCSVTags_c, &b_jetDeepCSVTags_c);
   //fChain->SetBranchAddress("jetDeepCSVTags_cc", &jetDeepCSVTags_cc, &b_jetDeepCSVTags_cc);
   fChain->SetBranchAddress("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg, &b_jetDeepCSVTags_udsg);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetHadFlvr", &jetHadFlvr, &b_jetHadFlvr);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("jetP4Smear", &jetP4Smear, &b_jetP4Smear);
   fChain->SetBranchAddress("jetP4SmearUp", &jetP4SmearUp, &b_jetP4SmearUp);
   fChain->SetBranchAddress("jetP4SmearDo", &jetP4SmearDo, &b_jetP4SmearDo);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetPUFullID", &jetPUFullID, &b_jetPUFullID);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetNNP", &jetNNP, &b_jetNNP);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
   fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
   fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);
   //fChain->SetBranchAddress("nAK8Jetpuppi", &nAK8Jetpuppi, &b_nAK8Jetpuppi);
   fChain->SetBranchAddress("nAK8Jet", &nAK8Jet, &b_nAK8Jet);
   fChain->SetBranchAddress("AK8JetPt", &AK8JetPt, &b_AK8JetPt);
   fChain->SetBranchAddress("AK8JetEn", &AK8JetEn, &b_AK8JetEn);
   fChain->SetBranchAddress("AK8JetRawPt", &AK8JetRawPt, &b_AK8JetRawPt);
   fChain->SetBranchAddress("AK8JetRawEn", &AK8JetRawEn, &b_AK8JetRawEn);
   fChain->SetBranchAddress("AK8JetEta", &AK8JetEta, &b_AK8JetEta);
   fChain->SetBranchAddress("AK8JetPhi", &AK8JetPhi, &b_AK8JetPhi);
   fChain->SetBranchAddress("AK8JetMass", &AK8JetMass, &b_AK8JetMass);
   fChain->SetBranchAddress("AK8Jet_tau1", &AK8Jet_tau1, &b_AK8Jet_tau1);
   fChain->SetBranchAddress("AK8Jet_tau2", &AK8Jet_tau2, &b_AK8Jet_tau2);
   fChain->SetBranchAddress("AK8Jet_tau3", &AK8Jet_tau3, &b_AK8Jet_tau3);
   //fChain->SetBranchAddress("AK8Jet_CHStau1", &AK8Jet_CHStau1, &b_AK8Jet_CHStau1);
   //fChain->SetBranchAddress("AK8Jet_CHStau2", &AK8Jet_CHStau2, &b_AK8Jet_CHStau2);
   //fChain->SetBranchAddress("AK8Jet_CHStau3", &AK8Jet_CHStau3, &b_AK8Jet_CHStau3);
   //fChain->SetBranchAddress("AK8Jet_CHStau4", &AK8Jet_CHStau4, &b_AK8Jet_CHStau4);
   fChain->SetBranchAddress("AK8JetCHF", &AK8JetCHF, &b_AK8JetCHF);
   fChain->SetBranchAddress("AK8JetNHF", &AK8JetNHF, &b_AK8JetNHF);
   fChain->SetBranchAddress("AK8JetCEF", &AK8JetCEF, &b_AK8JetCEF);
   fChain->SetBranchAddress("AK8JetNEF", &AK8JetNEF, &b_AK8JetNEF);
   fChain->SetBranchAddress("AK8JetNCH", &AK8JetNCH, &b_AK8JetNCH);
   fChain->SetBranchAddress("AK8JetNNP", &AK8JetNNP, &b_AK8JetNNP);
   fChain->SetBranchAddress("AK8JetMUF", &AK8JetMUF, &b_AK8JetMUF);
   fChain->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
   fChain->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
   fChain->SetBranchAddress("AK8JetPFTightLepVetoId", &AK8JetPFTightLepVetoId, &b_AK8JetPFTightLepVetoId);
   fChain->SetBranchAddress("AK8JetSoftDropMass", &AK8JetSoftDropMass, &b_AK8JetSoftDropMass);
   fChain->SetBranchAddress("AK8JetSoftDropMassCorr", &AK8JetSoftDropMassCorr, &b_AK8JetSoftDropMassCorr);
   fChain->SetBranchAddress("AK8JetPrunedMass", &AK8JetPrunedMass, &b_AK8JetPrunedMass);
   fChain->SetBranchAddress("AK8JetPrunedMassCorr", &AK8JetPrunedMassCorr, &b_AK8JetPrunedMassCorr);
   fChain->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
   fChain->SetBranchAddress("AK8JetDSVnewV4", &AK8JetDSVnewV4, &b_AK8JetDSVnewV4);
   fChain->SetBranchAddress("AK8JetCSV", &AK8JetCSV, &b_AK8JetCSV);
   fChain->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
   fChain->SetBranchAddress("AK8JetL2L3corr", &AK8JetL2L3corr, &b_AK8JetL2L3corr);
   fChain->SetBranchAddress("AK8puppiPt", &AK8puppiPt, &b_AK8puppiPt);
   fChain->SetBranchAddress("AK8puppiMass", &AK8puppiMass, &b_AK8puppiMass);
   fChain->SetBranchAddress("AK8puppiEta", &AK8puppiEta, &b_AK8puppiEta);
   fChain->SetBranchAddress("AK8puppiPhi", &AK8puppiPhi, &b_AK8puppiPhi);
   fChain->SetBranchAddress("AK8puppiTau1", &AK8puppiTau1, &b_AK8puppiTau1);
   fChain->SetBranchAddress("AK8puppiTau2", &AK8puppiTau2, &b_AK8puppiTau2);
   fChain->SetBranchAddress("AK8puppiTau3", &AK8puppiTau3, &b_AK8puppiTau3);
   //fChain->SetBranchAddress("AK8puppiTau4", &AK8puppiTau4, &b_AK8puppiTau4);
   fChain->SetBranchAddress("AK8puppiSDL2L3corr", &AK8puppiSDL2L3corr, &b_AK8puppiSDL2L3corr);
   fChain->SetBranchAddress("AK8puppiSDMass", &AK8puppiSDMass, &b_AK8puppiSDMass);
   fChain->SetBranchAddress("AK8puppiSDMassL2L3Corr", &AK8puppiSDMassL2L3Corr, &b_AK8puppiSDMassL2L3Corr);
   fChain->SetBranchAddress("AK8JetPartonID", &AK8JetPartonID, &b_AK8JetPartonID);
   fChain->SetBranchAddress("AK8JetHadFlvr", &AK8JetHadFlvr, &b_AK8JetHadFlvr);
   fChain->SetBranchAddress("AK8JetGenJetIndex", &AK8JetGenJetIndex, &b_AK8JetGenJetIndex);
   fChain->SetBranchAddress("AK8JetGenJetEn", &AK8JetGenJetEn, &b_AK8JetGenJetEn);
   fChain->SetBranchAddress("AK8JetGenJetPt", &AK8JetGenJetPt, &b_AK8JetGenJetPt);
   fChain->SetBranchAddress("AK8JetGenJetEta", &AK8JetGenJetEta, &b_AK8JetGenJetEta);
   fChain->SetBranchAddress("AK8JetGenJetPhi", &AK8JetGenJetPhi, &b_AK8JetGenJetPhi);
   fChain->SetBranchAddress("AK8JetGenPartonID", &AK8JetGenPartonID, &b_AK8JetGenPartonID);
   fChain->SetBranchAddress("AK8JetGenEn", &AK8JetGenEn, &b_AK8JetGenEn);
   fChain->SetBranchAddress("AK8JetGenPt", &AK8JetGenPt, &b_AK8JetGenPt);
   fChain->SetBranchAddress("AK8JetGenEta", &AK8JetGenEta, &b_AK8JetGenEta);
   fChain->SetBranchAddress("AK8JetGenPhi", &AK8JetGenPhi, &b_AK8JetGenPhi);
   fChain->SetBranchAddress("AK8JetGenPartonMomID", &AK8JetGenPartonMomID, &b_AK8JetGenPartonMomID);
   fChain->SetBranchAddress("AK8JetP4Smear", &AK8JetP4Smear, &b_AK8JetP4Smear);
   fChain->SetBranchAddress("AK8JetP4SmearUp", &AK8JetP4SmearUp, &b_AK8JetP4SmearUp);
   fChain->SetBranchAddress("AK8JetP4SmearDo", &AK8JetP4SmearDo, &b_AK8JetP4SmearDo);
   fChain->SetBranchAddress("nAK8SDSJ", &nAK8SDSJ, &b_nAK8SDSJ);
   fChain->SetBranchAddress("AK8SDSJPt", &AK8SDSJPt, &b_AK8SDSJPt);
   fChain->SetBranchAddress("AK8SDSJEta", &AK8SDSJEta, &b_AK8SDSJEta);
   fChain->SetBranchAddress("AK8SDSJPhi", &AK8SDSJPhi, &b_AK8SDSJPhi);
   fChain->SetBranchAddress("AK8SDSJMass", &AK8SDSJMass, &b_AK8SDSJMass);
   fChain->SetBranchAddress("AK8SDSJE", &AK8SDSJE, &b_AK8SDSJE);
   fChain->SetBranchAddress("AK8SDSJCharge", &AK8SDSJCharge, &b_AK8SDSJCharge);
   fChain->SetBranchAddress("AK8SDSJFlavour", &AK8SDSJFlavour, &b_AK8SDSJFlavour);
   fChain->SetBranchAddress("AK8SDSJCSV", &AK8SDSJCSV, &b_AK8SDSJCSV);
   fChain->SetBranchAddress("nAK8puppiSDSJ", &nAK8puppiSDSJ, &b_nAK8puppiSDSJ);
   fChain->SetBranchAddress("AK8puppiSDSJPt", &AK8puppiSDSJPt, &b_AK8puppiSDSJPt);
   fChain->SetBranchAddress("AK8puppiSDSJEta", &AK8puppiSDSJEta, &b_AK8puppiSDSJEta);
   fChain->SetBranchAddress("AK8puppiSDSJPhi", &AK8puppiSDSJPhi, &b_AK8puppiSDSJPhi);
   fChain->SetBranchAddress("AK8puppiSDSJMass", &AK8puppiSDSJMass, &b_AK8puppiSDSJMass);
   fChain->SetBranchAddress("AK8puppiSDSJE", &AK8puppiSDSJE, &b_AK8puppiSDSJE);
   fChain->SetBranchAddress("AK8puppiSDSJCharge", &AK8puppiSDSJCharge, &b_AK8puppiSDSJCharge);
   fChain->SetBranchAddress("AK8puppiSDSJFlavour", &AK8puppiSDSJFlavour, &b_AK8puppiSDSJFlavour);
   fChain->SetBranchAddress("AK8puppiSDSJCSV", &AK8puppiSDSJCSV, &b_AK8puppiSDSJCSV);
   Notify();
}

Bool_t TTbar_mc_bkg_PUPPI_030219::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TTbar_mc_bkg_PUPPI_030219::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TTbar_mc_bkg_PUPPI_030219::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TTbar_mc_bkg_PUPPI_030219_cxx
