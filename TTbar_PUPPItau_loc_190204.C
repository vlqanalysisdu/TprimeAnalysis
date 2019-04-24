#define TTbar_PUPPItau_loc_190204_cxx
#include "TTbar_PUPPItau_loc_190204.h"
#include <TH2.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include  <map>

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2){ 
		std::cout<<"pass the text file containing root file names "<<std::endl;
		exit(0);
	}
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <map>");
	TTbar_PUPPItau_loc_190204  a(argv[1]);
	TString InputTxtFile = argv[1];
	TString OutputFileName = InputTxtFile.ReplaceAll(".txt","_");
	a.Loop(OutputFileName.Data(), 1);
	a.Loop(OutputFileName.Data(), 2);
     	a.Loop(OutputFileName.Data(), 3);
     	a.Loop(OutputFileName.Data(), 4);
     	a.Loop(OutputFileName.Data(), 5);

}

void TTbar_PUPPItau_loc_190204::Loop(TString OutputFileName, int y)
{
	if (fChain == 0) return;
	//TString cut, root = "Tb_tH_1500_Muon_" ;
	int var1 = y; 
	
	//if(!(var1 ==1 || var1 == 2) ) exit(0) ;

	if( var1 == 1 ) OutputFileName = OutputFileName +  "METCut100_05-03-19.root" ;
	if( var1 == 2 ) OutputFileName = OutputFileName +  "METCut150_05-03-19.root" ;
	if( var1 == 3 ) OutputFileName = OutputFileName +  "METCut200_05-03-19.root" ;
	if( var1 == 4 ) OutputFileName = OutputFileName +  "METCut250_05-03-19.root" ;
	if( var1 == 5 ) OutputFileName = OutputFileName +  "METCut300_05-03-19.root" ;
//	if( var1 == 2 ) OutputFileName = OutputFileName +  "Group3_Jet+mu_55Wtag80Toptag_PUPPI70corr_290119.root" ;
//	if( var1 == 3 ) OutputFileName = OutputFileName +  "Group5_Jet+mu_55Wtag80Toptag_PUPPI70corr_290119.root" ;

	TFile* f2 = new TFile(OutputFileName.Data(),"recreate");  
	//TFile* f2 = new TFile( root.Data() ,"recreate");  
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<"\nTotal Eventt = " << nentries <<endl ;

	//===========Histogram Functions ==============================	
	/*	
		Define_Mt_Histo() ;  */
	Define_2DMass_Histo();
	dRHisto_MCRecoObject() ;   	
	DefineMC_NPtEta_Histo() ;	 		
	dRHisto_MCObject() ;     
	Define_NPtEta_Histo();
	dRHisto_RecoObject();
	Define_Tag_Jet_Histo() ;  
	Define_Reco_tagjetHisto() ;   
	dR_tagjetHisto() ;     
	Category_Object_Histo() ;    
	Category_Object_dRHisto() ;     
	Category_Object_MtHisto() ; 

	Ptbjet_dRW   = new TH2F ("Ptbjet_WdR", " Ptbjet_wrt_ dR(W,b)", 1500, 0.0 , 1500.0, 400, 0.0 ,400.0 );
	Ptbjet_dRmu  = new TH2F ("Ptbjet_mudR", " Ptbjet_wrt_ dR(mu,b)", 1500, 0.0 , 1500.0, 400, 0.0 ,400.0 );

	//==========================Event Study ===================

	Long64_t nbytes = 0, nb = 0;
	int j, st = 0 , muon = 0, bjet = 0 , AK8 = 0, mc = 0, ele= 0, bquark = 0, mu_bjet = 0, top = 0, lep = 0, event_muon = 0;
	int  q_forw, topbjet = -1;


	int b_match = 0;
	int W_match = 0;
	int qj = -1 ;
	float dR = 0.0 ;
	float dR1 = 0.0 ;
	double pz ;
	vr = var1 ;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		//  for (Long64_t jentry=0; jentry<5000;jentry++) { }
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//if (ientry == 2) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		Clear_Vector() ;
		j = -1;
		b_asso = -1 ;
		b_top = -1 ;
		q_forw = -1;

		if (jentry ==  nentries-1 ) {
		
		cout << " Total events = " << nentries << 
		" -> AK8jet Cut passed = " << AK8 <<
		" -> Lepton cut passed = " << lep <<
		" -> Muon Cut passed   = " << event_muon << endl ;

		cout << " 0Wtag Events = " << event_0W <<
		" -> Higgs Cut  passed = " << event_Higgs <<
		" -> bjet  Cut 	passed = " << event_bjet <<
		" -> dR Cut 	passed = " << event_dR_Hmu <<
		" -> LepIso Cut passed = " << event_LepIso_bmu <<
		" -> HiggsPt Cut pass  = " << event_HiggsPt << endl ;

		cout << " 1Wtag Events = " << event_1W <<
		" -> top   Cut  passed = " << event_top <<
		" -> dR Cut 	passed = " << event_dR_Wt_tmu <<
		" -> LepIso Cut passed = " << event_LepIso_Wmu_cat2 << 
		" -> TopPt Cut  passed = " << event_topPt << endl;

		cout << " 2Wtag Events = " << event_2W <<
		" -> bjet  Cut  passed = " << event_bjet_cat3 <<
		" -> dR Cut 	passed = " << event_dR_Wb_bmu <<
		" -> LepIso Cut passed = " << event_LepIso_Wmu_cat3 <<
		" -> WPt Cut    passed = " << event_W_Pt << endl;
 
		}

		//===========================MC info========================
		T_top = -1;
		b_top = -1 ;
		T_higgs = -1;
		Higgs_W = -1 ;
		Top_W = -1 ;
		int check =0;

	
 		for(int h =0; h< n_MC; h++){

                        //if ( abs((*mc_MomPID)[h]) == 24 ) cout << "\n , daughterID = " << (*mc_PID)[h] << " , parentID = " << (*mc_MomPID)[h] << " , Status = " << (*mc_Status)[h] ;
                        //continue;

       		        if( abs((*mc_PID)[h]) == 24 && abs((*mc_MomPID)[h]) ==  6 ) Top_W = h ;
                        if(abs((*mc_PID)[h]) <= 4 &&  (*mc_Status)[h] <= 25 && abs((*mc_MomPID)[h]) == 24 )  Higgs_Jet.push_back(h);

                }

		//=====================   Object selections  ================
		//Ptbjet_dRW->Fill((*AK8_JetPt)[0], (*AK8_JetPrunedMass)[0] ) ;
		//Ptbjet_dRmu->Fill((*AK8_JetPt)[1], (*AK8_JetPrunedMass)[1] ) ;
		//======================AK8 jet selection ==================
		for(int f = 0 ; f < N_AK8Jet ; f ++)
		{
			Cut_AK8jet(f, var1)  ;
		}
		if( n_AK8Jet.size() == 0)continue ;
		AK8 ++ ;
		//======================electron selection==========================
		for( int f = 0; f < n_Ele; f++)
		{
			if((*ele_Pt)[f] < 40.0 ) continue ;
			n_ele.push_back(f);
		}
		ele = n_ele.size();
		//======================muon selection==========================
		for( int f = 0; f < N_Mu ; f++)
		{
			if (Cut_Muon(f) ) n_Mu.push_back(f);
		}

		muon =  n_Mu.size() ;
		//if ( n_Mu.size() == 0    ) continue;		
		if((ele == 0 && muon == 0) ) continue ;
		lep ++ ;
		//=====================bjet selection=================================
		for( int f = 0 ; f < n_Jet ; f++ )
		{
			Cut_bjet(f,1) ;   // 1 for loose, 2 for medium, 3 for tight
		}

		//if ( n_forwjet.size() == 0) continue ;
		// ==================AK8 jetstudy==============================
		int size8  =  n_AK8Jet.size() ;
		int g2 = -1 ; 
		int lvl = -1 ;
		int b2 = -1;
		int t2 = -1 ;
		//if (size8 < 4) continue ;
		bjet ++ ;
	
		Higgs_selection(vr) ; 
		Top_selection(vr) ;		
		Wjet_selection(vr);

	
		if( vr == 1) {
		Higgsjet_Plots(1) ;			
		Topjet_Plots(vr);
		Wjet_Plots(vr);
		TagJets_dRPlots() ;
		}


			
		//========= Signal Category Study starts from here========================			
		if ( n_Mu.size() == 0    ) continue;
		event_muon ++ ;

	   	if (W_boson.size() == 0)  {
		   Wtag0_Category() ; 
		   event_0W ++;	
		   }

                if (W_boson.size() == 1 ) {
		   Wtag1_Category() ;
	           event_1W ++ ;				
		   }

		if (W_boson.size() != 2)   continue ; 
                   event_2W ++ ; 
		   if (b_jet.size() != 0 )  {
		   event_bjet_cat3 ++ ;
		   WW_lvbjet_Plots() ;    
		   }	
		//==============================

	}

	f2->Write();
	f2->Close();
}
