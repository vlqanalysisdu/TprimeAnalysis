#define TTbar_mc_bkg_PUPPI_030219_cxx
#include "TTbar_mc_bkg_PUPPI_030219.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


using namespace std;

int main(int argc, char **argv)
{
        if (argc < 2){ 
                std::cout<<"pass the text file containing root file names "<<std::endl;
                exit(0);
        }
        gROOT->ProcessLine("#include <vector>");
        gROOT->ProcessLine("#include <map>");
        TTbar_mc_bkg_PUPPI_030219 a(argv[1]);
        TString InputTxtFile = argv[1];
        TString OutputFileName = InputTxtFile.ReplaceAll(".txt",".root");
        a.Loop(OutputFileName.Data());

} 
void TTbar_mc_bkg_PUPPI_030219::Loop(TString OutputFileName)
{
        //---------------------------------------------------------------------------------
        if (fChain == 0) return;
        //..................................................................................
        TFile* f2 = new TFile(OutputFileName.Data(), "recreate");

        //..................................................................................

        Long64_t nentries = fChain->GetEntriesFast();
        Long64_t nbytes = 0, nb = 0;

        cout<<"Total number of events = "<<nentries<<endl;
        TTree *v_tree = new TTree("v_tree","Tree with vectors");

        Define_tree(v_tree);


        int event_lep = 0;
        int event_jet = 0;
        int i_tH ;
        int event_higg = 0;
        int event_Wb = 0;

        int event_bjet = 0;
        int iso_ele = -1;
        int iso_muon = -1;
        int pass_ID = 0;
        int pass_eID = 0;
        float j_ID = 0.0;
        int input = 0 ;
        int zth = 0;
        int eth = 0;
        int mth = -1;
        float dR =  99.0;
        float max = 99.0;
        int puppi_jet = -1 ;
        vector <int> particle_ID ;

//.....................PUPPI SD Correction Variables.........................
        fi = TFile::Open( "puppiCorr.root","READ");
        puppisd_corrGEN      = (TF1*)fi->Get("puppiJECcorr_gen");
        puppisd_corrRECO_cen = (TF1*)fi->Get("puppiJECcorr_reco_0eta1v3");
        puppisd_corrRECO_for = (TF1*)fi->Get("puppiJECcorr_reco_1v3eta2v5");
        //.................................................................

        for (Long64_t jentry=0; jentry<nentries;jentry++) {
                //      for (Long64_t jentry=0; jentry< 5000 ;jentry++) {}
                Long64_t ientry = LoadTree(jentry);
                if (ientry < 0) break;
                nb = fChain->GetEntry(jentry);
                nbytes += nb;
                //              Clear_Vector();

               // if (jentry == 39 ) break ;
                iso_ele = -1;
                iso_muon = -1;
                //              eth = -1;
                mth = -1;
                j_ID = 0.0;

                if (jentry ==  nentries-1 ) {
                }
                //---Filling of branches starts-------------------------------------------------------
                Clear_Branches();
                N_AK8Jet = 0 ;
                n_Jet = 0;
                N_Mu = 0 ;
                n_Ele = 0 ;
                n_MC = 0 ;
                i_tH = 0;

                for (int c_MC = 0; c_MC < nMC; c_MC++ )
	{
			Fill_MCevent(c_MC);

		}


		if (nAK8Jet < 1 && nJet < 1 && (nEle < 1 || nMu < 1))continue;
		// For AK8jet Variables 
		for(int c_jet8  = 0; c_jet8<nAK8Jet ;c_jet8 ++) {

			dR =  99.0;
			max = 99.0;				
			puppi_jet = -1 ;
			if (!(Cut_AK8Jet(c_jet8))) continue;
		        Fill_Histo_AK8Jets(c_jet8) ;
	
		}



		// For AK4jet Variables

		for(int c_jet  = 0; c_jet< nJet; c_jet++)
		{			
			if( Cut_Jet(c_jet) )  Fill_Histo_Jets(c_jet);
		}
		//		if(n_Jet < 1) continue;
		pass_ID ++;
		// For electron variables.
		for(int c_ele  = 0; c_ele  < nEle; c_ele ++)
		{
			if(Cut_Electron(c_ele )) Fill_Histo_Ele(c_ele) ;
		} 


		// For Muon Variables
		for (int c_muon = 0; c_muon< nMu; c_muon++)
		{
			if(Cut_Muon(c_muon)) Fill_Histo_Muon(c_muon);           
		} 



		event_lep ++ ;
		v_event++ ;  
		n_Vtx = nVtx;
		n_GoodVtx = nGoodVtx;
		//n_TrksPV  = nTrksPV;
		is_PVGood  = isPVGood;
		v_vtx      = vtx;
		v_vty      = vty; 
		v_vtz      = vtz;

		//2. MET Variables
		pf_MET     = pfMET ;
		//pf_METPhi  = pfMETPhi;
		pf_METsumEt= 0;//pfMETsumEt;
		pf_METmEtSig = 0;// pfMETmEtSig;
		//pf_METSig  = pfMETSig;

		gen_MET   = genMET ;
		gen_METPhi = genMETPhi;
	      //=============================================================================================

		event_bjet ++;
		v_tree->Fill();  // check 
	}
	f2->Write();
	fi->Close();
	f2->Close();
}



