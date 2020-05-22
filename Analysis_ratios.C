////////////////////////////////////////////////////////////////////
// Program description: Ratio Analysis of two ormore root files.  //
////////////////////////////////////////////////////////////////////


/////////////////
/////////////////
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TBranchElement.h"


//////////////////


#include <fstream>
#include <TFile.h>
#include <TBranch.h>
//#include <TMath.h>
//#include <TString.h>
//#include <TH1.h>
//#include <TH2.h>
#include <stdio.h>
//#include <TTreeReader.h>
//#include <TTreeReaderValue.h>
//#include <TTreeReaderArray.h>
//#include "TROOT.h"
//#include "TTree.h"

#include "TSystem.h"
#include "TKey.h"
#include "TDirectoryFile.h"

using namespace std;

void Analysis_ratios()
{
	string dirname = "rootFiles/";
	string PID[2] = {"eXe","eD"};
	string file[2];

	TFile *fin[2];

	//Reading ROOT files

	for (int i = 0; i < 2;i++)
	{
	
		file[i] = (dirname+"Analysis_ratio_"+PID[i]+"_490GeV_fixed_target_40k.root");
		fin[i] = new TFile((file[i]).c_str());
	
	}

	cout << "file 1: " << file[0] << endl;
	cout << "file 2: " << file[1] << endl;

	//Listing of histograms

	//fin[0]->GetListOfKeys()->Print();
	//fin[0]->ls();

	//fin[1]->GetListOfKeys()->Print();
	//fin[1]->ls();
	
	//################################ Creating and extracting histograms ###################################

	TH1F *h_1   = (TH1F*)fin[0]->Get("h_nTracks_p_gt");
	TH1F *h_2   = (TH1F*)fin[1]->Get("h_nTracks_proton");

	// ################################ RATIOS ############################################


	TH1F *Ratio_p_gt  = new TH1F("Ratio_p_gt","Ratio_p_gt",100,0,10);

	Ratio_p_gt->Sumw2();
	Ratio_p_gt->Divide(h_1,h_2,1,0.8566);

	TCanvas *c1 = new TCanvas("c1","Ratios protons");
	//c1->Divide(3,1);
	c1->cd();
	Ratio_p_gt->SetTitle("Ng(eXe) / N(eD) pdg:2212; Ng; R");
	Ratio_p_gt->Draw();



	return;
}
