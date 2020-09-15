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
	string dirname = "/Users/carorg/projects/eic/eic-centrality/root_files/analysis-ratios/";
	string PID[2] = {"eXe","eD"};
	string file[2];

	TFile *fin[2];

	//Reading ROOT files

	for (int i = 0; i < 2;i++)
	{
	
		file[i] = (dirname+"Analysis_ratio_"+PID[i]+"_490GeV_fixed_target_40k_TEST.root");
		fin[i] = new TFile((file[i]).c_str());
	
	}

	cout << "file 1: " << file[0] << endl;
	cout << "file 2: " << file[1] << endl;

	//Listing of histograms

	//fin[0]->GetListOfKeys()->Print();
	//fin[0]->ls();

	//fin[1]->GetListOfKeys()->Print();
	//fin[1]->ls();

	//h_nTracks_proton->GetMean()
//(double) 1.1788324
	
	//################################ Creating and extracting histograms ###################################

	TH1F *h_1   = (TH1F*)fin[0]->Get("h_nTracks_p_gt");
	TH1F *h_2   = (TH1F*)fin[1]->Get("h_nTracks_proton");
	TH1F *h_3 	= (TH1F*)fin[0]->Get("h_nTracks");

	// ################################ RATIOS ############################################


	TH1F *h_3  = new TH1F("h_3","h_3",10,0,10);
	TH1F *Ratio_p_gt  = new TH1F("Ratio_p_gt","Ratio_p_gt",11,-0.5,10.5);

	//TH1F *h_3 = *h_1->Scale(1/(h_2->GetMean()));
	//TH1F *h_3 = *h_1->Scale(1/5);
	h_1->Scale(1/h_2->GetMean());


	


	Ratio_p_gt->Sumw2();
	Ratio_p_gt->Divide(h_1,h_2,1,0.8566);

	

	TCanvas *c1 = new TCanvas("c1","Ratio Protons");
	//c1->Divide(3,1);
	c1->cd();
	Ratio_p_gt->SetTitle("N_{g}(eXe)/N(eD); N_{g}");
	Ratio_p_gt->SetTitle("Ng(eXe) / N(eD) Protons; Ng; R = Ng(eXe) / N(eD)");
	Ratio_p_gt->SetMarkerStyle(kStar);
	//Ratio_p_gt->SetMarkerColor(kBlue);
	Ratio_p_gt->Draw("P");

	TCanvas *c2 = new TCanvas("c2","c2");
	c2->cd();
	gPad->SetLogy();
	h_1->Draw();

	//Ratio_p_gt_test->Draw();



	return;
}
