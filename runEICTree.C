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

#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

using namespace std;
using namespace erhic;

void runEICTree(const TString filename="/Users/carorg/projects/Analysis_BeAGLE/rootFiles/eXe_490GeV_fixed_target_40k", const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);
	
	TFile *out = new TFile ("rootFiles/Analysis_ratio_eXe_490GeV_fixed_target_40k_TEST.root", "RECREATE");
	//####################  Creating Histograms ####################################

	TH1F *h_nTracks = new TH1F("h_nTracks", "h_nTracks", 400, 0, 400);
	TH1F *h_pdg = new TH1F("h_pdg", "h_pdg", 500,-2500, 2500);
	TH1F *h_Energy = new TH1F("h_Energy", "h_Energy", 100, 0,6);
	TH1F *h_dE_dN = new TH1F("h_idE_dN", "h_dE_dN", 100, 0, 200);
	TH1F *h_Status = new TH1F("h_Status", "h_Status", 50, 0, 25);
	TH1F *h_pt = new TH1F("h_pt","h_pt", 100,0, 2);
	TH1F *h_eta = new TH1F("h_eta", "h_eta", 100, -10, 10);
	TH1F *h_phi = new TH1F("h_phi", "h_phi", 100,-10, 10);
	TH1F *h_rapidity = new TH1F("h_rapidity","h_rapidity", 100, -7, 7);
	TH1F *h_mass = new TH1F("h_mass", "h_mass", 100, 0, 10);
	TH1F *h_b = new TH1F("h_b", "h_b", 100, 0, 10);
	TH1F *h_trueQ2 = new TH1F("h_trueQ2", "h_trueQ2", 100, 0, 20);
	TH1F *h_nTracks_cuts = new TH1F("h_nTracks_cuts","h_nTracks_cuts", 100, 0,50);
	TH1F *h_Energy_neutron = new TH1F("h_Energy_neutron", "h_Energy_neutron", 100, 0, 3);
	TH1F *h_nTracks_neutron = new TH1F("h_nTracks_neutron", "h_nTracks_neutron", 100, 0, 50);
	TH1F *h_pt_neutron = new TH1F("h_pt_neutron", "h_pt_neutron", 100, 0, 2);
	TH1F *h_Energy_proton = new TH1F("h_Energy_proton", "h_Energy_proton", 100, 0, 3);
	TH1F *h_Energy_proton_gt = new TH1F("h_Energy_proton_gt", "h_Energy_proton_gt", 100, 0, 3);
	TH1F *h_nTracks_proton = new TH1F("h_nTracks_proton", "h_nTracks_proton", 100, 0, 10);
	TH1F *h_pt_proton = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
	TH1F *h_pt_proton_gt= new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
	TH1F *h_nTracks_gt = new TH1F("h_nTracks_gt","h_nTracks_gt",100, 0, 50);
	TH1F *h_pt_gt = new TH1F("h_pt_gt", "h_pt_gt", 100, 0, 1);
	TH1F *h_Energy_gt = new TH1F("h_Energy_gt", "h_Energy_gt", 100,0, 6);
	TH1F *h_rapidity_gt = new TH1F("h_rapidity_gt","h_rapidity_gt",100,-10,10);
	TH1F *h_nTracks_p_gt= new TH1F("h_nTracks_p_gt","h_nTracks_p_gt", 100, 0, 10);
	TH1F *h_Q2_p= new TH1F("h_Q2_p","h_Q2_p", 100, 0, 50);
	TH1F *h_Q2_p_gt = new TH1F("h_Q2_p_gt", "h_Q2_p_gt", 100, 0, 50);
	TH1F *h_w2 = new TH1F("h_W2","h_W2", 100, 0, 1000);
	TH1F *h_nu = new TH1F("h_nu", "h_nu",100,0,1000);
	TH1F *h_bx = new TH1F("h_bx","h_bx",100,0,1);

	TH2F *h_rap_eta = new TH2F("h_rap_eta","h_rap_eta", 100, 0, 1, 100, 0, 1);
	TH2F *h_mass_mom = new TH2F("h_mass_mom","h_mass_mom", 100, 0, 50, 100, 0, 2);
	TH2F *h_mass_pt = new TH2F("h_mass_pt","h_mass_pt", 100, 0, 6, 100, 0, 2);

	TH1F *h_mom = new TH1F("h_momn", "h_mom", 100, 0, 50);
	
	int counter_total_particles = 0;
	int counter_event = 0;
	//int counter_neu = 0;
	string title_collision = "e - Xe";

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);
		
		double pzlep = event->pzlep;
		double pztarg = event->pztarg;
		int struck_nucleon = event->nucleon;
		double MASS_NUCLEON = MASS_PROTON;
		if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

		TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON));
		TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON));
		TLorentzVector e_scattered(0.,0.,0.,0.);

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->nu;
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		

		double impact_parameter = event->b;
		double Tb = event->Thickness;
		double distance = event->d1st;
		int N_nevap = event->Nnevap;
		int N_pevap = event->Npevap;
		
		//h_nTracks->Fill(nParticles);
		

		//event cuts
		if( event_process != 99 ) continue;	      //DIS process for pythia
		if( trueQ2 < 1 ) continue;
		if( trueX < 0.002 ) continue;
		if( trueNu < 50 || trueNu > 400 )continue;
		if( trueW2 < 16) continue;
		//if( trueQ2 < 1. || trueQ2 > 20. ) continue;
		//if( trueY > 0.95 || trueY < 0.1 ) continue;

		//do analysis, or fill historgrams for event levels
		
		counter_event++;

		//particle loop
		int counter = 0;
		int counter_neu = 0;
		int counter_p = 0;
		int gt_counter = 0;
		int counter_p_gt = 0;

		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			int index = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			double phi = particle->GetPhi();
			double rap = particle->GetRapidity();
			double mass = particle->GetM();
			double theta = particle->GetTheta(); 
			theta = theta*1000.0; //change to mrad;
			double mom = particle->GetP();
			double Energy = particle->GetE();
			
			//cout << "pdg ~ " << pdg << endl;
			
			//do analysis track-by-track
			
			//particle cuts
			if (status != 1) continue;
			//if (pt < 0.8) continue;
			if (pt > 0.2 && pt < 0.8)
			{
				gt_counter++;
				h_pt_gt->Fill(pt);
				h_Energy_gt->Fill(Energy);
				if (pdg == 2212)
				{
					counter_p_gt++;
					h_Energy_proton_gt->Fill(Energy);
					h_pt_proton_gt->Fill(pt);
					h_rapidity_gt->Fill(rap);
				}
			}	
			//Particles without the Pt Cut

			if (pdg == 2112)
			{

			counter_neu++;
			h_Energy_neutron->Fill(Energy);
			h_pt_neutron->Fill(pt);
			
			}

			if (pdg == 2212)
			{

				counter_p++;
				h_Energy_proton->Fill(Energy);
		    	h_pt_proton->Fill(pt);

			}
			

			//if (pdg != 2212 && pdg != 211 && pdg!=-211  && pdg!=321 && pdg!=-321 && pdg!=310 && pdg!= 331 && pdg!=1114 && pdg!=2214 ) continue;
			//if (pt > 0.2 || pt < 0.01) continue;
					
			counter++;
			counter_total_particles++;
			h_Energy->Fill(Energy);
						
			/*for (int i = 0; i < h_Energy->GetNbinsX(); ++i)
			 {
				 double binvalue = h_Energy->GetBinContent(i+1);
				 double binwidth = h_Energy->GetBinWidth(i+1);

				 h_Energy->SetBinContent (j+1, binvalue/binwidth);
			 }*/ 
																										
			//Fill histograms particle loop

			h_pdg->Fill(pdg);
			//h_Energy->Fill(Energy);
			h_Status->Fill(status);
			h_pt->Fill(pt);
			h_eta->Fill(eta);
			h_phi->Fill(phi);
			h_rapidity->Fill(rap);
			h_mass->Fill(mass);
			h_mom->Fill(mom);
			h_rap_eta->Fill(rap,eta,1);
			h_mass_mom->Fill(mom,mass,1);
			h_mass_pt->Fill(pt,mass,1);

		} // end of particle loop
	
		
		//fill histograms
		h_nTracks->Fill(nParticles);
		h_b->Fill(impact_parameter);
		h_trueQ2->Fill(trueQ2);
		h_nTracks_cuts->Fill(counter);
		h_nTracks_neutron->Fill(counter_neu);
		h_nTracks_proton->Fill(counter_p);//all tracks protons
		h_nTracks_gt->Fill(gt_counter);
		h_nTracks_p_gt->Fill(counter_p_gt);//grey tracks protons
		h_Q2_p->Fill(trueQ2);
		h_w2->Fill(trueW2);
		h_nu->Fill(trueNu);
		h_bx->Fill(trueX);
	}
	cout << "Number of particles: " << counter_total_particles << endl;
	cout << "Events after event cuts : " << counter_event << endl;
	
	//################################### Drawing for e - Xe #############################
	TCanvas *t15 = new TCanvas("t15","t15");
	t15->cd();
	gPad->SetLogy();
	h_nTracks_p_gt->SetTitle((title_collision+ " Log Scale  nTracks  pdg:2212 ; Grey Tracks ; nEvents").c_str());
	h_nTracks_p_gt->Draw();
	
	TCanvas *t23 = new TCanvas("t23","t23");
	t23->cd();
	gPad->SetLogy();
	h_nTracks_proton->SetTitle((title_collision+ " Log Scale nTracks pdg:2212 ; nTracks; nEvents").c_str());
	h_nTracks_proton->Draw();

	out->Write();

	return;
	
	TCanvas *t = new TCanvas("t","t");
	t->cd();
	h_pdg->SetTitle((title_collision+ " PDG All particles; PDG;nParticles").c_str());
	h_pdg->Draw();
	
	TCanvas *t10 = new TCanvas("t10","t10");
	t10->cd();
	h_nTracks_cuts->SetTitle((title_collision+ " nTracks All particles  ; nTracks ; nEvents").c_str());
	h_nTracks_cuts->Draw();

	TCanvas *t14 = new TCanvas("t14","t14");
	t14->cd();
	h_Energy_proton->SetTitle((title_collision+ " Energy pdg:2212 ; E(GeV) ; nParticles").c_str());
	h_Energy_proton->Draw();
				    

	TCanvas *t28 = new TCanvas("t28","t28");
	t28->cd();
	h_Energy_proton_gt->SetTitle((title_collision+ " Energy Gt pdg:2212  ; E(GeV) ; nParticles").c_str());
	h_Energy_proton_gt->Draw();

	
	TCanvas *t5 = new TCanvas("t5","t5");
	t5->cd();
	h_rapidity->SetTitle((title_collision+ " Rapidity All stable and final particles; Rapidity ;nParticles").c_str());	
	h_rapidity->Draw();
	
	TCanvas *t27 = new TCanvas("t27","t27");
	t27->cd();
	h_bx->SetTitle((title_collision+ " Bj x; x; nEvents").c_str());
	h_bx->Draw();
	
	TCanvas *t26 = new TCanvas("t26","t26");
	t26->cd();
	h_nu->SetTitle((title_collision+ " Nu ; Nu; nEvents").c_str());
	h_nu->Draw();
	
	TCanvas *t24 = new TCanvas("t24", "t24");
	t24->cd();
	h_Q2_p->SetTitle((title_collision+ " Q^2; Q^2; nEvents").c_str());
	h_Q2_p->Draw();
	
	TCanvas *t25 =new TCanvas("t25","t25");
	t25->cd();
	h_w2->SetTitle((title_collision+ " W^2 ; w^2; nEvents").c_str());
	h_w2->Draw();
   
	TCanvas *t22 = new TCanvas("t22","t22");
	t22->cd();
	h_rapidity_gt->SetTitle((title_collision+ " Rapidity Grey Tracks; Grey Tracks; nParticles").c_str());
	h_rapidity_gt->Draw();
		
	//TCanvas *t21 = new TCanvas("t21","t21");
	//t21->cd();
	//h_nTracks_gt->SetTitle((title_collision+ " All stable grey tracks;Grey Tracks; nEvents").c_str());
	//h_nTracks_gt->Draw();

	//TCanvas *t14 = new TCanvas("t14","t14");
	//t14->cd();
	//h_Energy_proton->SetTitle((title_collision+ " Energy pdg:2212 ; E(GeV) ; nParticles").c_str());
	//h_Energy_proton->Draw();
	//t14->SaveAs("plots/Energy_proton.png");

	//TCanvas *t15 = new TCanvas("t15","t15");
	//t15->cd();
	//gPad->SetLogy();
	//h_nTracks_p_gt->SetTitle((title_collision+ " nTracks  pdg:2212 ; Grey Tracks ; nEvents").c_str());
	//h_nTracks_p_gt->Draw();
	//t15->SaveAs("plots/nTracks_proton.png");
	
	
	//TCanvas *t23 = new TCanvas("t23","t23");
	//t23->cd();
	//h_nTracks_proton->SetTitle((title_collision+ " nTracks pdg:2212 ; nTracks; nEvents").c_str());
	//h_nTracks_proton->Draw();
	
	TCanvas *t16 = new TCanvas("t16","t16");
	t16->cd();
	h_pt_proton->SetTitle((title_collision+ " Pt  pdg:2212 ; Pt ; nParticles").c_str());
	h_pt_proton->Draw();
	
	TCanvas *t29 = new TCanvas("t29","t29");
	t29->cd();
	h_pt_proton_gt->SetTitle((title_collision+ " Pt Gt pdg:2212 ; Pt ; nParticles").c_str());
	h_pt_proton_gt->Draw();
	
	out->Write();
	
	return;

	
	
	TCanvas *t17 = new TCanvas("t17","t17");
	t17->cd();
	h_rap_eta->SetTitle((title_collision+ " All stable particles; rapidity; eta").c_str());
	h_rap_eta->Draw();
	
	TCanvas *t18 = new TCanvas("t18","t18");
	t18->cd();
	h_mass_mom->SetTitle((title_collision+ " All stable particles; p; M").c_str());
	h_mass_mom->Draw();
	

	TCanvas *t19 = new TCanvas("t19","t19");
	t19->cd();
	h_mom->SetTitle((title_collision+ " All stable particles; p(GeV); nParticles").c_str());
	h_mom->Draw();

	TCanvas *t20 = new TCanvas("t20","t20");
	t20->cd();
	h_mass_pt->SetTitle((title_collision+ " All stable particles; Pt(GeV); nParticles").c_str());
	h_mass_pt->Draw();

		
	
	//TCanvas *t = new TCanvas("t","t");
	//t->cd();
	//h_pdg->SetTitle((title_collision+ " PDG All particles; PDG;nParticles").c_str());
	//h_pdg->Draw();
	//t->SaveAs("plots/PDG_Charged_Particles.png");

	TCanvas *t1 = new TCanvas("t1","t1");
    t1->cd();
	h_Energy->SetMarkerStyle(kFullCircle);
    h_Energy->SetTitle((title_collision+ " Energy All  particles;E(GeV) ;nParticles").c_str());
	h_Energy->Draw();
	//t1->SaveAs("plots/Energy_charged_particles.png");
	//h_Energy->SetTitle("e - Pb Energy pdg:2112; E(GeV);dE/dN (1/GeV)");
	//h_Energy->Draw("PLC PMC");

	TCanvas *t2 = new TCanvas("t2","t2");
	t2->cd();
	h_pt->SetTitle((title_collision+ " Pt All particles; Pt(GeV);nParticles").c_str());
	h_pt->Draw();
	//t2->SaveAs("plots/Pt_charged_particles.png");

	TCanvas *t3 = new TCanvas("t3","t3");
	t3->cd();
	h_eta->SetTitle((title_collision+ " Eta All particles; Eta ;nParticles").c_str());
	h_eta->Draw();
	//t3->SaveAs("plots/Eta_charged_part.png");

	TCanvas *t4 = new TCanvas("t4","t4");
	t4->cd();
	h_phi->SetTitle((title_collision+ " Phi All  particles; Phi ;nParticles").c_str());
	h_phi->Draw();
	//t4->SaveAs("plots/Phi_charged_part.png");

	
	//return;

	TCanvas *t6 = new TCanvas("t6","t6");
	t6->cd();
	h_mass->SetTitle((title_collision+ " Mass All  particles; Mass(GeV);nParticles").c_str());
	h_mass->Draw();
	//t6->SaveAs("plots/mass_charged_part.png");
  	
	TCanvas *t7 = new TCanvas("t7","t7");
	t7->cd();
	h_nTracks->SetTitle((title_collision+ " nTracks All particles; nTracks;nEvents").c_str());
	h_nTracks->Draw();
	//t7->SaveAs("plots/nTracks_charged_part.png");

	TCanvas *t8 = new TCanvas("t8","t8");
	t8->cd();
	h_b->SetTitle((title_collision+ " Impact Parameter All particles; b ; nEvents").c_str());
	h_b->Draw();
	//t8->SaveAs("plots/b_charged_part.png");

	TCanvas *t9 = new TCanvas("t9","t9");
	t9->cd();
	h_trueQ2->SetTitle((title_collision+ " Q^2 ; Q^2 ; nEvents").c_str());
    h_trueQ2->Draw();
	//t9->SaveAs("plots/q2.png");

	//TCanvas *t10 = new TCanvas("t10","t10");
	//t10->cd();
	//h_nTracks_cuts->SetTitle((title_collision+ " nTracks All particles  ; nTracks ; nEvents").c_str());
	//h_nTracks_cuts->Draw();
	//t10->SaveAs("plots/nTracks_hadrons.png");
	
	TCanvas *t11 = new TCanvas("t11","t11");
	t11->cd();
	h_Energy_neutron->SetTitle((title_collision+ " Energy pdg:2112 ; E(GeV) ; nParticles").c_str());
	h_Energy_neutron->Draw();
	//t11->SaveAs("plots/Energy_neutron.png");

	TCanvas *t12 = new TCanvas("t12","t12");
	t12->cd();
	h_nTracks_neutron->SetTitle((title_collision+ " nTracks  pdg:2112 ; nTracks ; nEvents").c_str());
	h_nTracks_neutron->Draw();
	//t12->SaveAs("plots/nTracks_neutron.png");

	TCanvas *t13 = new TCanvas("t13","t13");
	t13->cd();
	h_pt_neutron->SetTitle((title_collision+ " Pt  pdg:2112 ; Pt ; nParticles").c_str());
	h_pt_neutron->Draw();
	//t13->SaveAs("plots/Pt_neutron.png");
	
	///
	///TCanvas *t14 = new TCanvas("t14","t14");
    ///t14->cd();
	///h_Energy_proton->SetTitle((title_collision+"Energy pdg:2212 ; E(GeV) ; nParticles").c_str());
	///h_Energy_proton->Draw();
	/////t14->SaveAs("plots/Energy_proton.png");

	///TCanvas *t15 = new TCanvas("t15","t15");
	///t15->cd();
	///h_nTracks_proton->SetTitle((title_collision+" nTracks  pdg:2212 ; nTracks ; nEvents").c_str());
	///h_nTracks_proton->Draw();
	//t15->SaveAs("plots/nTracks_proton.png");

	///TCanvas *t16 = new TCanvas("t16","t16");
	///t16->cd();
	///h_pt_proton->SetTitle((title_collision+" Pt  pdg:2212 ; Pt ; nParticles").c_str());
	///h_pt_proton->Draw();
	
	// out->Write();

}
