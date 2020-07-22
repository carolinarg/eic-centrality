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
#include "TProfile.h"

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
#define ENERGY_TARGET 490.0
//#define MASS_NUCLEON  0.93903412 //average (N*MASS_NEUTRON+Z*MASS_PROTON)/Z+N

using namespace std;
using namespace erhic;

//11631 number events D
//10805 number events Xe



void runEICTree(const TString filename="/Users/carorg/projects/eic/eic-centrality/root_files/root-files-2020/eXe_490GeV_fixed_target_40k", const int nEvents = 40000){

	TChain *tree = new TChain("EICTree");
	tree->Add( filename+".root" );
	
	EventBeagle* event(NULL);
	tree->SetBranchAddress("event", &event);
	
	TFile *out = new TFile ("/Users/carorg/projects/eic/eic-centrality/root_files/analysis-ratios/Analysis_ratio_eXe_490GeV_fixed_target_40k_TEST.root", "RECREATE");
	//####################  Creating Histograms ####################################

	TH1F *h_nTracks = new TH1F("h_nTracks", "h_nTracks", 100, 90, 300);
	TH1F *h_pdg = new TH1F("h_pdg", "h_pdg", 500,-2500, 2500);
	TH1F *h_Energy = new TH1F("h_Energy", "h_Energy", 100, 0,6);
	TH1F *h_dE_dN = new TH1F("h_idE_dN", "h_dE_dN", 100, 0, 200);
	TH1F *h_Status = new TH1F("h_Status", "h_Status", 50, 0, 25);
	TH1F *h_pt = new TH1F("h_pt","h_pt", 100,0, 2);
	TH1F *h_eta = new TH1F("h_eta", "h_eta", 100, -10, 10);
	TH1F *h_phi = new TH1F("h_phi", "h_phi", 100,-10, 10);
	TH1F *h_rapidity = new TH1F("h_rapidity","h_rapidity", 100, -8, 4);
	TH1F *h_mass = new TH1F("h_mass", "h_mass", 100, 0, 10);
	TH1F *h_b = new TH1F("h_b", "h_b", 100, 0, 10);
	TH1F *h_trueQ2 = new TH1F("h_trueQ2", "h_trueQ2", 100, 0, 20);
	TH1F *h_nTracks_cuts = new TH1F("h_nTracks_cuts","h_nTracks_cuts", 50, 0,100);
	TH1F *h_Energy_neutron = new TH1F("h_Energy_neutron", "h_Energy_neutron", 100, 0, 3);
	TH1F *h_nTracks_neutron = new TH1F("h_nTracks_neutron", "h_nTracks_neutron", 100, 0, 50);
	TH1F *h_pt_neutron = new TH1F("h_pt_neutron", "h_pt_neutron", 100, 0, 2);
	TH1F *h_Energy_proton = new TH1F("h_Energy_proton", "h_Energy_proton", 100, 0, 3);
	TH1F *h_Energy_proton_gt = new TH1F("h_Energy_proton_gt", "h_Energy_proton_gt", 100, 0, 3);
	TH1F *h_xf = new TH1F("h_xf", "h_xf", 20, -4, 4);
	
	TH1F *h_pt_proton = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
	TH1F *h_pt_proton_gt= new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
	TH1F *h_nTracks_gt = new TH1F("h_nTracks_gt","h_nTracks_gt",100, 0, 50);
	TH1F *h_pt_gt = new TH1F("h_pt_gt", "h_pt_gt", 100, 0, 1);
	TH1F *h_Energy_gt = new TH1F("h_Energy_gt", "h_Energy_gt", 100,0, 6);
	TH1F *h_rapidity_gt = new TH1F("h_rapidity_gt","h_rapidity_gt",100,-10,10);
	TH1F *h_nTracks_p_gt= new TH1F("h_nTracks_p_gt","h_nTracks_p_gt", 11, -0.5, 10.5); //****
	TH1F *h_nTracks_proton = new TH1F("h_nTracks_proton", "h_nTracks_proton", 11, -0.5, 10.5); //****

	TH1F *h_nTracks_charged = new TH1F("h_nTracks_charged","h_nTracks_charged", 21, -0.5, 20.5);
	TH1F *h_nTracks_charged_gt = new TH1F("h_nTracks_charged_gt","h_nTracks_charged_gt", 21, -0.5, 20.5);
	TH1F *h_nTracks_mesons = new TH1F("h_nTracks_mesons","h_nTracks_mesons",21,-0.5, 20.5);
	TH1F *h_nTracks_mesons_gt = new TH1F("h_nTracks_mesons_gt","h_nTracks_mesons_gt",21,-0.5, 20.5);

		
	
	TH1F *h_w2 = new TH1F("h_W2","h_W2", 100, 0, 1000);
	TH1F *h_nu = new TH1F("h_nu", "h_nu",100,50,400);
	TH1F *h_bx = new TH1F("h_bx","h_bx",100,0,1);

	TH2F *h_rap_eta = new TH2F("h_rap_eta","h_rap_eta", 100, 0, 1, 100, 0, 1);
	TH2F *h_mass_mom = new TH2F("h_mass_mom","h_mass_mom", 100, 0, 50, 100, 0, 2);
	TH2F *h_mass_pt = new TH2F("h_mass_pt","h_mass_pt", 100, 0, 6, 100, 0, 2);

	TH2F *h_rap_gap = new TH2F("h_rap_gap","rapidity gaps",100,-7,3,11,-0.5,10.5);

	TH1F *h_mom = new TH1F("h_mom", "h_mom", 100, 0, 500);
	TH1F *h_Energy_cms = new TH1F("h_Energy_cms","h_Energy_cms",35,0,35);
	TH1F *h_p_cms = new TH1F("h_p_cms","h_p_cms",20,0,20);
	TH1F *h_p_long_cms = new TH1F("h_p_long_cms","h_p_long_cms", 20,0,20);
	TH1F *h_x_f_cms = new TH1F("h_x_f_cms","h_x_f_cms", 20, -4,4);
	TH1F *h_theta = new TH1F("h_theta","h_theta",100,0,360);
	TH1F *h_y_cms = new TH1F("h_y_cms","h_y_cms", 100,-8,4);
	TH1F *h__y_cms = new TH1F("h__y_cms","h__y_cms", 100,-8,4);


	
	TProfile *prof = new TProfile("prof", "Rapidity;<ng>", 21, -7, 3);
	TH2F *hist = new TH2F("hist","Rapidity;<ng>", 21, -4, 4,11,-0.5,10.5);
	
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

			
		//TLorentzVector e_scattered(0.,0.,0.,0.); //final electron

		
		

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
		if( trueW2 < 16 || trueW2 > 900) continue;
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
		int counter_mesons = 0;
		int counter_mesons_gt = 0;
		int counter_charged_particles = 0;
		int counter_charged_particles_gt = 0;
		double rap_gt = 0;
		double Energy_cms = 0;
		double p_cms = 0;
		double x_f_cms = 0;
		double p_long_cms = 0;
		double beta_target_cms = 0.9980;
		double y_cms = 0;
		double y_cms_target = 0;
		double _y_cms = 0;


		double e_scattered_x = 0;
		double e_scattered_y = 0;
		double e_scattered_z = 0;
		double e_scattered_E = 0;

		double p_gamma = 0;
		double p_gamma_x = 0;
		double p_gamma_y = 0;
		double p_gamma_z = 0;
		double p_gamma_E = 0;


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
			double p = particle->GetP();
			double Energy = particle->GetE();
			double xf = particle->GetXFeynman();
			double px = particle->GetPx();
			double py = particle->GetPy();
			double pz = particle->GetPz();

			TLorentzVector e_beam(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON)); //initial electron
			TLorentzVector p_beam(0.,0.,pztarg,sqrt(pztarg*pztarg+MASS_NUCLEON*MASS_NUCLEON)); //target

			double e_beam_x = e_beam.X();
			double e_beam_y = e_beam.Y();
			double e_beam_z = e_beam.Z();
			double e_beam_E = e_beam.E();

			if (pdg==11 && status == 1){

				TLorentzVector e_scattered(px,py,pz,sqrt(pz*pz+MASS_ELECTRON*MASS_ELECTRON));

				e_scattered_x = e_scattered.X();
				e_scattered_y = e_scattered.Y();
				e_scattered_z = e_scattered.Z();
				e_scattered_E = e_scattered.E();
				

			}
			TVector3 e_u= e_scattered.Vect().Unit();
			//cout << "px = "<<e_scattered_x << "py = "<<e_scattered_y << "pz = "<<e_scattered_z <<endl;
			//Virtual photon

			TLorentzVector p_gamma(e_beam_x-e_scattered_x, e_beam_y-e_scattered_y, e_beam_z-e_scattered_z, sqrt((e_beam_x-e_scattered_x)*(e_beam_x-e_scattered_x)));

			p_gamma_x = p_gamma.X();
			p_gamma_y = p_gamma.Y();
			p_gamma_z = p_gamma.Z();
			p_gamma_E = p_gamma.E();

			//Unit vectors
			TVector3 uz = p_gamma.Vect().Unit();
			//TVector3 ux = (e_beam.Vect().Cross(e_scattered.Vect())).Unit();

			//Rotations
			//ux.Rotate(3.*TMath::Pi()/2,uz);

			
			// p_cms = p*MASS_NUCLEON/Energy_cms;
			// p_long_cms = p_cms*TMath::Cos(theta);	
			// x_f_cms = 2*p_long_cms/TMath::Sqrt(trueW2);

			// //y_cms_target=0.5*log(1+beta_target_cms/1-beta_target_cms);
			// //cout << "y_cms_target" << y_cms_target << endl;
			
			// y_cms = rap - 3.453; 
			// //cout << beta_target_cms<< endl;
			// _y_cms= 0.5*log(Energy_cms+p_long_cms/Energy_cms-p_long_cms);


			
			
			//do analysis track-by-track
			
			//particle cuts
			if (status != 1) continue; //final state particles
			//if (pt < 0.8) continue;
			if (p > 0.2 && p < 0.8) //Grey Tracks
			{
				gt_counter++;
				h_pt_gt->Fill(pt);
				h_Energy_gt->Fill(Energy);
				if (pdg == 2212)
				{
					counter_p_gt++;
					rap_gt = rap;
					h_Energy_proton_gt->Fill(Energy);
					h_pt_proton_gt->Fill(pt);
					h_rapidity_gt->Fill(rap);
					
				}
				if (pdg ==211)
				{
					counter_mesons_gt++;

				}

				counter_charged_particles_gt=counter_p_gt+counter_mesons_gt;
			}
			//############### Particles without the Pt Cut ###########

			if (pdg == 211)
			{
				counter_mesons++;

			}			

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
			

			counter_charged_particles = counter_p+counter_mesons;


					
			counter++;
			counter_total_particles++;
			h_Energy->Fill(Energy);
			h_Energy_cms->Fill(Energy_cms);
			h_p_cms->Fill(p_cms);
			h_p_long_cms->Fill(p_long_cms);
			h_x_f_cms->Fill(x_f_cms);
			h_theta->Fill(theta);
						
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
			h_mom->Fill(p);
			h_rap_eta->Fill(rap,eta,1);
			h_mass_mom->Fill(p,mass,1);
			h_mass_pt->Fill(pt,mass,1);
			h_rap_gap->Fill(counter_p_gt,rap);
			h_xf->Fill(xf);
			h_y_cms->Fill(y_cms);
			h__y_cms->Fill(_y_cms);


		} // end of particle loop

		// for(int j(0); j < nParticles; ++j ) {

		// prof->Fill(rap, counter_gt);
		// hist->Fill(rap, counter_gt);

		// }
	
		
		//fill histograms
		h_nTracks->Fill(nParticles);
		h_b->Fill(impact_parameter);
		h_trueQ2->Fill(trueQ2);
		h_nTracks_cuts->Fill(counter);
		h_nTracks_neutron->Fill(counter_neu);
		h_nTracks_proton->Fill(counter_p);//all tracks protons
		h_nTracks_gt->Fill(gt_counter);
		h_nTracks_p_gt->Fill(counter_p_gt);//grey tracks protons
		h_nTracks_charged->Fill(counter_charged_particles);
		h_nTracks_mesons->Fill(counter_mesons);
		h_nTracks_mesons_gt->Fill(counter_mesons_gt);
		h_nTracks_charged_gt->Fill(counter_charged_particles_gt);
		//cout << "number mesons: " << counter_mesons <<endl;
		
		
		h_w2->Fill(trueW2);
		h_nu->Fill(trueNu);
		h_bx->Fill(trueX);

		//cout << "e_beam" << e_beam[3] <<endl;
		//cout << "nucleon_target" << p_beam[3] <<endl;
		//cout << "Energy_cms: " << Energy_cms << endl;

	}// end event loop

	cout << "Number of particles: " << counter_total_particles << endl;
	cout << "Events after event cuts : " << counter_event << endl;
	cout << "mean value of All protons: " <<h_nTracks_proton->GetMean() << endl;

	
	//################################### Drawing for e - Xe #############################
	TCanvas *t15 = new TCanvas("t15","t15");
	t15->Divide(2,1);
	t15->cd(1);
	gPad->SetLogy();
	h_nTracks_p_gt->SetTitle((title_collision+ "Grey Tracks Protons ; Grey Tracks ; nEvents").c_str());
	h_nTracks_p_gt->SetMarkerStyle(kStar);
	h_nTracks_p_gt->SetMarkerColor(kBlue);
	h_nTracks_p_gt->Draw("P");
	
	t15->cd(2);
	gPad->SetLogy();
	h_nTracks_proton->SetTitle((title_collision+ " NTracks Protons ; nTracks; nEvents").c_str());
	h_nTracks_proton->SetMarkerStyle(kStar);
	h_nTracks_proton->SetMarkerColor(kBlue);
	h_nTracks_proton->Draw("P");

	

	TCanvas *c7 = new TCanvas("c7","c7");
	c7->Divide(2,1);
	c7->cd(1);
	gPad->SetLogy();
	h_nTracks_mesons_gt->SetTitle((title_collision+ " Grey Tracks Pions ;  Grey Tracks; nEvents").c_str());
	h_nTracks_mesons_gt->SetMarkerStyle(kStar);
	h_nTracks_mesons_gt->SetMarkerColor(kRed);
	h_nTracks_mesons_gt->Draw("P");

	c7->cd(2);
	gPad->SetLogy();
	h_nTracks_mesons->SetTitle((title_collision+ " NTracks Pions ; nTracks; nEvents").c_str());
	h_nTracks_mesons->SetMarkerStyle(kStar);
	h_nTracks_mesons->SetMarkerColor(kRed);
	h_nTracks_mesons->Draw("P");


	TCanvas *c8 = new TCanvas("c8","c8");
	c8->Divide(2,1);
	c8->cd(1);
	gPad->SetLogy();
	h_nTracks_charged_gt->SetTitle((title_collision+ " Grey Tracks Charged Particles ; Grey Tracks; nEvents").c_str());
	h_nTracks_charged_gt->SetMarkerStyle(kStar);
	h_nTracks_charged_gt->SetMarkerColor(kRed);
	h_nTracks_charged_gt->Draw("P");

	c8->cd(2);
	gPad->SetLogy();
	h_nTracks_charged->SetTitle((title_collision+ " NTracks Charged Particles ; nTracks; nEvents").c_str());
	h_nTracks_charged->SetMarkerStyle(kStar);
	h_nTracks_charged->SetMarkerColor(kRed);
	h_nTracks_charged->Draw("P");
	
	

	TCanvas *t14 = new TCanvas("t14","t14");
	t14->cd();
	h_Energy_proton->SetTitle((title_collision+ " Energy All Protons ; E(GeV) ; nParticles").c_str());
	h_Energy_proton->Draw();
				    

	TCanvas *t28 = new TCanvas("t28","t28");
	t28->cd();
	h_Energy_proton_gt->SetTitle((title_collision+ " Energy Grey Tracks Protons  ; E(GeV) ; nParticles").c_str());
	h_Energy_proton_gt->Draw();
	
	
	TCanvas *t22 = new TCanvas("t22","t22");
	t22->cd();
	h_rapidity_gt->SetTitle((title_collision+ " Rapidity Grey Tracks Protons; Grey Tracks; nParticles").c_str());
	h_rapidity_gt->Draw();


	//TCanvas *t21 = new TCanvas("t21","t21");
	//t21->cd();
	//h_nTracks_gt->SetTitle((title_collision+ " All stable Grey Tracks;Grey Tracks; nEvents").c_str());
	//h_nTracks_gt->Draw();

	
	
	TCanvas *t16 = new TCanvas("t16","t16");
	t16->cd();
	h_pt_proton->SetTitle((title_collision+ " Pt Protons ; Pt ; nParticles").c_str());
	h_pt_proton->Draw();
	
	TCanvas *t29 = new TCanvas("t29","t29");
	t29->cd();
	h_pt_proton_gt->SetTitle((title_collision+ " Pt Grey Tracks Protons ; Pt ; nParticles").c_str());
	h_pt_proton_gt->Draw();


	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd();
	h_rap_gap->SetTitle((title_collision+ " Grey Tracks Protons vs Rapidity ; y ; Grey Tracks").c_str());
	h_rap_gap->SetMarkerStyle(kFullCircle);
	h_rap_gap->Draw("P");
	
	TCanvas *c2 = new TCanvas("c2","c2");
	c2->cd();
	h_Energy_cms->SetTitle((title_collision+ " Energy CMS; Energy CMS ;nParticles").c_str());
	h_Energy_cms->Draw();

	TCanvas *c3 = new TCanvas("c3","c3");
	c3->cd();
	//gPad->SetLogy();
	h_p_cms->SetTitle((title_collision+ " Momentum CMS; p CMS ;nParticles").c_str());
	h_p_cms->Draw();

	TCanvas *c4 = new TCanvas("c4","c4");
	c4->cd();
	h_p_long_cms->SetTitle((title_collision+ " Longitudinal Momentum CMS; p longitudinal CMS ;nParticles").c_str());
	h_p_long_cms->Draw();

	TCanvas *c5 = new TCanvas("c5","c5");
	c5->Divide(2,1);
	c5->cd(1);
	//gPad->SetLogy();
	h_xf->SetTitle((title_collision+ " Feynman - X BeAGLE; Xf;nParticles").c_str());
	h_xf->Draw();

	c5->cd(2);
	h_x_f_cms->SetTitle((title_collision+ " Feynman - X ; Xf;nParticles").c_str());
	h_x_f_cms->Draw();

	TCanvas *c6 = new TCanvas("c6","c6");
	c6->cd();
	h_theta->SetTitle((title_collision+ " THETA; THETA ;nParticles").c_str());
	h_theta->Draw();

	
	//-------------------------------------- All stable Particles ----------------------------------------------

	// TCanvas *t5 = new TCanvas("t5","t5");
	// t5->Divide(2,1);
	// t5->cd(1);
	// h_rapidity->SetTitle((title_collision+ " Rapidity All stable and final particles; Rapidity ;nParticles").c_str());
	// h_rapidity->Draw();

	// t5->cd(2);
	// h_y_cms->SetMarkerStyle(kStar);
	// h_y_cms->SetTitle((title_collision+ " Rapidity All stable and final particles; Rapidity CMS ;nParticles").c_str());
	// h_y_cms->Draw();

	// t5->cd(3);
	// h_y_cms->SetTitle((title_collision+ " Rapidity All stable and final particles; Rapidity CMS 2 ;nParticles").c_str());	
	// h__y_cms->Draw();

	// t5->cd(2);
	// h_y_cms->SetTitle((title_collision+ " Rapidity CMS; Rapidity CMS ;nParticles").c_str());
	// h_y_cms->Draw();

	TCanvas *t26 = new TCanvas("t26","t26");
	t26->cd();
	h_nu->SetTitle((title_collision+ " Nu All Particles ; Nu; nEvents").c_str());
	h_nu->Draw();

	TCanvas *t7 = new TCanvas("t7","t7");
	t7->cd();
	h_nTracks->SetTitle((title_collision+ " nTracks All particles ; nTracks;nEvents").c_str());
	
	//h_nTracks->Scale(1/10805);
	h_nTracks->Draw();

	TCanvas *t10 = new TCanvas("t10","t10");
	t10->cd();
	h_nTracks_cuts->SetTitle((title_collision+ " nTracks All particles with events cuts  ; nTracks ; nEvents").c_str());
	h_nTracks_cuts->Scale(0.0000921);
	h_nTracks_cuts->Draw();

	TCanvas *t19 = new TCanvas("t19","t19");
	t19->cd();
	h_mom->SetTitle((title_collision+ " All stable particles P; p(GeV); nParticles").c_str());
	gPad->SetLogy();
	h_mom->Draw();

	TCanvas *t = new TCanvas("t","t");
	t->cd();
	h_pdg->SetTitle((title_collision+ " PDG All particles; PDG;nParticles").c_str());
	h_pdg->Draw();

	out->Write();

	return;
	

	

	

	
	TCanvas *t17 = new TCanvas("t17","t17");
	t17->cd();
	h_rap_eta->SetTitle((title_collision+ " All stable particles Rapidity vs Eta; rapidity; eta").c_str());
	h_rap_eta->Draw();
	
	// TCanvas *t18 = new TCanvas("t18","t18");
	// t18->cd();
	// h_mass_mom->SetTitle((title_collision+ " All stable particles; p; M").c_str());
	// h_mass_mom->Draw();
	

	

	TCanvas *t20 = new TCanvas("t20","t20");
	t20->cd();
	h_mass_pt->SetTitle((title_collision+ " All stable particles M vs Pt; Pt(GeV); Mass").c_str());
	h_mass_pt->Draw();
		
	
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
	

	TCanvas *t6 = new TCanvas("t6","t6");
	t6->cd();
	h_mass->SetTitle((title_collision+ " Mass All  particles; Mass(GeV);nParticles").c_str());
	h_mass->Draw();
	//t6->SaveAs("plots/mass_charged_part.png");
  	

	//t7->SaveAs("plots/nTracks_charged_part.png");

	TCanvas *t8 = new TCanvas("t8","t8");
	t8->cd();
	h_b->SetTitle((title_collision+ " Impact Parameter All particles; b ; nEvents").c_str());
	h_b->Draw();
	//t8->SaveAs("plots/b_charged_part.png");

	TCanvas *t9 = new TCanvas("t9","t9");
	t9->cd();
	h_trueQ2->SetTitle((title_collision+ " Q^2 All  ; Q^2 ; nEvents").c_str());
    h_trueQ2->Draw();
	//t9->SaveAs("plots/q2.png");

	TCanvas *t27 = new TCanvas("t27","t27");
	t27->cd();
	h_bx->SetTitle((title_collision+ " Bj x; x; nEvents").c_str());
	h_bx->Draw();
	
	
	

	
	TCanvas *t25 =new TCanvas("t25","t25");
	t25->cd();
	h_w2->SetTitle((title_collision+ " W^2 All; w^2; nEvents").c_str());
	h_w2->Draw();
	//-------------------------------------------- Neutrons ----------------------------------------------
	
	TCanvas *t11 = new TCanvas("t11","t11");
	t11->cd();
	h_Energy_neutron->SetTitle((title_collision+ " Energy Neutrons ; E(GeV) ; nParticles").c_str());
	h_Energy_neutron->Draw();
	//t11->SaveAs("plots/Energy_neutron.png");

	TCanvas *t12 = new TCanvas("t12","t12");
	t12->cd();
	h_nTracks_neutron->SetTitle((title_collision+ " nTracks Neutrons ; nTracks ; nEvents").c_str());
	h_nTracks_neutron->Draw();
	//t12->SaveAs("plots/nTracks_neutron.png");

	TCanvas *t13 = new TCanvas("t13","t13");
	t13->cd();
	h_pt_neutron->SetTitle((title_collision+ " Pt Neutrons ; Pt ; nParticles").c_str());
	h_pt_neutron->Draw();
	//t13->SaveAs("plots/Pt_neutron.png");	

	
	//out->Write();

}
