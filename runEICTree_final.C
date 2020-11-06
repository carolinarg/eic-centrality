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

#include "TLorentzVector.h"
#include "TBranchElement.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <TH1I.h>
#include "RootUtils.C"



#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI     3.09688
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

map<string, TH1I*> hist_distribution_tracks;

const vector<string> GreyTracks = {"1","2","3","4","5","6","7","8"};


void runEICTree_final(const TString filename="/Users/carorg/projects/eic/eic-centrality/root_files/root-files-2020/eXe_490GeV_fixed_target_40k", const int nEvents = 40000){

    TChain *tree = new TChain("EICTree");
    tree->Add( filename+".root" );
    
    EventBeagle* event(NULL);
    tree->SetBranchAddress("event", &event);
    
    TFile *out = new TFile ("/Users/carorg/projects/eic/eic-centrality/root_files/analysis-ratios/Analysis_ratio_eXe_490GeV_fixed_target_40k_TEST.root", "RECREATE");
    //####################  Creating Histograms ####################################

    //TTree *t = new TTree("t","t filter");

    for (const string& gt : GreyTracks) {
    hist_distribution_tracks[gt] = new TH1I(Form("hist_distribution_tracks%s", gt.c_str()), "title;Multiplicity Distributions;Events", 100, -0.5, 10.5);
    }
    
    TH1F *h_nTracks = new TH1F("h_nTracks", "h_nTracks", 100, 90, 300);
    TH1F *h_pdg = new TH1F("h_pdg", "h_pdg", 500,-2500, 2500);
    TH1F *h_Energy = new TH1F("h_Energy", "h_Energy", 100, 0,6);
    TH1F *h_dE_dN = new TH1F("h_idE_dN", "h_dE_dN", 100, 0, 200);
    TH1F *h_Status = new TH1F("h_Status", "h_Status", 50, 0, 25);
    TH1F *h_pt = new TH1F("h_pt","h_pt", 100,0, 2);
    TH1F *h_pz = new TH1F("h_pz","h_pz", 80,-40, 40);
    TH1F *h_eta = new TH1F("h_eta", "h_eta", 100, -10, 10);
    
    TH1F *h_rapidity = new TH1F("h_rapidity","h_rapidity", 100, -8, 5);
    TH1F *h_mass = new TH1F("h_mass", "h_mass", 100, 0, 10);
    TH1F *h_b = new TH1F("h_b", "h_b", 100, 0, 10);
    TH1F *h_trueQ2 = new TH1F("h_trueQ2", "h_trueQ2", 100, 0, 20);
    TH1F *h_nTracks_cuts = new TH1F("h_nTracks_cuts","h_nTracks_cuts", 50, 0,100);
    TH1F *h_Energy_proton = new TH1F("h_Energy_proton", "h_Energy_proton", 100, 0, 3);
    TH1F *h_Energy_proton_gt = new TH1F("h_Energy_proton_gt", "h_Energy_proton_gt", 100, 0, 3);
    TH1F *h_xf = new TH1F("h_xf", "h_xf", 100, -1.1, 1.1);
    
    TH1F *h_pt_proton       = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
    TH1F *h_pt_proton_gt    = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
    TH1F *h_nTracks_gt      = new TH1F("h_nTracks_gt","h_nTracks_gt",100, 0, 50);
    TH1F *h_pt_gt           = new TH1F("h_pt_gt", "h_pt_gt", 100, 0, 1);
    TH1F *h_Energy_gt       = new TH1F("h_Energy_gt", "h_Energy_gt", 100,0, 6);
    TH1F *h_rapidity_gt     = new TH1F("h_rapidity_gt","h_rapidity_gt",100,-10,10);
    TH1F *h_nTracks_p_gt    = new TH1F("h_nTracks_p_gt","h_nTracks_p_gt", 11, -0.5, 10.5); //****
    TH1F *h_nTracks_proton  = new TH1F("h_nTracks_proton", "h_nTracks_proton", 11, -0.5, 10.5); //****
    TH1F *h_pz_proton       = new TH1F("h_pz_proton","h_pz_proton", 40,-10,10);
    TH1F *h_pz_proton_cms  = new TH1F("h_pz__proton_cms","h_pz__proton_cms", 40,-10,10);

    
   
    TH1F *h_nTracks_pos_pions = new TH1F("h_nTracks_pos_pions","h_nTracks_pos_pions",21,-0.5, 20.5);
    TH2F *h_nTracks_Grey_tracks= new TH2F("h_nTracks_Grey_tracks","h_nTracks_Grey_tracks", 11,0,11,60,0,15);
    TProfile *prof_nTracks_Grey_tracks= new TProfile("prof_nTracks_Grey_tracks","prof_nTracks_Grey_tracks", 11,-0.5,10.5); 
    TProfile *prof_nTracks_target_positive = new TProfile("prof_nTracks_target_positive","prof_nTracks_target_positive", 11,-0.5, 10.5);
    TProfile *prof_nTracks_central_positive = new TProfile("prof_nTracks_central_positive","prof_nTracks_central_positive", 11,-0.5, 10.5);
    TProfile *prof_nTracks_projectile_positive = new TProfile("prof_nTracks_projectile_positive","prof_nTracks_projectile_positive", 11,-0.5, 10.5);
    TProfile *prof_nTracks_target_negative = new TProfile("prof_nTracks_target_negative","prof_nTracks_target_negative", 11,-0.5, 10.5);
    TProfile *prof_nTracks_central_negative = new TProfile("prof_nTracks_central_negative","prof_nTracks_central_negative", 11,-0.5, 10.5);
    TProfile *prof_nTracks_projectile_negative = new TProfile("prof_nTracks_projectile_negative","prof_nTracks_projectile_negative", 11,-0.5, 10.5);
    TProfile *prof_nTracks_target_charged = new TProfile("prof_nTracks_target_charged","prof_nTracks_target_charged", 11,-0.5, 10.5);
    TProfile *prof_nTracks_central_charged = new TProfile("prof_nTracks_central_charged","prof_nTracks_central_charged", 11,-0.5, 10.5);
    TProfile *prof_nTracks_projectile_charged = new TProfile("prof_nTracks_projectile_charged","prof_nTracks_projectile_charged", 11,-0.5, 10.5);

        
    
    TH1F *h_w2 = new TH1F("h_W2","h_W2", 100, 0, 1000);
    TH1F *h_w2_computed = new TH1F("h_w2_computed","h_w2_computed", 100, 0, 1000);
    TH1F *h_W2_BeAGLE_computed = new TH1F("h_W2_BeAGLE_computed","h_W2_BeAGLE_computed", 100,0,1000);
    TH1F *h_nu = new TH1F("h_nu", "h_nu",100,50,400);
    TH1F *h_bx = new TH1F("h_bx","h_bx",100,0,1);

    TH2F *h_rap_eta = new TH2F("h_rap_eta","h_rap_eta", 100, 0, 1, 100, 0, 1);
    TH2F *h_mass_mom = new TH2F("h_mass_mom","h_mass_mom", 100, 0, 50, 100, 0, 2);
    TH2F *h_mass_pt = new TH2F("h_mass_pt","h_mass_pt", 100, 0, 6, 100, 0, 2);

    TH2F *h_rap_gap = new TH2F("h_rap_gap","rapidity gaps",100,-7,3,11,-0.5,10.5);

    TH1F *h_mom = new TH1F("h_mom", "h_mom", 60, -30, 30);
    
    
    TH1F *h_theta = new TH1F("h_theta","h_theta",100,0,360);

    TH1F *h_p_cms       = new TH1F("h_p_cms", "h_p_cms", 60,0,30);
    TH1F *h_Energy_cms = new TH1F("h_Energy_cms","h_Energy_cms", 100,0,6);
    TH1F *h_pz_cms      = new TH1F("h_pz_cms","h_pz_cms",80, -40,40);
    TH1F *h_pt_cms      = new TH1F("h_pt_cms","h_pt_cms",100, 0,2);
    TH1F *h_rap_cms     = new TH1F("h_rap_cms","h_rap_cms",48,-4,4);
    TH1F *h_xf_cms      = new TH1F("h_xf_cms","h_xf_cms",100,-1.1,1.1);
    TH1F *h_xf_cms_s      = new TH1F("h_xf_cms_s","h_xf_cms_s",100,-1.1,1.1);


    TH1F *h_xf_2 = new TH1F("h_xf_2","h_xf_2",100,-1.1,1.1);
    TH1F *h_rap_lab = new TH1F("h_rap_lab","h_rap_lab", 48,-4,4);
    TH1F *h_rap_lab_computed = new TH1F("h_rap_lab","h_rap_lab", 48,-4,4);
    TH1F *h_rap_e_lab = new TH1F("h_rap_e_lab","h_rap_e_lab", 68,-6,4);

    
    TProfile *prof = new TProfile("prof", "Rapidity;<ng>", 21, -7, 3);
    TH2F *hist = new TH2F("hist","Rapidity;<ng>", 21, -4, 4,11,-0.5,10.5);
    
    int counter_total_particles = 0;
    int counter_event           = 0;
    
    
    string title_collision = "e - Xe";

    TLorentzVector e_scattered;
    TLorentzVector e_beam;
    TLorentzVector p_gamma;
    TLorentzVector p_final;
    TLorentzVector p_beam;
    TLorentzVector p_gamma_beagle;

    TLorentzVector p_final_rotated_x;
    TLorentzVector p_final_rotated_z;
    TLorentzVector p_final_boosted_z; 

    TVector3 ux;
    TVector3 uz;
    TVector3 uz_beagle;
    TVector3 ux_beagle;

    TVector3 p_final_rotated;

    TMatrixD rot_x(3,3);
    TMatrixD rot_z(3,3);
    TMatrixD boost_z(4,4);
    
    TMatrixD P_1;
    
    double rap_lab              = 0;
    double rap_cms              = 0;
    double rap_lab_computed     = 0;
    double xf_cms               = 0;
    double W2                   = 0;
    double W2_BeAGLE_computed   = 0;
    double W2_boosted           = 0;
    double rap_e_lab            = 0;
    double Beta_BeAGLE          = 0;
    double Beta                 = 0;
    int counter_final_e         = 0;
    int p_gamma_counter         = 0;
    int trys =0, trys1=0;
    int protons_lab = 0, neutrons_lab = 0;
    int struck_neutron          = 0;
    int struck_proton           = 0;
    // int counter_p               = 0;
    // int counter_p_gt            = 0;

    double sum_N_prot_gt_0=0;
    double average_N_prot = 0;
    
    
    int S_boosted               = 0;
    double S                    = 0;
    double gamma                = 0;
    int cms_protons_xf_cut      = 0;
    
    double mean_value_trks_D    = 1.17;

    auto graph_n_vs_ng = TH1TOTGraph(prof_nTracks_Grey_tracks);

    for(int i(0); i < nEvents; ++i ) {
      
        // Read the next entry from the tree.
        tree->GetEntry(i);

                
        double pzlep = event->pzlep;
        double pztarg = event->pztarg;
        int struck_nucleon = event->nucleon;
        double MASS_NUCLEON = MASS_PROTON;
        if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;

                

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
        if( event_process != 99 ) continue;       //DIS process for pythia
        if( trueQ2 < 1 ) continue;
        if( trueX < 0.002 ) continue;
        if( trueNu < 50 || trueNu > 400 )continue;
        if( trueW2 < 16 || trueW2 > 900) continue;
        //if( trueQ2 < 1. || trueQ2 > 20. ) continue;
        //if( trueY > 0.95 || trueY < 0.1 ) continue;

        //do analysis, or fill historgrams for event levels
        
        counter_event++;

        //particle loop
        int counter                      = 0;
               
        int gt_counter                   = 0;
        
        int counter_pos_pions               = 0;
        int counter_neg_pions             = 0;
        int counter_p_2loop               = 0;
        int counter_p_gt_2loop            = 0;
       
               
                       
        double pt_cms                    = 0;
        double pt_lab                    = 0;
        double xf_2                      = 0;
        double xf_beagle                 = 0;
        double xf_cms_s                  = 0;
        double p_final_pz                = 0;  

        trys = 0;
        trys1= 0;    
        int counter_p_gt = 0;
        int counter_p = 0;
        int gt_protons_forw=0, gt_protons_zero=0,gt_protons_back=0;
        int negative_nTracks_target= 0, negative_nTracks_central = 0, negative_nTracks_projectile = 0;
        int positive_nTracks_target=0,positive_nTracks_central=0,positive_nTracks_projectile=0;
        int charged_nTracks_target=0,charged_nTracks_central=0,charged_nTracks_projectile=0;
        
      
        for(int j(0); j < nParticles; ++j ) {


            const erhic::ParticleMC* particle = event->GetTrack(j);

            int pdg         = particle->GetPdgCode();
            int status      = particle->GetStatus();
            int index       = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
            double pt       = particle->GetPt();
            double eta      = particle->GetEta();
            double phi      = particle->GetPhi();
            double rap      = particle->GetRapidity();
            double mass     = particle->GetM();
            double theta    = particle->GetTheta(); 
            theta = theta*1000.0; //change to mrad;
            double p        = particle->GetP();
            double Energy   = particle->GetE();
            double xf       = particle->GetXFeynman();
            double px       = particle->GetPx();
            double py       = particle->GetPy();
            double pz       = particle->GetPz();
            int parent_index = particle->GetParentIndex();
            double parent_id =particle->GetParentId();
            //double event_particle = particle->GetEvent();

            p_beam.SetPxPyPzE(0.,0.,0.,sqrt(MASS_NUCLEON*MASS_NUCLEON)); //fixed target
            e_beam.SetPxPyPzE(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON)); //initial electron
            S = (e_beam+p_beam)*(e_beam+p_beam);
            
                 
            if (pdg==11 && status==1 &&trys==0 ){

                e_scattered.SetPxPyPzE(px,py,pz,Energy);
                p_gamma = e_beam-e_scattered;
                W2 = (p_beam+p_gamma)*(p_beam+p_gamma);
                trys++;
                rap_e_lab = rap;
                h_rap_e_lab->Fill(rap_e_lab);
                counter_final_e++; 

            }
            
            if (pdg == 22 && status==21 && trys1==0){

                p_gamma_beagle.SetPxPyPzE(px,py,pz,Energy); //virtual photon
                trys1++;
                W2_BeAGLE_computed = (p_beam+p_gamma_beagle)*(p_beam+p_gamma_beagle);
                Beta_BeAGLE = sqrt((p_gamma_beagle.E()*p_gamma_beagle.E())+trueQ2)/(p_gamma_beagle.E()+MASS_NUCLEON);
                p_gamma_counter++;
                        
            }

            if (pdg==2112 && status==21){struck_neutron++;}

            if (pdg==2212 && status==21){struck_proton++;}          

            //particle cuts
            if (status != 1) continue; //final state particles 

            //Initializing Particle 

            p_final.SetPxPyPzE(px,py,pz,Energy);
            pt_lab = sqrt(p_final.Px()*p_final.Px()+p_final.Py()*p_final.Py());            
            Beta = sqrt((p_gamma.E()*p_gamma.E())+trueQ2)/(p_gamma.E()+MASS_NUCLEON);
            gamma = 1/sqrt(1-(Beta*Beta));            
            
            //LAB FRAME 
           
            if (pdg==2212)
            {
                xf_2 = 2*p_final.Pz()/sqrt(S);
                rap_lab = rap;
                rap_lab_computed = 0.5*log((p_final.E()+p_final.Z())/(p_final.E()-p_final.Z()));
                p_final_pz = (p_final.Pz());
                h_pz_proton->Fill(p_final_pz);
                xf_beagle = xf;                
                protons_lab++;
                h_xf->Fill(xf_beagle);
                h_xf_2->Fill(xf_2);
                h_rap_lab->Fill(rap_lab);
                h_rap_lab_computed->Fill(rap_lab_computed);

                
            }

            if (pdg==2112) {neutrons_lab++; }                      


            h_Energy->Fill(p_final.E());
            h_pt->Fill(pt_lab);
            h_pz->Fill(p_final.Pz());
            h_mom->Fill(p_final.P());

            rot_x[0][0] = cos(e_scattered.Phi());
            rot_x[0][1] = -sin(e_scattered.Phi());
            rot_x[1][0] = sin(e_scattered.Phi());
            rot_x[1][1] = cos(e_scattered.Phi());
            rot_x[2][2] = 1;
            //rot_x[3][3] = 1;                        

            p_final_rotated = rot_x*(p_final.Vect());

            rot_z[0][0] = cos(p_gamma.Theta());
            rot_z[0][2] = -sin(p_gamma.Theta());
            rot_z[2][0] = sin(p_gamma.Theta());
            rot_z[2][2] = cos(p_gamma.Theta());
            rot_z[1][1] = 1;
           

            p_final_rotated = rot_z*(p_final_rotated);

            p_final.SetPxPyPzE(p_final_rotated.X(),p_final_rotated.Y(),p_final_rotated.Z(),Energy);


            TVectorD p_final_boost(4);
            TVectorD p_final_boost_after(4);

            p_final_boost(0) = p_final.Px();
            p_final_boost(1) = p_final.Py();
            p_final_boost(2) = p_final.Pz();
            p_final_boost(3) = p_final.E();           


            boost_z[0][0] = 1;
            boost_z[1][1] = 1;
            boost_z[2][2] = gamma;
            boost_z[2][3] = -gamma*Beta_BeAGLE;
            boost_z[3][2] = -gamma*Beta_BeAGLE;
            boost_z[3][3] = gamma;


            p_final_boost_after = boost_z*p_final_boost;
            p_final.SetPxPyPzE(p_final_boost_after(0),p_final_boost_after(1),p_final_boost_after(2),p_final_boost_after(3));
          
            //----------- CMS FRAME -----------------------------------------

            pt_cms = sqrt(p_final.Px()*p_final.Px()+p_final.Py()*p_final.Py());//***            
            h_p_cms->Fill(p_final.P());
            h_Energy_cms->Fill(p_final.E());
            h_pz_cms->Fill(p_final.Pz());
            h_pt_cms->Fill(pt_cms);
                        

            //Rapidity in the new frame
            
            rap_cms = 0.5*log((p_final.E()+p_final.Z())/(p_final.E()-p_final.Z()));//signos cambiados en z

            h_rap_cms->Fill(rap_cms);

            if (rap_cms < -1){gt_protons_back++;}                         
            
            if (rap_cms > -0.5 && rap_cms < 0.5 ){gt_protons_zero++;}
                            
            if (rap_cms >= 2){gt_protons_forw++;}

            //############ Rapidity cuts by charged, positive and negative particles###############

           if (pdg == 11 || pdg == -211)
           {

                if (rap_cms < -1){negative_nTracks_target++;}                       
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){negative_nTracks_central++;}
                            
                if (rap_cms >= 2){negative_nTracks_projectile++;}

           }
           
           if (pdg == 2212 || pdg == 211)
           {
                if (rap_cms < -1){positive_nTracks_target++;}                         
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){positive_nTracks_central++;}
                            
                if (rap_cms >= 2){positive_nTracks_projectile++;}
           }

           if (pdg == 2212 || pdg == 211 || pdg ==11 || pdg ==-211)
           {
                if (rap_cms < -1){charged_nTracks_target++;}                         
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){charged_nTracks_central++;}
                            
                if (rap_cms >= 2){charged_nTracks_projectile++;}
           }
           //###################################################################################

            if (pdg==2212)
            {
                counter_p++;
                counter_p_2loop++;
                xf_cms = 2*p_final.Pz()/sqrt(trueW2);
                xf_cms_s = 2*p_final.Pz()/sqrt(S);               
                //cout<<"A :" <<pow(2*p_final.Pz()/xf,2)<<endl;
                //rap_cms = 0.5*log((p_final.E()+p_final.Z())/(p_final.E()-p_final.Z()));

                h_nTracks_proton->Fill(counter_p_2loop);


                if(xf_cms > -0.2)continue;

                cms_protons_xf_cut++;


                    if (p>0.2 && p < 0.8 ){

                        counter_p_gt++;
                        counter_p_gt_2loop++;
                        //cout <<"counter_p_gt: "<< counter_p_gt_2loop << endl;
                        
                        h_Energy_proton_gt->Fill(Energy);
                        h_pt_proton_gt->Fill(pt);
                        h_rapidity_gt->Fill(rap);
                        h_nTracks_p_gt->Fill(counter_p_gt_2loop); 

                    }
                    
    
            h_Energy_proton->Fill(Energy);
            h_pt_proton->Fill(pt);
            h_pz_proton_cms->Fill(p_final.Pz());         
            h_xf_cms->Fill(xf_cms);
            h_xf_cms_s->Fill(xf_cms_s);
            //h_rap_cms->Fill(rap_cms);
            }

            if (pdg == 211){ counter_pos_pions++;} 
            if (pdg == -211) {counter_neg_pions++;}

            
            if (p > 0.2 && p < 0.8) //Grey Tracks All Particles
            {
                gt_counter++;
                h_pt_gt->Fill(pt);
                h_Energy_gt->Fill(Energy);
                        
            }

            //############### Particles without the Pt Cut ###########

            

            S_boosted = (e_beam+p_beam)*(e_beam+p_beam);

               
            if (pdg==11 && status==1  ){

            W2_boosted = (p_beam+p_gamma)*(p_beam+p_gamma);
           
           }
   

            //Fill histograms particle loop
                    
            counter++;
            counter_total_particles++;
            
            h_theta->Fill(theta);   
            h_pdg->Fill(pdg);
            h_Status->Fill(status);         
            h_eta->Fill(eta);
            
            h_rapidity->Fill(rap);
            h_mass->Fill(mass);         
            h_rap_eta->Fill(rap,eta,1);
            h_mass_mom->Fill(p,mass,1);
            h_mass_pt->Fill(pt,mass,1);
            h_rap_gap->Fill(counter_p_gt,rap);   

                
                        

        } // end of particle loop

        // if(counter_p_gt==0)
        // {
        //     sum_N_prot_gt_0 = counter_p+sum_N_prot_gt_0;

        //     // average_N_prot=sum_N_prot/counter_p;
        //     cout << "counter_p_gt: " << counter_p_gt << endl;
        //     cout << "counter_p: " << counter_p << endl;
        //     cout << "sum_N_prot : "<<sum_N_prot_gt_0 << endl;
        // }

        // cout << "sum_N_prot OUTSIDE LOOP IF: "<<sum_N_prot_gt_0 << endl;

        // cout << "average: " << average_N_prot << endl;
            
        
        //fill histograms
        h_nTracks->Fill(nParticles);
        h_b->Fill(impact_parameter);
        h_trueQ2->Fill(trueQ2);
        h_nTracks_cuts->Fill(counter);
       
        
        h_nTracks_gt->Fill(gt_counter);
        
        
        h_nTracks_pos_pions->Fill(counter_pos_pions);

        h_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);
        prof_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);

        prof_nTracks_target_positive->Fill(counter_p_gt,positive_nTracks_target);
        prof_nTracks_central_positive->Fill(counter_p_gt,positive_nTracks_central);
        prof_nTracks_projectile_positive->Fill(counter_p_gt,positive_nTracks_projectile);
        prof_nTracks_target_negative->Fill(counter_p_gt,negative_nTracks_target);
        prof_nTracks_central_negative->Fill(counter_p_gt,negative_nTracks_central);
        prof_nTracks_projectile_negative->Fill(counter_p_gt,negative_nTracks_projectile);
        prof_nTracks_target_charged->Fill(counter_p_gt,charged_nTracks_target);
        prof_nTracks_central_charged->Fill(counter_p_gt,charged_nTracks_central);
        prof_nTracks_projectile_charged->Fill(counter_p_gt,charged_nTracks_projectile);

        for (int p=0; p<graph_n_vs_ng->GetN(); ++p) {
        graph_n_vs_ng->GetY()[p] = graph_n_vs_ng->GetY()[p] * 1.0/1.17; // aquí le pides al eje y que te de entrada por entrada
        graph_n_vs_ng->GetEY()[p] = graph_n_vs_ng->GetY()[p] * 1.0/1.17; // tambien tienes que escalar el error
        }
        
             
        h_w2_computed->Fill(W2);
        h_W2_BeAGLE_computed->Fill(W2_BeAGLE_computed);        
        h_w2->Fill(trueW2);
        h_nu->Fill(trueNu);
        h_bx->Fill(trueX);

    //cout << "number of total protons CMS frame & final state: " << counter_p << endl;
        //cout << " ends of event loop: " << endl;
    }// end event loop

 

    cout << "Number of particles: " << counter_total_particles << endl;
    cout << "Events after event cuts : " << counter_event << endl;
    cout << "mean value of All protons: " <<h_nTracks_proton->GetMean() << endl;
    cout << "number scattered electrons (1 per event): " << counter_final_e << endl;
    cout << "number virtual photons (1 per event): " << p_gamma_counter << endl;
    
    cout << "number of total neutrons lab frame & final state: " << neutrons_lab << endl;
    cout << "number struck nucleons(1 per event) : " << struck_proton+struck_neutron << endl;
    cout << "number struck neutron : " << struck_neutron <<" percentage: " << (100*struck_neutron)/(struck_proton+struck_neutron)<< endl;
    cout << "number struck proton : " << struck_proton << " percentage: "<<(100*struck_proton)/(struck_proton+struck_neutron) << endl;
    
    //cout << "number of total protons CMS frame & final state: " << counter_p << endl;
    cout << "number of total protons lab frame & final state: " << protons_lab << endl;
    cout << "number protons xf cut: " << cms_protons_xf_cut <<endl;
    // cout << "number of protons rap < -1: " << gt_protons_back << endl;
    // cout << "number of protons -0.5 < rap < 0.5 : " << gt_protons_zero << endl;
    // cout << "number of protons rap > 2: " << gt_protons_forw << endl;


    cout << "W2 Before:  " << W2 << endl;
    cout << "W2 after:  " << W2_boosted << endl;

       // cout << "counter_p_gt OUTSIDE EVENT LOOP: " << counter_p_gt << endl;
       //  cout << "counter_p: " << counter_p << endl;

    
    
    //################################### Drawing for e - Xe #############################
    TCanvas *p4 = new TCanvas("p4","p4");
    p4->Divide(3,1);
    
    p4->cd(1);
    prof_nTracks_target_charged->SetTitle("Target negative;Target charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_charged->Scale(1.0/1.17);
    prof_nTracks_target_charged->Draw();
    p4->cd(2);
    prof_nTracks_central_charged->SetTitle("Central negative;Central Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_charged->Scale(1.0/1.17);
    prof_nTracks_central_charged->Draw();
    p4->cd(3);
    prof_nTracks_projectile_charged->SetTitle("Projectile negative;Projectile Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_charged->Scale(1.0/1.17);
    prof_nTracks_projectile_charged->Draw();

    TCanvas *p3 = new TCanvas("p3","p3");
    p3->Divide(3,1);
    
    p3->cd(1);
    prof_nTracks_target_negative->SetTitle("Target negative;Target - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_negative->Scale(1.0/1.17);
    prof_nTracks_target_negative->Draw();
    p3->cd(2);
    prof_nTracks_central_negative->SetTitle("Central negative;Central Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_negative->Scale(1.0/1.17);
    prof_nTracks_central_negative->Draw();
    p3->cd(3);
    prof_nTracks_projectile_negative->SetTitle("Projectile negative;Projectile Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_negative->Scale(1.0/1.17);
    prof_nTracks_projectile_negative->Draw();

    TCanvas *p2 = new TCanvas("p2","p2");
    p2->Divide(3,1);
    
    p2->cd(1);
    prof_nTracks_target_positive->SetTitle("Target negative;Target Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}"); 
    prof_nTracks_target_positive->Scale(1.0/1.17);
    prof_nTracks_target_positive->Draw();
    p2->cd(2);
    prof_nTracks_central_positive->SetTitle("Central positive;Central Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_positive->Scale(1.0/1.17);
    prof_nTracks_central_positive->Draw();
    p2->cd(3);
    prof_nTracks_projectile_positive->SetTitle("Projectile positive;Projectile + Region n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_positive->Scale(1.0/1.17);
    prof_nTracks_projectile_positive->Draw();
    


    TCanvas *tg1 = new TCanvas("tg1","tg1");
    tg1->cd();
    graph_n_vs_ng->SetTitle((title_collision+"Tracks; Grey Tracks n_{g}; All Tracks n").c_str());
    graph_n_vs_ng->Draw();


    
    TCanvas *p1 = new TCanvas("p1","p1");
    p1->cd();
    prof_nTracks_Grey_tracks->SetTitle((title_collision+"Tracks;Protons n_{g}; Protons <n>").c_str());
    prof_nTracks_Grey_tracks->Draw();



    TCanvas *s1 = new TCanvas("s1","tracks");
    s1->cd();
    //gStyle->SetOptStat(111111);
    h_nTracks_Grey_tracks->SetTitle((title_collision+"Tracks; Grey Tracks n_{g}; All Tracks n ").c_str());
    h_nTracks_Grey_tracks->Draw();


    TCanvas *t9 = new TCanvas("t9","t9");
    t9->cd();
    h_trueQ2->SetTitle((title_collision+ " Q^2 All  ; Q^2 ; nEvents").c_str());
    h_trueQ2->Draw();
    
    TCanvas *t25 =new TCanvas("t25","t25");
    t25->cd();
    h_w2_computed->SetTitle((title_collision+ " W^{2} All; w^{2}; nEvents").c_str());
    h_w2_computed->SetLineColor(kBlue);
    h_W2_BeAGLE_computed->SetLineColor(kRed);
    
    h_w2->Draw();
    
    h_w2_computed->Draw("same");
    h_W2_BeAGLE_computed->Draw("same");
    TLegend *lw2 = new TLegend(0.3, 0.75, 0.4, 0.9);
    lw2->SetBorderSize(0); // no border
    lw2->SetFillStyle(0);
    lw2->SetFillColor(0); // Legend background should be white
    lw2->SetTextFont(32);
    lw2->SetTextSize(0.04); // Increase entry font size!
    //lw2->SetHeader("Protons");
    lw2->AddEntry(h_w2, "W^{2} BeAGLE", "l");
    lw2->AddEntry(h_W2_BeAGLE_computed,"W^{2} BeAGLE Computed","l");
    lw2->AddEntry(h_w2_computed, "W^{2} Computed", "l");

    lw2->Draw();


    TCanvas *f7 = new TCanvas("f7","Y lab electron");
    f7->cd();
    h_rap_e_lab->SetTitle((title_collision+ "Y electrons  ; Y lab electrons  ; nParticles").c_str());
    h_rap_e_lab->SetLineColor(kBlue);
    //h_rap_e_lab->SetLineColor(kRed);
    h_rap_e_lab->Draw();
    //h_rap_e_lab->Draw("same");
    

    TCanvas *f6 = new TCanvas("f6","Xf cms");
    f6->cd();
    h_xf_cms_s->SetTitle((title_collision+ "Xf  ; Xf CMS  ; nParticles").c_str());
    h_xf_cms_s->SetLineColor(kBlue);
    h_xf_cms->SetLineColor(kRed);
    h_xf_cms_s->Draw();
    h_xf_cms->Draw("same");
    
    
    TLegend *l56 = new TLegend(0.2, 0.75, 0.25, 0.9);
    l56->SetBorderSize(0); // no border
    l56->SetFillStyle(0);
    l56->SetFillColor(0); // Legend background should be white
    l56->SetTextFont(32);
    l56->SetTextSize(0.04); // Increase entry font size!
    l56->AddEntry(h_xf_cms, "Xf sqrt{W^{2}} CMS Frame", "l");
    l56->AddEntry(h_xf_cms_s, "Xf sqrt{S} CMS Frame", "l");
    l56->Draw();


    TCanvas *f5 = new TCanvas("f5","Rap cms");
    f5->cd();
    h_rap_cms->SetTitle((title_collision+ "Rapidity  ; Rapidity CMS  ; nParticles").c_str());
    h_rap_cms->SetLineColor(kRed);
    h_rap_cms->Draw();

    TCanvas *f4 = new TCanvas("f4","Pz protons");
    f4->cd();
    h_pz_proton->SetTitle((title_collision+ "Pz  ; Pz  ; nParticles").c_str());
    h_pz_proton_cms->SetLineColor(kRed);
    h_pz_proton->Draw();
    h_pz_proton_cms->Draw("same");
    TLegend *l55 = new TLegend(0.2, 0.75, 0.25, 0.9);
    l55->SetBorderSize(0); // no border
    l55->SetFillStyle(0);
    l55->SetFillColor(0); // Legend background should be white
    l55->SetTextFont(32);
    l55->SetTextSize(0.04); // Increase entry font size!
    l55->SetHeader("Protons");
    l55->AddEntry(h_pz_proton, "P_{z} Lab Frame", "l");
    l55->AddEntry(h_pz_proton_cms, "P_{z} Momemtum Frame", "l");
    l55->Draw();


    TCanvas *f3 = new TCanvas("f3","f3");
    f3->cd();
    h_pt->SetTitle((title_collision+ "Pt  ; Pt ; nParticles").c_str());
    h_pt_cms->SetLineColor(kRed);   
    h_pt->Draw();
    h_pt_cms->Draw("same");
    TLegend *leg = new TLegend(0.35, 0.75, 0.4, 0.9);
    leg->SetBorderSize(0); // no border
    leg->SetFillStyle(0);
    leg->SetFillColor(0); // Legend background should be white
    leg->SetTextFont(32);
    leg->SetTextSize(0.04); // Increase entry font size!
    leg->AddEntry(h_pt, "P_{t} Lab Frame", "l");
    leg->AddEntry(h_pt_cms, "P_{t} Momemtum Frame", "l");
    leg->Draw();
    
    
    TCanvas *f2 = new TCanvas("f2","Pz");
    f2->cd();
    h_pz->SetTitle((title_collision+ "Pz  ; Pz  ; nParticles").c_str());
    h_pz_cms->SetLineColor(kRed);
    h_pz->Draw();
    h_pz_cms->Draw("same");
    TLegend *l4 = new TLegend(0.2, 0.75, 0.25, 0.9);
    l4->SetBorderSize(0); // no border
    l4->SetFillStyle(0);
    l4->SetFillColor(0); // Legend background should be white
    l4->SetTextFont(32);
    l4->SetTextSize(0.04); // Increase entry font size!
    l4->AddEntry(h_pz, "P_{z} Lab Frame", "l");
    l4->AddEntry(h_pz_cms, "P_{z} CMS Frame", "l");
    l4->Draw();

    TCanvas *f1 = new TCanvas("f1","f1");
    f1->cd();
    h_p_cms->SetLineColor(kRed);
    h_mom->SetTitle((title_collision+ "P  ; P  ; nParticles").c_str());
    h_mom->Draw();
    h_p_cms->Draw("same");
    TLegend *l = new TLegend(0.35, 0.75, 0.4, 0.92);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);
    l->SetFillColor(0); // Legend background should be white
    l->SetTextFont(32);
    l->SetTextSize(0.04); // Increase entry font size!
    l->AddEntry(h_mom, "P Lab Frame", "l");
    l->AddEntry(h_p_cms, "P Momemtum Frame", "l");
    l->Draw();


    

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
    
 
    TCanvas *t22 = new TCanvas("t22","t22");
    t22->cd();
    h_rapidity_gt->SetTitle((title_collision+ " Rapidity Grey Tracks Protons; Grey Tracks; nParticles").c_str());
    h_rapidity_gt->Draw();
    

    // c7->cd(2);
    // gPad->SetLogy();
    // h_nTracks_pos_pions->SetTitle((title_collision+ " NTracks positive Pions ; nTracks; nEvents").c_str());
    // h_nTracks_pos_pions->SetMarkerStyle(kStar);
    // h_nTracks_pos_pions->SetMarkerColor(kRed);
    // h_nTracks_pos_pions->Draw("P");   
   


    TCanvas *t14 = new TCanvas("t14","t14");
    t14->cd();
    h_Energy_proton->SetTitle((title_collision+ " Energy All Protons ; E(GeV) ; nParticles").c_str());
    h_Energy_proton->Draw();
                    

    TCanvas *t28 = new TCanvas("t28","t28");
    t28->cd();
    h_Energy_proton_gt->SetTitle((title_collision+ " Energy Grey Tracks Protons  ; E(GeV) ; nParticles").c_str());
    h_Energy_proton_gt->Draw();
    
   


    //TCanvas *t21 = new TCanvas("t21","t21");
    //t21->cd();
    //h_nTracks_gt->SetTitle((title_collision+ " All stable Grey Tracks;Grey Tracks; nEvents").c_str());
    //h_nTracks_gt->Draw();

    
    
    TCanvas *t16 = new TCanvas("t16","Pt protons");
    t16->cd();
    h_pt_proton->SetTitle((title_collision+ " Pt Protons ; Pt ; nParticles").c_str());
    h_pt_proton->Draw();
    
    TCanvas *t29 = new TCanvas("t29","Pt protons Gt");
    t29->cd();
    h_pt_proton_gt->SetTitle((title_collision+ " Pt Grey Tracks Protons ; Pt ; nParticles").c_str());
    h_pt_proton_gt->Draw();


    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    h_rap_gap->SetTitle((title_collision+ " Grey Tracks Protons vs Rapidity ; y ; Grey Tracks").c_str());
    h_rap_gap->SetMarkerStyle(kFullCircle);
    h_rap_gap->Draw("P");
        

    

    TCanvas *c5 = new TCanvas("c5","c5");
    c5->cd();
    //gPad->SetLogy();
    //gStyle->SetOptStat(111111);
    h_xf->SetTitle((title_collision+ " Feynman - X protons; X_{f};nParticles").c_str());    
    //h_xf_cms->SetTitle((title_collision+ " Feynman - X CMS; Xf;nParticles").c_str());
    h_xf_cms->SetLineColor(kRed);
    //h_xf_2->SetLineColor(kBlue);
    //cout<<h_xf_2->GetEntries()<<endl;
    cout<<h_xf->GetEntries()<<endl;
    cout<<h_xf_cms->GetEntries()<<endl;
    //h_xf_2->Draw(); 
    h_xf->Draw();
    h_xf_cms->Draw("same");    
    
     
    

    
    TLegend *l5 = new TLegend(0.35, 0.75, 0.45, 0.92);
    l5->SetBorderSize(0); // no border
    l5->SetFillStyle(0);
    l5->SetFillColor(0); // Legend background should be white
    l5->SetTextFont(32);
    l5->SetTextSize(0.035);
    l5->SetHeader("Protons"); // Increase entry font size!
    l5->AddEntry(h_xf, "X_{f} eicsmear ", "l");
    //l5->AddEntry(h_xf_2,"X_{f} Computed Lab Frame", "l");
    l5->AddEntry(h_xf_cms, "X_{f} CMS Frame", "l");
    l5->Draw(); 


    
    TCanvas *c6 = new TCanvas("c6","c6");
    c6->cd();
    h_theta->SetTitle((title_collision+ " THETA; THETA ;nParticles").c_str());
    h_theta->Draw();

    
    //-------------------------------------- All stable Particles ----------------------------------------------

    TCanvas *t5 = new TCanvas("t5","t5");
    t5->cd();
    h_rap_lab_computed->SetTitle((title_collision+ " Final stable Protons ; Rapidity ;nParticles").c_str());   
    h_rap_cms->SetLineColor(kRed);
    h_rap_lab_computed->SetLineColor(kBlue);
    h_rap_lab_computed->Draw();
    h_rap_cms->Draw("same");
    h_rap_lab->Draw("same");
    
    TLegend *l1 = new TLegend(0.56, 0.75, 0.6, 0.92);

    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);
    l1->SetFillColor(0); // Legend background should be white
    l1->SetTextFont(32);
    l1->SetTextSize(0.035);
    l1->SetHeader("Protons"); // Increase entry font size!
    l1->AddEntry(h_rap_lab, "Rapidity Lab Frame BeAGLE", "l");
    l1->AddEntry(h_rap_lab_computed, "Rapidity Lab Frame computed", "l");
    l1->AddEntry(h_rap_cms, "Rapidity CMS Frame", "l");
    l1->Draw(); 

    

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

    TCanvas *t1 = new TCanvas("t1","t1");
    t1->cd();
    h_Energy_cms->SetLineColor(kRed);
    h_Energy->SetTitle((title_collision+ " Energy All  particles;E(GeV) ;nParticles").c_str());    
    h_Energy->Draw();
    h_Energy_cms->Draw("same");
    TLegend *l2 = new TLegend(0.35, 0.75, 0.4, 0.92);
    l2->SetBorderSize(0); // no border
    l2->SetFillStyle(0);
    l2->SetFillColor(0); // Legend background should be white
    l2->SetTextFont(32);
    l2->SetTextSize(0.04); // Increase entry font size!
    l2->AddEntry(h_Energy, "Energy Lab Frame", "l");
    l2->AddEntry(h_Energy_cms, "Energy Momemtum Frame", "l");
    l2->Draw();
    

    

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

    

    TCanvas *t27 = new TCanvas("t27","t27");
    t27->cd();
    h_bx->SetTitle((title_collision+ " Bj x; x; nEvents").c_str());
    h_bx->Draw();
    
    
  

    
    //out->Write();

}
