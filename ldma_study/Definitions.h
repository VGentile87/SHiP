//*******************************
// Revision 1: 25/10/2015 (Valerio)
// - revison of recoil classification:
//   S included both in heavy and light nuclei
//
//*******************************

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <cmath>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"

using namespace std;

//// TEST DI CONTROLLO //////////////////////////////////////////////////////////
/*
const Int_t n_el = 6;
TString el[n_el]={"O","Pb","Xe","Ge","Ar","Si"};
Double_t M_riv = 1; //ton
//Double_t M_el[n_el]={14.8966,192.924,122.247,67.6257,37.1955,26.1499}; // uma
Double_t A_el[n_el]={15.999,207.2,131.293,72.63,39.948,28.085}; // uma
Double_t Z_el[n_el]={8,82,54,32,18,14}; // uma
Double_t atomic_fract_el[n_el]={1,1,1,1,1,1}; // uma
Double_t mass_fract[n_el]={1,1,1,1,1,1}; // uma
TString flavour[3]={"#nu_{e}","#bar{#nu}_{e}","#nu_{x}"};
//Double_t N_el[n_el]={7.999,125.2,77.293,40.63,21.948,14.085}; // uma
Double_t M_el[n_el]={};
Double_t N_el[n_el]={};
*/
//////////////////////////////////////////////////////////////////////////////


////// TEST NEWS   ////////////////////////////////////////////////////////////

const Int_t n_el = 7;
Double_t M_riv = 1; //ton
// Int_t index=0;
TString el[n_el]={"Ag","Br","I","C","O","N","S"};
TString flavour[3]={"#nu_{e}","#bar{#nu}_{e}","#nu_{x}"};
Int_t Z_el[n_el]={47,35,53,6,8,7,16};
Double_t A_el[n_el]={107.8682,79.904,126.90447,12.011,15.999,14.00674,32.066};
Double_t atomic_fract_el[n_el]={0.10,0.10,0.004,0.214,0.118,0.049,0.003}; // uma
Double_t mass_fract[n_el]={0.44,0.32,0.019,0.101,0.074,0.027,0.003}; // uma
//Double_t thr_el[n_el]={110,85,120,17,22,20,44}; //keV   cut 50 nm
//Double_t thr_el[n_el]={55,45,60,10,11.5,13,25}; //keV     cut 30 nm
Double_t thr_el[n_el]={250,180,275,35,45,41,90}; //keV   cut 100 nm   
Double_t M_el[n_el]={};
Double_t N_el[n_el]={};

//////////////////////////////////////////////////////////////////////////////

Double_t E_nu_max=100;   //100
Double_t num_ev=0;
Double_t num_ev_el=0;
Double_t num_ev_miss=0;
Double_t err_num_ev=0;
Double_t err_num_ev_el=0;
Double_t pi = TMath::Pi();
Double_t GF = 1.16637*TMath::Power(10,-11); // MeV^-2
Double_t sin2wnb = 0.239;
Double_t s=0.9; // fm 
Double_t res = 5; // degree    
Double_t uma = 0.9310986; //GeV
//Double_t dist_sn = 10; //kpc
Double_t int_const=0.;
Float_t L_true=0.;
Double_t L_proj=0.;

Double_t num_ev_bkg=0;
Double_t flux_boron=5.58*TMath::Power(10,6);
Double_t int_const_bkg=0.;
Double_t dN_vs_dE_bkg=0.;
Double_t flux_boron_vs_E[17]={0,74438,228550,412460,554100,661470,744380,701700,661470,554100,412460,272830,151180,66147,11253,587,15};
Double_t iflux_boron_vs_E=0.;

Int_t d_theta = 900000;//900;//360/5;

Double_t fm_in_GeV=1./0.1975; 
Double_t GeV_in_MeV=1000;
Double_t iMeV_in_fm=197.5;
Double_t kpc_in_cm =3.086*TMath::Power(10,21);
Double_t kpc_in_fm =3.086*TMath::Power(10,34);
Double_t uma_in_ton = 1.66054*TMath::Power(10,-30);
Double_t erg_in_MeV = 624150.65; 

Double_t    N_nu_e = 2.8*TMath::Power(10,57);
Double_t    err_N_nu_e = 0.3*TMath::Power(10,57);
Double_t    meanE_nu_e = 11; // MeV
Double_t    err_meanE_nu_e = 1; // MeV
Double_t    beta_nu_e = (3./meanE_nu_e); // MeV^-1
Double_t    intE_nu_e=0.0;
/////
Double_t    tot_ene_e = 0.5*TMath::Power(10,53); // erg
Double_t    mean_ene_e = 9.5; //MeV
Double_t    alpha_e = 2.5;
   
Double_t    N_antinu_e = 2.1*TMath::Power(10,57);
Double_t    err_N_antinu_e = 0.4*TMath::Power(10,57);
Double_t    meanE_antinu_e = 15; // MeV
Double_t    err_meanE_antinu_e = 3; // MeV
Double_t    beta_antinu_e = (3./meanE_antinu_e); // MeV^-1
Double_t    intE_antinu_e=0.0;
/////
Double_t    tot_ene_antie = 0.5*TMath::Power(10,53); // erg
Double_t    mean_ene_antie = 12; //MeV
Double_t    alpha_antie = 2.5;
   
Double_t    N_nu_x = 1.5*TMath::Power(10,57);
Double_t    err_N_nu_x = 0.4*TMath::Power(10,57);
Double_t    meanE_nu_x = 21; // MeV
Double_t    err_meanE_nu_x = 6; // MeV
Double_t    beta_nu_x = (3./meanE_nu_x);  // MeV^-1
Double_t    intE_nu_x=0.0;
/////
Double_t    tot_ene_x = 2*TMath::Power(10,53); // erg
Double_t    mean_ene_x = 15.6; //MeV
Double_t    alpha_x = 2.5;

Double_t spect_nu_e=0.0;
Double_t spect_antinu_e=0.0;
Double_t spect_nu_x=0.0;

Double_t    E_nu=0.;
Double_t    E_rec=0.;
Double_t    q=0.;
Double_t    qrn=0.;
Double_t    qs=0.;
Double_t    Fq=0.;
Double_t    crsect=0.;
Double_t    cos_theta=0.;
Double_t    cos_theta_norm=0.;
Double_t    C=0.;
Double_t    dN_vs_dE=0.;
Double_t    sum_E=0.;
Double_t    sum_cos_theta=0.;
Double_t    sum_crsect=0.;
//Int_t    index=0;
TString    iel = " ";
Int_t    iZ_el = 0;
Double_t    iA_el = 0.;
Double_t    ifract_el = 0.; // uma
Double_t imass_fract=0.;
Double_t    iM_el = 0.;
Double_t    iN_el = 0.;
Double_t    ithr_el = 0.;
Double_t iQw;
Double_t irn;  // fm
Double_t cos_theta_rec=0.;
Double_t domega=0.;
Double_t theta_rec=0.;
Double_t theta_nu=0.;
Double_t phi_rec=0.;

TGraph *gr_sigma[n_el] = {new TGraph()};
TMultiGraph *mg = new TMultiGraph();
TCanvas *cmg;
TCanvas *cc7;

TGraph *gr_spect[3] = {new TGraph()};
TMultiGraph *mg2 = new TMultiGraph();
TCanvas *cmg2;

TH1D *h1[7];
TH1D *h11[7];
TH1D *h7 = new TH1D("L_true","L_true",100,0.1,1);
TProfile *h6 = new TProfile("h6","h6",1000,0,2000,0,35);
TH1D *h2 = new TH1D("the","the",100,0,pi);
TH1D *h22 = new TH1D("thL","thL",100,0,pi);
TH1D *h3 = new TH1D("th","th",50,-pi/2.,pi/2.);
//TH1D *h3 = new TH1D("th","th",50,0,pi);
TH1D *h33 = new TH1D("costhSN","costhSN",50,0,1);
TH1D *h4 = new TH1D("ene_rec","ene_rec",400,0,10000);
TH1D *h44 = new TH1D("ene_rec_thr_1kev","ene_rec_thr_1kev",400,0,10000);
TH1D *h5 = new TH1D("phiSN","phiSN",18,-pi,pi);
TH1D *h55 = new TH1D("phib","phib",100,0,2*pi);

TH1D *hh0 = new TH1D("L_proj Ag","L_proj Ag",50,0.05,1);
TH1D *hh1 = new TH1D("L_proj Br","L_proj Br",50,0.05,1);
TH1D *hh2 = new TH1D("L_proj I","L_proj I",50,0.05,1);
TH1D *hh3 = new TH1D("L_proj C","L_proj C",50,0.05,1);
TH1D *hh4 = new TH1D("L_proj O","L_proj O",50,0.05,1);
TH1D *hh5 = new TH1D("L_proj N","L_proj N",50,0.05,1);
TH1D *hh6 = new TH1D("L_proj S","L_proj S",50,0.05,1);

TH2D *henelenAg = new TH2D("enelen Ag","enelen Ag",100,0,1,100,0,1);
TH2D *henelenBr = new TH2D("enelen Br","enelen Br",100,0,1,100,0,1);
TH2D *henelenC = new TH2D("enelen C","enelen C",100,0,1,100,0,1);
TH2D *henelenO = new TH2D("enelen O","enelen O",100,0,1,100,0,1);
TH2D *henelenN = new TH2D("enelen N","enelen N",100,0,1,100,0,1);


TH1D *h777 = new TH1D("the_sn","the_sn",100,-1,1);
TH1D *h77 = new TH1D("theta_rec_range","theta_rec_range",100,0,pi);
TH1D *hLen = new TH1D("L_true_range","L_true_range",1000,0,0.05);


TH1D *h_Enu = new TH1D("E_nu","",100,0,100);
TH1D *h_Enu0 = new TH1D("E_nu1","",100,0,100);
TH1D *h_Enu1 = new TH1D("E_nu2","",100,0,100);
TH1D *h_Enu2 = new TH1D("E_nu3","",100,0,100);
TH1D *h_Enu3 = new TH1D("E_nu4","",100,0,100);
TH1D *h_Enu4 = new TH1D("E_nu5","",100,0,100);
TH1D *h_Enu5 = new TH1D("E_nu6","",100,0,100);
TH1D *h_Enu6 = new TH1D("E_nu7","",100,0,100);
TH1D *hphi_thr = new TH1D("phi_thr","phi_thr",18,-pi,pi);
TH1D *hthe_thr = new TH1D("the_thr","the_thr",18,-pi,pi);
TH1D *hcosthe_thr = new TH1D("costhe_thr","costhe_thr",100,-1,1);

//TH1D *hene_thr = new TH1D("ene_thr","ene_thr",200,0,10000);
//TH1D *hlen_thr = new TH1D("lenSN","lenSN",200,0.,5000000);
//TH2D *hlen_ene = new TH2D("len_ene","len_ene",200,0,10000,200,0.,5000000);
//TH1D *hene_thr = new TH1D("ene_nu_rec","",1000,0,10000);
//TH1D *hlen_thr = new TH1D("len_nu_rec","",10000,0.,100000);
//TH2D *hlen_ene = new TH2D("len_ene_nu_rec","",1000,0,10000,10000,0.,100000);


TH3D *hxyzrec = new TH3D("hxyzrec","",100,-1,1,100,-1,1,100,-1,1);
TH1D *hphirec = new TH1D("hphirec","",60,-pi,pi);
TH1D *htherec = new TH1D("htherec","",60,-pi/2.,pi/2.);
TH1D *hcostherec = new TH1D("hcostherec","",100,-1.,1);
Double_t xrec=0;
Double_t yrec=0;
Double_t zrec=0;
Double_t phirec=0;
Double_t therec=0;

Double_t E_nu_cm =0;
Double_t E_rec_cm =0; 
Double_t rnd_cos_theta_cm=0; 
Double_t rnd_phi_cm=0; 
Double_t px_cm=0; 
Double_t py_cm=0; 
Double_t pz_cm=0;
Double_t p_mod_cm=0;
Double_t beta_cm=0;
Double_t gamma_cm=0;
Double_t px_lab=0;
Double_t py_lab=0;
Double_t E_rec_lab=0;
Double_t E_tot_lab=0;
Double_t p_mod_lab=0;
Double_t the_rec=0;

TH1D *h00 = new TH1D("costheSN","costheSN",50,0,1);
TH1D *h0 = new TH1D("h0","h0",100,-pi,pi);
TH1D *hpx_cm = new TH1D("hpx_cm"," ",200,-100,100);
TH1D *hpy_cm = new TH1D("hpy_cm"," ",200,-100,100);
TH1D *hpz_cm = new TH1D("hpz_cm"," ",200,-100,100);
TH1D *hpx_lab = new TH1D("hpx_lab"," ",200,-100,100);
TH1D *hphi = new TH1D("phi_rec","phi_rec",60,-pi,pi);
TH1D *hthe = new TH1D("the_rec","the_rec",100,-pi,pi);
TH1F *htthe_recx  = new TH1F("htthe_recx","",60,-pi,pi);
TH1F *htthe_recy  = new TH1F("htthe_recy","",60,-pi,pi);

Double_t px_sum = 0;
Double_t py_sum = 0;
Double_t pz_sum = 0;
Double_t phi_sum = 0;
Double_t theta_sum = 0;
Double_t cos_theta_sum = 0;
Double_t costheta_sc=0;
Double_t theta_sc=0;
Double_t phi_sc=0;
Double_t pphi=0;
Double_t pthe=0;

////// WIMP DEFINITIONS /////////
Double_t MN=0;
Double_t Mp=0;
Double_t xsec=0;
//Double_t Md=0;
Double_t r=0;
Double_t R0=0;
Double_t densdm=0.3; //0.4
Double_t v0=220.0;
Double_t numav=6.02*TMath::Power(10,26);
Double_t exp1=0;
//Double_t vE=220;
Double_t E0=0;
Double_t c=299792.458; //Km/s
//Double_t vE       = c*0.995;//(300000*0.995)/c;  //km/s (lab velocity)
Double_t vmin=0;
//Double_t p_wimp=5; // GeV
Double_t vesc=533;//533;
Double_t num_wimp_ev=0;
Double_t crsect_wimp=9.2*TMath::Power(10,-45);   // 5 per Xenon
Double_t C1=0.751;
Double_t C2=0.561;

////////////////////////
Float_t vE=0;
Float_t velocity=0;
Float_t theta_x, theta_y, the_rec_x, the_rec_y;
Float_t phi0, theta0, phi, theta, angle0, angle3D;
Float_t vDM_x, vDM_y, vDM_z, vDM;
Float_t vCM_x, vCM_y, vCM_z, vCM;
Float_t vT_CM_x, vT_CM_y, vT_CM_z, vT_CM;
Float_t vT_x, vT_y, vT_z, vT;
Float_t energy, length;
Float_t Edm_cm, pz_lab;

TRandom3 *rn = new TRandom3();
TRandom3 *rn2 = new TRandom3();
TH1F *htT0CUT  = new TH1F("htT0CUT","",60,-pi,pi);
TH1F *htT0CUT3D  = new TH1F("htT0CUT3D","",60,-pi/2.,pi/2.);
TH1F *htthex  = new TH1F("htthex","",60,-pi,pi);
TH1F *htthey  = new TH1F("htthey","",60,-pi,pi);
//TH1F *he       = new TH1F("he","",100,0,10000);
//TH1F *hl       = new TH1F("hl","",100,0,10000); 
