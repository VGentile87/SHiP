//  AUTHOR: V. GENTILE  2017 //
//
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>
#include <time.h>

#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLine.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveStats.h"

#include "Definitions.h"
//#include "Functions.h"
#include "EnergyRangeCorrelation.C"
#include "TRandom3.h"
#include <TSystem.h>
#include "TProfile2D.h"


using namespace std;

void ship_ldma() {

   Double_t threshold;
  Double_t max_range;
  Double_t MD;
  Double_t p_wimp;
  int choice;
  //cout << "Inserisci massa DM [GeV]" << endl;
  cout << "Masse disponibili [MeV]" << endl;
  cout << "1) 3.3 MeV \t 2) 10 MeV \t 3) 16.7 MeV \t 4) 26.7 MeV" << endl;
  cout << "5) 50 MeV \t 6) 103.3 MeV \t 7) 170 MeV \t 8) 333 MeV" << endl;
  cout << "Scegli opzione" << endl;
  cin  >> choice;

  if(choice==1)MD=3.3; // GeV
  if(choice==2)MD=10.0; // GeV
  if(choice==3)MD=16.7; // GeV
  if(choice==4)MD=26.7; // GeV
  if(choice==5)MD=50.0; // GeV
  if(choice==6)MD=103.3; // GeV
  if(choice==7)MD=170.0; // GeV
  if(choice==8)MD=333.0; // GeV
  
  //cout << "Inserisci energia nu [MeV]" << endl;
  //cin  >> E_nu;
  cout << "Inserisci soglia minima [nm]" << endl;
  cin  >> threshold;
  cout << "Inserisci range max [nm]" << endl;
  cin  >> max_range;
  //cout << "Inserisci impulso particella [GeV]" << endl;
  //cin >> p_wimp;

  int nwimp=0;
  int nnu=0;

  float rate_wimp=0;
  float rate_nu=0;
  
  TH1F *hpdm0      = new TH1F("hpdm0","",50,0,200);
  TH1F *hpdm      = new TH1F("hpdm","",50,0,200);
  TFile *f1 = new TFile(Form("./DM_spectrum/histo_dm_%.1fMeV.root",MD),"READ");
  hpdm = (TH1F*)f1->Get("h_p_dm_lab_all");
  hpdm->Scale(1./hpdm->Integral());
  //f1->Close();

  if(choice==1)MD/=1000.; // GeV
  if(choice==2)MD/=1000.; // GeV
  if(choice==3)MD/=1000.; // GeV
  if(choice==4)MD/=1000.; // GeV
  if(choice==5)MD/=1000.; // GeV
  if(choice==6)MD/=1000.; // GeV
  if(choice==7)MD/=1000.; // GeV
  if(choice==8)MD/=1000.; // GeV
  

  TH1F *henu0      = new TH1F("henu0","",100,0,400);
  TH1F *henu1      = new TH1F("henu1","",400,0,400);
  TH1F *henu2      = new TH1F("henu2","",400,0,400);
  TFile *f2 = new TFile("spectra_interacting_single.root","READ");
  henu1 = (TH1F*)f2->Get("hspectrum_nu_mu_intve_nc");
  henu2 = (TH1F*)f2->Get("hspectrum_nu_mu_bar_intve_nc");
  //henu2->Scale(henu1->Integral()/henu2->Integral());
  henu1->Add(henu2);
  henu1->Rebin(4); //4MeV a bin
  henu1->Scale(1./henu1->Integral());
  //f2->Close();

  const int npoints=100000;
  double p_rnd_wimp[npoints]={};
  double e_rnd_nu[npoints]={};

  for(int i=0;i<npoints;i++){
    p_rnd_wimp[i] = hpdm->GetRandom();
    e_rnd_nu[i] = henu1->GetRandom();
  }

  TH1F *hv0      = new TH1F("hv0","",100,0,600);
  TH1F *hEr      = new TH1F("Er","",100,0,1000);

  TH1F* hlen_nu[n_el];
  TH1F* hlen_dm[n_el];
  

  //Float_t Edm = TMath::Sqrt(p_wimp*p_wimp + MD*MD); // GeV
  Float_t Edm =0;  

  TRandom3 *tt = new TRandom3();
    
  gStyle->SetOptStat(1111);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);                
  gStyle->SetStatW(0.3);                
  gStyle->SetStatH(0.15);                

  TRandom3 *r = new TRandom3();
  Double_t t = 0;
  TRandom3 *rr = new TRandom3();
  Double_t a = 0;


  ///// BIN LOGARITMICO
  const int nbins=100;
  Double_t xbins_len[nbins+1]={};
  double dx_len = 8./nbins;  // il num per gli ordini di grandezza
  Double_t xbins_ene[nbins+1]={};
  double dx_ene = 6./nbins;  // il num per gli ordini di grandezza
  double l10 = TMath::Log(10);
  for (int ibin=0;ibin<=nbins;ibin++) {
    xbins_len[ibin] = TMath::Exp(l10*ibin*dx_len-l10*1);
    xbins_ene[ibin] = TMath::Exp(l10*ibin*dx_ene/*-l10*3*/);
    //printf("i %d, xbins_ene %4.2f, xbins_len %4.2f\n ",ibin, xbins_ene[ibin], xbins_len[ibin]);
  }  // il num per l'origine (3 sta per 10e-3)
  TH1F *hlen_thr = new TH1F("len_thr","",nbins,xbins_len);
  TH1F *hene_thr = new TH1F("ene_thr","",nbins,xbins_ene);
  TH1F *he = new TH1F("ene_dm_rec","",nbins,xbins_ene);
  TH1F *hl = new TH1F("len_dm_rec","",nbins,xbins_len);
  TH2F *hle = new TH2F("len_ene_dm_rec","",nbins,xbins_ene,nbins,xbins_len);
  TProfile2D *hlen_ene[7];
  TProfile2D *hlen_ene_dm[7];
  for (int i=0;i<7;i++) {
    hlen_ene[i] = new TProfile2D(Form("len_ene_%d",i),"",nbins,xbins_ene,nbins,xbins_len);
    hlen_ene_dm[i] = new TProfile2D(Form("len_ene_dm_%d",i),"",nbins,xbins_ene,nbins,xbins_len);
  }


  //// STARTING POINT
  
  for(int i=0; i<n_el;i++){
    h1[i] = new TH1D(Form("cos_el_%d",i),"",100,-pi,pi);
    iel = el[i];
    iZ_el = Z_el[i];
    iA_el = A_el[i];
    ifract_el = atomic_fract_el[i]; // uma
    imass_fract = mass_fract[i];
    iM_el=iA_el*uma;
    iN_el = N_el[i];

    hlen_nu[i] = new TH1F(Form("nu_len_"+el[i]),"",nbins,xbins_len);
    hlen_dm[i] = new TH1F(Form("dm_len_"+el[i]),"",nbins,xbins_len);
        
    //////////////////// SPLINES //////////////////////////////////
    std::ifstream inputFile;
    if(iZ_el==6)inputFile.open("DATA/C.txt", std::ifstream::in);
    if(iZ_el==7)inputFile.open("DATA/N.txt", std::ifstream::in);
    if(iZ_el==8)inputFile.open("DATA/O.txt", std::ifstream::in);
    if(iZ_el==16)inputFile.open("DATA/S.txt", std::ifstream::in);
    if(iZ_el==35)inputFile.open("DATA/Br.txt", std::ifstream::in);
    if(iZ_el==47)inputFile.open("DATA/Ag.txt", std::ifstream::in);
    if(iZ_el==53)inputFile.open("DATA/I.txt", std::ifstream::in);
    
    int num_line=0;
    std::string line;
    
    Double_t x[156]; // vedi num line in data
    Double_t y[156];
    
    Int_t Narray = 0;
    
    while(getline(inputFile, line)) {
      num_line++;
    if (!line.length() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    double ene = 0.;
    char ene_unit[3] = {};
    double fake=0.;
    double fake2=0;
    double len=0.;
    char len_unit[3] = {};

    if(num_line<=156){
      iss>>ene>>ene_unit>>fake>>fake2>>len>>len_unit;
      //cout << ene_unit << endl;
      if(ene_unit[0]=='k')ene=ene; //keV
      if(ene_unit[0]=='M')ene*=1000;  //keV
      if(ene_unit[0]=='G')ene*=1000000;  //keV
      if(len_unit[0]=='A')len/=10; // nm
      if(len_unit[0]=='u')len*=1000; // nm
      if(len_unit[0]=='m')len*=1000000; // nm
      //std::cout<<"point: "<<ene<< " keV " << len<<' '<< " nm " <<std::endl;
      //printf("i %d -> ENE %4.4f LEN %4.4f\n",Narray,ene,len);
      
      y[Narray] = len;
      x[Narray] = ene;
      Narray++;
    }
  }

  gr = new TGraph(Narray,x,y);
 
  
    for(int k=0;k<npoints;k++){
      
      //E_nu = 0.01; // cross-check
      E_nu = e_rnd_nu[k];      
      cos_theta = rn2->Uniform(-1,1);
      E_rec = ((E_nu*E_nu*(1-cos_theta))/(E_nu*(1-cos_theta)+iM_el))*TMath::Power(10,3); // energia nucluear recoil (MeV)
      
      costheta_sc = ((E_nu + iM_el)/(E_nu))*TMath::Sqrt(E_rec/(2*iM_el*TMath::Power(10,3)+E_rec));  // angolo nuclear recoil rispetto neutrino incidente
      theta_sc = TMath::ACos(costheta_sc);
      L_true = gr->Eval(E_rec*1000); // nm

      henu0->Fill(E_nu);
      if(i==1 && E_rec>10000 && theta_sc==0)cout << "neutrino  " << " " << iZ_el << " " << E_rec*1000 << " " << L_true << " " << iM_el << " " << TMath::ACos(costheta_sc) << " " << costheta_sc << endl;    
      
      
      if(L_true>=threshold && L_true<=max_range){	
	h1[i]->Fill(TMath::ACos(costheta_sc));
	//h1[i]->Fill(theta_sc);
	hlen_thr->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]); //nm
	hene_thr->Fill(E_rec,atomic_fract_el[i]*A_el[i]*A_el[i]); //MeV
	hlen_ene[i]->Fill(E_rec,L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	///
	a = r->Uniform(-pi,pi);
	zrec = TMath::Cos(theta_sc);
	xrec = TMath::Sin(theta_sc)*TMath::Cos(a);
	yrec = TMath::Sin(theta_sc)*TMath::Sin(a);
	therec = TMath::ATan(TMath::Sqrt(TMath::Power(xrec,2)+TMath::Power(yrec,2))/zrec);
	phirec = TMath::ATan(yrec/xrec);
	the_rec_x = TMath::ATan(xrec/zrec);
	the_rec_y = TMath::ATan(yrec/zrec);
	htthe_recx->Fill(the_rec_x,atomic_fract_el[i]*A_el[i]*A_el[i]);
	htthe_recy->Fill(the_rec_y,atomic_fract_el[i]*A_el[i]*A_el[i]);
	hphirec->Fill(phirec,atomic_fract_el[i]*A_el[i]*A_el[i]);
	htherec->Fill(therec,atomic_fract_el[i]*A_el[i]*A_el[i]);

	hlen_nu[i]->Fill(L_true/1000.,atomic_fract_el[i]*A_el[i]*A_el[i]);
	///
	nnu++;
      }
    }

    ////////////////// WIMP ////////////////
    for(int j=0;j<npoints;j++){
      //////////////////////////////////////
      p_wimp = p_rnd_wimp[j]; //hpdm->GetRandom();
      hpdm0->Fill(p_wimp);
      //p_wimp=10; // GeV
      Edm = TMath::Sqrt(p_wimp*p_wimp + MD*MD); // GeV
      beta_cm = p_wimp / (Edm + iM_el); // adim
      gamma_cm = (Edm + iM_el)/(TMath::Sqrt(MD*MD + iM_el*iM_el + 2*Edm*iM_el)); //adim

      Edm_cm = gamma_cm*(Edm - beta_cm*p_wimp);  // GeV          
      E_rec_cm = TMath::Sqrt(Edm_cm*Edm_cm - MD*MD + iM_el*iM_el);  //GeV
      
      rnd_cos_theta_cm = tt->Uniform(-1,1);
      rnd_phi_cm = tt->Uniform(0,2*TMath::Pi());
      
      px_cm = beta_cm*E_rec_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Sin(rnd_phi_cm);  // punta al centro galattico
      py_cm = beta_cm*E_rec_cm*TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2)))*TMath::Cos(rnd_phi_cm);  // punta alla direzione del Cigno
      pz_cm = beta_cm*E_rec_cm*rnd_cos_theta_cm;                                                          // normale al piano galattico
      p_mod_cm = TMath::Sqrt(px_cm*px_cm + py_cm*py_cm + pz_cm*pz_cm);            
      pz_lab = beta_cm*gamma_cm*E_rec_cm + gamma_cm*pz_cm;
            
      E_rec_lab = gamma_cm*E_rec_cm + beta_cm*gamma_cm*pz_cm - iM_el; //GeV
      p_mod_lab = TMath::Sqrt(E_tot_lab*E_tot_lab - iM_el*iM_el); // GeV
      the_rec = TMath::ATan2(TMath::Sqrt((1-TMath::Power(rnd_cos_theta_cm,2))),(gamma_cm*(beta_cm*(E_rec_cm/p_mod_cm)+rnd_cos_theta_cm)));  //buono
      //phi_rec = TMath::ATan2(px_lab,px_cm);
      phi_rec = TMath::ATan2(py_cm,px_lab);

      angle3D=TMath::ATan(TMath::Sqrt(px_cm*px_cm+py_cm*py_cm)/pz_lab);
      angle0=TMath::ATan(py_cm/px_cm);
      theta_x = TMath::ATan(px_cm/pz_lab);
      theta_y = TMath::ATan(py_cm/pz_lab);
      energy = E_rec_lab*1E6;
      /////////////////////////////////////
      
      /*
      angle3D=TMath::ATan(TMath::Sqrt(vT_x*vT_x+vT_y*vT_y)/vT_z);
      angle0=TMath::ATan(vT_y/vT_x);
      theta_x = TMath::ATan(vT_x/vT_z);
      theta_y = TMath::ATan(vT_y/vT_z);
      if(vT_x<0) angle0=TMath::ATan(vT_y/-vT_x);   
      energy = 0.5*iM_el*vT*vT;
      energy = energy*1E6; // keV

      printf("%lf %lf %lf %lf \n",vCM,vT_CM, vT, energy);
      */
      length = gr->Eval(energy);// nm
      
      //cout << "wimp  " << " " << iZ_el << " " << energy << " " << length << endl;
      
    if(length>=threshold && length<=max_range){      
      hl->Fill(length*0.001,atomic_fract_el[i]*A_el[i]*A_el[i]); //um
      he->Fill(energy*0.001,atomic_fract_el[i]*A_el[i]*A_el[i]); //MeV
      //hle->Fill(length*0.001,energy*0.001,atomic_fract_el[i]*A_el[i]*A_el[i]);
      hlen_ene_dm[i]->Fill(energy*0.001,length*0.001,atomic_fract_el[i]*A_el[i]*A_el[i]);
      htthex->Fill(theta_x,atomic_fract_el[i]*A_el[i]*A_el[i]);
      htthey->Fill(theta_y,atomic_fract_el[i]*A_el[i]*A_el[i]);
      htT0CUT->Fill(angle0,atomic_fract_el[i]*A_el[i]*A_el[i]);
      htT0CUT3D->Fill(angle3D,atomic_fract_el[i]*A_el[i]*A_el[i]);

      hlen_dm[i]->Fill(length*0.001,atomic_fract_el[i]*A_el[i]*A_el[i]);
      nwimp++;
    }
    }
    gr->Clear();
  }

  rate_nu = ((float)nnu)/(float (npoints*n_el));
  rate_wimp = ((float)nwimp)/(float (npoints*n_el));

  cout << "Frazione di neutrini selezionati: " << rate_nu*100 << endl;
  cout << "Frazione di wimps selezionate: " << rate_wimp*100 << endl;

  ofstream log(Form("./Logs/log_dm_%.2fMeV_range_%.2fum_%.2fum.txt",MD*1000,threshold/1000.,max_range/1000.));
  log.is_open();
  log << "Frazione di neutrini selezionati: " << rate_nu*100 << endl;
  log << "Frazione di wimps selezionate: " << rate_wimp*100 << endl;	       
  
  
  TFile *fout = new TFile(Form("./Plots/histos_dm_%.1fMeV_range_%.1fum_%.1fum.root",MD*1000,threshold/1000.,max_range/1000.),"RECREATE");


  TCanvas *c01 = new TCanvas("c01","c01",1200,600);
  c01->Divide(2,1);
  c01->cd(1);
  hpdm->Draw("HIST");
  c01->cd(2);
  hpdm0->Draw("HIST");

  TCanvas *c02 = new TCanvas("c02","c02",1200,600);
  c02->Divide(2,1);
  c02->cd(1);
  henu1->Draw("HIST");
  c02->cd(2);
  henu0->Draw("HIST");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(4,2);
  //for(int s=0;s<7;s++){
  //c1->cd(s+1);
  //h1[s]->Draw("");
  //}
  c1->cd(1);
  h1[2]->Draw("");
  h1[2]->SetLineWidth(2);
  h1[2]->SetName("I");
  h1[2]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[2]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(2);
  h1[0]->Draw("");
  h1[0]->SetLineWidth(2);
  h1[0]->SetName("Ag");
  h1[0]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[0]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(3);
  h1[1]->Draw("");
  h1[1]->SetLineWidth(2);
  h1[1]->SetName("Br");
  h1[1]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[1]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(4);
  h1[6]->Draw("");
  h1[6]->SetLineWidth(2);
  h1[6]->SetName("S");
  h1[6]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[6]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(5);
  h1[4]->Draw("");
  h1[4]->SetLineWidth(2);
  h1[4]->SetName("O");
  h1[4]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[4]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(6);
  h1[5]->Draw("");
  h1[5]->SetLineWidth(2);
  h1[5]->SetName("N");
  h1[5]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[5]->GetXaxis()->SetTitleOffset(1.10);
  c1->cd(7);
  h1[3]->Draw("");
  h1[3]->SetLineWidth(2);
  h1[3]->SetName("C");
  h1[3]->GetXaxis()->SetTitle("#theta_{sc} [rad]");
  h1[3]->GetXaxis()->SetTitleOffset(1.10);

  
  
  TCanvas *clen = new TCanvas("clen","clen",800,800);
  hlen_thr->Draw("HIST");
  clen->SetLogy();
  clen->SetLogx();
  hlen_thr->SetLineWidth(2);
  hlen_thr->GetXaxis()->SetTitle("L [#mum]");
  hlen_thr->GetXaxis()->SetTitleOffset(1.20);
  //hlen_thr->Scale(num_ev/hlen_thr->Integral());
  
  TCanvas *cene = new TCanvas("cene","cene",800,800);
  hene_thr->Draw("HIST");
  cene->SetLogy();
  cene->SetLogx();
  hene_thr->SetLineWidth(2);
  hene_thr->GetXaxis()->SetTitle("E [MeV]");
  hene_thr->GetXaxis()->SetTitleOffset(1.20);
  //hene_thr->Scale(num_ev/hene_thr->Integral());

  TCanvas *clen_ene = new TCanvas("clen_ene","clen_ene",800,800);
  hlen_ene[3]->Draw("prof");
  hlen_ene[3]->SetStats(0);
  hlen_ene[3]->GetXaxis()->SetTitle("E [MeV]");
  hlen_ene[3]->GetYaxis()->SetTitle("L [#mum]");
  hlen_ene[3]->GetXaxis()->SetTitleOffset(1.20);
  hlen_ene[3]->GetYaxis()->SetTitleOffset(1.20);
  clen_ene->SetLogy();
  clen_ene->SetLogx();
  hlen_ene[3]->SetMarkerColor(kBlue);
  hlen_ene[2]->Draw("same");
  hlen_ene[2]->SetMarkerColor(kPink+1);
  hlen_ene[6]->Draw("same");
  hlen_ene[6]->SetMarkerColor(kGreen);
  hlen_ene[4]->Draw("same");
  hlen_ene[4]->SetMarkerColor(kViolet);
  hlen_ene[5]->Draw("same");
  hlen_ene[5]->SetMarkerColor(kCyan);
  hlen_ene[1]->Draw("same");
  hlen_ene[1]->SetMarkerColor(kOrange);
  hlen_ene[0]->Draw("same");
  hlen_ene[0]->SetMarkerColor(kRed);
  
  TCanvas *cphi = new TCanvas("cphi","cphi",1500,500);
  cphi->Divide(3,1);
  cphi->cd(1);
  hphirec->Draw("HIST");
  //hphirec->Rebin(3);
  hphirec->SetLineWidth(2);
  hphirec->GetXaxis()->SetTitle("tz [rad]");
  hphirec->SetTitleOffset(1.10);
  //hphirec->Scale(num_ev/hphirec->Integral());
  cphi->cd(2);
  htthe_recx->Draw("HIST");
  htthe_recx->SetLineWidth(2);
  htthe_recx->GetXaxis()->SetTitle("tx [rad]");
  cphi->cd(3);
  htthe_recy->Draw("HIST");
  htthe_recy->SetLineWidth(2);
  htthe_recy->GetXaxis()->SetTitle("ty [rad]");
  
  TCanvas *cthe = new TCanvas("cthe","cthe",800,800);
  htherec->Draw("HIST");
  htherec->SetLineWidth(2);
  htherec->GetXaxis()->SetTitle("#theta [rad]");
  htherec->SetTitleOffset(1.10);
  //htherec->Scale(num_ev/htherec->Integral());

  ////////// WIMP ////////////////////

  TCanvas *cenergy = new TCanvas("cenergy","",800,800);
  he->Draw("HIST");
  cenergy->SetLogy();
  cenergy->SetLogx();
  he->GetXaxis()->SetTitle("E [MeV]");
  he->SetLineWidth(2);
  he->GetXaxis()->SetTitleOffset(1.20);
  //he->SetLineColor(kBlack);

  TCanvas *clength = new TCanvas("clength","",800,800);
  hl->Draw("HIST");
  clength->SetLogy();
  clength->SetLogx();
  hl->GetXaxis()->SetTitle("L [#mum]");
  //hl->SetLineColor(kBlack);
  hl->SetLineWidth(2);
  hl->GetXaxis()->SetTitleOffset(1.20);

  TCanvas *cee = new TCanvas("cee","",800,800);
  hene_thr->Draw("HIST");
  cee->SetLogy();
  cee->SetLogx();
  he->Draw("sames&&HIST");
  hene_thr->GetXaxis()->SetTitle("E [MeV]");
  hene_thr->SetLineWidth(2);
  hene_thr->GetXaxis()->SetTitleOffset(1.20);
  hene_thr->SetTitle("E_{#nu} = GeV");
  he->SetLineWidth(2);
  he->SetLineColor(kRed);
  he->Scale(1./he->Integral());
  hene_thr->Scale(1./hene_thr->Integral());
  he->SetTitle(Form("M_{dm} = %.2f MeV",MD*1000));
  cee->BuildLegend();

  TCanvas *cll = new TCanvas("cll","",800,800);
  hlen_thr->Draw("HIST");
  cll->SetLogy();
  cll->SetLogx();
  hl->Draw("sames && HIST");
  hlen_thr->GetXaxis()->SetTitle("L [#mum]");
  hlen_thr->SetLineWidth(2);
  hlen_thr->GetXaxis()->SetTitleOffset(1.20);
  hlen_thr->SetTitle("E_{#nu} = GeV");
  hl->SetLineWidth(2);
  hl->SetLineColor(kRed);
  hl->SetTitle(Form("M_{dm} = %.2f MeV",MD*1000));
  hl->Scale(1./hl->Integral());
  hlen_thr->Scale(1./hlen_thr->Integral());
  cll->BuildLegend();


  TCanvas *cle = new TCanvas("clen_ene_dm","clen_ene_dm",800,800);
  hlen_ene_dm[3]->Draw("prof");
  hlen_ene_dm[3]->SetStats(0);
  hlen_ene_dm[3]->GetXaxis()->SetTitle("E [MeV]");
  hlen_ene_dm[3]->GetYaxis()->SetTitle("L [#mum]");
  hlen_ene_dm[3]->GetXaxis()->SetTitleOffset(1.20);
  hlen_ene_dm[3]->GetYaxis()->SetTitleOffset(1.20);
  cle->SetLogy();
  cle->SetLogx();
  hlen_ene_dm[3]->SetMarkerColor(kBlue);
  hlen_ene_dm[2]->Draw("same");
  hlen_ene_dm[2]->SetMarkerColor(kPink+1);
  hlen_ene_dm[6]->Draw("same");
  hlen_ene_dm[6]->SetMarkerColor(kGreen);
  hlen_ene_dm[4]->Draw("same");
  hlen_ene_dm[4]->SetMarkerColor(kViolet);
  hlen_ene_dm[5]->Draw("same");
  hlen_ene_dm[5]->SetMarkerColor(kCyan);
  hlen_ene_dm[1]->Draw("same");
  hlen_ene_dm[1]->SetMarkerColor(kOrange);
  hlen_ene_dm[0]->Draw("same");
  hlen_ene_dm[0]->SetMarkerColor(kRed);

  TCanvas *cangleTCUT = new TCanvas("cangleTCUT","",1500,500);
  cangleTCUT->Divide(3,1);
  cangleTCUT->cd(1);
  //htT0CUT->SetTitle("Recoil Angle 2D");
  htT0CUT->GetXaxis()->SetTitle("tz [rad]");
  htT0CUT->SetLineWidth(2);
  htT0CUT->Draw("HIST");
  cangleTCUT->cd(2);
  htthex->GetXaxis()->SetTitle("tx [rad]");
  htthex->SetLineWidth(2);
  htthex->Draw("HIST");
  cangleTCUT->cd(3);
  htthey->GetXaxis()->SetTitle("ty [rad]");
  htthey->SetLineWidth(2);
  htthey->Draw("HIST");

  TCanvas *cangle3DTCUT = new TCanvas("cangle3DTCUT","",800,800);
  //htT0CUT3D->SetTitle("Recoil Angle 3D");
  htT0CUT3D->GetXaxis()->SetTitle("#theta [rad]");
  htT0CUT3D->Draw("HIST");
  htT0CUT3D->SetLineWidth(2);

  TCanvas *ctt = new TCanvas("ctt","",800,800);
  htherec->Draw("HIST");
  htT0CUT3D->Draw("sames&&HIST");
  htherec->GetXaxis()->SetTitle("#theta [rad]");
  htherec->SetLineWidth(2);
  htherec->GetXaxis()->SetTitleOffset(1.20);
  htherec->SetTitle("E_{#nu} = GeV");
  htT0CUT3D->SetLineWidth(2);
  htT0CUT3D->SetLineColor(kRed);
  htT0CUT3D->SetTitle(Form("M_{dm} = %.2f MeV",MD*1000));
  htherec->Scale(1./htherec->Integral());
  htT0CUT3D->Scale(1./htT0CUT3D->Integral());
  ctt->BuildLegend();

  TCanvas *cff = new TCanvas("cff","",800,800);
  htthe_recx->Draw("HIST");
  htthex->Draw("sames&&HIST");
  htthe_recx->GetXaxis()->SetTitle("#phi [rad]");
  htthe_recx->SetLineWidth(2);
  htthe_recx->GetXaxis()->SetTitleOffset(1.20);
  htthe_recx->SetTitle("E_{#nu} = GeV");
  htthex->SetLineWidth(2);
  htthex->SetLineColor(kRed);
  htthex->SetTitle(Form("M_{dm} = %.2f MeV",MD*1000));
  htthe_recx->Scale(1./htthe_recx->Integral());
  htthex->Scale(1./htthex->Integral());
  cff->BuildLegend();


  TCanvas *clen_nu = new TCanvas("clen_nu","clen_nu",800,800);
  hlen_nu[3]->Draw("HIST");
  hlen_nu[3]->SetStats(0);
  hlen_nu[3]->GetXaxis()->SetTitle("L [#mum]");
  hlen_nu[3]->GetXaxis()->SetTitleOffset(1.20);
  clen_nu->SetLogy();
  clen_nu->SetLogx();
  hlen_nu[3]->SetLineColor(kBlue);
  hlen_nu[2]->Draw("HIST same");
  hlen_nu[2]->SetLineColor(kPink+1);
  hlen_nu[6]->Draw("HIST same");
  hlen_nu[6]->SetLineColor(kGreen);
  hlen_nu[4]->Draw("HIST same");
  hlen_nu[4]->SetLineColor(kViolet);
  hlen_nu[5]->Draw("HIST same");
  hlen_nu[5]->SetLineColor(kCyan);
  hlen_nu[1]->Draw("HIST same");
  hlen_nu[1]->SetLineColor(kOrange);
  hlen_nu[0]->Draw("HIST same");
  hlen_nu[0]->SetLineColor(kRed);
  clen_nu->BuildLegend();

  TCanvas *clen_dm = new TCanvas("clen_dm","clen_dm",800,800);
  hlen_dm[3]->Draw("HIST");
  hlen_dm[3]->SetStats(0);
  hlen_dm[3]->GetXaxis()->SetTitle("L [#mum]");
  hlen_dm[3]->GetXaxis()->SetTitleOffset(1.20);
  clen_dm->SetLogy();
  clen_dm->SetLogx();
  hlen_dm[3]->SetLineColor(kBlue);
  hlen_dm[2]->Draw("HIST same");
  hlen_dm[2]->SetLineColor(kPink+1);
  hlen_dm[6]->Draw("HIST same");
  hlen_dm[6]->SetLineColor(kGreen);
  hlen_dm[4]->Draw("HIST same");
  hlen_dm[4]->SetLineColor(kViolet);
  hlen_dm[5]->Draw("HIST same");
  hlen_dm[5]->SetLineColor(kCyan);
  hlen_dm[1]->Draw("HIST same");
  hlen_dm[1]->SetLineColor(kOrange);
  hlen_dm[0]->Draw("HIST same");
  hlen_dm[0]->SetLineColor(kRed);
  clen_dm->BuildLegend();
  
  
  c01->Write();
  c02->Write();
  c1->Write();
  clen->Write();
  cene->Write();
  clen_ene->Write();
  cphi->Write();
  cthe->Write();
  cenergy->Write();
  clength->Write();
  cangle3DTCUT->Write();
  cangleTCUT->Write();
  cle->Write();
  cll->Write();
  cee->Write();
  ctt->Write();
  cff->Write();
  clen_nu->Write();
  clen_dm->Write();
  fout->Close();
  
}
