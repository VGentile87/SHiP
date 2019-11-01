// Classificazione dei vertici simulati ricostruiti in FEDRA
// Training per analisi multivariata
// V. Gentile 2019

#include "Definitions.h"
#include <TH2.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <TF1.h>
#include <TStyle.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TObjArray.h>

#include <TParticlePDG.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <EdbVertex.h>

#include <algorithm>    // std::equal
#include <vector>       // std::vector

using namespace std;


// OUTPUT FILE STRUCTURE

struct Vertex
{
  Int_t event;   // index of the tree 
  Int_t fe_id;  // fedra_id
  Int_t mc_id;  // monte carlo id
  Int_t ntrk;   // number of tracks attached to the vertex
  Int_t mc_nd[2] ; // number montecarlo charged daugthers
  Double_t mc_dz[2];   // delta z between charm and primary vertex
  Double_t mc_dl[2];   // decay length between charm and primary vertex
  Double_t x;
  Double_t y;
  Double_t z;
  };
  
struct Track
{
  vector<int> icharm;
  vector<int> fe_itrk;              // number fedra charged daugthers 
  vector<int>  fe_id;
  vector<int>  mc_id;
  vector<int>  pdg;                // track pdg
  vector<double>  fe_dz;           // starting coordinates of the track
  vector<double> fe_dl;           // 2d angle projections in xz and yz planes
  vector<double> ip;              // impact parameter
  vector<double> ka;              // kink angle
  vector<int> nseg;               // number of segments
  vector<int> type;               // 0 - track alone, 1 - in primary vtx, 2 - in other vtx same event, 3 - in other vtx, 4 - fake vertex
  vector<double> x;              // xpos
  vector<double> y;              // xpos
  vector<double> z;              // xpos
  vector<double> min_xy;
  vector<double> min_xyz;
  vector<int> frd_xy;
  vector<int> frd_xyz;
};

Vertex vtx;
Track trk;

void FillVertex(int ievent, int ife_id, int imc_id, int intrk, int imc_nd1, int imc_nd2, double imc_dz1, double imc_dz2, double imc_dl1, double imc_dl2, double xv, double yv, double zv){
  vtx.event = ievent;
  vtx.fe_id = ife_id;
  vtx.mc_id = imc_id;
  vtx.ntrk = intrk;
  vtx.mc_nd[0] = imc_nd1;
  vtx.mc_nd[1] = imc_nd2;
  vtx.mc_dz[0] = imc_dz1;
  vtx.mc_dz[1] = imc_dz2;
  vtx.mc_dl[0] = imc_dl1;
  vtx.mc_dl[1] = imc_dl2;
}
void FillTrack(int ife_itrk, int iicharm, int ife_id, int imc_id, int ipdg, int inseg,  int itype, double ife_dz, double ife_dl, double iip, double ika, double xt, double yt, double zt){
  trk.fe_itrk.push_back(ife_itrk);
  trk.icharm.push_back(iicharm);
  trk.fe_id.push_back(ife_id);
  trk.mc_id.push_back(imc_id);
  trk.pdg.push_back(ipdg);
  trk.fe_dz.push_back(ife_dz);
  trk.fe_dl.push_back(ife_dl);
  trk.ip.push_back(iip);
  trk.ka.push_back(ika);
  trk.nseg.push_back(inseg);
  trk.type.push_back(itype);
  trk.x.push_back(xt);
  trk.y.push_back(yt);
  trk.z.push_back(zt);
}

void Set0() {
  vtx.event=0, vtx.fe_id=0, vtx.mc_id=0, vtx.ntrk=0, vtx.mc_nd[0]=0, vtx.mc_nd[1]=0, vtx.mc_dz[0]=0, vtx.mc_dz[1]=0, vtx.mc_dl[0]=0, vtx.mc_dl[1]=0;
  vtx.x=0, vtx.y=0, vtx.z=0;
  trk.fe_itrk.clear(), trk.icharm.clear(), trk.fe_id.clear(), trk.mc_id.clear(), trk.pdg.clear(), trk.fe_dz.clear(), trk.fe_dl.clear(), trk.ip.clear(), trk.ka.clear(), trk.nseg.clear(), trk.type.clear(), trk.x.clear(), trk.y.clear(), trk.z.clear(), trk.min_xy.clear(), trk.min_xyz.clear(), trk.frd_xy.clear(), trk.frd_xyz.clear();
}

// FUNCTIONS

// Impact parameter of a track with respect to a vertex
float IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
 
 float dz = -vertexpos(2) + trackstartpos(2);
 float ipx = -tracktx * dz + trackstartpos(0) - vertexpos(0);
 float ipy = -trackty * dz + trackstartpos(1) - vertexpos(1);

 float ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

 return ip;
}

// IPmax over IPrms calculated on segments of a track with respet to a vertex
void SegIPtoVertex2(TVector3 vertexpos, int nseg, double* seg_x, double* seg_tx, double* seg_y, double* seg_ty, double* seg_z, float IpSeg[2]){
 
  float ipseg[nseg];
  float delta_ipseg[nseg-1];
  float ipseg_nomax[nseg-2];
  for (int i = 0; i < nseg; i++){   
    float dz = -vertexpos(2) + seg_z[i];
    float ipx = -seg_tx[i] * dz + seg_x[i] - vertexpos(0);
    float ipy = -seg_ty[i] * dz + seg_y[i] - vertexpos(1);
    ipseg[i] = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
  }

  for (int i = 0; i < nseg-1; i++){
    delta_ipseg[i] = ipseg[i+1] - ipseg[i];
  }

  float ipsegmax = TMath::MaxElement(nseg-1, delta_ipseg);
  int index=0;
  IpSeg[0]=0;
  for (int i = 0; i < nseg-2; i++){
    if(delta_ipseg[i]!=ipsegmax){
      ipseg_nomax[index]=delta_ipseg[i];
      index++;
    }
    else IpSeg[0]=i;
  }
  IpSeg[0]++;
  float ipsegrms = TMath::RMS(index, ipseg_nomax);
  
  //IpSeg[1] = ipsegmax/ipsegrms;
  if(ipsegrms!=0) IpSeg[1] = ipsegmax/ipsegrms;
  else IpSeg[1] = -1;
}

//// Kinkmax over Kinkrms calculated on segments of a track with respet to a vertex
void FedraTrackKink2(int nseg, double* seg_tx, double* seg_ty, float *dtheta){

  // cout << "num seg " << nseg << endl;
  
  float kinkangles[nseg-1];
  float kinkangles_nomax[nseg-2];
  //loop on subsequent segment pairs of the tracks
 
  
  for (int i = 0; i < nseg-1; i++){
    kinkangles[i]=TMath::Sqrt(pow(seg_tx[i+1]-seg_tx[i],2)+pow(seg_ty[i+1]-seg_ty[i],2));
  }
  //getting maximum and rms
  float deltathetamax = TMath::MaxElement(nseg-1, kinkangles);
  int index=0;
  for (int i = 0; i < nseg-2; i++){
    if(kinkangles[i]!=deltathetamax){
      kinkangles_nomax[index]=kinkangles[i];
      index++;
    }
    else dtheta[0]=i;
  }
  dtheta[0]++;
  float deltathetarms = TMath::RMS(index, kinkangles_nomax);
  //dtheta[1] = deltathetamax/deltathetarms;
  if(deltathetarms!=0)dtheta[1] = deltathetamax/deltathetarms;
  else dtheta[1]= -1;
}

// vector elements comparison
static bool compareVectors(vector<int> a, vector<int> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return (a == b);
}


// Function to find daughters linked to the primary vertex

tuple <vector<vector<int>>, vector<double>, vector<double>, vector<double>, vector<int>, vector<int>,vector<vector<int>>, vector<double>, vector<double>> dbscan(double epsilon, double epsilon_phi, int min_points, int ntrack, vector<double> xpos, vector<double> ypos, vector<double> zpos, vector<double> tx, vector<double> ty){

  //// CLUSTERING PER EVENTO
     int icl=0;
     int icl_rel=0;
     int sphere_points=0;
     int sphere_points2=0;

     std::vector<std::vector<bool>> visited;
     std::vector<bool> visited2;
     std::vector<int> n_el;
     std::vector<int> n_el_rel;
     std::vector<double> mean_dist;
     std::vector<double> mean_dist_phi;
     std::vector<double> mean_dist_vtx;
     std::vector<double> mean_dist_vtx_rel;
     std::vector<double> first_seg_dist;
     std::vector <std::vector<int> > cl_el;
     std::vector <std::vector<int> > cl_el_rel;
     float dist=0;
     float dist_phi=0;

     float min_dist=10000000000;
     int tmp_1=-1;
     int tmp_2=-1;

     vector<bool> taken_track;
     

     vector<int> same_z(ntrack);
     same_z.clear();

     vector<double> x_cross;
     vector<double> y_cross;
     vector<double> xy_cross;
     vector<double> dist_rel;
     vector<double> phi_cross;
     vector<int> cross_id1;
     vector<int> cross_id2;
     int index_cross=0;

     int ncross = ntrack*(ntrack-1)/2.;

     vector<int> count_far;
     vector<int> count_close;
     count_far.clear();
     count_close.clear();
     count_far.resize(ntrack);
     count_close.resize(ntrack);

     icl=0;
     icl_rel=0;
     sphere_points=0;
     sphere_points2=0;
     dist=0;
     dist_phi=0;

     taken_track.clear();
     taken_track.resize(ntrack);
     visited.clear();
     visited2.clear();
     n_el.clear();
     mean_dist.clear();
     mean_dist_phi.clear();
     first_seg_dist.clear();

     cl_el.clear();
     visited2.resize(ncross);
     n_el.resize(ncross);
     mean_dist.resize(ncross);
     mean_dist_phi.resize(ncross);
     mean_dist_vtx.resize(ncross);
     cl_el.resize(ncross, std::vector<int>(ncross));
     visited.resize(ncross,std::vector<bool>(ncross));
     

     /// CONTA LE OCCORRENZE PER Z 
     bool zflag=false;;
     for(int h=0; h<ntrack; h++){
       //cout << "dbscan " << zpos.at(h) << " " << xpos.at(h) << " " << ypos.at(h) << endl;
       zflag=0;
       for (int i=0;i<h;i++){
	 if (zpos[h]==zpos[i]){ // controlla che quel numero non sia già stato considerato in precedenza
	   zflag=true;
	   break;
	 }
       }
	 if (zflag==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	   same_z[h]=1; 
	   for (int i=h+1;i<ntrack;i++){		 
	     if(zpos[h] == zpos[i]) same_z[h]++; 
	   }
	 }
	 //cout << "samez " << same_z[h] << endl;
     }

     // aggiusto le posizioni xy sullo stesso piano z
     vector<int>::iterator max_z = max_element(same_z.begin(), same_z.end()); // c++11
     auto max_z_occurrence = *max_element(same_z.begin(), same_z.end()); // c++11    
     int iz = distance(same_z.begin(), max_z);
     
     // calcola le posizioni nel piano trasverso xy
     for(int j=0;j<ntrack;j++){
       //cout << "dbscan2 " << j << " " << ntrack << " " <<  zpos.at(j) << " " << xpos.at(j) << " " << ypos.at(j) << endl;
       if(zpos.at(iz)==zpos.at(j))continue;
       else{
	 
	 xpos[j] = xpos[j] - tx[j]*(zpos.at(j)-zpos.at(iz));
	 ypos[j] = ypos[j] - ty[j]*(zpos.at(j)-zpos.at(iz));
	 zpos[j] = zpos[iz];
       }
       //cout << "dbscan2 " << j << " " << ntrack << " " <<  zpos.at(j) << " " << xpos.at(j) << " " << ypos.at(j) << endl;
     }

 
     // calcolo dei punti di intersezione in xy
     for(int j=0;j<ntrack;j++){
       taken_track.push_back(false);
       for(int l=(j+1);l<ntrack;l++){

	 x_cross.push_back((ypos[l]-ypos[j] + (ty[j]/tx[j])*xpos[j] - (ty[l]/tx[l])*xpos[l])/(ty[j]/tx[j] - ty[l]/tx[l])) ;
	 y_cross.push_back(ypos[j] + (ty[j]/tx[j])*(x_cross[index_cross] - xpos[j]));
	 cross_id1.push_back(j);
	 cross_id2.push_back(l);
	 
	 xy_cross.push_back(TMath::Sqrt(TMath::Power(x_cross[index_cross]-vx,2)+TMath::Power(y_cross[index_cross]-vy,2)));
	 dist_rel.push_back(TMath::Sqrt(TMath::Power(xpos[l]-xpos[j],2)+TMath::Power(ypos[l]-ypos[j],2)));
	 phi_cross.push_back(TMath::ATan2(y_cross[index_cross]-vy, x_cross[index_cross]-vx));
	 
	 if(l!=j){
	   if(xy_cross[index_cross]>10)count_far[l]++;
	   if(xy_cross[index_cross]<10)count_close[j]++;
	   //cout << j << " " << l << " " << xy_cross[index_cross] << " " << phi_cross[index_cross] << " " << x_cross[index_cross] << " " << y_cross[index_cross] << " " << dist_rel[index_cross] << endl;
	 }
	 index_cross++;
       }       
     }


     // do-while per un dbscan dinamico     
     do{

       icl=0;
       icl_rel=0;
       sphere_points=0;
       sphere_points2=0;
       dist=0;
       dist_phi=0;
       
       taken_track.clear();
       taken_track.resize(ntrack);
       visited.clear();
       visited2.clear();
       n_el.clear();
       n_el_rel.clear();
       mean_dist.clear();
       mean_dist_phi.clear();
       mean_dist_vtx.clear();
       mean_dist_vtx_rel.clear();
       first_seg_dist.clear();

       cl_el.clear();
       cl_el_rel.clear();
       
       visited.resize(ncross,std::vector<bool>(ncross));
       visited2.resize(ncross);
       n_el.resize(ncross);
       mean_dist.resize(ncross);
       mean_dist_phi.resize(ncross);
       mean_dist_vtx.resize(ncross);
       cl_el.resize(ncross, std::vector<int>(ncross));

        for(int j=0;j<ncross;j++){
	  for(int k=0;k<ncross;k++){
	    visited[j][k]=false;
	  }
	}
       
       
       for(int j=0;j<ncross;j++){
	 //cout << "dbscan3 " <<j << " " <<  zpos.at(j) << " " << xpos.at(j) << " " << ypos.at(j) << endl;
	 //cout << "dbscan3 " <<j << " " <<  xy_cross.at(j) << " " << phi_cross.at(j)  << endl;

	 
	 int npoints=0;
	 for(int k=0;k<taken_track.size();k++){
	   taken_track[k]=false;
	 }
	 
	 double sum_dist=0;
	 double sum_dist_phi=0;
	 double sum_dist_vtx=0;
	   
	 for(int k=0;k<ncross;k++){ // distanza relativa
	   dist = xy_cross.at(k)-xy_cross.at(j);
	   dist_phi = phi_cross.at(k)-phi_cross.at(j);
	   //cout << "dist " << dist << " " << dist_phi << endl;
	     if(dist < min_dist && k!=j){
	       min_dist = dist;
	       tmp_1=k;
	       tmp_2=j;
	     }
	     if(abs(dist)<epsilon && abs(dist_phi)<epsilon_phi && k!=j && !visited[j][k] && !visited[k][j]/*&& !visited2.at(j)*/){
	       //cout << "dist " << cross_id1.at(j) << " " << cross_id2.at(j) << " " << cross_id1.at(k) << " " << cross_id2.at(k) << " " <<  dist << " " << dist_phi << endl;

	       //visited2[k]=true;
	       //cout <<"a "<< icl << " " <<  j << " " << k << " " << sphere_points << " " << dist <<  endl;
	       //cl_el[icl][sphere_points]=k;

	       if(!taken_track[cross_id1.at(j)]){
		 cl_el[icl][npoints]=cross_id1.at(j);
		 taken_track[cross_id1.at(j)]=true;
		 npoints++;
	       }

	       if(!taken_track[cross_id2.at(j)]){
		 cl_el[icl][npoints]=cross_id2.at(j);
		 taken_track[cross_id2.at(j)]=true;
		 npoints++;
	       }

	       if(!taken_track[cross_id1.at(k)]){
		 cl_el[icl][npoints]=cross_id1.at(k);
		 taken_track[cross_id1.at(k)]=true;
		 npoints++;
	       }

	       if(!taken_track[cross_id2.at(k)]){
		 cl_el[icl][npoints]=cross_id2.at(k);
		 taken_track[cross_id2.at(k)]=true;
		 npoints++;
	       }

	       sum_dist += dist;
	       sum_dist_phi += dist_phi;
	       
	       sphere_points++; // crash code

	       visited[j][k]=true;
	       visited[k][j]=true;
	     }
	 }
	 // cluster di 3 o più elementi
	   if(sphere_points!=0){
	     
	     n_el[icl]=npoints;
	     
	       mean_dist[icl] = sum_dist/sphere_points;
	       mean_dist_phi[icl] = sum_dist_phi/sphere_points;
	       
	       for(int k=0;k<npoints;k++){
		 //cout << "element " << k << " " << cl_el[icl][k] << endl;
		 sum_dist_vtx = xy_cross.at(cl_el[icl][k]);
	       }	     
	       mean_dist_vtx[icl] = sum_dist_vtx / npoints;
	       //if(npoints!=0)cout << "cluster_info " << icl << " " << sphere_points << " " << npoints << " " << n_el[icl] << " " << mean_dist[icl] << " " << mean_dist_phi[icl] << " " << mean_dist_vtx[icl] << endl;
	       
	        icl++;

	   }

	   // per un cluster di soli due elementi
	   if(dist_rel.at(j)<10){
	     
	     n_el_rel.push_back(2);
	     vector<int> tmp_el;
	     tmp_el.push_back(cross_id1.at(j));
	     tmp_el.push_back(cross_id2.at(j));	       
	     cl_el_rel.push_back(tmp_el);
	     mean_dist_vtx_rel.push_back(xy_cross.at(j));
	     first_seg_dist.push_back(dist_rel.at(j));
	     
	     //cout << "cluster_info " << icl_rel << " " << cross_id1.at(j) << " " << cross_id2.at(j) << " " << sphere_points << " " << npoints << " " << n_el_rel[icl_rel] << " " << first_seg_dist[icl_rel] << endl;
	     icl_rel++;
	     
	   }
	   
	   }
       //cout << "epsilon " << epsilon << " " << " max_el_cl " << *max_element(n_el.begin(), n_el.end()) << endl;
       epsilon -= 0.1;
       epsilon_phi -=0.01;
       }while(*max_element(n_el.begin(), n_el.end())>(ntrack*0.5));
     
 
     return make_tuple(cl_el,mean_dist,mean_dist_phi,mean_dist_vtx,n_el,n_el_rel,cl_el_rel,mean_dist_vtx_rel,first_seg_dist);
     
}


void Loop()
{

   ofstream log_vtx;
  log_vtx.open("log_vtx_search.txt",ios::out);
  log_vtx << "//----V. GENTILE 10/03/2019 ------// \n VERTEX SEARCH ON MC SIMULATION" << endl;
  
  ofstream log_decay;
  log_decay.open("log_data_search.txt",ios::out);
  log_decay << "//----V. GENTILE ------// \n VERTEX SEARCH ON DATA" << endl;

  TFile * fedrafile = TFile::Open("vertextree_newformat.root");  // full
  //TFile * trackfile = TFile::Open("verticesandtracks.root");  // full

  //TFile * bdtfile = TFile::Open("vtx_BDT_first_quarter_15_04_19_evaluated.root");  // full

  TFile * simulationfile = TFile::Open("ship.conical.Pythia8CharmOnly-TGeant4.root");
 

  TTree *simtree = (TTree*) simulationfile->Get("cbmsim");     
  TTree *vtx_fedratree = (TTree*) fedrafile->Get("vtx");
  //TTree *trk_fedratree = (TTree*) trackfile->Get("tracks");
  //TTree *bdt_tree = (TTree*) bdtfile->Get("bdt");

  if(vtx_fedratree == 0) return;
      
  vtx_reader_Fedra(vtx_fedratree);
  //tracks_reader_Fedra(trk_fedratree);  
  //bdt_reader(bdt_tree);
  cbmsim_reader_MC(simtree);

  //trk_fedratree->BuildIndex("trid");

  
  /* BDT
  Long64_t nentries_bdt = bdt_tree->GetEntriesFast();
  log_decay<<"tot entries bdt "<< nentries_bdt << endl;

  vector<double> vtx_bdt;
  vtx_bdt.clear();

  for (Long64_t ientry=0; ientry<nentries_bdt;ientry++) { // Vertex index
    
    if(ientry%100==0)cout << "BDT " << ientry << endl;
    
    bdt_tree->GetEntry(ientry);

    vtx_bdt.push_back(bdt_value);
    //cout << ientry << " " << bdt_value <<" " << vtx_bdt.at(ientry) << endl;

  }
  */
  
  Long64_t nentries_vtx_fedra = vtx_fedratree->GetEntriesFast();
  log_decay<<"tot entries vtx fedra "<< nentries_vtx_fedra << endl;

  //Long64_t nentries_trk_fedra = trk_fedratree->GetEntriesFast();
  //log_decay<<"tot entries trk fedra "<< nentries_trk_fedra << endl;
   
  vtx_fedratree->BuildIndex("vID");

  // INPUT FILE
     ifstream input;
   int fa,fb,fc,fd,fe,ff,fg,fh,fi,fl,fm,fn;
   const int nsim=1000;
   int iline=0;
   
   std::vector<std::vector<int> > fe_trk_id;
   std::vector<std::vector<int> > fe_vtx_trk;
   std::vector<std::vector<int> > mc_chr_id;
   std::vector<std::vector<int> > vtx_type;
   fe_trk_id.resize(nsim);
   fe_vtx_trk.resize(nsim);
   mc_chr_id.resize(nsim);
   vtx_type.resize(nsim);
   // find primary and charm vertex

   // (ANTONIO FILES)

   input.open("MC_vertexlist.txt",ios::in);
   while(1){
     input >> fa >> fb >> fc >> fd >> fe >> ff  >> fg >> fh >> fi >> fl >> fm >> fn;
     //cout << fd << " " << ff << " " << fb << " " << fa << endl;
     fe_trk_id[fd].push_back(fc);
     fe_vtx_trk[fd].push_back(fb);
     mc_chr_id[fd].push_back(ff);
     vtx_type[fd].push_back(fn);
     if(input.eof())break;
     if(!input.good())break;
   }
   
   // end of primary and charm vertex


  //bool ev_with_primary[nsim];
  
  for(int i=0;i<nsim;i++){
    //ev_with_primary[nsim]=false;
    for(int j=0;j<fe_trk_id[i].size();j++){
      //if(vtx_type[i][j]==)
      //cout << i << " " << j << " " << fe_trk_id[i][j] << " " << mc_chr_id[i][j] <<  endl;
    }
  }

  ////////

  int vbinx=13;//52;
  int vbiny=6;//24;
  int vbinz=4;//16;

  double epsilon=1; //micron
  double epsilon_phi=0.1; //rad
  int min_points=1;
  
  TH3D * hvertex = new TH3D("hv","hv",vbinx,0,130000,vbiny,20000,80000,vbinz,-37000,3000);

  int width=hvertex->GetNbinsX();
  int height=hvertex->GetNbinsY();
  int depth=hvertex->GetNbinsZ();

  int xmin = hvertex->GetXaxis()->GetXmin();
  int xmax = hvertex->GetXaxis()->GetXmax();
  int ymin = hvertex->GetYaxis()->GetXmin();
  int ymax = hvertex->GetYaxis()->GetXmax();
  int zmin = hvertex->GetZaxis()->GetXmin();
  int zmax = hvertex->GetZaxis()->GetXmax();
  float xrange = (xmax-xmin)/width;
  float yrange = (ymax-ymin)/height;
  float zrange = (zmax-zmin)/depth; 
  //int counts=3;
  
  //cout << width << " " << height << " " << depth << endl;
  //cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << endl;
  //cout << xrange << " " << yrange << " " << zrange << endl;

  log_decay << "\nxrange [" << xmin << " , " << xmax << "] ; yrange [" << ymin << " , " << ymax << "] ; zrange [" << zmin << " , " << zmax << "]" << endl;
  log_decay << "Cube dimensions XYZ (" << xrange << " , " << yrange << " , " << zrange << ")" << endl;
  log_decay << "N cubes in XYZ "<< width << " " << height << " " << depth << endl;
  

  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_list;
    std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_list_nearby;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_list;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_entry;   
  vtx_list.resize(width);
  vtx_list_nearby.resize(width);
  trk_list.resize(width);
  trk_entry.resize(width);
  for(int i=0;i<width;i++){
    //y axis size
    vtx_list[i].resize(height);
    vtx_list_nearby[i].resize(height);
    trk_list[i].resize(height);
    trk_entry[i].resize(height);
    for(int j=0;j<height;j++){
      //z axis size
      vtx_list[i][j].resize(depth);
      vtx_list_nearby[i][j].resize(depth);
      trk_list[i][j].resize(depth);
      trk_entry[i][j].resize(depth);
    }
  }
  

  std::vector<bool> event_primary_proton;
  std::vector<int> vtx_trackmother;
  std::vector<int> vtx_trackprocess;
  std::vector<int> vtx_trackpdgcode;
  std::vector<int> vtx_motherpdgcode;
  std::vector<int> vtx_trackid;
  std::vector<double> vtx_trackstartX;
  std::vector<double> vtx_trackstartY;
  std::vector<double> vtx_trackstartZ;
  std::vector<double> vtx_trackeventId;
  std::vector<double> vtx_trackmom;
  
  float mean_freq=0;
  int max_freq=0;
  float mean_ev_freq=0;
  int max_ev_freq=0;
  int mp_motherID=0;
  int mp_eventID=0;
  int mp_procID=0;
  int mp_pdgID=0;
  int mp_ev_pdgID=0;
  int n_electrons=0;

  bool flag_trk=0;
  bool flag_ev=0;
  //bool evt_good=false;
  //bool daughter_proton=false;
  bool vtx_good=false;
  

  float mp_vx=0;
  float mp_vy=0;
  float mp_vz=0;
  float mp_vdx=-1;
  float mp_vdy=-1;
  float mp_vdz=-1;
  TH1F *hoccurrence = new TH1F("occurrence / ntracks","",50,0,1);
  TH1F *hdeltaz = new TH1F("dz sec vtx","",500,-100,20000);
  TH1F *hntrk = new TH1F("ntrk sec vtx","vtx2_ntrk",50,0,50);
  TH1F *hnseg = new TH1F("nseg trk sec vtx","vtx2_nseg",50,0,50);
  TH1F *hflag = new TH1F("flag sec vtx","vtx2_flag",20,-10,10);
  TH1F *hka2 = new TH1F("ka trk sec vtx to prim vtx","vtx1_ka2",100,0,1);
  TH1F *hip2 = new TH1F("ip trk sec vtx to prim vtx","vtx1_ip2",100,0,1000);
  TH2F *hkaip2 = new TH2F("ka vs ip trk sec vtx to prim vtx","vtx1_kaip2",100,0,1000,100,0,1);

  TH1F *hnseg1st = new TH1F("nseg trk primary vtx","vtx1_nseg",50,0,50);
  TH1F *hka = new TH1F("ka trk prim vtx","vtx1_ka",100,0,1);
  TH1F *hip = new TH1F("ip trk prim vtx","vtx1_ip",100,0,1000);
  TH2F *hkaip = new TH2F("ka vs ip trk prim vtx","vtx1_kaip",100,0,1000,100,0,1);
  TH1F *hkaMC = new TH1F("kaMC trk prim vtx","vtx1_kaMC",100,0,1);
  TH1F *hipMC = new TH1F("ipMC trk prim vtx","vtx1_ipMC",100,0,1000);
  TH2F *hkaipMC = new TH2F("kaMC vs ipMC trk prim vtx","vtx1_kaipMC",100,0,1000,100,0,1);

  
  TH1F *hrmax = new TH1F("rmax trk prim vtx","vtx1_rmax",100,0,100);
  TH1F *hipmax = new TH1F("ipmax trk prim vtx","vtx1_ipmax",100,0,100);
  TH2F *hrmaxipmax = new TH2F("rmax vs ipmax trk prim vtx","vtx1_rmaxipmax",100,0,100,100,0,100);
  TH1F *hrmaxsame = new TH1F("same rmax trk prim vtx","vtx1_rmaxsame",100,0,100);
  TH1F *hipmaxsame = new TH1F("same ipmax trk prim vtx","vtx1_ipmaxsame",100,0,100);
  TH2F *hrmaxipmaxsame = new TH2F("same rmax vs ipmax trk prim vtx","vtx1_rmaxipmaxsame",100,0,100,100,0,100);
  TH1F *hsameseg = new TH1F("same seg ipmax rmax","same_seg",10,-1,2);
  TH1F *hrmaxMC = new TH1F("rmaxMC trk prim vtx","vtx1_rmaxMC",100,0,100);
  TH1F *hipmaxMC = new TH1F("ipmaxMC trk prim vtx","vtx1_ipmaxMC",100,0,100);
  TH2F *hrmaxipmaxMC = new TH2F("rmaxMC vs ipmaxMC trk prim vtx","vtx1_rmaxipmaxMC",100,0,100,100,0,100);
  TH1F *hrmaxsameMC = new TH1F("same rmaxMC trk prim vtx","vtx1_rmaxMCsame",100,0,100);
  TH1F *hipmaxsameMC = new TH1F("same ipmaxMC trk prim vtx","vtx1_ipmaxMCsame",100,0,100);
  TH2F *hrmaxipmaxsameMC = new TH2F("same rmaxMC vs ipmaxMC trk prim vtx","vtx1_rmaxipmaxMCsame",100,0,100,100,0,100);
  TH1F *hsamesegMC = new TH1F("same seg ipmaxMC rmaxMC","same_segMC",10,-1,2);
  
   std::vector<int> vtx_good_vec(nentries_vtx_fedra);
   std::vector<int> vtx_mp_eventID(nentries_vtx_fedra);
   std::vector<int> vtx_mp_motherID(nentries_vtx_fedra);

   // VERTEX STUDY

   //int first_entry=0;
   //int last_entry=nentries_vtx_fedra;
   
   for (Long64_t ientry=0; ientry<nentries_vtx_fedra;ientry++) { // Vertex index

     if(ientry>=0 /*|| ientry==1183*/){ // for debugging
     
     if(ientry%100==0)cout << "VID " << ientry << endl;
     log_vtx << "\nFedra Vertex number "<< ientry << endl;
     log_vtx << "MC_EventID / MCtrkID / MCMotherID / MCMotherpdgCode / MCpdgCode / MCProcessId" << endl;
     vtx_fedratree->GetEntry(ientry);
     vtx_good_vec[ientry]=-1;
     vtx_mp_eventID[ientry]=-1;
     vtx_mp_motherID[ientry]=-1;
     
     if(vx>xmin && vx<xmax && vy>ymin && vy<ymax && vz>zmin && vz<zmax ){

     vtx_trackmother.clear();
     vtx_trackprocess.clear();
     vtx_trackpdgcode.clear();
     vtx_motherpdgcode.clear();
     vtx_trackid.clear();
     vtx_trackstartX.clear();
     vtx_trackstartY.clear();
     vtx_trackstartZ.clear();
     vtx_trackeventId.clear();
    


     int ix = floor((vx-xmin)/xrange);
     int iy = floor((vy-ymin)/yrange);
     int iz = floor((vz-zmin)/zrange);
     //cout << vx << " " << vy <<" " << vz<< endl;
     //cout << ix << " " << iy << " " << iz << endl;
     
     vtx_list.at(ix).at(iy).at(iz).push_back(vID);
     //    if(n<=5)vtx_list_nearby.at(ix).at(iy).at(iz).push_back(vID);
     //vtx_ntrk.at(ix).at(iy).at(iz).push_back(n);

     //cout <<  "size " << vtx_list.at(ix).at(iy).at(iz).size() << endl;
     
     for (int itrk=0; itrk<n;itrk++) {  // Track index
       
       //vtx_trackeventId.push_back(MCEventID[itrk]);
       simtree->GetEntry(MCEventID[itrk]);
       //cout << "tracks " << ientry << " " << n << " " << MCTrackID[itrk] <<  endl;
       
       
       //std::vector<int> vtx_motherid_event(MCTrack_);
       //cout << MCTrack_ << endl;
       
       for (Long64_t jn=0; jn<MCTrack_;jn++) {
	 
	 if(jn==MCTrackID[itrk]){
	   vtx_trackeventId.push_back(MCEventID[itrk]);
	   vtx_trackpdgcode.push_back(MCTrack_fPdgCode[jn]);
	   vtx_trackmother.push_back(MCTrack_fMotherId[jn]);
	   vtx_trackprocess.push_back(MCTrack_fProcID[jn]);
	   vtx_trackid.push_back(jn);
	   vtx_trackstartX.push_back(MCTrack_fStartX[jn]);
	   vtx_trackstartY.push_back(MCTrack_fStartY[jn]);
	   vtx_trackstartZ.push_back(MCTrack_fStartZ[jn]);
	 }
       }
     }
     
     
     for(int h=0; h<vtx_trackmother.size(); h++){
       int tmp_jn=-1;
       if(vtx_trackmother.at(h) == -1) vtx_motherpdgcode.push_back(2212);
       else{
	 bool found_mother=false;
	 simtree->GetEntry(vtx_trackeventId.at(h));
	 for (Long64_t jn=0; jn<MCTrack_;jn++) {
	   if(vtx_trackmother.at(h) == jn && !found_mother){
	     found_mother=true;
	     tmp_jn = jn;
	     vtx_motherpdgcode.push_back(MCTrack_fPdgCode[jn]);
	   }
	 }
       }
       //cout <<"LOG "<< ientry << " " << vtx_trackeventId.at(h)  << " " << h << " " << vtx_trackid.at(h) << " " << vtx_trackmother.at(h) << " " << vtx_trackpdgcode.at(h) << " " << vtx_motherpdgcode.at(h) << " " << tmp_jn << endl;
       if(abs(vtx_motherpdgcode.at(h))==411)log_vtx << "Figlio del D+ " << endl;
       if(abs(vtx_motherpdgcode.at(h))==421)log_vtx << "Figlio del D0 " << endl;
       if(abs(vtx_motherpdgcode.at(h))==431)log_vtx << "Figlio del Ds " << endl;
       if(abs(vtx_motherpdgcode.at(h))==441)log_vtx << "Figlio dell'Eta_c " << endl;
       if(abs(vtx_motherpdgcode.at(h))==4122)log_vtx << "Figlio del Lambda_c+ " << endl;
       if(abs(vtx_motherpdgcode.at(h))==4132)log_vtx << "Figlio della Xsi_c0 " << endl;
       if(abs(vtx_motherpdgcode.at(h))==4232)log_vtx << "Figlio della Xsi_c+ " << endl;
       if(abs(vtx_motherpdgcode.at(h))==4332)log_vtx << "Figlio dell'Omega_c0 " << endl;
       log_vtx << vtx_trackeventId.at(h) <<  " "  << vtx_trackid.at(h) << " " << vtx_trackmother.at(h) << " " << vtx_trackpdgcode.at(h) << " " << vtx_motherpdgcode.at(h) << endl;
       
     }
     
     
     std::vector<int> vtx_ntracks_same_mother(n);
     //vtx_ntracks_same_mother.clear();
     std::vector<int> vtx_ntracks_same_event(n);
     //vtx_ntracks_same_event.clear();
     
     /// CONTA LE OCCORRENZE PER MOTHER
     for(int h=0; h<vtx_trackmother.size(); h++){
       flag_trk=0;
       for (int i=0;i<h;i++){
	 if (vtx_trackmother[h]==vtx_trackmother[i]){// controlla che quel numero non sia già stato considerato in precedenza
	   flag_trk=1;
	   break;
	 }
	 //cout << "before "<< ientry << " " << h << " " << i << " " << flag << " " << vtx_trackmother[h] << " " << vtx_ntracks_same_mother[h] << endl;
       }
	 if (flag_trk==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	   vtx_ntracks_same_mother[h]=1; // assegna 1 alla mother della prima traccia che incontra
	   for (int i=h+1;i<vtx_trackmother.size();i++){		 
	     if(vtx_trackmother[h] == vtx_trackmother[i] && vtx_trackid[h]!=vtx_trackid[i] && vtx_trackeventId[h]==vtx_trackeventId[i]) vtx_ntracks_same_mother[h]++;
	     //  cout << "after "<< ientry << " " <<  h << " " << i << " " << flag << " " <<" " << vtx_trackmother[i] << " " <<  vtx_trackmother[h] << " " << vtx_ntracks_same_mother[h] << endl;
	   }
	   //cout << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	   //log << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	   }
	 //cout << "ntracks "<< ientry << " " <<  h << " " <<  vtx_trackmother[h] << endl;
     }

     
     /// CONTA LE OCCORRENZE PER EVENTO
     for(int h=0; h<vtx_trackmother.size(); h++){
       flag_ev=0;
       for (int i=0;i<h;i++){
	 if (vtx_trackeventId[h]==vtx_trackeventId[i]){// controlla che quel numero non sia già stato considerato in precedenza
	   flag_ev=1;
	   break;
	 }
       }
       if (flag_ev==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	 vtx_ntracks_same_event[h]=1;
	 for (int i=h+1;i<vtx_trackmother.size();i++){		 
	   if(vtx_trackeventId[h] == vtx_trackeventId[i] && vtx_trackid[h]!=vtx_trackid[i]) vtx_ntracks_same_event[h]++;
	 }
	 //cout << ientry << " " << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackid[h] << "trk Frequency " << vtx_ntracks_same_event[h] << " times" << endl;
       }
     }

     
     
     vector<int>::iterator max_mID = max_element(vtx_ntracks_same_mother.begin(), vtx_ntracks_same_mother.end()); // c++11
     auto max_occurrence = *max_element(vtx_ntracks_same_mother.begin(), vtx_ntracks_same_mother.end()); // c++11    
     int iel = distance(vtx_ntracks_same_mother.begin(), max_mID);
     
     vector<int>::iterator max_evID = max_element(vtx_ntracks_same_event.begin(), vtx_ntracks_same_event.end()); // c++11
     auto max_ev_occurrence = *max_element(vtx_ntracks_same_event.begin(), vtx_ntracks_same_event.end()); // c++11
     int iel_ev = distance(vtx_ntracks_same_event.begin(), max_evID);
     
     vector <float> vtx_x;
     vector <float> vtx_y;
     vector <float> vtx_z;
     vtx_x.clear();
     vtx_y.clear();
     vtx_z.clear();

     int max_motherId = *max_mID;
     //auto max_eventId = *max_evID;

     //cout << "max occ "<<ientry << " " << iel << " " << max_occurrence << endl;

     n_electrons=0;
     //daughter_proton=false;
     
     if(max_ev_occurrence<=1){
       mp_eventID=-1;
       //primary_proton=false;
     }
     else {
       //cout << "iel_ev " << max_occurrence << " " << iel_ev << endl;
       mp_eventID = vtx_trackeventId.at(iel_ev);
       //cout << ientry << " " << mp_eventID << endl;
       //primary_proton=true;
       for(int i=0; i<vtx_trackpdgcode.size(); i++){
	 //if(abs(vtx_trackpdgcode.at(i))==2212 && vtx_motherpdgcode.at(i)==-1)daughter_proton=true;
	 if(abs(vtx_trackpdgcode.at(i))==11 && vtx_trackeventId.at(i)==mp_eventID)n_electrons++;
	 }
     }
     
     if(max_occurrence==1){
       mp_motherID=-2;
       mp_procID=-2;
       mp_pdgID=-2;
       mp_vx=0;
       mp_vy=0;
       mp_vz=0;
       mp_vdx=-1;
       mp_vdy=-1;
       mp_vdz=-1;
     }
     if(max_occurrence>1){
              
       int trk_index=0;

       
       for(unsigned int itrk=0;itrk<vtx_trackmother.size();itrk++){
	 //cout <<"crash "<< ientry << " " << vtx_trackmother.size() << " " <<  max_occurrence << " " << itrk << " " << iel  << endl;
	 
	 if(vtx_trackmother.at(itrk)==vtx_trackmother.at(iel)){
	   vtx_x.push_back(vtx_trackstartX.at(itrk));
	   vtx_y.push_back(vtx_trackstartY.at(itrk));
	   vtx_z.push_back(vtx_trackstartZ.at(itrk));	   
	   //cout << ientry << " " << trk_index << " " << vtx_x.at(trk_index) << " " << vtx_y.at(trk_index) << " " << vtx_z.at(trk_index) << endl;
	   trk_index++;
	 }
	 
       }
       
       mp_vx = accumulate(vtx_x.begin(), vtx_x.end(), 0.0)/vtx_x.size();
       mp_vy = accumulate(vtx_y.begin(), vtx_y.end(), 0.0)/vtx_y.size();
       mp_vz = accumulate(vtx_z.begin(), vtx_z.end(), 0.0)/vtx_z.size();
      
       mp_vx = mp_vx*pow(10,4) + 62500;
       mp_vy = mp_vy*pow(10,4) + 49500;
       mp_vz = (mp_vz-125.6)*pow(10,4);       

       auto max_vtx_x = *max_element(vtx_x.begin(), vtx_x.end()); 
       auto max_vtx_y = *max_element(vtx_y.begin(), vtx_y.end()); 
       auto max_vtx_z = *max_element(vtx_z.begin(), vtx_z.end());
       auto min_vtx_x = *min_element(vtx_x.begin(), vtx_x.end()); 
       auto min_vtx_y = *min_element(vtx_y.begin(), vtx_y.end()); 
       auto min_vtx_z = *min_element(vtx_z.begin(), vtx_z.end());        

       mp_vdx = max_vtx_x - min_vtx_x;
       mp_vdy = max_vtx_y - min_vtx_y;
       mp_vdz = max_vtx_z - min_vtx_z;
                    
       mp_motherID = vtx_trackmother.at(iel);
       mp_procID = vtx_trackprocess.at(iel);
       mp_pdgID = vtx_motherpdgcode.at(iel);

       //cout << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
       //cout << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
        
       //cout << ientry << " " << vtx_trackmother.at(iel) << " " << endl;
       
     }
     mean_freq=(1.0*max_occurrence)/n;
     max_freq=max_occurrence;

     mean_ev_freq=(1.0*max_ev_occurrence)/n;
     max_ev_freq=max_ev_occurrence;

     log_vtx << "The max number of tracks with same mother ID "<< mp_motherID << " in the Vertex " << ientry << " is " << max_occurrence << endl;
     log_vtx << "The max number of tracks with same event ID "<< mp_eventID << " in the Vertex " << ientry << " is " << max_ev_occurrence << endl;


     //------ DESCRIPTION  -----------//
     
  
       if(mp_pdgID!=-2){
	 if(mp_motherID==-1){
	   log_vtx << "It's a primary proton vertex with pdg code " << mp_pdgID << endl;
	   vtx_good=true;
	   vtx_good_vec[ientry]=0;
	   vtx_mp_eventID[ientry]=mp_eventID;
	   vtx_mp_motherID[ientry]=mp_motherID;
	 }
	 if(mp_motherID==0){
	   log_vtx << "It's a charm vertex with pdg code " << mp_pdgID << endl;
	   vtx_good=true;
	   vtx_good_vec[ientry]=1;
	   vtx_mp_eventID[ientry]=mp_eventID;
	   vtx_mp_motherID[ientry]=mp_motherID;
	 }
	 if( mp_motherID>0){
	   log_vtx << "It's a  secondary vertex with pdg code " << mp_pdgID << endl;
	   vtx_good=true;
	   vtx_good_vec[ientry]=2;
	   vtx_mp_eventID[ientry]=mp_eventID;
	   vtx_mp_motherID[ientry]=mp_motherID;
	 }
	 
       }
     else {
       log_vtx << "It's a fake\n";
       vtx_good=false;
     }

     //------------------------- //
     
       //if(daughter_proton)log << "Warning! A primary proton is a daughter!\n";
     //if(mp_motherID==0 && mp_pdgID==2212)log << "It's a primary proton\n";
     //if(mp_motherID>0 && mp_pdgID==2212)log << "It's a secondary proton\n";
     //if(mp_motherID==-1 && mp_pdgID==2212)log << "Warning! A primary proton is a daughter!\n";

     log_vtx << "VID / MostProbable Vertex X (Y,Z) from MC / Vertex X (Y,Z) from Fedra\n";
     log_vtx << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
     log_vtx << "Gap between MC and Fedra Vertex position (X,Y,Z)\n";
     log_vtx << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
     
     hoccurrence->Fill(mean_freq);
     //Tree_out->Fill();
     
     //cout << vx << " " << vy << " " << vz << endl;
     hvertex->Fill(vx,vy,vz);
     }
   } // only for special vertex
     }
   log_vtx << "VERTEXING LOG TERMINED" << endl;
   
   
   TFile * f_ds = new TFile("ds_data_result.root","RECREATE");
   TTree *TreeDS = new TTree("ds","decay search");
   TreeDS->Branch("vtx.event",&vtx.event,"vtx.event/I");
   TreeDS->Branch("vtx.fe_id",&vtx.fe_id,"vtx.fe_id/I");
   TreeDS->Branch("vtx.mc_id",&vtx.mc_id,"vtx.mc_id/I");
   TreeDS->Branch("vtx.ntrk",&vtx.ntrk,"vtx.ntrk/I");
   TreeDS->Branch("vtx.mc_nd",&vtx.mc_nd,"vtx.mc_nd[2]/I");
   TreeDS->Branch("vtx.mc_dz",&vtx.mc_dz,"vtx.mc_dz[2]/D");
   TreeDS->Branch("vtx.mc_dl",&vtx.mc_dl,"vtx.mc_dl[2]/D");
   TreeDS->Branch("vtx.x",&vtx.x,"vtx.x/D");
   TreeDS->Branch("vtx.y",&vtx.y,"vtx.y/D");
   TreeDS->Branch("vtx.z",&vtx.z,"vtx.z/D");
   TreeDS->Branch("trk.fe_itrk",&trk.fe_itrk);
   TreeDS->Branch("trk.icharm",&trk.icharm);
   TreeDS->Branch("trk.fe_id",&trk.fe_id);
   TreeDS->Branch("trk.mc_id",&trk.mc_id);
   TreeDS->Branch("trk.pdg",&trk.pdg);
   TreeDS->Branch("trk.nseg",&trk.nseg);
   TreeDS->Branch("trk.type",&trk.type);
   TreeDS->Branch("trk.fe_dz",&trk.fe_dz);
   TreeDS->Branch("trk.fe_dl",&trk.fe_dl);
   TreeDS->Branch("trk.ip",&trk.ip);
   TreeDS->Branch("trk.ka",&trk.ka);
   TreeDS->Branch("trk.x",&trk.x);
   TreeDS->Branch("trk.y",&trk.y);
   TreeDS->Branch("trk.z",&trk.z);
   TreeDS->Branch("trk.min_xy",&trk.min_xy);
   TreeDS->Branch("trk.min_xyz",&trk.min_xyz);
   TreeDS->Branch("trk.frd_xy",&trk.frd_xy);
   TreeDS->Branch("trk.frd_xyz",&trk.frd_xyz);
   
  
   
   // DECAY SEARCH

   int first_entry=0;//26000;
   int last_entry=nentries_vtx_fedra;
   int event = 0;

   
   for (Long64_t ientry=first_entry; ientry<last_entry;ientry++) {

     // Vertex index
     
     if(ientry%1==0)cout << "VID DS " << ientry << endl;    

       Set0();
       int nfound=0;
       
       int tmp_i=-1;
       int tmp_j=-1;
       int tmp_k=-1;
       int tmp_l=-1;
       int tmp_i_p1=-1;
       int tmp_j_p1=-1;
       int tmp_k_p1=-1;
       int tmp_i_m1=-1;
       int tmp_j_m1=-1;
       int tmp_k_m1=-1;
       int tmp_vID=-1;
       double tmp_vx=-1;
       double tmp_vy=-1;
       double tmp_vz=-1;

       double mc_tx=0;
       double mc_ty=0;

       TVector3 vtx_pos;
       TVector3 vtx_pos_2nd;
       TVector3 vtx_pos_nearby;
       TVector3 trk_pos;

       vector<int> trk_mc_id;
       vector<int> trk_z;
       trk_mc_id.clear();
       trk_z.clear();

       bool splitted_track=false;
       bool found_secondary=false;
       
       double fdz = 0;
       double fdl = 0;
       double ip = 0;
       double ka = 0;
       int ty = 0;
       int mc_nd1=0;
       int mc_nd2=0;
       int fnd=0;
       int icharm=0;

       vector <int> nvx;
       vector <int> nvy;
       vector <int> nvz;

       //float *dtheta;
       //float *ipseg;
       bool no_secondary_search=false;
       
       vtx_fedratree->GetEntry(ientry);
       //cout << n << endl;
       // RICERCA DEL CUBO IESIMO DEL VERTICE IESIMO
       
       int mean_vtx_seg = s__/n; // mi assicuro che non ci siano troppe tracce con pochi segmenti
       //cout << "mean-Seg "<<mean_vtx_seg << " " << s_ << " " << n << endl;
       // cerco solo nei vertici buoni (bdt >0) con più di 10 tracce (probabili primary proton vertexes)
       if(n>10 && probability>0.8 /*&& mean_vtx_seg>=8*/ /*&& vtx_bdt[ientry]>0*/){ // da scommentare solo se il file in input è stato testato con la bdt
	 //cout <<"vtx pos " <<  vx << " " << vy << " " << vz << endl;
	 for(int i=0;i<width;i++){
	 for(int j=0;j<height;j++){
	   for(int k=0;k<depth;k++){
	     //cout << "size " << vtx_list.at(i).at(j).at(k).size() << endl;
	     if(vtx_list.at(i).at(j).at(k).size()>0){
	       for(int l=0;l<vtx_list.at(i).at(j).at(k).size();l++){	       
		 //cout << i << " " << j << " " << k << " " << l <<  " " << vtx_list.at(i).at(j).at(k).at(l) <<  endl;
		 //log2 << i << " " << j << " " << k << " " << l <<  " " << vtx_list.at(i).at(j).at(k).at(l) << " " << "0" <<  endl;
		 if(vtx_list.at(i).at(j).at(k).at(l)==vID){ // trovo il cubetto di appartenenza del vertice
		   log_decay << "\nVtx " << vtx_list.at(i).at(j).at(k).at(l) << " found in coord " << i << " " << j << " " << k << " " << l << "  with Fedra xyz " << vx << " " << vy << " " << vz << " and ntracks " << n <<  endl;

		  
		   //cout <<"vtx pos " <<  vx << " " << vy << " " << vz << endl;
		   vtx_pos.SetXYZ(vx,vy,vz);
		   tmp_vID=vID;
		   tmp_vx=vx;
		   tmp_vy=vy;
		   tmp_vz=vz;
		   tmp_i=i;
		   tmp_j=j;
		   tmp_k=k;
		   tmp_l=l;
		   vtx.fe_id=vID;
		   vtx.event=event;
		   vtx.ntrk=n;
		   vtx.x=vx;
		   vtx.y=vy;
		   vtx.z=vz;


		    ////
		   int mpevent = vtx_mp_eventID.at(ientry);
		   vector<int> vtx_same_ID;
		   for(int j=0;j<fe_trk_id[mpevent].size();j++){
		     if(mc_chr_id[mpevent][j]==1 || mc_chr_id[mpevent][j]==2){
		       //cout <<"vtx_id " <<  fe_vtx_trk[mpevent][j] << " " << tmp_vID << endl;
		       if(fe_vtx_trk[mpevent][j]==tmp_vID){
			 //log_decay << "MC charm daugther in primary vertex with TrackID " << fe_trk_id[mpevent][j] << "and " << nseg[itrk] << " segments" << endl;
			 for(int itrk=0;itrk<n;itrk++){
			   //cout << itrk << endl;
			   vtx_fedratree->GetEntry(ientry);
			   if(TrackID[itrk]==fe_trk_id[mpevent][j]){
			     log_decay << "MC charm daugther in primary vertex with TrackID " << fe_trk_id[mpevent][j] << " and " << nseg[itrk] << " segments" << endl;
			      hnseg1st->Fill(nseg[itrk]);
			   }
			 }
			 mpevent++;
		       }
		       else{
			 int tmp_index = vtx_fedratree->GetEntryNumberWithIndex(fe_vtx_trk[mpevent][j]); // indice sul vertice secondario
			 vtx_fedratree->GetEntry(tmp_index);
			 float deltaz = vz-tmp_vz;
			 vtx_same_ID.push_back(fe_vtx_trk[mpevent][j]);
			 log_decay << "MC charm daugther not in primary vertex with TrackID " << fe_trk_id[mpevent][j] << " and linked to the secondary vertex ID " << fe_vtx_trk[mpevent][j] << " at a distance from the primary vtx " << deltaz << " with ntrack " << n << " and flag ID " << flag << endl;
			 if(count(vtx_same_ID.begin(), vtx_same_ID.end(), fe_vtx_trk[mpevent][j])>1) log_decay << "Secondary vertex with at least two daughters with ID " << fe_vtx_trk[mpevent][j] << endl;
			 for(int itrk=0;itrk<n;itrk++){
			   //cout << itrk << endl;
			   if(TrackID[itrk]==fe_trk_id[mpevent][j]){
			     log_decay << "2vtx Track " << TrackID[itrk] << " has " << nseg[itrk] << " segments" << endl;
			     hnseg->Fill(nseg[itrk]);

			     trk_pos.SetXYZ(t__eX[itrk],t__eY[itrk],t__eZ[itrk]);
			     fdz = t__eZ[itrk] - vz;
			     ip = IPtoVertex(vtx_pos, trk_pos, t__eTX[itrk], t__eTY[itrk]);
			     mc_tx = (t__eX[itrk]- vx)/fdz;
			     mc_ty = (t__eY[itrk]- vy)/fdz;
			     ka = sqrt(pow(mc_tx - t__eTX[itrk],2) + pow(mc_ty - t__eTY[itrk],2));
			     log_decay << "2vtx mc info - Ip to primary vertex " << ip << " and kink angle " << ka << endl;
			     log_decay << "2vtx mc info - TX " << t__eTX[itrk] << " and TY " << t__eTY[itrk] << endl;

			     hka2->Fill(ka);
			     hip2->Fill(ip);
			     hkaip2->Fill(ip,ka);
			     
			   }
			 }			 
			 hdeltaz->Fill(deltaz);
			 hntrk->Fill(n);
			 hflag->Fill(flag);
			 found_secondary=true;
		       }
		     }
		     //cout << i << " " << j << " " << fe_trk_id[mpevent][j] << " " << mc_chr_id[mpevent][j] <<  endl;
		   }
		   /*if(found_secondary){
		     vector<int>::iterator max_vID = max_element(vtx_same_ID.begin(), vtx_same_ID.end()); // c++11
		     auto max_vtx_occurrence = *max_element(vtx_same_ID.begin(), vtx_same_ID.end()); // c++11    
		     int i_vtx = distance(vtx_same_ID.begin(), max_vID);
		     cout << "ivtx " << i_vtx << " " << max_vID << endl;
		     if(max_vtx_occurrence>1)log_decay << "Secondary vertex with at least two daughters with ID " << vtx_same_ID.at(i_vtx) << endl;
		       }*/
		   ////

		   vtx_fedratree->GetEntry(ientry);
		   
		   float dtheta[2]={};
		   float ipseg[2]={};
		   float dthetaMC[2]={};
		   float ipsegMC[2]={};

		   int sum_seg=0;
		   int sum_segMC=0;

		   
		   //cout << "t__" << t__ << endl;

		   //cout << tmp_vz << " " << tmp_vx << " " << tmp_vy << endl;

		   vector <double> xpos;
		   vector <double> ypos;
		   vector <double> zpos;
		   vector <double> tx;
		   vector <double> ty;
		   vector <int> trk_good_id;

		   double ip_vec[n];
		   double ka_vec[n];

		   vector <double> ip_vec2;
		   vector <double> ka_vec2;

		   // loop su tracce del vertice in esame
		   for(int itrk=0;itrk<t__;itrk++){		    

		     double seg_x[nseg[itrk]];
		     double seg_y[nseg[itrk]];
		     double seg_z[nseg[itrk]];
		     double seg_tx[nseg[itrk]];
		     double seg_ty[nseg[itrk]];
		     double seg_mc_trk[nseg[itrk]];
		     double seg_mc_evt[nseg[itrk]];

		     bool found_daughter= false;

		     //xpos.push_back(t__eX[itrk]);
		     //ypos.push_back(t__eY[itrk]);
		     //zpos.push_back(t__eZ[itrk]);

		     trk_pos.SetXYZ(t__eX[itrk],t__eY[itrk],t__eZ[itrk]);
		     fdz = t__eZ[itrk] - tmp_vz;
		     ip = IPtoVertex(vtx_pos, trk_pos, t__eTX[itrk], t__eTY[itrk]);
		     mc_tx = (t__eX[itrk]- tmp_vx)/fdz;
		     mc_ty = (t__eY[itrk]- tmp_vy)/fdz;
		     ka = sqrt(pow(mc_tx - t__eTX[itrk],2) + pow(mc_ty - t__eTY[itrk],2));
		     ip_vec[itrk]=ip;
		     ka_vec[itrk]=ka;

		     //cout << itrk << " " << TrackID[itrk] << " " << ka << " " << ip << " " << fdz << " " << t__eX[itrk] << " " << t__eY[itrk] <<  endl;

		     // MC INFO
		     
		     for(int j=0;j<fe_trk_id[mpevent].size();j++){

		       //if(j==0) cout  << "Primary vertex track info " << endl;

		       if(fe_trk_id[mpevent][j]==TrackID[itrk]){
			
			 
			 if(mc_chr_id[mpevent][j]==1 || mc_chr_id[mpevent][j]==2){			

			   double segMC_x[nseg[itrk]];
			   double segMC_y[nseg[itrk]];
			   double segMC_z[nseg[itrk]];
			   double segMC_tx[nseg[itrk]];
			   double segMC_ty[nseg[itrk]];
			   double segMC_mc_trk[nseg[itrk]];
			   double segMC_mc_evt[nseg[itrk]];

			   //cout <<"mctrack " << fe_trk_id[mpevent][j] << " " << mc_chr_id[mpevent][j] << endl;
			   //cout << itrk << " " << TrackID[itrk] << " " << vID << " " << sum_seg+nseg[itrk] << endl;
			   
			   for(int u=sum_seg;u<(sum_seg+nseg[itrk]);u++){
			     //cout << u << " " << itrk << " " << TrackID[itrk] << " " << s__eTrack[u] << " " <<  vID << endl;
			     if(s__eTrack[u]==TrackID[itrk]){
			       //cout << "segment mctrack" << endl;
			       //cout << "seg " << u << u-sum_seg << " " << itrk << " " << s__eTrack[u] << endl;
			       segMC_x[u-sum_seg]=s__eX[u];
			       segMC_y[u-sum_seg]=s__eY[u];
			       segMC_z[u-sum_seg]=s__eZ[u];
			       segMC_tx[u-sum_seg]=s__eTX[u];
			       segMC_ty[u-sum_seg]=s__eTY[u];
			       }
			   }
			   
			   // calcolo di rmax e ipmax per trovare eventuali 1 prong
			   FedraTrackKink2(nseg[itrk],segMC_tx,segMC_ty,dthetaMC);
			   SegIPtoVertex2(vtx_pos, nseg[itrk], segMC_x, segMC_tx, segMC_y, segMC_ty, segMC_z, ipsegMC);
			   
			   hkaMC->Fill(ka);
			   hipMC->Fill(ip);
			   hkaipMC->Fill(ip,ka);

			   hrmaxMC->Fill(dthetaMC[1]);
			   hipmaxMC->Fill(ipsegMC[1]);
			   hrmaxipmaxMC->Fill(ipsegMC[1],dthetaMC[1]);
			   if(dthetaMC[0]==ipsegMC[0]){
			     hsamesegMC->Fill(1);
			     hrmaxsameMC->Fill(dthetaMC[1]);
			     hipmaxsameMC->Fill(ipsegMC[1]);
			     hrmaxipmaxsameMC->Fill(ipsegMC[1],dthetaMC[1]);
			   }
			   else hsamesegMC->Fill(0);
			   
			   log_decay << "1vtx mc info - TrackID " << TrackID[itrk] << " with nseg " << nseg[itrk] << endl;
			   log_decay << "1vtx mc info - Rmax seg " << dthetaMC[0] << " and value " << dthetaMC[1] << endl;
			   log_decay << "1vtx mc info - Ipmax seg " << ipsegMC[0] << " and value " << ipsegMC[1] << endl;
		     	   log_decay << "1vtx mc info - Ip to primary vertex " << ip << " and kink angle " << ka << endl;
			   log_decay << "1vtx mc info - TX " << t__eTX[itrk] << " and TY " << t__eTY[itrk] << endl;

			   found_daughter=true;
			 }		        
		       }
		     }

		     ////
		     
		     
		     // loop su segmenti delle tracce del vertice in esame che hanno nseg >= 8
		     if(nseg[itrk]>=8){
		       for(int u=sum_seg;u<(sum_seg+nseg[itrk]);u++){ // perché i segmenti non tornano a zero al cambiare della nuova traccia
			 seg_x[u-sum_seg]=s__eX[u];
			 seg_y[u-sum_seg]=s__eY[u];
			 seg_z[u-sum_seg]=s__eZ[u];
			 seg_tx[u-sum_seg]=s__eTX[u];
			 seg_ty[u-sum_seg]=s__eTY[u];
			 //cout <<"a "<< i << " " << j << " " << k <<" " <<  l << " " << n << " " << itrk << " " << s__eTrack[u] << " " << nseg[itrk] <<  " " <<  (u-sum_seg)  << " " << seg_x[u-sum_seg] << " " <<  seg_y[u-sum_seg] << endl;
		       }
		       
		       // calcolo di rmax e ipmax per trovare eventuali 1 prong
		       FedraTrackKink2(nseg[itrk],seg_tx,seg_ty,dtheta);
		       SegIPtoVertex2(vtx_pos, nseg[itrk], seg_x, seg_tx, seg_y, seg_ty, seg_z, ipseg);

		       if(!found_daughter){
			 hrmax->Fill(dtheta[1]);
			 hipmax->Fill(ipseg[1]);
			 hrmaxipmax->Fill(ipseg[1],dtheta[1]);
			 
			 if(dtheta[0]==ipseg[0]){
			   hsameseg->Fill(1);
			   hrmaxsame->Fill(dtheta[1]);
			   hipmaxsame->Fill(ipseg[1]);
			   hrmaxipmaxsame->Fill(ipseg[1],dtheta[1]);
			 }
			 else hsameseg->Fill(0);
		       }
		       
		       // solo se è lo stesso segmento che ha ipmax>10 e rmax>10
		       //if(dtheta[0]==ipseg[0] && dtheta[1]>10 && ipseg[1]>10){
		       if(dtheta[1]>5 || ipseg[1]>5){
			 log_decay << "ds info Possible charm 1 prong with TrackID " << TrackID[itrk] <<  endl;
			 log_decay << "ds info Rmax seg " << dtheta[0] << " and value " << dtheta[1] << endl;
			 log_decay << "ds info Ipmax seg " << ipseg[0] << " and value " << ipseg[1] << endl;
		       }

		       if(!found_daughter){
			 hka->Fill(ka);
			 hip->Fill(ip);
			 hkaip->Fill(ip,ka);
		       }
		      
		       
		     }

		     //cout << itrk << " " << TrackID[itrk] << " " << ka_vec[itrk] << " " << ip_vec[itrk] << " " << fdz << " " << t__eTX[itrk] << " " << t__eTY[itrk] <<  " " <<  t__eX[itrk] << " " << t__eY[itrk] << t__eZ[itrk] << endl;
		     
		     // if(ip_vec[itrk]>5 && ka_vec[itrk]>0.01 && abs(t__eTX[itrk])>0.01 && abs(t__eTY[itrk])>0.01){
		     if(nseg[itrk]>=8){
		       xpos.push_back(t__eX[itrk]);
		       ypos.push_back(t__eY[itrk]);
		       zpos.push_back(t__eZ[itrk]);
		       tx.push_back(t__eTX[itrk]);
		       ty.push_back(t__eTY[itrk]);
		       ip_vec2.push_back(ip);
		       ka_vec2.push_back(ka);
		       trk_good_id.push_back(TrackID[itrk]);
		     }
		     
		       /*if(nseg[itrk]>=8){
			log_decay << "Interesting track linked to the primary vertex " << endl;
			log_decay << "1vtx ds info - TrackID " << TrackID[itrk] << " with nseg " << nseg[itrk] << endl;
			log_decay << "1vtx ds info - Ip to primary vertex " << ip_vec[itrk] << " and kink angle " << ka_vec[itrk] << endl;
			log_decay << "1vtx ds info - TX " << t__eTX[itrk] << " and TY " << t__eTY[itrk] << endl;
			}*/
		       // }
		      
		     sum_seg += nseg[itrk];
		   }

		   vector<vector<int>> cl_el;
		   vector<vector<int>> cl_el_rel;
		   vector<double> mean_dist;
		   vector<double> mean_dist_phi;
		   vector<double> mean_dist_vtx;
		   vector<double> mean_dist_vtx_rel;
		   vector<double> dist_rel;
		   vector<int> n_el;
		   vector<int> n_el_rel;
		   float min_dist=-1;
		   //cout << xpos.size() << endl;
		   if(xpos.size()>1){
		     
		     tie(cl_el,mean_dist,mean_dist_phi,mean_dist_vtx,n_el,n_el_rel,cl_el_rel,mean_dist_vtx_rel,dist_rel) = dbscan(epsilon,epsilon_phi,min_points,xpos.size(),xpos,ypos,zpos,tx,ty);
		     int cl_index=0;

		     for(int icl=0;icl<n_el.size();icl++){

		       float sum_ka=0;
		       float sum_ip=0;
		       bool not_new_2vtx=false;
		       if(n_el[icl]!=0){
			 
			 for(int p=0;p<icl;p++){
			   not_new_2vtx = compareVectors(cl_el[icl],cl_el[p]);
			   if(not_new_2vtx)break;
			 }

			 if(!not_new_2vtx){
			   for(int jcl=0;jcl<n_el[icl];jcl++){
			     
			     sum_ka += ka_vec2[cl_el[icl][jcl]];
			     sum_ip += ip_vec2[cl_el[icl][jcl]];
			     if(jcl==(n_el[icl]-1)){
			       sum_ka /= n_el[icl];
			       sum_ip /= n_el[icl];
			       
			     if(sum_ip>50 && sum_ka>0.01 && n_el[icl]<=5 && abs(mean_dist.at(icl))<0.1 && abs(mean_dist_phi.at(icl))<0.01 && abs(mean_dist_vtx.at(icl))>10){
			       log_decay << "Possible cluster " << icl << " with nTracks " << n_el[icl] << " mean cross dist " << mean_dist.at(icl) << " mean phi dist " << mean_dist_phi.at(icl) << " mean vtx dist " << mean_dist_vtx.at(icl) << endl;
			       for(int s=0;s<n_el[icl];s++){
				 log_decay << "Element " << s << " with Track ID " << trk_good_id[cl_el[icl][s]] << endl;
				 log_decay << "ip " << ip_vec2[cl_el[icl][s]] << " ka " << ka_vec2[cl_el[icl][s]] << " tx " << t__eTX[cl_el[icl][s]] << " ty " << t__eTY[cl_el[icl][s]] << endl;
			       }
			       log_decay << "Average ip " << sum_ip << " and ka " << sum_ka << endl;
			     }
			     }
			   }
		       }
		       }
		     }
		     
		     for(int icl=0;icl<n_el_rel.size();icl++){
		       float sum_ka=0;
		       float sum_ip=0;
		       bool not_good=false;
		       for(int jcl=0;jcl<n_el_rel[icl];jcl++){

			 sum_ka += ka_vec2[cl_el_rel[icl][jcl]];
			 sum_ip += ip_vec2[cl_el_rel[icl][jcl]];

			 //if(abs(t__eTX[cl_el_rel[icl][jcl]])<0.01 || t__eTY[cl_el_rel[icl][jcl]]<0.01)not_good=true;
			 
			 if(jcl==(n_el_rel[icl]-1)){
			   sum_ka /= n_el_rel[icl];
			   sum_ip /= n_el_rel[icl];
			   if(sum_ip>50 && sum_ka>0.01 && abs(mean_dist_vtx_rel.at(icl))>10 && !not_good){
			     log_decay << "Possible  2vtx with mean dist vtx " << mean_dist_vtx_rel.at(icl) << " and a relative distance first seg " << dist_rel.at(icl) <<  endl;
			     for(int s=0;s<n_el_rel[icl];s++){
			       log_decay << "Element " << s << " with Track ID " << trk_good_id[cl_el_rel[icl][s]] << endl;
			       log_decay << "ip " << ip_vec2[cl_el_rel[icl][s]] << " ka " << ka_vec2[cl_el_rel[icl][s]] << " tx " << t__eTX[cl_el_rel[icl][s]] << " ty " << t__eTY[cl_el_rel[icl][s]] << endl;
			     }
			     log_decay << "Average ip " << sum_ip << " and ka " << sum_ka << endl;
			   }
			 }
		       }
		     }
		   }
		   
		   event++;
		   break;
		 }
	       }
	     }
	     //cout << "size " << vtx_list.at(i).at(j).at(k).size() << endl;
	     if(vtx_list.at(i).at(j).at(k).size()==0) no_secondary_search=true;
	   }
	 }
       }

	 
	 
       // RICERCA SUI PRIMI VICINI   (otto cubetti intorno + il cubetto medesimo)

	 if(!no_secondary_search){ 
       tmp_i_m1 = tmp_i - 1;
       tmp_j_m1 = tmp_j - 1;
       tmp_k_m1 = tmp_k;// - 2;

       tmp_i_p1 = tmp_i + 2;
       tmp_j_p1 = tmp_j + 2;
       tmp_k_p1 = tmp_k + 1;// + 2;

       if(tmp_i_m1<0)tmp_i_m1=0;
       if(tmp_j_m1<0)tmp_j_m1=0;
       if(tmp_k_m1<0)tmp_k_m1=0;

       if(tmp_i_p1>vbinx)tmp_i_p1=vbinx;
       if(tmp_j_p1>vbiny)tmp_j_p1=vbiny;
       if(tmp_k_p1>vbinz)tmp_k_p1=vbinz;
       
       cout << tmp_i << " " << tmp_i_m1 << " " << tmp_i_p1 << endl;
       cout << tmp_j << " " << tmp_j_m1 << " " << tmp_j_p1 << endl;
       cout << tmp_k << " " << tmp_k_m1 << " " << tmp_k_p1 << endl;

       
	   // loop sui primi vicini
	   for(int i=tmp_i_m1;i<tmp_i_p1;i++){
	     for(int j=tmp_j_m1;j<tmp_j_p1;j++){
	       for(int k=tmp_k_m1;k<tmp_k_p1;k++){
		 //cout << "cubs " << i << " " << j << " " << k << endl;
	     
		 if(!vtx_list.at(i).at(j).at(k).empty()){
		   // cout <<"size "<< vtx_list.at(i).at(j).at(k).size() << endl;
		   // SEARCH IN VERTEX NEARBY
	     
		   for(int l=0;l< vtx_list.at(i).at(j).at(k).size();l++){
		     int tmp_index = vtx_fedratree->GetEntryNumberWithIndex(vtx_list.at(i).at(j).at(k).at(l)); // indice sul vertice secondario
		     vtx_fedratree->GetEntry(tmp_index);

		     //cout << "2nd vtx " << tmp_index << endl;
	       
		     vtx_pos_nearby.SetXYZ(vx,vy,vz);
		     fdz = vz - tmp_vz;

		     bool trk_found=false;


		     //cout << vID << endl;
	       
		     // condizioni sul vertice secondario 
		     if(fdz < 10000 && fdz > 500 && tmp_vID!=vID && (flag==4 || flag==0 || flag==3) && n<=5){
		 
		       // mi assicuro che le traccie abbiano almeno 8 seg escluso il parent se è attaccato al vertice
		       int low_seg=0;
		       for(int itrk=0;itrk<n;itrk++){		  
			 if(nseg[itrk]<6 && incoming[itrk]!=0)low_seg++;
		       }

		       // per i carichi ne posso avere solo 1 traccia (il parent) con un numero di seg minore di 8
		       if(low_seg==0){
			 // cerco solo vertici secondari linkati che hanno meno di 4 tracce
			 if(flag==4){
			   //log_decay <<"Nearby linked vertex found " << vtx_list.at(i).at(j).at(k).at(l) << " at distance " << fdz << " with " << n << " tracks " << endl; 
			   for(int itrk=0;itrk<n;itrk++){
			     trk_pos.SetXYZ(t__eX[itrk],t__eY[itrk],t__eZ[itrk]);
			     fdz = t__eZ[itrk] - tmp_vz;
			     ip = IPtoVertex(vtx_pos, trk_pos, t__eTX[itrk], t__eTY[itrk]);
			     mc_tx = (t__eX[itrk]-tmp_vx)/fdz;
			     mc_ty = (t__eY[itrk]-tmp_vy)/fdz;
			     ka = sqrt(pow(mc_tx - t__eTX[itrk],2) + pow(mc_ty - t__eTY[itrk],2));
			     if(nseg[itrk]>=8 && ip>50 && ip<1000 && fdz<10000 && fdz>500  && ka<0.2){ // stabilire se necessario
			       if(!trk_found){
				 log_decay << "DS info on secondary vertex " << endl;
				 log_decay <<"Nearby linked vertex found " << vtx_list.at(i).at(j).at(k).at(l) << " at distance " << (vz-tmp_vz) << " with " << n << " tracks " << endl;
				 trk_found=true;
			       } 
			       log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << TrackID[itrk] << " and " << nseg[itrk] << " segments" << " and zpos " << incoming[itrk] << endl;
			       log_decay << "track_pos (XYZ) "<< t__eX[itrk] << " " << t__eY[itrk] << " " << t__eZ[itrk]  << endl;
			       log_decay << "track_dir (TX - TY) " << t__eTX[itrk] << " " << t__eTY[itrk] << endl;
			       log_decay << "ip: " << ip << " - ka: " << ka << endl;
			     }
			   }
		     
			 }		   
		   
			 if((flag==0 || flag==3) && n<4){ // vertici non linkati (magari neutral charmed daughters)
			   //if(vertexobject_tmp->Zpos(itrk))
			   //log_decay <<"Nearby neutral vertex found " << vtx_list.at(i).at(j).at(k).at(l) << " at distance " << fdz << " with " << n << " tracks " << endl;
			   for(int itrk=0;itrk<n;itrk++){
		       
			     trk_pos.SetXYZ(t__eX[itrk],t__eY[itrk],t__eZ[itrk]);
			     fdz = t__eZ[itrk] - tmp_vz;
			     ip = IPtoVertex(vtx_pos, trk_pos, t__eTX[itrk], t__eTY[itrk]);
			     mc_tx = (t__eX[itrk]-tmp_vx)/fdz;
			     mc_ty = (t__eY[itrk]-tmp_vy)/fdz;
			     ka = sqrt(pow(mc_tx - t__eTX[itrk],2) + pow(mc_ty - t__eTY[itrk],2));

			     //if(TrackID[itrk]==24090)cout << ip << " " << ka << " " << t__eTX[itrk] << " " << t__eTY[itrk] << endl;
			     // condizioni per le tracce del vertice secondario
			     //if(nseg[itrk]>=8 && ip>50 && ip<400 && fdz<10000 && fdz>-100 && (t__eTX[itrk]>0.1 || t__eTY[itrk]>0.1) && ka<0.2){
			     if(nseg[itrk]>=8 && ip>50 && ip<1000 && fdz<10000 && fdz>500  && ka<0.2){
			       if(!trk_found){
				 log_decay << "DS info on secondary vertex " << endl;
				 log_decay <<"Nearby neutral vertex found " << vtx_list.at(i).at(j).at(k).at(l) << " at distance " << (vz-tmp_vz) << " with " << n << " tracks " << endl;
				 trk_found=true;
			       }
			       log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << TrackID[itrk] << " and " << nseg[itrk] << " segments" << " and zpos " << incoming[itrk] << endl;
			       log_decay << "track_pos (XYZ) "<< t__eX[itrk] << " " << t__eY[itrk] << " " << t__eZ[itrk]  << endl;
			       log_decay << "track_dir (TX - TY) " << t__eTX[itrk] << " " << t__eTY[itrk] << endl;
			       log_decay << "ip: " << ip << " - ka: " << ka << endl;

			       // assegno le coordinate del neutral nearby vertex (mi serviranno per la ricerca di tracce staccate che possano essere il parent, quindi un charged)
			       nvx.push_back(vx);
			       nvy.push_back(vy);
			       nvz.push_back(vz);
			     }
			   }
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	   }
       
	   TreeDS->Fill();
	 }
       }
   }
   
   TreeDS->Write();
   f_ds->Close();
   log_decay << "DECAY SEARCH LOG TERMINED" << endl;
   log_decay.close();


   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
   c1->Divide(3,2);
   c1->cd(1);
   hip->Draw("");
   c1->cd(2);
   hka->Draw("");
   c1->cd(3);
   hkaip->Draw("");
   c1->cd(4);
   hipMC->Draw("");
   c1->cd(5);
   hkaMC->Draw("");
   c1->cd(6);
   hkaipMC->Draw("");

   TCanvas *c1d = new TCanvas("c1d","c1d",500,500);
   hnseg1st->Draw("");

   TCanvas *c1a = new TCanvas("c1a","c1a",1000,500);
   c1a->Divide(2,1);
   c1a->cd(1);
   hsameseg->Draw("");
   c1a->cd(2);
   hsamesegMC->Draw("");
   

   TCanvas *c1b = new TCanvas("c1b","c1b",1500,1000);
   c1b->Divide(3,2);
   c1b->cd(1);
   hipmax->Draw("");
   c1b->cd(2);
   hrmax->Draw("");
   c1b->cd(3);
   hrmaxipmax->Draw("");
   c1b->cd(4);
   hipmaxMC->Draw("");
   c1b->cd(5);
   hrmaxMC->Draw("");
   c1b->cd(6);
   hrmaxipmaxMC->Draw("");

   TCanvas *c1c = new TCanvas("c1c","c1c",1500,1000);
   c1c->Divide(3,2);
   c1c->cd(1);
   hipmaxsame->Draw("");
   c1c->cd(2);
   hrmaxsame->Draw("");
   c1c->cd(3);
   hrmaxipmaxsame->Draw("");
   c1c->cd(4);
   hipmaxsameMC->Draw("");
   c1c->cd(5);
   hrmaxsameMC->Draw("");
   c1c->cd(6);
   hrmaxipmaxsameMC->Draw("");
   
   TCanvas *c2 = new TCanvas("c2","c2",1500,1000);
   c2->Divide(2,3);
   c2->cd(1);
   hdeltaz->Draw("");
   c2->cd(2);
   hntrk->Draw("");
   c2->cd(3);
   hnseg->Draw("");
   c2->cd(4);
   hflag->Draw("");

   TCanvas *c2a = new TCanvas("c2a","c2a",1500,500);
   c2a->Divide(3,1);
   c2a->cd(1);
   hka2->Draw("");
   c2a->cd(2);
   hip2->Draw("");
   c2a->cd(3);
   hkaip2->Draw("");
 
  
}

int myrun(){
  Loop();
  return 0;
}




     /*
     DBSCAN(D, epsilon, min_points):
      C = 0
      for each unvisited point P in dataset
            mark P as visited
            sphere_points = regionQuery(P, epsilon)
            if sizeof(sphere_points) < min_points
                  ignore P
            else
                  C = next cluster
                  expandCluster(P, sphere_points, C, epsilon, min_points)

expandCluster(P, sphere_points, C, epsilon, min_points):
      add P to cluster C
      for each point P’ in sphere_points
            if P’ is not visited
                  mark P’ as visited
                  sphere_points’ = regionQuery(P’, epsilon)
                  if sizeof(sphere_points’) >= min_points
                        sphere_points = sphere_points joined with sphere_points’
                  if P’ is not yet member of any cluster
                        add P’ to cluster C

regionQuery(P, epsilon):
      return all points within the n-dimensional sphere centered at P with radius epsilon (including P)

     */




