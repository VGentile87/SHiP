void reader_evaluation_vtx(){

   TFile fileIn("vtx_BDT_DS_evaluated.root");   
   TTreeReader theReaderBDT("bdt",&fileIn);
   TTreeReader theReaderVTX("vertices",&fileIn);
   TTreeReaderValue<int> Vtx_id(theReaderVTX, "vtx_id");
   TTreeReaderValue<int> Vtx2_id(theReaderVTX, "vtx2_id");
   TTreeReaderValue<int> Vtx_mc_ev(theReaderVTX, "vtx_mc_ev");
   TTreeReaderValue<int> Vtx_ntrk(theReaderVTX, "vtx_ntrk");
   TTreeReaderValue<int> zpos(theReaderVTX, "zpositive");
   TTreeReaderValue<int> samev(theReaderVTX, "samevent");
   TTreeReaderValue<int> charmdaug(theReaderVTX, "charmdaughter");
   TTreeReaderValue<int> Goodvtx(theReaderVTX, "goodvtx");
   TTreeReaderValue<float> Pointing(theReaderVTX, "pointing");
   TTreeReaderValue<float> imp(theReaderVTX, "impact");
   TTreeReaderValue<float> Kink(theReaderVTX, "kink");
   TTreeReaderValue<int> chr1(theReaderVTX, "charm1");
   TTreeReaderValue<int> chr2(theReaderVTX, "charm2");
   //TTreeReaderValue<int> Nseg(theReaderVTX, "nseg");
   TTreeReaderValue<float> Dist_tr(theReaderVTX, "dist_tr");
   TTreeReaderValue<float> Bdt(theReaderBDT, "bdt_value");
   TTreeReaderValue<float> Vtx_z(theReaderVTX, "vtx_z");
   TTreeReaderValue<float> Vtx2_z(theReaderVTX, "vtx2_z");
   // TTreeReaderValue<float> Trk_z(theReaderVTX, "trk_z");
   TTreeReaderValue<float> Decaylength(theReaderVTX, "decaylength");
   //TTreeReaderValue<float> Trk_fill(theReaderVTX, "trk_fill");
   TTreeReaderValue<float> Trk_pms(theReaderVTX, "trk_pms");
   TTreeReaderValue<float> Tka_rms(theReaderVTX, "tka_rms");


   ofstream log_dat("vtx_list.txt");
   
   // dare in input
   const int nevents = 10000;
   const int nvtx = 200000;
   
   TH1I *hch1 = new TH1I("hch1","hch1",nevents,0,nevents);
   TH1I *hch2 = new TH1I("hch2","hch2",nevents,0,nevents);
   TH1I *hch1_nocut = new TH1I("hch1_nocut","hch1_nocut",nevents,0,nevents);
   TH1I *hch2_nocut = new TH1I("hch2_nocut","hch2_nocut",nevents,0,nevents);
   TH1I *hch1_ipmax = new TH1I("hch1_ipmax","hch1_ipmax",nevents,0,nevents);
   TH1I *hch2_ipmax = new TH1I("hch2_ipmax","hch2_ipmax",nevents,0,nevents);
   TH1I *hch1_ipmin = new TH1I("hch1_ipmin","hch1_ipmin",nevents,0,nevents);
   TH1I *hch2_ipmin = new TH1I("hch2_ipmin","hch2_ipmin",nevents,0,nevents);
   TH1I *hch1_kink = new TH1I("hch1_kink","hch1_kink",nevents,0,nevents);
   TH1I *hch2_kink = new TH1I("hch2_kink","hch2_kink",nevents,0,nevents);
   TH1I *hch1_rms = new TH1I("hch1_rms","hch1_rms",nevents,0,nevents);
   TH1I *hch2_rms = new TH1I("hch2_rms","hch2_rms",nevents,0,nevents);
   TH1I *hch1_decaylength = new TH1I("hch1_decaylength","hch1_decaylength",nevents,0,nevents);
   TH1I *hch2_decaylength = new TH1I("hch2_decaylength","hch2_decaylength",nevents,0,nevents);
   TH1I *hch1_dxy = new TH1I("hch1_dxy","hch1_dxy",nevents,0,nevents);
   TH1I *hch2_dxy = new TH1I("hch2_dxy","hch2_dxy",nevents,0,nevents);
   TH1I *hch1_bdt = new TH1I("hch1_bdt","hch1_bdt",nevents,0,nevents);
   TH1I *hch2_bdt = new TH1I("hch2_bdt","hch2_bdt",nevents,0,nevents);
   TH1I *hch1_dzvtx = new TH1I("hch1_dzvtx","hch1_dzvtx",nevents,0,nevents);
   TH1I *hch2_dzvtx = new TH1I("hch2_dzvtx","hch2_dzvtx",nevents,0,nevents);
   TH1I *hch1_dztrk = new TH1I("hch1_dztrk","hch1_dztrk",nevents,0,nevents);
   TH1I *hch2_dztrk = new TH1I("hch2_dztrk","hch2_dztrk",nevents,0,nevents);
   TH1I *hch1_fill = new TH1I("hch1_fill","hch1_fill",nevents,0,nevents);
   TH1I *hch2_fill = new TH1I("hch2_fill","hch2_fill",nevents,0,nevents);
   TH1I *hch1_mom = new TH1I("hch1_mom","hch1_mom",nevents,0,nevents);
   TH1I *hch2_mom = new TH1I("hch2_mom","hch2_mom",nevents,0,nevents);
   TH1I *hch1_pointing = new TH1I("hch1_pointing","hch1_pointing",nevents,0,nevents);
   TH1I *hch2_pointing = new TH1I("hch2_pointing","hch2_pointing",nevents,0,nevents);
   
   TH1I *hsig = new TH1I("hsig","hsig",nevents,0,nevents);
   TH1I *hsig_nocut = new TH1I("hsig_nocut","hsig_nocut",nevents,0,nevents);
   TH1I *hbkg = new TH1I("hbkg","hbkg",nevents,0,nevents);
   TH1I *hdat = new TH1I("hdat","hdat",nevents,0,nevents);

   TH1I *hsig_vtx = new TH1I("hsig_vtx","hsig_vtx",nvtx,0,nvtx);
   TH1I *hbkg_vtx = new TH1I("hbkg_vtx","hbkg_vtx",nvtx,0,nvtx);
   TH1I *hdat_vtx = new TH1I("hdat_vtx","hdat_vtx",nvtx,0,nvtx);

   vector<vector<int>> couples_sig(nvtx);
   vector<vector<int>> couples_bkg(nvtx);
   vector<vector<int>> couples_bkg2(nvtx);


   int vtx_id, vtx_mc_ev, zpositive, samevent, charmdaughter, nseg, charm1, charm2, vtx2_id, goodvtx, vtx_ntrk;
   float impact, kink, dist_tr, bdt_value, vtx_z, vtx2_z, trk_z, decaylength, trk_fill, trk_pms, tka_rms, pointing;

   // counters
   int silver=0;
   int golden=0;
   int nsig=0;
   int silver_nocut=0;
   int golden_nocut=0;
   int nsig_nocut=0;
   int nbkg=0;
   int nsig_vtx=0;
   int nbkg_vtx=0;
   int mixed=0;
   int mixed2=0;
   int ndat=0;
   int ndat_vtx=0;

   int silver_ipmax=0;
   int silver_ipmin=0;
   int silver_kink=0;  
   int silver_rms=0;
   int silver_decaylength=0;
   int silver_dxy=0;  
   int silver_bdt=0;
   int silver_dzvtx=0;
   int silver_dztrk=0;
   int silver_fill=0;
   int silver_mom=0;
   int silver_pointing=0;

   int golden_ipmax=0;
   int golden_ipmin=0;
   int golden_kink=0;  
   int golden_rms=0;
   int golden_decaylength=0;
   int golden_dxy=0;  
   int golden_bdt=0;
   int golden_dzvtx=0;
   int golden_dztrk=0;
   int golden_fill=0;
   int golden_mom=0;
   int golden_pointing=0;
   
   //cuts
   float cut_dxy=2000;//200;
   float cut_bdt=0;
   float cut_ipmax=10000;//500;
   float cut_ipmin=0;//10
   float cut_decaylength=30000; //9000
   float cut_ka=0.010;//0.4; 
   float cut_dz=500;//500;
   float cut_dz_trk=0;
   float cut_trk_fill=0;
   float cut_momentum=0; //1.5
   float cut_rms=0.075;
   float cut_pointing=0.2;
   float cut_ctau =30;//*TMath::Power(10,-12);
   int cut_vtx1_ntrk =6;
   
   const Int_t ntrk = theReaderVTX.GetEntries();
   //cout << nvertices << endl;
   for (int ivtx=0;ivtx<ntrk;ivtx++){

     theReaderVTX.Next();
     theReaderBDT.Next();

     vtx_id = *Vtx_id;
     vtx2_id = *Vtx2_id;
     vtx_mc_ev = *Vtx_mc_ev;
     vtx_z = *Vtx_z;
     vtx2_z = *Vtx2_z;
     vtx_ntrk = *Vtx_ntrk;
     //cout << vtx_ntrk << endl;
     zpositive = *zpos;
     samevent = *samev;
     charmdaughter = *charmdaug;
     goodvtx = *Goodvtx;
     pointing = *Pointing;
     impact = *imp;
     kink = *Kink;
     charm1 = *chr1;
     charm2 = *chr2;
     //nseg = *Nseg;
     dist_tr = *Dist_tr;
     bdt_value = *Bdt;
     decaylength = *Decaylength;
     //trk_fill =*Trk_fill;
     trk_pms = *Trk_pms;
     tka_rms = *Tka_rms;
     // cout << ivtx << " " << vtx_id << endl;

     if(zpositive  && goodvtx && pointing!=10 && zpositive && samevent && charmdaughter && vtx_ntrk>=cut_vtx1_ntrk){
       hsig_nocut->Fill(vtx_mc_ev);
       if(charm1)hch1_nocut->Fill(vtx_mc_ev);
       if(charm2)hch2_nocut->Fill(vtx_mc_ev);
     }


     
     if(zpositive &&  pointing!=10 && zpositive && samevent && charmdaughter && vtx_ntrk>=cut_vtx1_ntrk){
       if(impact<cut_ipmax){                              if(charm1)hch1_ipmax->Fill(vtx_mc_ev);       if(charm2)hch2_ipmax->Fill(vtx_mc_ev);
	 if(impact>cut_ipmin){                            if(charm1)hch1_ipmin->Fill(vtx_mc_ev);       if(charm2)hch2_ipmin->Fill(vtx_mc_ev);
	   if(kink>cut_ka){                               if(charm1)hch1_kink->Fill(vtx_mc_ev);        if(charm2)hch2_kink->Fill(vtx_mc_ev);
	     if(tka_rms<=cut_rms){                        if(charm1)hch1_rms->Fill(vtx_mc_ev);         if(charm2)hch2_rms->Fill(vtx_mc_ev); 
	       if(decaylength<cut_decaylength){           if(charm1)hch1_decaylength->Fill(vtx_mc_ev); if(charm2)hch2_decaylength->Fill(vtx_mc_ev); 
		 if(dist_tr<cut_dxy){                     if(charm1)hch1_dxy->Fill(vtx_mc_ev);         if(charm2)hch2_dxy->Fill(vtx_mc_ev); 
		   if(bdt_value>cut_bdt){                 if(charm1)hch1_bdt->Fill(vtx_mc_ev);         if(charm2)hch2_bdt->Fill(vtx_mc_ev); 
		     if((vtx2_z-vtx_z)>cut_dz){           if(charm1)hch1_dzvtx->Fill(vtx_mc_ev);       if(charm2)hch2_dzvtx->Fill(vtx_mc_ev); 
		       if((trk_z-vtx_z)>cut_dz_trk){      if(charm1)hch1_dztrk->Fill(vtx_mc_ev);       if(charm2)hch2_dztrk->Fill(vtx_mc_ev); 
			 if(trk_fill>cut_trk_fill){       if(charm1)hch1_fill->Fill(vtx_mc_ev);        if(charm2)hch2_fill->Fill(vtx_mc_ev);
			   if(pointing<cut_pointing){     if(charm1)hch1_pointing->Fill(vtx_mc_ev);    if(charm2)hch2_pointing->Fill(vtx_mc_ev); 
			     if(abs(trk_pms)>=cut_momentum){ if(charm1)hch1_mom->Fill(vtx_mc_ev);      if(charm2)hch2_mom->Fill(vtx_mc_ev); 
			       
			       //if(((1.3/trk_pms)*(decaylength/(3*TMath::Power(10,14))))<cut_ctau){
				 if(charm1)hch1->Fill(vtx_mc_ev);
				 if(charm2)hch2->Fill(vtx_mc_ev);
				 
				 hsig->Fill(vtx_mc_ev);
				 hsig_vtx->Fill(vtx_id);
				 couples_sig.at(vtx_id).push_back(vtx2_id);
				 //}
			     }// if(abs(trk_pms)>cut_momentum){
			   }
			 }//if(trk_fill>cut_trk_fill){
		       }//if((trk_z-vtx_z)>cut_dz_trk){
		     }//if((vtx2_z-vtx_z)>cut_dz){
		   }//if(bdt_value>cut_bdt){
		 }//if(dist_tr<cut_dxy){
	       }//if(decaylength<cut_decaylength){
	     }//if(nseg>=cut_nseg){
	   }//if(kink<cut_ka){
	 }//if(impact>cut_ipmin){
       }//if(impact<cut_ipmax){
     }//if(zpositive && samevent && charmdaughter){
       
       
     if(zpositive && pointing!=10 && vtx_ntrk>=cut_vtx1_ntrk && !(samevent && charmdaughter)){
       if(impact<cut_ipmax){
	 if(impact>cut_ipmin){
	   if(kink>cut_ka){
	     if(tka_rms<=cut_rms){
	       if(decaylength<cut_decaylength){
		 if(dist_tr<cut_dxy){
		   if(bdt_value>cut_bdt){
		     if((vtx2_z-vtx_z)>cut_dz){
		       if((trk_z-vtx_z)>cut_dz_trk){
			 if(trk_fill>cut_trk_fill){
			   if(pointing<cut_pointing){ 
			     if(abs(trk_pms)>=cut_momentum){
			       //if(((1.3/trk_pms)*(decaylength/(3*TMath::Power(10,14))))<cut_ctau){
				 hbkg->Fill(vtx_mc_ev);
				 hbkg_vtx->Fill(vtx_id);
				 couples_bkg.at(vtx2_id).push_back(vtx_id);
				 couples_bkg2.at(vtx_id).push_back(vtx2_id);
				 //}
			     }// if(abs(trk_pms)>cut_momentum){
			   }
			 }//if(trk_fill>cut_trk_fill){
		       }//if((trk_z-vtx_z)>cut_dz_trk){
		     }//if((vtx2_z-vtx_z)>cut_dz){
		   }//if(bdt_value>cut_bdt){
		 }//if(dist_tr<cut_dxy){
	       }//if(decaylength<cut_decaylength){
	     }//if(nseg>=cut_nseg){
	   }//if(kink<cut_ka){
	 }//if(impact>cut_ipmin){
       }//if(impact<cut_ipmax){
     }//if(!(zpositive && samevent && charmdaughter)){ 

			   

     if(zpositive && pointing!=10 && vtx_ntrk>=cut_vtx1_ntrk){
       //cout << vtx_ntrk << " " << cut_vtx1_ntrk << endl;
        if(impact<cut_ipmax){
	 if(impact>cut_ipmin){
	   if(kink>cut_ka){
	     if(tka_rms<=cut_rms){
	       if(decaylength<cut_decaylength){
		 if(dist_tr<cut_dxy){
		   if(bdt_value>cut_bdt){
		     if((vtx2_z-vtx_z)>cut_dz){
		       if((trk_z-vtx_z)>cut_dz_trk){
			 if(trk_fill>cut_trk_fill){
			   if(pointing<cut_pointing){
			     if(abs(trk_pms)>=cut_momentum){
			       //if(((1.3/trk_pms)*(decaylength/(3*TMath::Power(10,14))))<cut_ctau){
				 hdat->Fill(vtx_mc_ev);
				 hdat_vtx->Fill(vtx_id);
				 log_dat << vtx_id << " " << vtx2_id << " " << bdt_value << endl;
				 //}
			     }// if(abs(trk_pms)>cut_momentum){
			   }
			 }//if(trk_fill>cut_trk_fill){
		       }//if((trk_z-vtx_z)>cut_dz_trk){
		     }//if((vtx2_z-vtx_z)>cut_dz){
		   }//if(bdt_value>cut_bdt){
		 }//if(dist_tr<cut_dxy){
	       }//if(decaylength<cut_decaylength){
	     }//if(nseg>=cut_nseg){
	   }//if(kink<cut_ka){
	 }//if(impact>cut_ipmin){
	}//if(impact<cut_ipmax){
     }//if(zpositive){ 
     

   } // end trk loop


   // contatori
   for(int i=0;i<hsig->GetNbinsX();i++){
     if(hch1->GetBinContent(i+1)>0 && hch2->GetBinContent(i+1)>0)golden++;
     if(hch1->GetBinContent(i+1)>0 || hch2->GetBinContent(i+1)>0)silver++;
     if(hsig->GetBinContent(i+1)>0)nsig++;
     if(hbkg->GetBinContent(i+1)>0)nbkg++;
     if(hdat->GetBinContent(i+1)>0)ndat++;

     if(hch1_nocut->GetBinContent(i+1)>0 && hch2_nocut->GetBinContent(i+1)>0)golden_nocut++;
     if(hch1_nocut->GetBinContent(i+1)>0 || hch2_nocut->GetBinContent(i+1)>0)silver_nocut++;
     if(hsig_nocut->GetBinContent(i+1)>0)nsig_nocut++;

     //ipmax
     if(hch1_ipmax->GetBinContent(i+1)>0 && hch2_ipmax->GetBinContent(i+1)>0)golden_ipmax++;
     if(hch1_ipmax->GetBinContent(i+1)>0 || hch2_ipmax->GetBinContent(i+1)>0)silver_ipmax++;
     //ipmin
     if(hch1_ipmin->GetBinContent(i+1)>0 && hch2_ipmin->GetBinContent(i+1)>0)golden_ipmin++;
     if(hch1_ipmin->GetBinContent(i+1)>0 || hch2_ipmin->GetBinContent(i+1)>0)silver_ipmin++;
     //kink
     if(hch1_kink->GetBinContent(i+1)>0 && hch2_kink->GetBinContent(i+1)>0)golden_kink++;
     if(hch1_kink->GetBinContent(i+1)>0 || hch2_kink->GetBinContent(i+1)>0)silver_kink++;
     //nseg
     if(hch1_rms->GetBinContent(i+1)>0 && hch2_rms->GetBinContent(i+1)>0)golden_rms++;
     if(hch1_rms->GetBinContent(i+1)>0 || hch2_rms->GetBinContent(i+1)>0)silver_rms++;
     //decaylenth
     if(hch1_decaylength->GetBinContent(i+1)>0 && hch2_decaylength->GetBinContent(i+1)>0)golden_decaylength++;
     if(hch1_decaylength->GetBinContent(i+1)>0 || hch2_decaylength->GetBinContent(i+1)>0)silver_decaylength++;
     //dxy
     if(hch1_dxy->GetBinContent(i+1)>0 && hch2_dxy->GetBinContent(i+1)>0)golden_dxy++;
     if(hch1_dxy->GetBinContent(i+1)>0 || hch2_dxy->GetBinContent(i+1)>0)silver_dxy++;
     //bdt
     if(hch1_bdt->GetBinContent(i+1)>0 && hch2_bdt->GetBinContent(i+1)>0)golden_bdt++;
     if(hch1_bdt->GetBinContent(i+1)>0 || hch2_bdt->GetBinContent(i+1)>0)silver_bdt++;
     //dzvtx
     if(hch1_dzvtx->GetBinContent(i+1)>0 && hch2_dzvtx->GetBinContent(i+1)>0)golden_dzvtx++;
     if(hch1_dzvtx->GetBinContent(i+1)>0 || hch2_dzvtx->GetBinContent(i+1)>0)silver_dzvtx++;
     //dztrk
     if(hch1_dztrk->GetBinContent(i+1)>0 && hch2_dztrk->GetBinContent(i+1)>0)golden_dztrk++;
     if(hch1_dztrk->GetBinContent(i+1)>0 || hch2_dztrk->GetBinContent(i+1)>0)silver_dztrk++;
     //fill
     if(hch1_fill->GetBinContent(i+1)>0 && hch2_fill->GetBinContent(i+1)>0)golden_fill++;
     if(hch1_fill->GetBinContent(i+1)>0 || hch2_fill->GetBinContent(i+1)>0)silver_fill++;
     //mom
     if(hch1_mom->GetBinContent(i+1)>0 && hch2_mom->GetBinContent(i+1)>0)golden_mom++;
     if(hch1_mom->GetBinContent(i+1)>0 || hch2_mom->GetBinContent(i+1)>0)silver_mom++;
     //pointing
     if(hch1_pointing->GetBinContent(i+1)>0 && hch2_pointing->GetBinContent(i+1)>0)golden_pointing++;
     if(hch1_pointing->GetBinContent(i+1)>0 || hch2_pointing->GetBinContent(i+1)>0)silver_pointing++;
     
   }

   for(int i=0;i<hsig_vtx->GetNbinsX();i++){
     if(hsig_vtx->GetBinContent(i+1)>0)nsig_vtx++;
     if(hbkg_vtx->GetBinContent(i+1)>0)nbkg_vtx++;
     if(hdat_vtx->GetBinContent(i+1)>0)ndat_vtx++;
   }

   
   for(int ivtx=0;ivtx<nvtx;ivtx++){
     bool mix_found=false;
     bool mix_found2=false;
     int idvtx=0;
     for(int i=0;i<couples_sig.at(ivtx).size();i++){
       idvtx = couples_sig.at(ivtx).at(i);
       for(int j=0;j<couples_bkg.at(idvtx).size();j++){
	 if(couples_bkg.at(idvtx).at(j)==ivtx){
	   mix_found=true;
	 }
       }
       for(int j=0;j<couples_bkg2.at(idvtx).size();j++){
	 if(couples_bkg2.at(idvtx).at(j)==ivtx){
	   mix_found2=true;
	 }
       }
     }
     if(mix_found){
       //cout << mixed << " " <<  idvtx << " " << ivtx << endl;
       mixed++;
     }
     if(mix_found2){
       //cout << mixed2 << " " <<  ivtx << " " << idvtx << endl;
       mixed2++;
     }
   }
   
   silver -= golden;
   silver_nocut -= golden_nocut;
   silver_ipmax -= golden_ipmax;
   silver_ipmin -= golden_ipmin;
   silver_kink -= golden_kink;
   silver_rms -= golden_rms;
   silver_decaylength -= golden_decaylength;
   silver_dxy -= golden_dxy;
   silver_bdt -= golden_bdt;
   silver_dzvtx -= golden_dzvtx;
   silver_dztrk -= golden_dztrk;
   silver_fill -= golden_fill;
   silver_mom -= golden_mom;
   silver_pointing -= golden_pointing;
   
   nbkg_vtx -= mixed;

   cout << "NO CUTS" << endl;
   cout << "Number of signal events " << nsig_nocut << endl;
   cout << "Number of silver events " << silver_nocut << endl;
   cout << "Number of golden events " << golden_nocut << endl;
   cout << "AFTER CUTS" << endl;
   cout << "Number of signal events " << nsig << endl;
   cout << "Number of silver events " << silver << endl;
   cout << "Number of golden events " << golden << endl;
   cout << "Number of signal primary vertices " << nsig_vtx << endl;
   cout << "Number of background primary vertices " << nbkg_vtx << endl;
   cout << "Number of dirty secondary vertices " << mixed << endl;
   cout << "Number of events in data " << ndat << endl;
   cout << "Number of primary vertices in data " << ndat_vtx << endl;
   cout << "BREAKDOWN EFFICIENCIES" << endl;
   cout << "        SILVER        " << endl;
   cout << "starting sample       " << silver_nocut       << endl;
   cout << "ipmax                 " << silver_ipmax       << " -> " <<  silver_ipmax      *1.0/silver_nocut << endl;
   cout << "ipmin                 " << silver_ipmin       << " -> " <<  silver_ipmin      *1.0/silver_nocut << endl;
   cout << "kink angle            " << silver_kink        << " -> " <<  silver_kink       *1.0/silver_nocut << endl;
   cout << "trk rms               " << silver_rms         << " -> " <<  silver_rms        *1.0/silver_nocut << endl;
   cout << "decay length          " << silver_decaylength << " -> " <<  silver_decaylength*1.0/silver_nocut << endl;
   cout << "dist xy               " << silver_dxy         << " -> " <<  silver_dxy        *1.0/silver_nocut << endl;
   cout << "bdt                   " << silver_bdt         << " -> " <<  silver_bdt        *1.0/silver_nocut << endl;
   cout << "dz vtx                " << silver_dzvtx       << " -> " <<  silver_dzvtx      *1.0/silver_nocut << endl;
   cout << "dz trk                " << silver_dztrk       << " -> " <<  silver_dztrk      *1.0/silver_nocut << endl;
   cout << "fill factor           " << silver_fill        << " -> " <<  silver_fill       *1.0/silver_nocut << endl;
   cout << "pointing              " << silver_pointing    << " -> " <<  silver_pointing   *1.0/silver_nocut << endl;
   cout << "momentum              " << silver_mom         << " -> " <<  silver_mom        *1.0/silver_nocut << endl;
   cout << "        ******        " << endl;
   cout << "BREAKDOWN EFFICIENCIES" << endl;
   cout << "        GOLDEN        " << endl;
   cout << "starting sample       " << golden_nocut       << endl;
   cout << "ipmax                 " << golden_ipmax       << " -> " <<  golden_ipmax      *1.0/golden_nocut << endl;
   cout << "ipmin                 " << golden_ipmin       << " -> " <<  golden_ipmin      *1.0/golden_nocut << endl;
   cout << "kink angle            " << golden_kink        << " -> " <<  golden_kink       *1.0/golden_nocut << endl;
   cout << "rms                  "  << golden_rms         << " -> " <<  golden_rms        *1.0/golden_nocut << endl;
   cout << "decay length          " << golden_decaylength << " -> " <<  golden_decaylength*1.0/golden_nocut << endl;
   cout << "dist xy               " << golden_dxy         << " -> " <<  golden_dxy        *1.0/golden_nocut << endl;
   cout << "bdt                   " << golden_bdt         << " -> " <<  golden_bdt        *1.0/golden_nocut << endl;
   cout << "dz vtx                " << golden_dzvtx       << " -> " <<  golden_dzvtx      *1.0/golden_nocut << endl;
   cout << "dz trk                " << golden_dztrk       << " -> " <<  golden_dztrk      *1.0/golden_nocut << endl;
   cout << "fill factor           " << golden_fill        << " -> " <<  golden_fill       *1.0/golden_nocut << endl;
   cout << "pointing              " << golden_pointing    << " -> " <<  golden_pointing   *1.0/golden_nocut << endl;
   cout << "momentum              " << golden_mom         << " -> " <<  golden_mom        *1.0/golden_nocut << endl;


   log_dat.close();
   
   ofstream log("log_ds_results.txt");
   log << "NO CUTS" << endl;
   log << "Number of signal events " << nsig_nocut << endl;
   log << "Number of silver events " << silver_nocut << endl;
   log << "Number of golden events " << golden_nocut << endl;
   log << "AFTER CUTS" << endl;
   log << "Number of signal events " << nsig << endl;
   log << "Number of silver events " << silver << endl;
   log << "Number of golden events " << golden << endl;
   log << "Number of signal primary vertices " << nsig_vtx << endl;
   log << "Number of background primary vertices " << nbkg_vtx << endl;
   log << "Number of dirty secondary vertices " << mixed << endl;
   log << "Number of events in data " << ndat << endl;
   log << "Number of primary vertices in data " << ndat_vtx << endl;
   log << endl;
   log << "Cut used:" << endl;
   log << "Cut impact min " << cut_ipmin << endl;
   log << "Cut impact max " << cut_ipmax << endl;
   log << "Cut kink " << cut_ka << endl;
   log << "Cut trk rms " << cut_rms << endl;
   log << "Cut transverse distance " << cut_dxy << endl;
   log << "Cut bdt DS " << cut_bdt << endl;
   log << "Cut dz between decay and primary vertices " << cut_dz << endl;
   log << "Cut dz between decay tracks and primary vertex " << cut_dz_trk << endl;
   log << "Cut track fill factor " << cut_trk_fill << endl;
   log << "Cut track pointing " << cut_pointing << endl;
   log << "Cut track momentum " << cut_momentum << endl;
   log << "Additional info: "<< endl;
   log << "Number of reverse couples " << mixed << endl;
   log.close();
}
