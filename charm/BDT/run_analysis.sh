#!bin/bash     # CREATED BY V.GENTILE
#
#
if [ ! -f all.root ]; then
    echo "Merging all plots"
    hadd all.root vz* n_*
fi

if [ ! -f all_tree.root ]; then
    echo "Merging all tree"
    hadd all_tree.root ../firstquarter/TMVA/vtx_BDT* ../secondquarter/TMVA/vtx_BDT* ../thirdquarter/TMVA/vtx_BDT* ../fourthquarter/TMVA/vtx_BDT*
fi


echo "Type the configuration (e.g. CH1R6)"
read name;

if [ ! -f vz_${name}_data.root ]; then
    root -l all.root <<EOC
vz_1st_data->Add(vz_2nd_data);
vz_1st_data->Add(vz_3rd_data);	
vz_1st_data->Add(vz_4th_data);	
n_1st_data->Add(n_2nd_data);
n_1st_data->Add(n_3rd_data);
n_1st_data->Add(n_4th_data);
vz_1st_data->GetXaxis()->SetTitle("z [mm]");
n_1st_data->GetXaxis()->SetTitle("multiplicity");
vz_1st_data->SaveAs("vz.root")
n_1st_data->SaveAs("n.root")
EOC

    mv vz.root vz_${name}_data.root
    mv n.root n_${name}_data.root
fi
