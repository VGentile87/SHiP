#!/bin/bash    # CREATED BY V.GENTILE
#
#

if [ ! -d ../../analysis ]; then
    mkdir -p ../../analysis;
fi

echo -e "Type the number of the  quarter \n"
read choice
if [ $choice -eq 1 ]; then
    root -l vtx_BDT_data_evaluated.root <<EOC
bdt->Draw("vz/1000.>>vz_1st_data(40,-40,0)","bdt_value>0.15 && ntracks>5 && vz<-4000")
vz_1st_data->SaveAs("vz_1st_data.root")
bdt->Draw("ntracks>>n_1st_data(70,0,70)","bdt_value>0.15 && ntracks>5 && vz<-4000")
n_1st_data->SaveAs("n_1st_data.root");
.q
EOC
    cp n_1st_data.root  vz_1st_data.root ../../analysis/
elif [ $choice -eq 2 ] ; then
     root -l vtx_BDT_data_evaluated.root <<EOC
bdt->Draw("vz/1000.>>vz_2nd_data(40,-40,0)","bdt_value>0.15 && ntracks>5 && vz<-4000")
vz_2nd_data->SaveAs("vz_2nd_data.root")
bdt->Draw("ntracks>>n_2nd_data(70,0,70)","bdt_value>0.15 && ntracks>5 && vz<-4000")
n_2nd_data->SaveAs("n_2nd_data.root");
.q
EOC
     cp n_2nd_data.root  vz_2nd_data.root ../../analysis/
elif [ $choice -eq 3 ] ; then
     root -l vtx_BDT_data_evaluated.root <<EOC
bdt->Draw("vz/1000.>>vz_3rd_data(40,-40,0)","bdt_value>0.15 && ntracks>5 && vz<-4000")
vz_3rd_data->SaveAs("vz_3rd_data.root")
bdt->Draw("ntracks>>n_3rd_data(70,0,70)","bdt_value>0.15 && ntracks>5 && vz<-4000")
n_3rd_data->SaveAs("n_3rd_data.root");
.q
EOC
     cp n_3rd_data.root  vz_3rd_data.root ../../analysis/
elif [ $choice -eq 4 ] ; then
     root -l vtx_BDT_data_evaluated.root <<EOC
bdt->Draw("vz/1000.>>vz_4th_data(40,-40,0)","bdt_value>0.15 && ntracks>5 && vz<-4000")
vz_4th_data->SaveAs("vz_4th_data.root")
bdt->Draw("ntracks>>n_4th_data(70,0,70)","bdt_value>0.15 && ntracks>5 && vz<-4000")
n_4th_data->SaveAs("n_4th_data.root");
.q
EOC
     cp n_4th_data.root  vz_4th_data.root ../../analysis/
else break;
fi
