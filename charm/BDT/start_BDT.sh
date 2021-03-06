#!/bin/bash   # CREATED BY V.GENTILE
#
#
echo -e "Type: \n 100 for BDT1 \n 200 for BDT2 \n 300 for BDTds \n"
echo -e "Type: \n 120 for BDT1 and BDT2 \n 123 for all BDTs\n"
echo -e "Type: \n 230 for BDT2 and BDTds \n"
read choice

if [ $choice -ge 100 ] && [ $choice -lt 200 ]; then
    root -l <<EOC
gROOT->ProcessLine(".L TMVAClassification.C");
prepareTMVAtree();
.q
EOC
    root -l TMVAClassificationApplication.C\(1\) <<EOC
.q    
EOC

    if [ $choice -eq 120 ]; then
	root -l <<EOC
gROOT->ProcessLine(".L TMVAClassification.C");
prepareTMVAtree2nd(0.15);
.q
EOC
	root -l TMVAClassificationApplication2nd.C\(1\) <<EOC
.q
EOC
    fi
fi

cp tmva_* vtx_BDT* ../
