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
prepareTMVAtree()
.q
EOC
    root -l TMVAClassificationApplication.C\(1\)
fi
