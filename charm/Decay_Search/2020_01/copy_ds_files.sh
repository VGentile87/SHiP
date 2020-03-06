#!/bin/bash   # CREATED BY V.GENTILE
#
#
echo -e "Type: 1 for MC study mode or 2 for data study mode \n"
read choice

if [ $choice -ge 1 ]; then
    ln -s $CHARM/Decay_Search/2020_01/vtx_data_study_MC.C .
fi

if [ $choice -ge 2 ]; then
    ln -s $CHARM/Decay_Search/2020_01/vtx_data_study.C .
fi

ln -s $CHARM/Decay_Search/2020_01/Definitions.h .
ln -s $CHARM/Decay_Search/2020_01/DsCuts.h .
