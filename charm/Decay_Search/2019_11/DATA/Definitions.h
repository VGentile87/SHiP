#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"
#include "EdbSegP.h"
#include "EdbID.h"


// FEDRA VARIABLES

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

// VTX TREE

   // Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t vkMaxt_ = 79;
   static constexpr Int_t vkMaxs = 700;
   static constexpr Int_t vkMaxsf = 700;

   // Declaration of leaf types
   Int_t           vID;
   Int_t           flag;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         maxaperture;
   Float_t         probability;
   Int_t           n;
   Int_t           t__;
   UInt_t          t__fUniqueID[vkMaxt_];   //[t._]
   UInt_t          t__fBits[vkMaxt_];   //[t._]
   Int_t           t__ePID[vkMaxt_];   //[t._]
   Int_t           t__eID[vkMaxt_];   //[t._]
   Int_t           t__eVid[vkMaxt_][2];   //[t._]
   Int_t           t__eAid[vkMaxt_][2];   //[t._]
   Int_t           t__eFlag[vkMaxt_];   //[t._]
   Int_t           t__eTrack[vkMaxt_];   //[t._]
   Float_t         t__eX[vkMaxt_];   //[t._]
   Float_t         t__eY[vkMaxt_];   //[t._]
   Float_t         t__eZ[vkMaxt_];   //[t._]
   Float_t         t__eTX[vkMaxt_];   //[t._]
   Float_t         t__eTY[vkMaxt_];   //[t._]
   Float_t         t__eSZ[vkMaxt_];   //[t._]
   Float_t         t__eChi2[vkMaxt_];   //[t._]
   Float_t         t__eProb[vkMaxt_];   //[t._]
   Float_t         t__eW[vkMaxt_];   //[t._]
   Float_t         t__eVolume[vkMaxt_];   //[t._]
   Float_t         t__eDZ[vkMaxt_];   //[t._]
   Float_t         t__eDZem[vkMaxt_];   //[t._]
   Float_t         t__eP[vkMaxt_];   //[t._]
   Int_t           t__eMCTrack[vkMaxt_];   //[t._]
   Int_t           t__eMCEvt[vkMaxt_];   //[t._]
   UInt_t          t__eScanID_fUniqueID[vkMaxt_];   //[t._]
   UInt_t          t__eScanID_fBits[vkMaxt_];   //[t._]
   Int_t           t__eScanID_eBrick[vkMaxt_];   //[t._]
   Int_t           t__eScanID_ePlate[vkMaxt_];   //[t._]
   Int_t           t__eScanID_eMajor[vkMaxt_];   //[t._]
   Int_t           t__eScanID_eMinor[vkMaxt_];   //[t._]
   Int_t           s__;
   UInt_t          s__fUniqueID[vkMaxs];   //[s__]
   UInt_t          s__fBits[vkMaxs];   //[s__]
   Int_t           s__ePID[vkMaxs];   //[s__]
   Int_t           s__eID[vkMaxs];   //[s__]
   Int_t           s__eVid[vkMaxs][2];   //[s__]
   Int_t           s__eAid[vkMaxs][2];   //[s__]
   Int_t           s__eFlag[vkMaxs];   //[s__]
   Int_t           s__eTrack[vkMaxs];   //[s__]
   Float_t         s__eX[vkMaxs];   //[s__]
   Float_t         s__eY[vkMaxs];   //[s__]
   Float_t         s__eZ[vkMaxs];   //[s__]
   Float_t         s__eTX[vkMaxs];   //[s__]
   Float_t         s__eTY[vkMaxs];   //[s__]
   Float_t         s__eSZ[vkMaxs];   //[s__]
   Float_t         s__eChi2[vkMaxs];   //[s__]
   Float_t         s__eProb[vkMaxs];   //[s__]
   Float_t         s__eW[vkMaxs];   //[s__]
   Float_t         s__eVolume[vkMaxs];   //[s__]
   Float_t         s__eDZ[vkMaxs];   //[s__]
   Float_t         s__eDZem[vkMaxs];   //[s__]
   Float_t         s__eP[vkMaxs];   //[s__]
   Int_t           s__eMCTrack[vkMaxs];   //[s__]
   Int_t           s__eMCEvt[vkMaxs];   //[s__]
   UInt_t          s__eScanID_fUniqueID[vkMaxs];   //[s__]
   UInt_t          s__eScanID_fBits[vkMaxs];   //[s__]
   Int_t           s__eScanID_eBrick[vkMaxs];   //[s__]
   Int_t           s__eScanID_ePlate[vkMaxs];   //[s__]
   Int_t           s__eScanID_eMajor[vkMaxs];   //[s__]
   Int_t           s__eScanID_eMinor[vkMaxs];   //[s__]
   Int_t           sf__;
   UInt_t          sf__fUniqueID[vkMaxsf];   //[sf__]
   UInt_t          sf__fBits[vkMaxsf];   //[sf__]
   Int_t           sf__ePID[vkMaxsf];   //[sf__]
   Int_t           sf__eID[vkMaxsf];   //[sf__]
   Int_t           sf__eVid[vkMaxsf][2];   //[sf__]
   Int_t           sf__eAid[vkMaxsf][2];   //[sf__]
   Int_t           sf__eFlag[vkMaxsf];   //[sf__]
   Int_t           sf__eTrack[vkMaxsf];   //[sf__]
   Float_t         sf__eX[vkMaxsf];   //[sf__]
   Float_t         sf__eY[vkMaxsf];   //[sf__]
   Float_t         sf__eZ[vkMaxsf];   //[sf__]
   Float_t         sf__eTX[vkMaxsf];   //[sf__]
   Float_t         sf__eTY[vkMaxsf];   //[sf__]
   Float_t         sf__eSZ[vkMaxsf];   //[sf__]
   Float_t         sf__eChi2[vkMaxsf];   //[sf__]
   Float_t         sf__eProb[vkMaxsf];   //[sf__]
   Float_t         sf__eW[vkMaxsf];   //[sf__]
   Float_t         sf__eVolume[vkMaxsf];   //[sf__]
   Float_t         sf__eDZ[vkMaxsf];   //[sf__]
   Float_t         sf__eDZem[vkMaxsf];   //[sf__]
   Float_t         sf__eP[vkMaxsf];   //[sf__]
   Int_t           sf__eMCTrack[vkMaxsf];   //[sf__]
   Int_t           sf__eMCEvt[vkMaxsf];   //[sf__]
   UInt_t          sf__eScanID_fUniqueID[vkMaxsf];   //[sf__]
   UInt_t          sf__eScanID_fBits[vkMaxsf];   //[sf__]
   Int_t           sf__eScanID_eBrick[vkMaxsf];   //[sf__]
   Int_t           sf__eScanID_ePlate[vkMaxsf];   //[sf__]
   Int_t           sf__eScanID_eMajor[vkMaxsf];   //[sf__]
   Int_t           sf__eScanID_eMinor[vkMaxsf];   //[sf__]
   Int_t           TrackID[79];   //[n]
   Int_t           nseg[79];   //[n]
   Int_t           nholes[79];   //[n]
   Int_t           maxgap[79];   //[n]
   Int_t           incoming[79];   //[n]
   Float_t         impactparameter[79];   //[n]
   Int_t           MCEventID[79];   //[n]
   Int_t           MCTrackID[79];   //[n]
   Int_t           MCMotherID[79];   //[n]

   // List of branches
   TBranch        *b_vID;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_maxaperture;   //!
   TBranch        *b_probability;   //!
   TBranch        *b_n;   //!
   TBranch        *b_t__;   //!
   TBranch        *b_t__fUniqueID;   //!
   TBranch        *b_t__fBits;   //!
   TBranch        *b_t__ePID;   //!
   TBranch        *b_t__eID;   //!
   TBranch        *b_t__eVid;   //!
   TBranch        *b_t__eAid;   //!
   TBranch        *b_t__eFlag;   //!
   TBranch        *b_t__eTrack;   //!
   TBranch        *b_t__eX;   //!
   TBranch        *b_t__eY;   //!
   TBranch        *b_t__eZ;   //!
   TBranch        *b_t__eTX;   //!
   TBranch        *b_t__eTY;   //!
   TBranch        *b_t__eSZ;   //!
   TBranch        *b_t__eChi2;   //!
   TBranch        *b_t__eProb;   //!
   TBranch        *b_t__eW;   //!
   TBranch        *b_t__eVolume;   //!
   TBranch        *b_t__eDZ;   //!
   TBranch        *b_t__eDZem;   //!
   TBranch        *b_t__eP;   //!
   TBranch        *b_t__eMCTrack;   //!
   TBranch        *b_t__eMCEvt;   //!
   TBranch        *b_t__eScanID_fUniqueID;   //!
   TBranch        *b_t__eScanID_fBits;   //!
   TBranch        *b_t__eScanID_eBrick;   //!
   TBranch        *b_t__eScanID_ePlate;   //!
   TBranch        *b_t__eScanID_eMajor;   //!
   TBranch        *b_t__eScanID_eMinor;   //!
   TBranch        *b_s__;   //!
   TBranch        *b_s__fUniqueID;   //!
   TBranch        *b_s__fBits;   //!
   TBranch        *b_s__ePID;   //!
   TBranch        *b_s__eID;   //!
   TBranch        *b_s__eVid;   //!
   TBranch        *b_s__eAid;   //!
   TBranch        *b_s__eFlag;   //!
   TBranch        *b_s__eTrack;   //!
   TBranch        *b_s__eX;   //!
   TBranch        *b_s__eY;   //!
   TBranch        *b_s__eZ;   //!
   TBranch        *b_s__eTX;   //!
   TBranch        *b_s__eTY;   //!
   TBranch        *b_s__eSZ;   //!
   TBranch        *b_s__eChi2;   //!
   TBranch        *b_s__eProb;   //!
   TBranch        *b_s__eW;   //!
   TBranch        *b_s__eVolume;   //!
   TBranch        *b_s__eDZ;   //!
   TBranch        *b_s__eDZem;   //!
   TBranch        *b_s__eP;   //!
   TBranch        *b_s__eMCTrack;   //!
   TBranch        *b_s__eMCEvt;   //!
   TBranch        *b_s__eScanID_fUniqueID;   //!
   TBranch        *b_s__eScanID_fBits;   //!
   TBranch        *b_s__eScanID_eBrick;   //!
   TBranch        *b_s__eScanID_ePlate;   //!
   TBranch        *b_s__eScanID_eMajor;   //!
   TBranch        *b_s__eScanID_eMinor;   //!
   TBranch        *b_sf__;   //!
   TBranch        *b_sf__fUniqueID;   //!
   TBranch        *b_sf__fBits;   //!
   TBranch        *b_sf__ePID;   //!
   TBranch        *b_sf__eID;   //!
   TBranch        *b_sf__eVid;   //!
   TBranch        *b_sf__eAid;   //!
   TBranch        *b_sf__eFlag;   //!
   TBranch        *b_sf__eTrack;   //!
   TBranch        *b_sf__eX;   //!
   TBranch        *b_sf__eY;   //!
   TBranch        *b_sf__eZ;   //!
   TBranch        *b_sf__eTX;   //!
   TBranch        *b_sf__eTY;   //!
   TBranch        *b_sf__eSZ;   //!
   TBranch        *b_sf__eChi2;   //!
   TBranch        *b_sf__eProb;   //!
   TBranch        *b_sf__eW;   //!
   TBranch        *b_sf__eVolume;   //!
   TBranch        *b_sf__eDZ;   //!
   TBranch        *b_sf__eDZem;   //!
   TBranch        *b_sf__eP;   //!
   TBranch        *b_sf__eMCTrack;   //!
   TBranch        *b_sf__eMCEvt;   //!
   TBranch        *b_sf__eScanID_fUniqueID;   //!
   TBranch        *b_sf__eScanID_fBits;   //!
   TBranch        *b_sf__eScanID_eBrick;   //!
   TBranch        *b_sf__eScanID_ePlate;   //!
   TBranch        *b_sf__eScanID_eMajor;   //!
   TBranch        *b_sf__eScanID_eMinor;   //!
   TBranch        *b_TrackID;   //!
   TBranch        *b_nseg;   //!
   TBranch        *b_nholes;   //!
   TBranch        *b_maxgap;   //!
   TBranch        *b_incoming;   //!
   TBranch        *b_impactparameter;   //!
   TBranch        *b_MCEventID;   //!
   TBranch        *b_MCTrackID;   //!
   TBranch        *b_MCMotherID;   //!

// TRACKS VARIABLES
/*
// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxt = 1;
   static constexpr Int_t kMaxs = 29;
   static constexpr Int_t kMaxsf = 29;

   // Declaration of leaf types
   Int_t           trid;
   Int_t           tnseg;
   Int_t           npl;
   Int_t           n0;
   Float_t         xv;
   Float_t         yv;
   Float_t         w;
   EdbSegP         *t_;
   UInt_t          t_TObject_fUniqueID;
   UInt_t          t_TObject_fBits;
   Int_t           t_ePID;
   Int_t           t_eID;
   Int_t           t_eVid[2];
   Int_t           t_eAid[2];
   Int_t           t_eFlag;
   Int_t           t_eTrack;
   Float_t         t_eX;
   Float_t         t_eY;
   Float_t         t_eZ;
   Float_t         t_eTX;
   Float_t         t_eTY;
   Float_t         t_eSZ;
   Float_t         t_eChi2;
   Float_t         t_eProb;
   Float_t         t_eW;
   Float_t         t_eVolume;
   Float_t         t_eDZ;
   Float_t         t_eDZem;
   Float_t         t_eP;
   Int_t           t_eMCTrack;
   Int_t           t_eMCEvt;
   UInt_t          t_eScanID_fUniqueID;
   UInt_t          t_eScanID_fBits;
   Int_t           t_eScanID_eBrick;
   Int_t           t_eScanID_ePlate;
   Int_t           t_eScanID_eMajor;
   Int_t           t_eScanID_eMinor;
   Int_t           s_;
   UInt_t          s_fUniqueID[kMaxs];   //[s_]
   UInt_t          s_fBits[kMaxs];   //[s_]
   Int_t           s_ePID[kMaxs];   //[s_]
   Int_t           s_eID[kMaxs];   //[s_]
   Int_t           s_eVid[kMaxs][2];   //[s_]
   Int_t           s_eAid[kMaxs][2];   //[s_]
   Int_t           s_eFlag[kMaxs];   //[s_]
   Int_t           s_eTrack[kMaxs];   //[s_]
   Float_t         s_eX[kMaxs];   //[s_]
   Float_t         s_eY[kMaxs];   //[s_]
   Float_t         s_eZ[kMaxs];   //[s_]
   Float_t         s_eTX[kMaxs];   //[s_]
   Float_t         s_eTY[kMaxs];   //[s_]
   Float_t         s_eSZ[kMaxs];   //[s_]
   Float_t         s_eChi2[kMaxs];   //[s_]
   Float_t         s_eProb[kMaxs];   //[s_]
   Float_t         s_eW[kMaxs];   //[s_]
   Float_t         s_eVolume[kMaxs];   //[s_]
   Float_t         s_eDZ[kMaxs];   //[s_]
   Float_t         s_eDZem[kMaxs];   //[s_]
   Float_t         s_eP[kMaxs];   //[s_]
   Int_t           s_eMCTrack[kMaxs];   //[s_]
   Int_t           s_eMCEvt[kMaxs];   //[s_]
   UInt_t          s_eScanID_fUniqueID[kMaxs];   //[s_]
   UInt_t          s_eScanID_fBits[kMaxs];   //[s_]
   Int_t           s_eScanID_eBrick[kMaxs];   //[s_]
   Int_t           s_eScanID_ePlate[kMaxs];   //[s_]
   Int_t           s_eScanID_eMajor[kMaxs];   //[s_]
   Int_t           s_eScanID_eMinor[kMaxs];   //[s_]
Int_t           sf_;
   UInt_t          sf_fUniqueID[kMaxsf];   //[sf_]
   UInt_t          sf_fBits[kMaxsf];   //[sf_]
   Int_t           sf_ePID[kMaxsf];   //[sf_]
   Int_t           sf_eID[kMaxsf];   //[sf_]
   Int_t           sf_eVid[kMaxsf][2];   //[sf_]
   Int_t           sf_eAid[kMaxsf][2];   //[sf_]
   Int_t           sf_eFlag[kMaxsf];   //[sf_]
   Int_t           sf_eTrack[kMaxsf];   //[sf_]
   Float_t         sf_eX[kMaxsf];   //[sf_]
   Float_t         sf_eY[kMaxsf];   //[sf_]
   Float_t         sf_eZ[kMaxsf];   //[sf_]
   Float_t         sf_eTX[kMaxsf];   //[sf_]
   Float_t         sf_eTY[kMaxsf];   //[sf_]
   Float_t         sf_eSZ[kMaxsf];   //[sf_]
   Float_t         sf_eChi2[kMaxsf];   //[sf_]
   Float_t         sf_eProb[kMaxsf];   //[sf_]
   Float_t         sf_eW[kMaxsf];   //[sf_]
   Float_t         sf_eVolume[kMaxsf];   //[sf_]
   Float_t         sf_eDZ[kMaxsf];   //[sf_]
   Float_t         sf_eDZem[kMaxsf];   //[sf_]
   Float_t         sf_eP[kMaxsf];   //[sf_]
   Int_t           sf_eMCTrack[kMaxsf];   //[sf_]
   Int_t           sf_eMCEvt[kMaxsf];   //[sf_]
   UInt_t          sf_eScanID_fUniqueID[kMaxsf];   //[sf_]
   UInt_t          sf_eScanID_fBits[kMaxsf];   //[sf_]
   Int_t           sf_eScanID_eBrick[kMaxsf];   //[sf_]
   Int_t           sf_eScanID_ePlate[kMaxsf];   //[sf_]
   Int_t           sf_eScanID_eMajor[kMaxsf];   //[sf_]
   Int_t           sf_eScanID_eMinor[kMaxsf];   //[sf_]

   // List of branches
   TBranch        *b_trid;   //!
   TBranch        *b_tnseg;   //!
   TBranch        *b_npl;   //!
   TBranch        *b_n0;   //!
   TBranch        *b_xv;   //!
   TBranch        *b_yv;   //!
   TBranch        *b_w;   //!
   TBranch        *b_t_TObject_fUniqueID;   //!
   TBranch        *b_t_TObject_fBits;   //!
   TBranch        *b_t_ePID;   //!
   TBranch        *b_t_eID;   //!
   TBranch        *b_t_eVid;   //!
   TBranch        *b_t_eAid;   //!
   TBranch        *b_t_eFlag;   //!
   TBranch        *b_t_eTrack;   //!
   TBranch        *b_t_eX;   //!
   TBranch        *b_t_eY;   //!
   TBranch        *b_t_eZ;   //!
   TBranch        *b_t_eTX;   //!
   TBranch        *b_t_eTY;   //!
   TBranch        *b_t_eSZ;   //!
   TBranch        *b_t_eChi2;   //!
   TBranch        *b_t_eProb;   //!
   TBranch        *b_t_eW;   //!
   TBranch        *b_t_eVolume;   //!
   TBranch        *b_t_eDZ;   //!
   TBranch        *b_t_eDZem;   //!
   TBranch        *b_t_eP;   //!
   TBranch        *b_t_eMCTrack;   //!
   TBranch        *b_t_eMCEvt;   //!
   TBranch        *b_t_eScanID_fUniqueID;   //!
   TBranch        *b_t_eScanID_fBits;   //!
   TBranch        *b_t_eScanID_eBrick;   //!
   TBranch        *b_t_eScanID_ePlate;   //!
   TBranch        *b_t_eScanID_eMajor;   //!
   TBranch        *b_t_eScanID_eMinor;   //!
   TBranch        *b_s_;   //!
   TBranch        *b_s_fUniqueID;   //!
   TBranch        *b_s_fBits;   //!
   TBranch        *b_s_ePID;   //!
   TBranch        *b_s_eID;   //!
   TBranch        *b_s_eVid;   //!
   TBranch        *b_s_eAid;   //!
   TBranch        *b_s_eFlag;   //!
   TBranch        *b_s_eTrack;   //!
   TBranch        *b_s_eX;   //!
   TBranch        *b_s_eY;   //!
   TBranch        *b_s_eZ;   //!
   TBranch        *b_s_eTX;   //!
   TBranch        *b_s_eTY;   //!
   TBranch        *b_s_eSZ;   //!
   TBranch        *b_s_eChi2;   //!
   TBranch        *b_s_eProb;   //!
   TBranch        *b_s_eW;   //!
   TBranch        *b_s_eVolume;   //!
   TBranch        *b_s_eDZ;   //!
   TBranch        *b_s_eDZem;   //!
   TBranch        *b_s_eP;   //!
   TBranch        *b_s_eMCTrack;   //!
   TBranch        *b_s_eMCEvt;   //!
   TBranch        *b_s_eScanID_fUniqueID;   //!
   TBranch        *b_s_eScanID_fBits;   //!
   TBranch        *b_s_eScanID_eBrick;   //!
   TBranch        *b_s_eScanID_ePlate;   //!
   TBranch        *b_s_eScanID_eMajor;   //!
   TBranch        *b_s_eScanID_eMinor;   //!
TBranch        *b_sf_;   //!
   TBranch        *b_sf_fUniqueID;   //!
   TBranch        *b_sf_fBits;   //!
   TBranch        *b_sf_ePID;   //!
   TBranch        *b_sf_eID;   //!
   TBranch        *b_sf_eVid;   //!
   TBranch        *b_sf_eAid;   //!
   TBranch        *b_sf_eFlag;   //!
   TBranch        *b_sf_eTrack;   //!
   TBranch        *b_sf_eX;   //!
   TBranch        *b_sf_eY;   //!
   TBranch        *b_sf_eZ;   //!
   TBranch        *b_sf_eTX;   //!
   TBranch        *b_sf_eTY;   //!
   TBranch        *b_sf_eSZ;   //!
   TBranch        *b_sf_eChi2;   //!
   TBranch        *b_sf_eProb;   //!
   TBranch        *b_sf_eW;   //!
   TBranch        *b_sf_eVolume;   //!
   TBranch        *b_sf_eDZ;   //!
   TBranch        *b_sf_eDZem;   //!
   TBranch        *b_sf_eP;   //!
   TBranch        *b_sf_eMCTrack;   //!
   TBranch        *b_sf_eMCEvt;   //!
   TBranch        *b_sf_eScanID_fUniqueID;   //!
   TBranch        *b_sf_eScanID_fBits;   //!
   TBranch        *b_sf_eScanID_eBrick;   //!
   TBranch        *b_sf_eScanID_ePlate;   //!
   TBranch        *b_sf_eScanID_eMajor;   //!
   TBranch        *b_sf_eScanID_eMinor;   //!
*/
// BDT VARIABLES

// Declaration of leaf types
   Float_t         bdt_value;

   // List of branches
   TBranch        *b_bdt_value;   //!



void  vtx_reader_Fedra(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("vID", &vID, &b_vID);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("maxaperture", &maxaperture, &b_maxaperture);
   fChain->SetBranchAddress("probability", &probability, &b_probability);
   fChain->SetBranchAddress("n", &n, &b_n);
   fChain->SetBranchAddress("t.", &t__, &b_t__);
   fChain->SetBranchAddress("t..fUniqueID", t__fUniqueID, &b_t__fUniqueID);
   fChain->SetBranchAddress("t..fBits", t__fBits, &b_t__fBits);
   fChain->SetBranchAddress("t..ePID", t__ePID, &b_t__ePID);
   fChain->SetBranchAddress("t..eID", t__eID, &b_t__eID);
   fChain->SetBranchAddress("t..eVid[2]", t__eVid, &b_t__eVid);
   fChain->SetBranchAddress("t..eAid[2]", t__eAid, &b_t__eAid);
   fChain->SetBranchAddress("t..eFlag", t__eFlag, &b_t__eFlag);
   fChain->SetBranchAddress("t..eTrack", t__eTrack, &b_t__eTrack);
   fChain->SetBranchAddress("t..eX", t__eX, &b_t__eX);
   fChain->SetBranchAddress("t..eY", t__eY, &b_t__eY);
   fChain->SetBranchAddress("t..eZ", t__eZ, &b_t__eZ);
   fChain->SetBranchAddress("t..eTX", t__eTX, &b_t__eTX);
   fChain->SetBranchAddress("t..eTY", t__eTY, &b_t__eTY);
   fChain->SetBranchAddress("t..eSZ", t__eSZ, &b_t__eSZ);
   fChain->SetBranchAddress("t..eChi2", t__eChi2, &b_t__eChi2);
   fChain->SetBranchAddress("t..eProb", t__eProb, &b_t__eProb);
   fChain->SetBranchAddress("t..eW", t__eW, &b_t__eW);
   fChain->SetBranchAddress("t..eVolume", t__eVolume, &b_t__eVolume);
   fChain->SetBranchAddress("t..eDZ", t__eDZ, &b_t__eDZ);
   fChain->SetBranchAddress("t..eDZem", t__eDZem, &b_t__eDZem);
   fChain->SetBranchAddress("t..eP", t__eP, &b_t__eP);
   fChain->SetBranchAddress("t..eMCTrack", t__eMCTrack, &b_t__eMCTrack);
   fChain->SetBranchAddress("t..eMCEvt", t__eMCEvt, &b_t__eMCEvt);
   fChain->SetBranchAddress("t..eScanID.fUniqueID", t__eScanID_fUniqueID, &b_t__eScanID_fUniqueID);
   fChain->SetBranchAddress("t..eScanID.fBits", t__eScanID_fBits, &b_t__eScanID_fBits);
   fChain->SetBranchAddress("t..eScanID.eBrick", t__eScanID_eBrick, &b_t__eScanID_eBrick);
   fChain->SetBranchAddress("t..eScanID.ePlate", t__eScanID_ePlate, &b_t__eScanID_ePlate);
   fChain->SetBranchAddress("t..eScanID.eMajor", t__eScanID_eMajor, &b_t__eScanID_eMajor);
   fChain->SetBranchAddress("t..eScanID.eMinor", t__eScanID_eMinor, &b_t__eScanID_eMinor);
   fChain->SetBranchAddress("s", &s__, &b_s__);
   fChain->SetBranchAddress("s.fUniqueID", s__fUniqueID, &b_s__fUniqueID);
   fChain->SetBranchAddress("s.fBits", s__fBits, &b_s__fBits);
   fChain->SetBranchAddress("s.ePID", s__ePID, &b_s__ePID);
   fChain->SetBranchAddress("s.eID", s__eID, &b_s__eID);
   fChain->SetBranchAddress("s.eVid[2]", s__eVid, &b_s__eVid);
   fChain->SetBranchAddress("s.eAid[2]", s__eAid, &b_s__eAid);
   fChain->SetBranchAddress("s.eFlag", s__eFlag, &b_s__eFlag);
   fChain->SetBranchAddress("s.eTrack", s__eTrack, &b_s__eTrack);
   fChain->SetBranchAddress("s.eX", s__eX, &b_s__eX);
   fChain->SetBranchAddress("s.eY", s__eY, &b_s__eY);
   fChain->SetBranchAddress("s.eZ", s__eZ, &b_s__eZ);
   fChain->SetBranchAddress("s.eTX", s__eTX, &b_s__eTX);
   fChain->SetBranchAddress("s.eTY", s__eTY, &b_s__eTY);
   fChain->SetBranchAddress("s.eSZ", s__eSZ, &b_s__eSZ);
   fChain->SetBranchAddress("s.eChi2", s__eChi2, &b_s__eChi2);
   fChain->SetBranchAddress("s.eProb", s__eProb, &b_s__eProb);
   fChain->SetBranchAddress("s.eW", s__eW, &b_s__eW);
   fChain->SetBranchAddress("s.eVolume", s__eVolume, &b_s__eVolume);
   fChain->SetBranchAddress("s.eDZ", s__eDZ, &b_s__eDZ);
   fChain->SetBranchAddress("s.eDZem", s__eDZem, &b_s__eDZem);
   fChain->SetBranchAddress("s.eP", s__eP, &b_s__eP);
   fChain->SetBranchAddress("s.eMCTrack", s__eMCTrack, &b_s__eMCTrack);
   fChain->SetBranchAddress("s.eMCEvt", s__eMCEvt, &b_s__eMCEvt);
   fChain->SetBranchAddress("s.eScanID.fUniqueID", s__eScanID_fUniqueID, &b_s__eScanID_fUniqueID);
   fChain->SetBranchAddress("s.eScanID.fBits", s__eScanID_fBits, &b_s__eScanID_fBits);
   fChain->SetBranchAddress("s.eScanID.eBrick", s__eScanID_eBrick, &b_s__eScanID_eBrick);
   fChain->SetBranchAddress("s.eScanID.ePlate", s__eScanID_ePlate, &b_s__eScanID_ePlate);
   fChain->SetBranchAddress("s.eScanID.eMajor", s__eScanID_eMajor, &b_s__eScanID_eMajor);
   fChain->SetBranchAddress("s.eScanID.eMinor", s__eScanID_eMinor, &b_s__eScanID_eMinor);
   fChain->SetBranchAddress("sf", &sf__, &b_sf__);
   fChain->SetBranchAddress("sf.fUniqueID", sf__fUniqueID, &b_sf__fUniqueID);
   fChain->SetBranchAddress("sf.fBits", sf__fBits, &b_sf__fBits);
   fChain->SetBranchAddress("sf.ePID", sf__ePID, &b_sf__ePID);
   fChain->SetBranchAddress("sf.eID", sf__eID, &b_sf__eID);
   fChain->SetBranchAddress("sf.eVid[2]", sf__eVid, &b_sf__eVid);
   fChain->SetBranchAddress("sf.eAid[2]", sf__eAid, &b_sf__eAid);
   fChain->SetBranchAddress("sf.eFlag", sf__eFlag, &b_sf__eFlag);
   fChain->SetBranchAddress("sf.eTrack", sf__eTrack, &b_sf__eTrack);
   fChain->SetBranchAddress("sf.eX", sf__eX, &b_sf__eX);
   fChain->SetBranchAddress("sf.eY", sf__eY, &b_sf__eY);
   fChain->SetBranchAddress("sf.eZ", sf__eZ, &b_sf__eZ);
   fChain->SetBranchAddress("sf.eTX", sf__eTX, &b_sf__eTX);
   fChain->SetBranchAddress("sf.eTY", sf__eTY, &b_sf__eTY);
   fChain->SetBranchAddress("sf.eSZ", sf__eSZ, &b_sf__eSZ);
   fChain->SetBranchAddress("sf.eChi2", sf__eChi2, &b_sf__eChi2);
   fChain->SetBranchAddress("sf.eProb", sf__eProb, &b_sf__eProb);
   fChain->SetBranchAddress("sf.eW", sf__eW, &b_sf__eW);
   fChain->SetBranchAddress("sf.eVolume", sf__eVolume, &b_sf__eVolume);
   fChain->SetBranchAddress("sf.eDZ", sf__eDZ, &b_sf__eDZ);
   fChain->SetBranchAddress("sf.eDZem", sf__eDZem, &b_sf__eDZem);
   fChain->SetBranchAddress("sf.eP", sf__eP, &b_sf__eP);
   fChain->SetBranchAddress("sf.eMCTrack", sf__eMCTrack, &b_sf__eMCTrack);
   fChain->SetBranchAddress("sf.eMCEvt", sf__eMCEvt, &b_sf__eMCEvt);
   fChain->SetBranchAddress("sf.eScanID.fUniqueID", sf__eScanID_fUniqueID, &b_sf__eScanID_fUniqueID);
   fChain->SetBranchAddress("sf.eScanID.fBits", sf__eScanID_fBits, &b_sf__eScanID_fBits);
   fChain->SetBranchAddress("sf.eScanID.eBrick", sf__eScanID_eBrick, &b_sf__eScanID_eBrick);
   fChain->SetBranchAddress("sf.eScanID.ePlate", sf__eScanID_ePlate, &b_sf__eScanID_ePlate);
   fChain->SetBranchAddress("sf.eScanID.eMajor", sf__eScanID_eMajor, &b_sf__eScanID_eMajor);
   fChain->SetBranchAddress("sf.eScanID.eMinor", sf__eScanID_eMinor, &b_sf__eScanID_eMinor);
   fChain->SetBranchAddress("TrackID", TrackID, &b_TrackID);
   fChain->SetBranchAddress("nseg", nseg, &b_nseg);
   fChain->SetBranchAddress("nholes", nholes, &b_nholes);
   fChain->SetBranchAddress("maxgap", maxgap, &b_maxgap);
   fChain->SetBranchAddress("incoming", incoming, &b_incoming);
   fChain->SetBranchAddress("impactparameter", impactparameter, &b_impactparameter);
   fChain->SetBranchAddress("MCEventID", MCEventID, &b_MCEventID);
   fChain->SetBranchAddress("MCTrackID", MCTrackID, &b_MCTrackID);
   fChain->SetBranchAddress("MCMotherID", MCMotherID, &b_MCMotherID);
  
}

/*
void tracks_reader_Fedra(TTree *tree)
{

   // Set object pointer
   //t_ = 0;
   //s_EdbTrack2D = 0;
   //sf_EdbTrack2D = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trid", &trid, &b_trid);
   fChain->SetBranchAddress("nseg", &tnseg, &b_tnseg);
   fChain->SetBranchAddress("npl", &npl, &b_npl);
   fChain->SetBranchAddress("n0", &n0, &b_n0);
   fChain->SetBranchAddress("xv", &xv, &b_xv);
   fChain->SetBranchAddress("yv", &yv, &b_yv);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("t.TObject.fUniqueID", &t_TObject_fUniqueID, &b_t_TObject_fUniqueID);
   fChain->SetBranchAddress("t.TObject.fBits", &t_TObject_fBits, &b_t_TObject_fBits);
   fChain->SetBranchAddress("t.ePID", &t_ePID, &b_t_ePID);
   fChain->SetBranchAddress("t.eID", &t_eID, &b_t_eID);
   fChain->SetBranchAddress("t.eVid[2]", t_eVid, &b_t_eVid);
   fChain->SetBranchAddress("t.eAid[2]", t_eAid, &b_t_eAid);
   fChain->SetBranchAddress("t.eFlag", &t_eFlag, &b_t_eFlag);
   fChain->SetBranchAddress("t.eTrack", &t_eTrack, &b_t_eTrack);
   fChain->SetBranchAddress("t.eX", &t_eX, &b_t_eX);
   fChain->SetBranchAddress("t.eY", &t_eY, &b_t_eY);
   fChain->SetBranchAddress("t.eZ", &t_eZ, &b_t_eZ);
   fChain->SetBranchAddress("t.eTX", &t_eTX, &b_t_eTX);
   fChain->SetBranchAddress("t.eTY", &t_eTY, &b_t_eTY);
   fChain->SetBranchAddress("t.eSZ", &t_eSZ, &b_t_eSZ);
   fChain->SetBranchAddress("t.eChi2", &t_eChi2, &b_t_eChi2);
   fChain->SetBranchAddress("t.eProb", &t_eProb, &b_t_eProb);
   fChain->SetBranchAddress("t.eW", &t_eW, &b_t_eW);
   fChain->SetBranchAddress("t.eVolume", &t_eVolume, &b_t_eVolume);
   fChain->SetBranchAddress("t.eDZ", &t_eDZ, &b_t_eDZ);
   fChain->SetBranchAddress("t.eDZem", &t_eDZem, &b_t_eDZem);
   fChain->SetBranchAddress("t.eP", &t_eP, &b_t_eP);
   fChain->SetBranchAddress("t.eMCTrack", &t_eMCTrack, &b_t_eMCTrack);
   fChain->SetBranchAddress("t.eMCEvt", &t_eMCEvt, &b_t_eMCEvt);
   fChain->SetBranchAddress("t.eScanID.fUniqueID", &t_eScanID_fUniqueID, &b_t_eScanID_fUniqueID);
   fChain->SetBranchAddress("t.eScanID.fBits", &t_eScanID_fBits, &b_t_eScanID_fBits);
   fChain->SetBranchAddress("t.eScanID.eBrick", &t_eScanID_eBrick, &b_t_eScanID_eBrick);
   fChain->SetBranchAddress("t.eScanID.ePlate", &t_eScanID_ePlate, &b_t_eScanID_ePlate);
   fChain->SetBranchAddress("t.eScanID.eMajor", &t_eScanID_eMajor, &b_t_eScanID_eMajor);
   fChain->SetBranchAddress("t.eScanID.eMinor", &t_eScanID_eMinor, &b_t_eScanID_eMinor);
   fChain->SetBranchAddress("s", &s_, &b_s_);
   fChain->SetBranchAddress("s.fUniqueID", s_fUniqueID, &b_s_fUniqueID);
   fChain->SetBranchAddress("s.fBits", s_fBits, &b_s_fBits);
   fChain->SetBranchAddress("s.ePID", s_ePID, &b_s_ePID);
   fChain->SetBranchAddress("s.eID", s_eID, &b_s_eID);
   fChain->SetBranchAddress("s.eVid[2]", s_eVid, &b_s_eVid);
   fChain->SetBranchAddress("s.eAid[2]", s_eAid, &b_s_eAid);
   fChain->SetBranchAddress("s.eFlag", s_eFlag, &b_s_eFlag);
   fChain->SetBranchAddress("s.eTrack", s_eTrack, &b_s_eTrack);
   fChain->SetBranchAddress("s.eX", s_eX, &b_s_eX);
   fChain->SetBranchAddress("s.eY", s_eY, &b_s_eY);
   fChain->SetBranchAddress("s.eZ", s_eZ, &b_s_eZ);
   fChain->SetBranchAddress("s.eTX", s_eTX, &b_s_eTX);
   fChain->SetBranchAddress("s.eTY", s_eTY, &b_s_eTY);
   fChain->SetBranchAddress("s.eSZ", s_eSZ, &b_s_eSZ);
   fChain->SetBranchAddress("s.eChi2", s_eChi2, &b_s_eChi2);
   fChain->SetBranchAddress("s.eProb", s_eProb, &b_s_eProb);
   fChain->SetBranchAddress("s.eW", s_eW, &b_s_eW);
   fChain->SetBranchAddress("s.eVolume", s_eVolume, &b_s_eVolume);
   fChain->SetBranchAddress("s.eDZ", s_eDZ, &b_s_eDZ);
   fChain->SetBranchAddress("s.eDZem", s_eDZem, &b_s_eDZem);
   fChain->SetBranchAddress("s.eP", s_eP, &b_s_eP);
   fChain->SetBranchAddress("s.eMCTrack", s_eMCTrack, &b_s_eMCTrack);
   fChain->SetBranchAddress("s.eMCEvt", s_eMCEvt, &b_s_eMCEvt);
   fChain->SetBranchAddress("s.eScanID.fUniqueID", s_eScanID_fUniqueID, &b_s_eScanID_fUniqueID);
   fChain->SetBranchAddress("s.eScanID.fBits", s_eScanID_fBits, &b_s_eScanID_fBits);
   fChain->SetBranchAddress("s.eScanID.eBrick", s_eScanID_eBrick, &b_s_eScanID_eBrick);
   fChain->SetBranchAddress("s.eScanID.ePlate", s_eScanID_ePlate, &b_s_eScanID_ePlate);
   fChain->SetBranchAddress("s.eScanID.eMajor", s_eScanID_eMajor, &b_s_eScanID_eMajor);
   fChain->SetBranchAddress("s.eScanID.eMinor", s_eScanID_eMinor, &b_s_eScanID_eMinor);
   fChain->SetBranchAddress("sf", &sf_, &b_sf_);
   fChain->SetBranchAddress("sf.fUniqueID", sf_fUniqueID, &b_sf_fUniqueID);
   fChain->SetBranchAddress("sf.fBits", sf_fBits, &b_sf_fBits);
   fChain->SetBranchAddress("sf.ePID", sf_ePID, &b_sf_ePID);
   fChain->SetBranchAddress("sf.eID", sf_eID, &b_sf_eID);
   fChain->SetBranchAddress("sf.eVid[2]", sf_eVid, &b_sf_eVid);
   fChain->SetBranchAddress("sf.eAid[2]", sf_eAid, &b_sf_eAid);
   fChain->SetBranchAddress("sf.eFlag", sf_eFlag, &b_sf_eFlag);
   fChain->SetBranchAddress("sf.eTrack", sf_eTrack, &b_sf_eTrack);
   fChain->SetBranchAddress("sf.eX", sf_eX, &b_sf_eX);
   fChain->SetBranchAddress("sf.eY", sf_eY, &b_sf_eY);
   fChain->SetBranchAddress("sf.eZ", sf_eZ, &b_sf_eZ);
   fChain->SetBranchAddress("sf.eTX", sf_eTX, &b_sf_eTX);
   fChain->SetBranchAddress("sf.eTY", sf_eTY, &b_sf_eTY);
   fChain->SetBranchAddress("sf.eSZ", sf_eSZ, &b_sf_eSZ);
   fChain->SetBranchAddress("sf.eChi2", sf_eChi2, &b_sf_eChi2);
   fChain->SetBranchAddress("sf.eProb", sf_eProb, &b_sf_eProb);
   fChain->SetBranchAddress("sf.eW", sf_eW, &b_sf_eW);
   fChain->SetBranchAddress("sf.eVolume", sf_eVolume, &b_sf_eVolume);
   fChain->SetBranchAddress("sf.eDZ", sf_eDZ, &b_sf_eDZ);
   fChain->SetBranchAddress("sf.eDZem", sf_eDZem, &b_sf_eDZem);
   fChain->SetBranchAddress("sf.eP", sf_eP, &b_sf_eP);
   fChain->SetBranchAddress("sf.eMCTrack", sf_eMCTrack, &b_sf_eMCTrack);
   fChain->SetBranchAddress("sf.eMCEvt", sf_eMCEvt, &b_sf_eMCEvt);
   fChain->SetBranchAddress("sf.eScanID.fUniqueID", sf_eScanID_fUniqueID, &b_sf_eScanID_fUniqueID);
   fChain->SetBranchAddress("sf.eScanID.fBits", sf_eScanID_fBits, &b_sf_eScanID_fBits);
   fChain->SetBranchAddress("sf.eScanID.eBrick", sf_eScanID_eBrick, &b_sf_eScanID_eBrick);
   fChain->SetBranchAddress("sf.eScanID.ePlate", sf_eScanID_ePlate, &b_sf_eScanID_ePlate);
   fChain->SetBranchAddress("sf.eScanID.eMajor", sf_eScanID_eMajor, &b_sf_eScanID_eMajor);
   fChain->SetBranchAddress("sf.eScanID.eMinor", sf_eScanID_eMinor, &b_sf_eScanID_eMinor);
}
*/
void bdt_reader(TTree *tree)
{
   
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bdt_value", &bdt_value, &b_bdt_value);
}
