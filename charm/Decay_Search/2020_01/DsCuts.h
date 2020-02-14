// Cuts used in the Decay Search //
//
// 
//
//-------------------------------//

int cut_nseg=6;                // min number of segments
float cut_prg1=5;              // min normalized sigma 
float bdt_cut1= 0.15;          // bdt cut on primaries
float bdt_cut2=-0.22;          // bdt cut on secondaries
float cut_lip=50;              // min cut on large ip
float cut_cip=0;               // min cut on ip for tracks in clusters
int cut_cls_ntrk=5;            // max number of tracks in a cluster
float cut_max_dz=10000;        // max distance between secondaries and primaries
float cut_min_dz=0;            // min distance between secondaries and primaries 
int cut_vtx2_ntrk=5;           // max number of tracks for a secondary vertex
float cut_vtx2_ka=0.2;         // max kink angle of charm from primary
float cut_max_ip=1000;         // max impact parameter
float cut_extrk_max_dz=10000;  // max distance between secondaries and primaries
float cut_extrk_min_dz=0;      // min distance between secondaries and primaries


int cut_ntrk=5;                // min number of tracks for interaction vertices
float cut_vz=-4000;            // vz upper cut for interaction vertices
float cut2_vz=-3000;            // vz upper cut for decay vertices
  
// in clustering function
// dist_rel = 20;              // distance between (x,y) proj of first segment of two tracks
// epsilon = 10;               // distance between intersection points;
// epsilon_phi = 0.1;          // angular distance between intersection points;
