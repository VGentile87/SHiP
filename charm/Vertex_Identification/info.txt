CREATED BY V.GENTILE 2020/10/02

Script name: vtx_MC_analysis.C 
Description: Code for the identification of vertices reconstructed by FEDRA using a MC simulation.
	     The code provides most probable event and mother ID for each vertex, the number of 
	     electrons, primaries and charm daughters. A flag, charm_daug tells about the presence
   	     of a charm daughter in a vertex.

List of variables:
	1) mp_eventID  --> Most probable montecarlo id of the vertex
	2) mp_motherID  --> Most probable montecarlo mother id of the vertex
	2) mp_pdgID  --> Most probable pdg id of the vertex
	3) charm_daug  --> If 0 no charm daughter tracks in the vertex, if 1 a least one charm daughter in the vertex
	4) n_electrons  --> Number of electrons in the vertex (even if the montecarlo id is different from mp_eventID)
	5) n_primaries  --> Number of primary daughter tracks (mother id of the track is -1)
	6) n_charmdaug  --> Number of charm daughter tracks (mother id of the track is 1 or 2)
	7) n_events --> Number of montecarlo events linked to the vertex (ideal case == 1)	
	8) n_mothers --> Number of mother IDs linked to the vertex (ideal case == 1)
