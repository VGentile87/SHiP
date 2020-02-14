Created by V.Gentile (2020/11/02)

Script name: TMVAClassification_MC.C  (only for training and test of simulation data)
Description: This code provides several functions to create input files for the training and test of the BDTs;
	     and a main function to perform the BDTs.

List of functions:
	1) prepareTMVAtree() --> It requires a vertex_file with MC information (mp_motherID, mp_eventID, charm_daug, etc...)
				 The file in input is usually called "vtx_MC_analysis.root"
				 It provides in output the file "tmva_input_vertices.root" that is used as input for the BDT on
 				 interaction vertices.

	2) prepareTMVAtree2nd() --> It requires a vertex_file with MC information (mp_motherID, mp_eventID, charm_daug, etc...) and a file
                                    with the BDT values evaluated on interaction vertices.
				    The file in input are usually called "vtx_MC_analysis.root" and "vtx_BDT_data_evaluated.root"
				    It provides in output the file "tmva_input_vertices2nd.root" that is used as input for the BDT on
			            decay vertices.

	3) prepareTMVAtreeDS() -->  It requires a file where the decay search is already performed.
				    The file in input is usually called "annotated_data_result.root".
				    It provides in output the file "tmva_input_verticesDS.root" that is used as input for the decay search 
				    BDT.

	4) TMVAClassification_MC() -->  The script allows to perform the BDT on interaction vertices, decay vertices and the BDT for the
				        decay search.
					Different cases are possible:
					Case 1: BDT on interaction vertices --> Input files like "tmva_input_vertices.root".
					Case 2: BDT on decay vertices --> 2 Input files like "tmva_input_vertices2nd.root" for the signal 
					(charm vertices) and the background (the remaining part).
					Case 3: BDT for the decay search --> Input files like "tmva_input_verticesDS.root".
					It provides in output the directories with the datasets with the weigths of the tested network and 
					the "TMVA.root" to check the BDT results and distributions.


Script name: TMVAClassification.C   (only for data and simulation evaluation)
Description: This code provides several functions to create input files for the evaluation of the BDTs;

List of functions:
	1) prepareTMVAtree() --> It requires a vertex file usually called "vertextree_test.root"
				 It provides in output the file "tmva_input_vertices.root" that is used as input for the BDT evaluation of
 				 interaction vertices usually called "vtx_BDT_data_evaluated.root".

	2) prepareTMVAtree2nd() --> It requires a vertex_file and a file with the BDT values evaluated on interaction vertices 
				    The file in input are usually called "vertextree_test.root" and "vtx_BDT_data_evaluated.root"
				    It provides in output the file "tmva_input_vertices2nd.root" that is used as input for the BDT
				    evaluation on decay vertices.

	3) prepareTMVAtreeDS() -->  It requires a file where the decay search is already performed.
				    The file in input is usually called "annotated_data_result.root".
				    It provides in output the file "tmva_input_verticesDS.root" that is used as input for the decay search 
				    BDT evaluation.


Usage:

# input files
root -l
.L TMVAClassification_MC.C or (.L TMVAClassification.C)
prepareTMVAtree()
prepareTMVAtree2nd()
prepareTMVAtreeDS()

# BDTs training and test
root -l ./TMVAClassification_MC.C\(\"BDT,Likelihood\"\)

# BDT results
root -l
TMVA::TMVAGui("TMVA.root")
	     
