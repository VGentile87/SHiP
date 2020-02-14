Created by V.Gentile (2020/02/14)

The decay search "DS" is performed on vertex and track files processed by FEDRA.
Usually these files are called "vertextree_test.root" and "linked_tracks.root"

Script name: Definitions.h 
Description: Header file where track and vertex variables and trees are declared.

Script name: DsCuts.h
Description: Datacard where cuts for the DS are defined (some cuts depends on dataset and should be set here)

Script name: vtx_data_study.C 
Description: Code for the DS which requires in input the following files:
		1) Definitions.h
		2) DsCuts.h
		3) vertextree_test.root
		4) linked_tracks.root
		5) vtx_BDT_data_evaluated.root (see README_BDT.txt)
		6) vtx_BDT_data_evaluated2nd.root (see README_BDT.txt)

	      It provides in output the following files:
		1) log_info_ds.txt  (preliminary info before the DS)
		2) log_ds_search.txt (report of the decay search)
		3) ds_data_result.root (output tree with the DS results)

Script name: selection_decaysearch_sim.C (created by A.Iuliano)
Description: The script provides additional selections to the DS tree.
	     It requires as input file "ds_data_result.root" and provides
	     in output the file "annotated_data_result.root"



Usage:

# DS run
root -l
.L vtx_data_study.C
myrun()
(type 2 for the decay search)

# new selections
root -l selection_decaysearch_sim.C
	     
