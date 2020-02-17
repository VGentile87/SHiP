// DECAY SEARCH 

	1) Creare directory di lavoro
	
	2) Copiare i seguenti file:
		- linked_tracks.root
		- vertextree.root
		- DSCuts.h
		- Definitions.h
		- vtx_data_study.C
		- selection_decaysearch_sim.C 	
		- VertexTrackDisplay.py   (opzionale)
	
	3) Verificare in Definitions.h:
		- kMaxt e tkMaxt
		- kMaxs e tkMaxs
		- kMaxsf e tkMaxsf
		(N.b. puoi usare MakeClass su linked_tracks e vertextree_test per verificare i numeri)

	4) Verificare i tagli in DSCuts.h

	5) Creare una sottocartella TMVA e copiare i seguenti file:

		- dataset1st/
		- dataset2nd/
		- datasetDS/
		- TMVAClassification.C
		- TMVAClassificationApplication.C
		- TMVAClassificationApplication2nd.C
		- TMVAClassificationApplicationDS.C
		- reader_evaluation_vtx.C


	6) Comandi per la BDT vtx interazione:
		
		- root -l
		[] .L TMVAClassification.C
		[] prepareTMVAtree()

		fornisce in output "tmva_input_vertices.root"

		- root -l TMVAClassificationApplication.C
		fornisce in output "vtx_BDT_data_evaluated.root"

	7) Comandi per la BDT vtx decadimento:

		- root -l
		[] .L TMVAClassification.C
		[] prepareTMVAtree2nd()

		fornisce in output "tmva_input_vertices2nd.root"

		- root -l TMVAClassificationApplication2nd.C
		fornisce in output "vtx_BDT_data_evaluated2nd.root"

	8) Copiare i file "evaluated" nella cartella principale

		- cp vtx_BDT_data_evaluated* ../

	9) Lanciare lo script per la Decay Search:
	
		- verificare in "vertextree_test.root" il valore massimo di MCEvent che serve come input a numero di simulazioni
		- root -l
		[] .L vtx_data_study.C+
		[] myrun()
		[] opzione 2
		[] inserire il max di MCEvent come "eventi nella simulazione"

		fornisce in output:
		- "ds_data_result.root"
		- "log_info_ds.txt"
		- "log_ds_search.txt"

	10) Lanciare lo script che aggiunge nuovi branch al file di DS:

		- root -l selection_decaysearch_sim.C 
		fornisce in output "annotated_data_result.root"

	11) Creare il file di input per la BDT di Decay search:

		cp annotated_data_result.root TMVA/
		- root -l 
		[] .L TMVAClassification.C
		[] prepareTMVAtreeDS()

		fornisce in output "tmva_input_verticesDS.root"

	12) Lanciare script BDT per la Decay Search:
		
		- root -l MVAClassificationApplicationDS.C
		fornisce in output "vtx_BDT_DS_evaluated.root"

	13) Lanciare script per efficienze e stime numeriche:

		- root -l reader_evaluation_vtx.C
		- settare numero vertici e/o eventi simulati

		fornisce in uscita:
			
		"log_DS_results.txt"  (efficienze e conteggi)
		"vtx_list.txt"  (lista dei vertici sopravvissuti)
				


