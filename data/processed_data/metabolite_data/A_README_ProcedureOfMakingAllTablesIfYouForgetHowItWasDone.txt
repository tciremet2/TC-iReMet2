1) calculate ratios from standart tables in excel (mainfile: 'metabolitedata_ratios_mean.ods') and save them by hand from mainfile in 
   seperate files called
	MTWT_ratio_mets_0d1d7d21d_#.csv 

2) use the ratios (from 1) to calculate mean and sd via matlab (script: 'CP_metabolite_ratio_calc_meanAndSd.m') and sort them so they
   match the dictionary
	MTWT_ratios_mets_experiment#_mean_sdev_sorted.txt
	MTWT_ratio_mets_dictionary_experiment#_sorted.txt

3) by hand produce the final tables (used tables are the ones from 2.) that can be used to put in the iReMet-Flux - THIS IS FINAL TABLE 
	WTMT_ratios_mets_experiment#_mean_sdev_sorted_fitforreadin.txt
