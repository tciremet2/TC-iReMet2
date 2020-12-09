%%%% CP run eflux on Transcriptomic data table 

%%% read data 
tdata = readtable('/home/mpimp-golm.mpg.de/pries1347/Desktop/winhome/main/MT/data_final_tables/transcription_data/Transcriptomic_data_aggregated.txt','Delimiter',',');

%%% init
atgcodelist = table2array(tdata(:,1)); % atg codes
