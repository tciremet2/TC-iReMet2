# Table of contents
* [TC-iReMet2](#TC-iReMet2)
* [Technologies](#Technologies)
* [Setup](#Setup)
* [Metabolite-Ratios](#Metabolite-Ratios)
* [Transcriptomic-Ratios](#Transcriptomic-Ratios)

## TC-iReMet2
Code-Repository for [TC-iReMet2](https://rdcu.be/ctIYt).

## Technologies
Project is created with:
* Matlab2018
* R (partly used for analysis)

## Setup
To run this project, go to the TC-iReMet2_algorithm folder and run nMAIN() in the Matlab console.
* If you wish to run TC-iReMet2 with different parameters, do so by changing the 'AnalysisParameters.init' file
* Paths might need to be adjusted based on your chosen folder structure
* Files neccessary to run TC-iReMet2 are located in the 'Data' folder
* Results are saved to the 'Results' folder with a following folder structure depending on the chosen run parameters
* [Metabolite](#Metabolite-Ratios)- and [Transcriptomic](#Transcriptomic-Ratios)-Ratios were calculated seperately

### Metabolite-Ratios
For the creation of the metabolite ratios see '/data/processed_data/metabolite_data'.
* See 'A_README_ProcedureOfMakingAllTables.txt' for further details

### Transcriptomic-Ratios
For the creation of the transcriptomic ratios see '/data/processed_data/transcription_data'
* See '01_README.txt' for further details.
* 'createefluxtable.m' creates and calculates the transcriptomic ratio constraints