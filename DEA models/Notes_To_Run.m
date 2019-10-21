%% Notes to the Grader/Reader

% Writers:

% Daghan Kendirli
% Ruitong Li
% Monica Le
% Christopher Loynes
% Olivier Kraaijeveld
% Giorgios Makridakis

% for any queries, email daglink@hotmail.com


% Use Files starting with T!!!
% A_DEA_Original is the original file found online
% B_BBC_IO_OO and C_CCR_IO_OO are the practice files

% T1 stands for Task 1, T2 for Task 2 and so on.

%%%%%%%% All Models are separated in two files: IO and OO %%%%%%%%

     %%%% Run T files separately => Total files to run = 8 %%%%
% e.g. 2 files for T1, 2 files for T2, 2 files for T3 and 2 files for T4

% Each file will output different tables in Excel or Table formats:

%% TASK 1: CCR-IO and OO

% Visualization Table will show EVERYTHING EXCEPT WEIGHTS:
% inputs, outputs, slacks, projections, percentage increase and decrease 
% and DEA score

% Other tables just separate all measures separately:

% Statistical Measures Table shows the min, max, mean, std and interquartile
% of measures

% Slacks and DEA scores Table shows the slacks and DEA scores

% Projection Table shows the projections

% Percentage increase and decrease Table

% FOR WEIGHTS: LOOK AT THE DEA_T1_IO_Output_Data_File.table file (txt)
   %   For OO:             DEA_T1_OO_Output_Data_File.table file (txt)

% Note that the data in workspace will be cleared each time a new file is run!

%% TASK 2: DUAL CCR-IO and OO

% Change epsilon to wished value (e.g. 1E-6 or 1E-9)

% For Table showing inputs, outputs, weights and DEA score:
% LOOK AT      DEA_T2_IO_Output_Data_File.table file Excel
   %  For OO:  DEA_T1_OO_Output_Data_File.table file Excel

% For statistics table of the weights: DEA_T2_IO_Statistics_Weights.xlsx
   %  For OO:                          DEA_T2_OO_Statistics_Weights.xlsx

%% TASK 3: BCC-IO and OO

% Visualization Table will show EVERYTHING EXCEPT WEIGHTS:
% inputs, outputs, slacks, projections, percentage increase and decrease 
% and DEA score

% Other tables just separate all measures separately:

% Statistical Measures Table shows the min, max, mean, std and interquartile
% of measures

% Slacks and DEA scores Table shows the slacks and DEA scores

% Projection Table shows the projections

% Percentage increase and decrease Table

% FOR WEIGHTS: LOOK AT THE DEA_T3_IO_Output_Data_File.table file (txt)
   %   For OO:             DEA_T3_OO_Output_Data_File.table file (txt)

% Note that the data in workspace will be cleared each time a new file is run!

%% TASK 4: DUAL BCC-IO and OO

% No code


%% T-extra

% Can run file for correlations and T-test after running all T1 AND T3
% files
% Can run a multi dimensional scaling efficiency frontier graph
