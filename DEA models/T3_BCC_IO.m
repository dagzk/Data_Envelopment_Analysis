%% PAMP DEA BCC-IO File;
clear all;

%% Task 3 BCC-IO

%% Model selection;
model='BCC';               % Choice of model    model= 'BCC' ;
orientation='io';          % orientation = 'io' or 'oo' (for input or output oriented)    ;

epsilon= 1E-6;              % epsilon (non archimedian number for the BCC & CCR models)    ;

%% Enters the inputs and outputs of all DMUs ;
%  Alternatively, use the "putmatrix" function in Excel to generate X and Y;


X1 = xlsread('DEA_Input_Data_File_2018.xls', 'B:B');

X = X1(:, 2:3);
Y = X1(:, 4:5);

% extracts the number of DMUs, inputs and outputs;
[n,m] = size(X);
[n,s] = size(Y);

%% Computes the results from the selected model;
%  The results from the selected model are in matrix Z;

switch model;
    
    % BCC model;
    case ('BCC')
    switch orientation;
        % Input oriented; THIS IS USED IN THIS CASE
        case ('io')
            Z = zeros(n,n+m+s+1);
            
            % Objective function of the BCC model: min(0*lambda - epsilon*(s+ + s-) + theta);
            f = [zeros(1,n) -epsilon*ones(1,s+m) 1];
            
            lblambda = zeros(n,1);                % Lower bounds for (n) lambdas;
            lboutput = zeros(s,1);                % Lower bounds for (s) outputs;
            lbinput = zeros(m,1);                 % Lower bounds for (m) inputs ;
            lb = [lblambda; lboutput; lbinput];   % Lower bounds for lambdas, outputs (s+) and inputs (s-);
            for j=1:n
                Aeq = [Y', -eye(s,s), zeros(s,m+1);
                      -X', zeros(m,s), -eye(m,m) X(j,:)';
                      ones(1,n), zeros(1,s), zeros(1,m+1)];
                beq = [Y(j,:)';zeros(m,1);1];
                z = linprog(f,[],[],Aeq,beq,lb);
                Z(j,:) = z;
            end
            Z

        % Output oriented;
        case ('oo')
            Z = zeros(n,n+m+s+1);

            % Objective function of the BCC_oo model: max(0*lambda + epsilon*(s+ + s-) + theta);
            f = -[zeros(1,n), epsilon*ones(1,s+m), 1];

            lblambda = zeros(n,1);                % Lower bounds for (n) lambdas;
            lboutput = zeros(s,1);                % Lower bounds for (s) outputs;
            lbinput = zeros(m,1);                 % Lower bounds for (m) inputs ;
            lb = [lblambda; lboutput; lbinput];   % Lower bounds for lambdas, outputs (s+) and inputs (s-);
            for j=1:n
                Aeq = [-Y', eye(s,s), zeros(s,m), Y(j,:)'; 
                        X', zeros(m,s), eye(m,m), zeros(m,1); 
                        ones(1,n), zeros(1,s+m+1)];
                beq = [zeros(s,1);X(j,:)';1];
                z = linprog(f,[],[],Aeq,beq,lb);
                Z(j,:) = z;
            end
            Z
            
    end
end

%% Analysis 

% Once we have the solution matrix Z
% Create a table in excel with the main findings of the analysis
% Create individual vectors that contain: DEA score(theta), slack x1,slack
% x2, slack y1, slack y2 and one matrix with all the densities for each
% DMUs
 
DEA_score_vector_T3_BBCIO = Z(:,n+s+m+1);
x1_input=X(:,1);
x2_input=X(:,2);
y1_input=Y(:,1);
y2_input=Y(:,2);
slack_x1_vector_BBCIO=Z(:,51);
slack_x2_vector_BBCIO=Z(:,52);
slack_y1_vector_BBCIO=Z(:,53);
slack_y2_vector_BBCIO=Z(:,54);
intensity_vector_BBCIO=Z(:,1:50);
 
% Proceed with the necessary calculation to calculate the sum of slacks, 
% projections of the inefficient DMUs, the percentage of input and output
% that need to be decreased or increased.
 
sum_of_x_slacks_BBCIO=slack_x1_vector_BBCIO+slack_x2_vector_BBCIO;
sum_of_y_slacks_BBCIO=slack_y1_vector_BBCIO+slack_y2_vector_BBCIO;
projections_x1_BBCIO= DEA_score_vector_T3_BBCIO.*x1_input-(slack_x1_vector_BBCIO);
projections_x2_BBCIO= DEA_score_vector_T3_BBCIO.*x2_input-(slack_x2_vector_BBCIO);           
projections_y1_BBCIO= y1_input+(slack_y1_vector_BBCIO);
projections_y2_BBCIO= y2_input+(slack_y2_vector_BBCIO);
 
% percentage_to_decrease/increase: calculate how far in percentage the
% current input is from the projection in the efficiency frontier.
 
percentage_to_decrease_x1_BBCIO=(x1_input-projections_x1_BBCIO)./x1_input;
percentage_to_decrease_x2_BBCIO=(x2_input-projections_x2_BBCIO)./x2_input;
percentage_to_decrease_x2_BBCIO(isnan(percentage_to_decrease_x2_BBCIO))=0
percentage_to_increase_y1_BBCIO=(projections_y1_BBCIO-y1_input)./y1_input;
percentage_to_increase_y1_BBCIO(isnan(percentage_to_increase_y1_BBCIO))=0
percentage_to_increase_y2_BBCIO=(projections_y2_BBCIO-y2_input)./y2_input;
percentage_to_increase_y2_BBCIO(isnan(percentage_to_increase_y2_BBCIO))=0

% Create tables

% Creation of a vector with the names of the DMUs 
temp1 = 'DMU_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
nameDMU = [temp3, temp4];
clear temp1 temp2 temp3 temp4
 
% Creation of the table for Slacks and DEA scores
table_analysis_BBCIO=table(nameDMU,x1_input,x2_input,y1_input,y2_input,slack_x1_vector_BBCIO,... 
    slack_x2_vector_BBCIO,slack_y1_vector_BBCIO,slack_y2_vector_BBCIO,...
    sum_of_x_slacks_BBCIO,sum_of_y_slacks_BBCIO,DEA_score_vector_T3_BBCIO);
filename = 'BBCIO_Slacks_and_DEA_scores.xlsx';
writetable(table_analysis_BBCIO,filename)
 
% Creation of the table Percentage inc and dec
table_analysis_percentages_BBCIO=table(nameDMU,percentage_to_decrease_x1_BBCIO, ...
    percentage_to_decrease_x2_BBCIO,percentage_to_increase_y1_BBCIO,...
    percentage_to_increase_y2_BBCIO,DEA_score_vector_T3_BBCIO)
filename = 'BBCIO_Percentages_increase_decrease.xlsx';
writetable(table_analysis_percentages_BBCIO,filename)
 
% Creation of Projections table
table_analysis_projections_efficient_frontier_BBCIO=table(nameDMU,...
    x1_input,x2_input,y1_input,y2_input,...
    DEA_score_vector_T3_BBCIO,projections_x1_BBCIO,projections_x2_BBCIO...
    ,projections_y1_BBCIO,projections_y2_BBCIO);
filename = 'BBCIO_Projection_efficient_frontier_BBCIO.xlsx';
writetable(table_analysis_projections_efficient_frontier_BBCIO,filename)
 
 
% The outputs are calculated, we use statistical measures to find 
% the meaning and patterns in the outcomes of the DEA
 
% statistical measures for slack_x1
Min_slack_x1_BBCIO=min(slack_x1_vector_BBCIO(slack_x1_vector_BBCIO>0)); 
Max_slack_x1_BBCIO=max(slack_x1_vector_BBCIO);
Mean_slack_x1_BBCIO=mean(slack_x1_vector_BBCIO);
Standard_deviation_slack_x1_BBCIO=std(slack_x1_vector_BBCIO);
interquartile_range_slack_x1_BBCIO=iqr(slack_x1_vector_BBCIO);
 
Stats_slack_x1_BBCIO=[Min_slack_x1_BBCIO,...
Max_slack_x1_BBCIO,Mean_slack_x1_BBCIO,...
    Standard_deviation_slack_x1_BBCIO,...
    interquartile_range_slack_x1_BBCIO];
Stats_slack_x1_BBCIO=Stats_slack_x1_BBCIO.'
 
% statistical measures for slack_x2
Min_slack_x2_BBCIO=min(slack_x2_vector_BBCIO(slack_x2_vector_BBCIO>0)); 
Max_slack_x2_BBCIO=max(slack_x2_vector_BBCIO);
Mean_slack_x2_BBCIO=mean(slack_x2_vector_BBCIO);
Standard_deviation_slack_x2_BBCIO=std(slack_x2_vector_BBCIO);
interquartile_range_slack_x2_BBCIO=iqr(slack_x2_vector_BBCIO);
 
Stats_slack_x2_BBCIO=[Min_slack_x2_BBCIO,...
Max_slack_x2_BBCIO,Mean_slack_x2_BBCIO,...
    Standard_deviation_slack_x2_BBCIO,...
    interquartile_range_slack_x2_BBCIO];
Stats_slack_x2_BBCIO=Stats_slack_x2_BBCIO.'
 
% statistical measures for slack_y1
Min_slack_y1_BBCIO=min(slack_y1_vector_BBCIO); 
Max_slack_y1_BBCIO=max(slack_y1_vector_BBCIO);
Mean_slack_y1_BBCIO=mean(slack_y1_vector_BBCIO);
Standard_deviation_slack_y1_BBCIO=std(slack_y1_vector_BBCIO);
interquartile_range_slack_y1_BBCIO=iqr(slack_y1_vector_BBCIO);
 
Stats_slack_y1_BBCIO=[Min_slack_y1_BBCIO,...
Max_slack_y1_BBCIO,...
Mean_slack_y1_BBCIO,Standard_deviation_slack_y1_BBCIO,...
interquartile_range_slack_y1_BBCIO];
 
Stats_slack_y1_BBCIO=Stats_slack_y1_BBCIO.'
 
% statistical measures for slack_y2
Min_slack_y2_BBCIO=min(slack_y2_vector_BBCIO(slack_y2_vector_BBCIO>0)); 
Max_slack_y2_BBCIO=max(slack_y2_vector_BBCIO);
Mean_slack_y2_BBCIO=mean(slack_y2_vector_BBCIO);
Standard_deviation_slack_y2_BBCIO=std(slack_y2_vector_BBCIO);
interquartile_range_slack_y2_BBCIO=iqr(slack_y2_vector_BBCIO);
 
Stats_slack_y2_BBCIO=[Min_slack_y2_BBCIO,...
Max_slack_y2_BBCIO,Mean_slack_y2_BBCIO,...
    Standard_deviation_slack_y2_BBCIO,...
    interquartile_range_slack_y2_BBCIO];
 
Stats_slack_y2_BBCIO=Stats_slack_y2_BBCIO.'
 
% statistical measures for DEA score 
Min_DEA_Score_BBCIO=min(DEA_score_vector_T3_BBCIO(DEA_score_vector_T3_BBCIO>0)); 
Max_DEA_Score_BBCIO=max(DEA_score_vector_T3_BBCIO);
Mean_DEA_Score_BBCIO=mean(DEA_score_vector_T3_BBCIO);
Standard_deviation_DEA_Score_BBCIO=std(DEA_score_vector_T3_BBCIO);
interquartile_range_DEA_Score_BBCIO=iqr(DEA_score_vector_T3_BBCIO);
 
Stats_DEA_Score_BBCIO=[Min_DEA_Score_BBCIO,...
Max_DEA_Score_BBCIO,Mean_DEA_Score_BBCIO,...
    Standard_deviation_DEA_Score_BBCIO,...
    interquartile_range_DEA_Score_BBCIO];
Stats_DEA_Score_BBCIO=Stats_DEA_Score_BBCIO.'
 
% statistical measures for percentage_to_decrease_x1_BCCIO
 
Min_percentage_to_decrease_x1_BBCIO=min(percentage_to_decrease_x1_BBCIO);
Max_percentage_to_decrease_x1_BBCIO=max(percentage_to_decrease_x1_BBCIO);
Mean_percentage_to_decrease_x1_BBCIO=mean(percentage_to_decrease_x1_BBCIO);
Standard_deviation_percentage_to_decrease_x1_BBCIO=std(percentage_to_decrease_x1_BBCIO);
interquartile_range_percentage_to_decrease_x1_BBCIO=iqr(percentage_to_decrease_x1_BBCIO);
 
Stats_percentage_to_decrease_x1_BBCIO=[Min_percentage_to_decrease_x1_BBCIO,...
   Max_percentage_to_decrease_x1_BBCIO,Mean_percentage_to_decrease_x1_BBCIO,...
    Standard_deviation_percentage_to_decrease_x1_BBCIO,...
    interquartile_range_percentage_to_decrease_x1_BBCIO];
Stats_percentage_to_decrease_x1_BBCIO=Stats_percentage_to_decrease_x1_BBCIO.'
 
% statistical measures for percentage_to_decrease_x2_BCCIO
Min_percentage_to_decrease_x2_BBCIO=min(percentage_to_decrease_x2_BBCIO); 
Max_percentage_to_decrease_x2_BBCIO=max(percentage_to_decrease_x2_BBCIO); 
Mean_percentage_to_decrease_x2_BBCIO=mean(percentage_to_decrease_x2_BBCIO); 
Standard_deviation_percentage_to_decrease_x2_BBCIO=std(percentage_to_decrease_x2_BBCIO); 
interquartile_range_percentage_to_decrease_x2_BBCIO=iqr(percentage_to_decrease_x2_BBCIO); 
 
Stats_percentage_to_decrease_x2_BBCIO=[Min_percentage_to_decrease_x2_BBCIO,...
    Max_percentage_to_decrease_x2_BBCIO,Mean_percentage_to_decrease_x2_BBCIO,...
    Standard_deviation_percentage_to_decrease_x2_BBCIO,...
   interquartile_range_percentage_to_decrease_x2_BBCIO];
Stats_percentage_to_decrease_x2_BBCIO=Stats_percentage_to_decrease_x2_BBCIO.'
% statistical measures for percentage_to_increase_y1_BBCIO 
Min_percentage_to_increase_y1_BBCIO=min(percentage_to_increase_y1_BBCIO); 
Max_percentage_to_increase_y1_BBCIO=max(percentage_to_increase_y1_BBCIO);
Mean_percentage_to_increase_y1_BBCIO=mean(percentage_to_increase_y1_BBCIO);
Standard_deviation_percentage_to_increase_y1_BBCIO=std(percentage_to_increase_y1_BBCIO);
interquartile_range_percentage_to_increase_y1_BBCIO=iqr(percentage_to_increase_y1_BBCIO);
 
Stats_percentage_to_increase_y1_BBCIO=[Min_percentage_to_increase_y1_BBCIO,...
    Max_percentage_to_increase_y1_BBCIO,Mean_percentage_to_increase_y1_BBCIO,...
    Standard_deviation_percentage_to_increase_y1_BBCIO,...
    interquartile_range_percentage_to_increase_y1_BBCIO];
Stats_percentage_to_increase_y1_BBCIO=Stats_percentage_to_increase_y1_BBCIO.'
 
% statistical measures for percentage_to_increase_y2_BBCIO 
Min_percentage_to_increase_y2_BBCIO=min(percentage_to_increase_y2_BBCIO); 
Max_percentage_to_increase_y2_BBCIO=max(percentage_to_increase_y2_BBCIO);
Mean_percentage_to_increase_y2_BBCIO=mean(percentage_to_increase_y2_BBCIO);
Standard_deviation_percentage_to_increase_y2_BBCIO=std(percentage_to_increase_y2_BBCIO);
interquartile_range_percentage_to_increase_y2_BBCIO=iqr(percentage_to_increase_y2_BBCIO);
 
Stats_percentage_to_increase_y2_BBCIO=[Min_percentage_to_increase_y2_BBCIO,...
    Max_percentage_to_increase_y2_BBCIO,Mean_percentage_to_increase_y2_BBCIO,...
    Standard_deviation_percentage_to_increase_y2_BBCIO,...
    interquartile_range_percentage_to_increase_y2_BBCIO];
 
Stats_percentage_to_increase_y2_BBCIO=Stats_percentage_to_increase_y2_BBCIO.'
    
%table with all the statistical measures
statistical_measures=["MIN","MAX","MEAN",...
    "STD","INTER_QUARTILE"]
statistical_measures=statistical_measures.'
  
table_statistical_measures_BBCIO=table(statistical_measures,...
    Stats_DEA_Score_BBCIO,...
    Stats_percentage_to_decrease_x1_BBCIO,...
    Stats_percentage_to_decrease_x2_BBCIO,...
    Stats_percentage_to_increase_y1_BBCIO,...
    Stats_percentage_to_increase_y2_BBCIO,...
    Stats_slack_x1_BBCIO,...
    Stats_slack_x2_BBCIO,...
    Stats_slack_y1_BBCIO,...
    Stats_slack_y2_BBCIO)
table_statistical_measures_BBCIO.Properties.RowNames={'MIN','MAX','MEAN',...
    'STD','INTER_QUARTILE'}
 
filename = 'BBCIO_Table_Statistical_Measures.xlsx';
writetable(table_statistical_measures_BBCIO,filename);
 
% visualization table: in order to visualize all the important parts in one
% table 
table_visualization_table_BBCIO=table(nameDMU,x1_input,x2_input,y1_input,...
    y2_input,DEA_score_vector_T3_BBCIO,...
    slack_x1_vector_BBCIO,... 
    slack_x2_vector_BBCIO,slack_y1_vector_BBCIO,slack_y2_vector_BBCIO,...
    sum_of_x_slacks_BBCIO,sum_of_y_slacks_BBCIO,projections_x1_BBCIO,...
    projections_x2_BBCIO,projections_y1_BBCIO,projections_y2_BBCIO,...
    percentage_to_decrease_x1_BBCIO,percentage_to_decrease_x2_BBCIO,...
    percentage_to_increase_y1_BBCIO,percentage_to_increase_y2_BBCIO);
 
filename = 'BBCIO_Visualization_Table.xlsx';
writetable(table_visualization_table_BBCIO,filename)


%% Generates the output file "DEA_T3_IO_Output_Data_File.table";

% THIS IS THE TABLE SHOWCASING THE LAMBDAS (WEIGHTS)

% nameDMU;
temp1 = 'DMU_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
nameDMU = [temp3, temp4];
clear temp1 temp2 temp3 temp4

% labelX;
temp1 = '     Input_';
temp2 = ones(m,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:m)');
temp5 = [temp3, temp4];
labelX = [];
for ind = 1:m
    labelX = [labelX, temp5(ind,:)];
end
labelX;
clear temp1 temp2 temp3 temp4 temp5

% labelY;
temp1 = '    output_';
temp2 = ones(s,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:s)');
temp5 = [temp3, temp4];
labelY = [];
for ind = 1:s
    labelY = [labelY,temp5(ind,:)];
end
labelY;
clear temp1 temp2 temp3 temp4 temp5

% labellambda;
temp1 = '   lambda_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
temp5 = [temp3, temp4];
labellambda = [];
for ind = 1:n
    labellambda = [labellambda,temp5(ind,:)];
end
labellambda;
clear temp1 temp2 temp3 temp4 temp5

% labelx_slack;
temp1 = '   x_slack_';
temp2 = ones(m,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:m)');
temp5 = [temp3, temp4];
labelx_slack = [];
for ind = 1:m
    labelx_slack = [labelx_slack,temp5(ind,:)];
end
labelx_slack;
clear temp1 temp2 temp3 temp4 temp5

% labely_slack;
temp1 = '   y_slack_';
temp2 = ones(s,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:s)');
temp5 = [temp3, temp4];
labely_slack = [];
for ind = 1:s
    labely_slack = [labely_slack,temp5(ind,:)];
end
labely_slack;
clear temp1 temp2 temp3 temp4 temp5

% formatX;
temp1 = '%12.2f';
formatX = [];
for ind = 1:m
    formatX = [formatX,temp1];
end
clear temp1;

% formatY;
temp1 =  '%12.2f';
formatY = [];
for ind = 1:s
    formatY = [formatY,temp1];
end
clear temp1;

% formatlambda;
temp1 = '%12.4f';
formatlambda = [];
for ind = 1:n
    formatlambda = [formatlambda,temp1];
end
clear temp1;

%% Generates the file 'DEA_T3_IO_Output_Data_File.table' with the results; 

fid = fopen('DEA_T3_IO_Output_Data_File.table','wt');
fprintf(fid,'DEA results for the ');

switch model

    case 'BCC'
    switch orientation
        case 'io'
        fprintf(fid,'input oriented BCC model\n\n');
        case 'oo'
        fprintf(fid,'output oriented BCC model\n\n');
    end
    fmtW = [formatX, formatY, formatlambda, formatX, formatY, '%8.2f \n'];
    W = [X, Y, Z ];
    fprintf(fid,['DMU name ', labelX, labelY, labellambda, labelx_slack, labely_slack '   Theta\n']);

end

% fprintf(fid,fmtW,W');
for ind = 1:n;
    fprintf(fid,nameDMU(ind,:));
    fprintf(fid,fmtW,W(ind,:)');
end
fclose(fid);
