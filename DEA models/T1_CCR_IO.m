% PAMP: DEA PROJECT;

%% Task 1 CCR-IO

% Daghan Kendirli
% Ruitong Li
% Monica Le
% Christopher Loynes
% Olivier Kraaijeveld
% Giorgios Makridakis
clear all;

%% Model selection; CCR IO and OO
model='CCR';               % Choice of model    model= 'CCR' ;
orientation='io';          % orientation = 'io' for input oriented ;

epsilon= 1E-6;              % epsilon (non archimedian number for the BCC & CCR models) ;

%% Enters the inputs and outputs of all DMUs

X1 = xlsread('DEA_Input_Data_File_2018.xls', 'B:B');

X = X1(:, 2:3);
Y = X1(:, 4:5);

% extracts the number of DMUs, inputs and outputs;
[n,m] = size(X);
[n,s] = size(Y);

%% Computes the results from the selected model;
%  The results from the selected model are put in matrix Z;

switch model; % allows switch between IO and OO

    % CCR model;
    case ('CCR')
    switch orientation;
        % Input oriented; THIS PART IS BEING RUN FOR IO
        case ('io') 
            Z = zeros(n,n+m+s+1);
            
            % Objective function of the CCR model: min(0*lambda - epsilon*(s+ + s-) + theta);
            f = [zeros(1,n) -epsilon*ones(1,s) -epsilon*ones(1,m) 1];

            lblambda = zeros(n,1);                % Lower bounds for (n) lambdas;
            lboutput = zeros(s,1);                % Lower bounds for (s) outputs;
            lbinput  = zeros(m,1);                % Lower bounds for (m) inputs ;
            lb = [lblambda; lboutput; lbinput];   % Lower bounds for lambdas, outputs (s+) and inputs (s-);
            for j=1:n
                Aeq = [Y', -eye(s,s), zeros(s,m+1);
                      -X', zeros(m,s), -eye(m,m), X(j,:)'];
                beq = [Y(j,:)';zeros(m,1)];
                z = linprog(f,[],[],Aeq,beq,lb);
                Z(j,:) = z;
            end
            Z

        % Output oriented;
        case ('oo')
            Z = zeros(n,n+m+s+1);
            
            % Objective function of the CCR_oo model: max(0*lambda + epsilon*(s+ + s-) + theta);
            f = -[zeros(1,n), epsilon*ones(1,s+m), 1];

            lblambda = zeros(n,1);                % Lower bounds for (n) lambdas;
            lboutput = zeros(s,1);                % Lower bounds for (s) outputs;
            lbinput  = zeros(m,1);                % Lower bounds for (m) inputs ;
            lb = [lblambda; lboutput; lbinput];   % Lower bounds for lambdas, outputs (s+) and inputs (s-);
            for j=1:n
                Aeq = [-Y', eye(s,s), zeros(s,m), Y(j,:)'; ...
                        X', zeros(m,s), eye(m,m), zeros(m,1)];
                beq = [zeros(s,1);X(j,:)'];
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
 
DEA_score_vector_CCRIO=Z(:,n+s+m+1);
x1_input=X(:,1);
x2_input=X(:,2);
y1_input=Y(:,1);
y2_input=Y(:,2);
slack_x1_vector_CCRIO=Z(:,51);
slack_x2_vector_CCRIO=Z(:,52);
slack_y1_vector_CCRIO=Z(:,53);
slack_y2_vector_CCRIO=Z(:,54);
intensity_vector_CCRIO=Z(:,1:50);
 
% Proceed with the necessary calculation to calculate the sum of slacks, 
% projections of the inefficient DMUs, the percentage of input and output
% that need to be decreased or increased.
 
sum_of_x_slacks_CCRIO=slack_x1_vector_CCRIO+slack_x2_vector_CCRIO;
sum_of_y_slacks_CCRIO=slack_y1_vector_CCRIO+slack_y2_vector_CCRIO;
projections_x1_CCRIO= DEA_score_vector_CCRIO.*x1_input-(slack_x1_vector_CCRIO);
projections_x2_CCRIO= DEA_score_vector_CCRIO.*x2_input-(slack_x2_vector_CCRIO);           
projections_y1_CCRIO= y1_input+(slack_y1_vector_CCRIO);
projections_y2_CCRIO=y2_input+(slack_y2_vector_CCRIO);
 
% percentage_to_decrease/increase: calculate how far in percentage the
% current input is from the projection in the efficiency frontier.
 
percentage_to_decrease_x1_CCRIO=(x1_input-projections_x1_CCRIO)./x1_input;
 
percentage_to_decrease_x2_CCRIO=(x2_input-projections_x2_CCRIO)./x2_input;
percentage_to_decrease_x2_CCRIO(isnan(percentage_to_decrease_x2_CCRIO))=0
percentage_to_increase_y1_CCRIO=(projections_y1_CCRIO-y1_input)./y1_input;
percentage_to_increase_y1_CCRIO(isnan(percentage_to_increase_y1_CCRIO))=0
percentage_to_increase_y2_CCRIO=(projections_y2_CCRIO-y2_input)./y2_input;
percentage_to_increase_y2_CCRIO(isnan(percentage_to_increase_y2_CCRIO))=0

% Create tables

% Creation of a vector with the names of the DMUs
temp1 = 'DMU_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
nameDMU = [temp3, temp4];
clear temp1 temp2 temp3 temp4
 
% Creation of the table for Slacks and DEA scores
table_analysis_CCRIO=table(nameDMU,x1_input,x2_input,y1_input,y2_input,slack_x1_vector_CCRIO,... 
    slack_x2_vector_CCRIO,slack_y1_vector_CCRIO,slack_y2_vector_CCRIO,...
    sum_of_x_slacks_CCRIO,sum_of_y_slacks_CCRIO,DEA_score_vector_CCRIO);
filename = 'CCRIO_Slacks_and_DEA_scores.xlsx';
writetable(table_analysis_CCRIO,filename)
 
% Creation of the table Percentage inc and dec
table_analysis_percentages_CCRIO=table(nameDMU,percentage_to_decrease_x1_CCRIO, ...
    percentage_to_decrease_x2_CCRIO,percentage_to_increase_y1_CCRIO,...
    percentage_to_increase_y2_CCRIO,DEA_score_vector_CCRIO)
filename = 'CCRIO_Percentages_increase_decrease.xlsx';
writetable(table_analysis_percentages_CCRIO,filename)

% Creation of Projections table
table_analysis_projections_in_efficient_frontier_CCRIO=table(nameDMU,...
    x1_input,x2_input,y1_input,y2_input,...
    DEA_score_vector_CCRIO,projections_x1_CCRIO,projections_x2_CCRIO...
    ,projections_y1_CCRIO,projections_y2_CCRIO);
filename = 'CCRIO_Projection_efficient_frontier_CCRIO.xlsx';
writetable(table_analysis_projections_in_efficient_frontier_CCRIO,filename)
 
 
% The outputs are calculated, we use statistical measures to find 
% the meaning and patterns in the outcomes of the DEA
 
% statistical measures for slack_x1
Min_slack_x1_CCRIO=min(slack_x1_vector_CCRIO(slack_x1_vector_CCRIO>0)); 
Max_slack_x1_CCRIO=max(slack_x1_vector_CCRIO);
Mean_slack_x1_CCRIO=mean(slack_x1_vector_CCRIO);
Standard_deviation_slack_x1_CCRIO=std(slack_x1_vector_CCRIO);
interquartile_range_slack_x1_CCRIO=iqr(slack_x1_vector_CCRIO);
 
Stats_slack_x1_CCRIO=[Min_slack_x1_CCRIO,...
Max_slack_x1_CCRIO,Mean_slack_x1_CCRIO,...
    Standard_deviation_slack_x1_CCRIO,...
    interquartile_range_slack_x1_CCRIO];
Stats_slack_x1_CCRIO=Stats_slack_x1_CCRIO.'
 
% statistical measures for slack_x2
Min_slack_x2_CCRIO=min(slack_x2_vector_CCRIO(slack_x2_vector_CCRIO>0)); 
Max_slack_x2_CCRIO=max(slack_x2_vector_CCRIO);
Mean_slack_x2_CCRIO=mean(slack_x2_vector_CCRIO);
Standard_deviation_slack_x2_CCRIO=std(slack_x2_vector_CCRIO);
interquartile_range_slack_x2_CCRIO=iqr(slack_x2_vector_CCRIO);
 
Stats_slack_x2_CCRIO=[Min_slack_x2_CCRIO,...
Max_slack_x2_CCRIO,Mean_slack_x2_CCRIO,...
    Standard_deviation_slack_x2_CCRIO,...
    interquartile_range_slack_x2_CCRIO];
Stats_slack_x2_CCRIO=Stats_slack_x2_CCRIO.'
 
% statistical measures for slack_y1
Min_slack_y1_CCRIO=min(slack_y1_vector_CCRIO); 
Max_slack_y1_CCRIO=max(slack_y1_vector_CCRIO);
Mean_slack_y1_CCRIO=mean(slack_y1_vector_CCRIO);
Standard_deviation_slack_y1_CCRIO=std(slack_y1_vector_CCRIO);
interquartile_range_slack_y1_CCRIO=iqr(slack_y1_vector_CCRIO);
 
Stats_slack_y1_CCRIO=[Min_slack_y1_CCRIO,...
Max_slack_y1_CCRIO,...
Mean_slack_y1_CCRIO,Standard_deviation_slack_y1_CCRIO,...
interquartile_range_slack_y1_CCRIO];
 
Stats_slack_y1_CCRIO=Stats_slack_y1_CCRIO.'
 
% statistical measures for slack_y2
Min_slack_y2_CCRIO=min(slack_y2_vector_CCRIO(slack_y2_vector_CCRIO>0)); 
Max_slack_y2_CCRIO=max(slack_y2_vector_CCRIO);
Mean_slack_y2_CCRIO=mean(slack_y2_vector_CCRIO);
Standard_deviation_slack_y2_CCRIO=std(slack_y2_vector_CCRIO);
interquartile_range_slack_y2_CCRIO=iqr(slack_y2_vector_CCRIO);
 
Stats_slack_y2_CCRIO=[Min_slack_y2_CCRIO,...
Max_slack_y2_CCRIO,Mean_slack_y2_CCRIO,...
    Standard_deviation_slack_y2_CCRIO,...
    interquartile_range_slack_y2_CCRIO];
 
Stats_slack_y2_CCRIO=Stats_slack_y2_CCRIO.'
 
% statistical measures for DEA score 
Min_DEA_Score_CCRIO=min(DEA_score_vector_CCRIO(DEA_score_vector_CCRIO>0)); 
Max_DEA_Score_CCRIO=max(DEA_score_vector_CCRIO);
Mean_DEA_Score_CCRIO=mean(DEA_score_vector_CCRIO);
Standard_deviation_DEA_Score_CCRIO=std(DEA_score_vector_CCRIO);
interquartile_range_DEA_Score_CCRIO=iqr(DEA_score_vector_CCRIO);
 
Stats_DEA_Score_CCRIO=[Min_DEA_Score_CCRIO,...
Max_DEA_Score_CCRIO,Mean_DEA_Score_CCRIO,...
    Standard_deviation_DEA_Score_CCRIO,...
    interquartile_range_DEA_Score_CCRIO];
Stats_DEA_Score_CCRIO=Stats_DEA_Score_CCRIO.'
 
% statistical measures for percentage_to_decrease_x1_CRRIO
 
Min_percentage_to_decrease_x1_CCRIO=min(percentage_to_decrease_x1_CCRIO);
Max_percentage_to_decrease_x1_CCRIO=max(percentage_to_decrease_x1_CCRIO);
Mean_percentage_to_decrease_x1_CCRIO=mean(percentage_to_decrease_x1_CCRIO);
Standard_deviation_percentage_to_decrease_x1_CCRIO=std(percentage_to_decrease_x1_CCRIO);
interquartile_range_percentage_to_decrease_x1_CCRIO=iqr(percentage_to_decrease_x1_CCRIO);
 
Stats_percentage_to_decrease_x1_CCRIO=[Min_percentage_to_decrease_x1_CCRIO,...
   Max_percentage_to_decrease_x1_CCRIO,Mean_percentage_to_decrease_x1_CCRIO,...
    Standard_deviation_percentage_to_decrease_x1_CCRIO,...
    interquartile_range_percentage_to_decrease_x1_CCRIO];
Stats_percentage_to_decrease_x1_CCRIO=Stats_percentage_to_decrease_x1_CCRIO.'
 
% statistical measures for percentage_to_decrease_x2_CRRIO
Min_percentage_to_decrease_x2_CCRIO=min(percentage_to_decrease_x2_CCRIO); 
Max_percentage_to_decrease_x2_CCRIO=max(percentage_to_decrease_x2_CCRIO); 
Mean_percentage_to_decrease_x2_CCRIO=mean(percentage_to_decrease_x2_CCRIO); 
Standard_deviation_percentage_to_decrease_x2_CCRIO=std(percentage_to_decrease_x2_CCRIO); 
interquartile_range_percentage_to_decrease_x2_CCRIO=iqr(percentage_to_decrease_x2_CCRIO); 
 
Stats_percentage_to_decrease_x2_CCRIO=[Min_percentage_to_decrease_x2_CCRIO,...
    Max_percentage_to_decrease_x2_CCRIO,Mean_percentage_to_decrease_x2_CCRIO,...
    Standard_deviation_percentage_to_decrease_x2_CCRIO,...
   interquartile_range_percentage_to_decrease_x2_CCRIO];
Stats_percentage_to_decrease_x2_CCRIO=Stats_percentage_to_decrease_x2_CCRIO.'
% statistical measures for percentage_to_increase_y1_CCRIO 
Min_percentage_to_increase_y1_CCRIO=min(percentage_to_increase_y1_CCRIO); 
Max_percentage_to_increase_y1_CCRIO=max(percentage_to_increase_y1_CCRIO);
Mean_percentage_to_increase_y1_CCRIO=mean(percentage_to_increase_y1_CCRIO);
Standard_deviation_percentage_to_increase_y1_CCRIO=std(percentage_to_increase_y1_CCRIO);
interquartile_range_percentage_to_increase_y1_CCRIO=iqr(percentage_to_increase_y1_CCRIO);
 
Stats_percentage_to_increase_y1_CCRIO=[Min_percentage_to_increase_y1_CCRIO,...
    Max_percentage_to_increase_y1_CCRIO,Mean_percentage_to_increase_y1_CCRIO,...
    Standard_deviation_percentage_to_increase_y1_CCRIO,...
    interquartile_range_percentage_to_increase_y1_CCRIO];
Stats_percentage_to_increase_y1_CCRIO=Stats_percentage_to_increase_y1_CCRIO.'
 
% statistical measures for percentage_to_increase_y2_CCRIO 
Min_percentage_to_increase_y2_CCRIO=min(percentage_to_increase_y2_CCRIO); 
Max_percentage_to_increase_y2_CCRIO=max(percentage_to_increase_y2_CCRIO);
Mean_percentage_to_increase_y2_CCRIO=mean(percentage_to_increase_y2_CCRIO);
Standard_deviation_percentage_to_increase_y2_CCRIO=std(percentage_to_increase_y2_CCRIO);
interquartile_range_percentage_to_increase_y2_CCRIO=iqr(percentage_to_increase_y2_CCRIO);
 
Stats_percentage_to_increase_y2_CCRIO=[Min_percentage_to_increase_y2_CCRIO,...
    Max_percentage_to_increase_y2_CCRIO,Mean_percentage_to_increase_y2_CCRIO,...
    Standard_deviation_percentage_to_increase_y2_CCRIO,...
    interquartile_range_percentage_to_increase_y2_CCRIO];
 
Stats_percentage_to_increase_y2_CCRIO=Stats_percentage_to_increase_y2_CCRIO.'
    
%table with all the statistical measures
statistical_measures=["MIN","MAX","MEAN",...
    "STD","INTER_QUARTILE"]
statistical_measures=statistical_measures.'
table_statistical_measures_CCRIO=table(statistical_measures,...
    Stats_DEA_Score_CCRIO,...
    Stats_percentage_to_decrease_x1_CCRIO,...
    Stats_percentage_to_decrease_x2_CCRIO,...
    Stats_percentage_to_increase_y1_CCRIO,...
    Stats_percentage_to_increase_y2_CCRIO,...
    Stats_slack_x1_CCRIO,...
    Stats_slack_x2_CCRIO,...
    Stats_slack_y1_CCRIO,...
    Stats_slack_y2_CCRIO)
 
 
filename = 'CCRIO_Table_Statistical_Measures.xlsx';
writetable(table_statistical_measures_CCRIO,filename);
    

% visualization table: in order to visualize all the important parts in one
% table 
table_visualization_table_CCRIO=table(nameDMU,x1_input,x2_input,y1_input,...
    y2_input,DEA_score_vector_CCRIO,...
    slack_x1_vector_CCRIO,... 
    slack_x2_vector_CCRIO,slack_y1_vector_CCRIO,slack_y2_vector_CCRIO,...
    sum_of_x_slacks_CCRIO,sum_of_y_slacks_CCRIO,projections_x1_CCRIO,...
    projections_x2_CCRIO,projections_y1_CCRIO,projections_y2_CCRIO,...
    percentage_to_decrease_x1_CCRIO,percentage_to_decrease_x2_CCRIO,...
    percentage_to_increase_y1_CCRIO,percentage_to_increase_y2_CCRIO);
 
filename = 'CCRIO_Visualization_Table.xlsx';
writetable(table_visualization_table_CCRIO,filename)


%% Generate the output file "DEA_T1_IO_Output_Data_File.table";

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

%% Generates the file 'DEA_T1_IO_Output_Data_File.table' with the results'; 

fid = fopen('DEA_T1_IO_Output_Data_File.table','wt');
fprintf(fid,'DEA results for the ');
    
 switch model;
    
    case 'CCR'
    switch orientation
        case 'io'
        fprintf(fid,'input oriented CCR model\n\n');
        case 'oo'
        fprintf(fid,'output oriented CCR model\n\n');
    end
    fprintf(fid,'CCR model\n\n');
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
    
    