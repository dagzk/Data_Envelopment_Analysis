%% PAMP DEA BCC-OO File;
clear all;

%% Task 3 BCC-OO

%% Model selection;
model='BCC';               % Choice of model    model= 'BCC' ;
orientation='oo';          % orientation = 'io' or 'oo' (for input or output oriented)    ;

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
        % Input oriented;
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

        % Output oriented; THIS IS USED IN THIS CASE
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
 
DEA_score_vector_BBCOO=Z(:,n+s+m+1);
x1_input=X(:,1);
x2_input=X(:,2);
y1_input=Y(:,1);
y2_input=Y(:,2);
slack_x1_vector_BBCOO=Z(:,51);
slack_x2_vector_BBCOO=Z(:,52);
slack_y1_vector_BBCOO=Z(:,53);
slack_y2_vector_BBCOO=Z(:,54);
intensity_vector_BBCOO=Z(:,1:50);
 
% Proceed with the necessary calculation to calculate the sum of slacks, 
% projections of the inefficient DMUs, the percentage of input and output
% that need to be decreased or increased.
 
sum_of_x_slacks_BBCOO=slack_x1_vector_BBCOO+slack_x2_vector_BBCOO;
sum_of_y_slacks_BBCOO=slack_y1_vector_BBCOO+slack_y2_vector_BBCOO;
projections_x1_BBCOO= x1_input-slack_x1_vector_BBCOO;
projections_x2_BBCOO= x2_input-slack_x2_vector_BBCOO;           
projections_y1_BBCOO= DEA_score_vector_BBCOO.*y1_input+(slack_y1_vector_BBCOO);
projections_y2_BBCOO=DEA_score_vector_BBCOO.*y2_input+(slack_y2_vector_BBCOO);
 
% percentage_to_decrease/increase: calculate how far in percentage the
% current input is from the projection in the efficiency frontier.
 
percentage_to_decrease_x1_BBCOO=(x1_input-projections_x1_BBCOO)./x1_input;
percentage_to_decrease_x1_BBCOO(isnan(percentage_to_decrease_x1_BBCOO))=0
percentage_to_decrease_x2_BBCOO=(x2_input-projections_x2_BBCOO)./x2_input;
percentage_to_decrease_x2_BBCOO(isnan(percentage_to_decrease_x2_BBCOO))=0
percentage_to_increase_y1_BBCOO=(projections_y1_BBCOO-y1_input)./y1_input;
percentage_to_increase_y1_BBCOO(isnan(percentage_to_increase_y1_BBCOO))=0
percentage_to_increase_y2_BBCOO=(projections_y2_BBCOO-y2_input)./y2_input;
percentage_to_increase_y2_BBCOO(isnan(percentage_to_increase_y2_BBCOO))=0
 
% Create tables

% Creation of a vector with the names of the DMUs 

temp1 = 'DMU_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
nameDMU = [temp3, temp4];
clear temp1 temp2 temp3 temp4
 
% Creation of the table for Slacks and DEA scores
table_analysis_BBCOO=table(nameDMU,x1_input,x2_input,y1_input,y2_input,slack_x1_vector_BBCOO,... 
    slack_x2_vector_BBCOO,slack_y1_vector_BBCOO,slack_y2_vector_BBCOO,...
    sum_of_x_slacks_BBCOO,sum_of_y_slacks_BBCOO,DEA_score_vector_BBCOO);
filename = 'BBCOO_Slacks_and_DEA_scores.xlsx';
writetable(table_analysis_BBCOO,filename)
 
% Creation of the table Percentage inc and dec
table_analysis_percentages_BBCOO=table(nameDMU,percentage_to_decrease_x1_BBCOO, ...
    percentage_to_decrease_x2_BBCOO,percentage_to_increase_y1_BBCOO,...
    percentage_to_increase_y2_BBCOO,DEA_score_vector_BBCOO)
filename = 'BBCOO_Percentage_increase_decrease.xlsx';
writetable(table_analysis_percentages_BBCOO,filename)
 
% Creation of Projections table
table_analysis_projections_efficient_frontier_BBCOO=table(nameDMU,...
    x1_input,x2_input,y1_input,y2_input,...
    DEA_score_vector_BBCOO,projections_x1_BBCOO,projections_x2_BBCOO...
    ,projections_y1_BBCOO,projections_y2_BBCOO);
filename = 'BBCOO_Projection_efficient_frontier_BBCOO.xlsx';
writetable(table_analysis_projections_efficient_frontier_BBCOO,filename)
 
 
% The outputs are calculated, we use statistical measures to find 
% the meaning and patterns in the outcomes of the DEA
 
% statistical measures for slack_x1
Min_slack_x1_BBCOO=min(slack_x1_vector_BBCOO(slack_x1_vector_BBCOO>0)); 
Max_slack_x1_BBCOO=max(slack_x1_vector_BBCOO);
Mean_slack_x1_BBCOO=mean(slack_x1_vector_BBCOO);
Standard_deviation_slack_x1_BBCOO=std(slack_x1_vector_BBCOO);
interquartile_range_slack_x1_BBCOO=iqr(slack_x1_vector_BBCOO);
 
Stats_slack_x1_BBCOO=[Min_slack_x1_BBCOO,...
Max_slack_x1_BBCOO,Mean_slack_x1_BBCOO,...
    Standard_deviation_slack_x1_BBCOO,...
    interquartile_range_slack_x1_BBCOO];
Stats_slack_x1_BBCOO=Stats_slack_x1_BBCOO.'
 
% statistical measures for slack_x2
Min_slack_x2_BBCOO=min(slack_x2_vector_BBCOO(slack_x2_vector_BBCOO>0)); 
Max_slack_x2_BBCOO=max(slack_x2_vector_BBCOO);
Mean_slack_x2_BBCOO=mean(slack_x2_vector_BBCOO);
Standard_deviation_slack_x2_BBCOO=std(slack_x2_vector_BBCOO);
interquartile_range_slack_x2_BBCOO=iqr(slack_x2_vector_BBCOO);
 
Stats_slack_x2_BBCOO=[Min_slack_x2_BBCOO,...
Max_slack_x2_BBCOO,Mean_slack_x2_BBCOO,...
    Standard_deviation_slack_x2_BBCOO,...
    interquartile_range_slack_x2_BBCOO];
Stats_slack_x2_BBCOO=Stats_slack_x2_BBCOO.'
 
% statistical measures for slack_y1
Min_slack_y1_BBCOO=min(slack_y1_vector_BBCOO); 
Max_slack_y1_BBCOO=max(slack_y1_vector_BBCOO);
Mean_slack_y1_BBCOO=mean(slack_y1_vector_BBCOO);
Standard_deviation_slack_y1_BBCOO=std(slack_y1_vector_BBCOO);
interquartile_range_slack_y1_BBCOO=iqr(slack_y1_vector_BBCOO);
 
Stats_slack_y1_BBCOO=[Min_slack_y1_BBCOO,...
Max_slack_y1_BBCOO,...
Mean_slack_y1_BBCOO,Standard_deviation_slack_y1_BBCOO,...
interquartile_range_slack_y1_BBCOO];
 
Stats_slack_y1_BBCOO=Stats_slack_y1_BBCOO.'
 
% statistical measures for slack_y2
Min_slack_y2_BBCOO=min(slack_y2_vector_BBCOO(slack_y2_vector_BBCOO>0)); 
Max_slack_y2_BBCOO=max(slack_y2_vector_BBCOO);
Mean_slack_y2_BBCOO=mean(slack_y2_vector_BBCOO);
Standard_deviation_slack_y2_BBCOO=std(slack_y2_vector_BBCOO);
interquartile_range_slack_y2_BBCOO=iqr(slack_y2_vector_BBCOO);
 
Stats_slack_y2_BBCOO=[Min_slack_y2_BBCOO,...
Max_slack_y2_BBCOO,Mean_slack_y2_BBCOO,...
    Standard_deviation_slack_y2_BBCOO,...
    interquartile_range_slack_y2_BBCOO];
 
Stats_slack_y2_BBCOO=Stats_slack_y2_BBCOO.'
 
% statistical measures for DEA score 
Min_DEA_Score_BBCOO=min(DEA_score_vector_BBCOO(DEA_score_vector_BBCOO>0)); 
Max_DEA_Score_BBCOO=max(DEA_score_vector_BBCOO);
Mean_DEA_Score_BBCOO=mean(DEA_score_vector_BBCOO);
Standard_deviation_DEA_Score_BBCOO=std(DEA_score_vector_BBCOO);
interquartile_range_DEA_Score_BBCOO=iqr(DEA_score_vector_BBCOO);
 
Stats_DEA_Score_BBCOO=[Min_DEA_Score_BBCOO,...
Max_DEA_Score_BBCOO,Mean_DEA_Score_BBCOO,...
    Standard_deviation_DEA_Score_BBCOO,...
    interquartile_range_DEA_Score_BBCOO];
Stats_DEA_Score_BBCOO=Stats_DEA_Score_BBCOO.'
 
% statistical measures for percentage_to_decrease_x1_BCCOO
 
Min_percentage_to_decrease_x1_BBCOO=min(percentage_to_decrease_x1_BBCOO);
Max_percentage_to_decrease_x1_BBCOO=max(percentage_to_decrease_x1_BBCOO);
Mean_percentage_to_decrease_x1_BBCOO=mean(percentage_to_decrease_x1_BBCOO);
Standard_deviation_percentage_to_decrease_x1_BBCOO=std(percentage_to_decrease_x1_BBCOO);
interquartile_range_percentage_to_decrease_x1_BBCOO=iqr(percentage_to_decrease_x1_BBCOO);
 
Stats_percentage_to_decrease_x1_BBCOO=[Min_percentage_to_decrease_x1_BBCOO,...
   Max_percentage_to_decrease_x1_BBCOO,Mean_percentage_to_decrease_x1_BBCOO,...
    Standard_deviation_percentage_to_decrease_x1_BBCOO,...
    interquartile_range_percentage_to_decrease_x1_BBCOO];
Stats_percentage_to_decrease_x1_BBCOO=Stats_percentage_to_decrease_x1_BBCOO.'
 
% statistical measures for percentage_to_decrease_x2_BCCIO
Min_percentage_to_decrease_x2_BBCOO=min(percentage_to_decrease_x2_BBCOO); 
Max_percentage_to_decrease_x2_BBCOO=max(percentage_to_decrease_x2_BBCOO); 
Mean_percentage_to_decrease_x2_BBCOO=mean(percentage_to_decrease_x2_BBCOO); 
Standard_deviation_percentage_to_decrease_x2_BBCOO=std(percentage_to_decrease_x2_BBCOO); 
interquartile_range_percentage_to_decrease_x2_BBCOO=iqr(percentage_to_decrease_x2_BBCOO); 
 
Stats_percentage_to_decrease_x2_BBCOO=[Min_percentage_to_decrease_x2_BBCOO,...
    Max_percentage_to_decrease_x2_BBCOO,Mean_percentage_to_decrease_x2_BBCOO,...
    Standard_deviation_percentage_to_decrease_x2_BBCOO,...
   interquartile_range_percentage_to_decrease_x2_BBCOO];
Stats_percentage_to_decrease_x2_BBCOO=Stats_percentage_to_decrease_x2_BBCOO.'
 
% statistical measures for percentage_to_increase_y1_BBCOO 
Min_percentage_to_increase_y1_BBCOO=min(percentage_to_increase_y1_BBCOO); 
Max_percentage_to_increase_y1_BBCOO=max(percentage_to_increase_y1_BBCOO);
Mean_percentage_to_increase_y1_BBCOO=mean(percentage_to_increase_y1_BBCOO);
Standard_deviation_percentage_to_increase_y1_BBCOO=std(percentage_to_increase_y1_BBCOO);
interquartile_range_percentage_to_increase_y1_BBCOO=iqr(percentage_to_increase_y1_BBCOO);
 
Stats_percentage_to_increase_y1_BBCOO=[Min_percentage_to_increase_y1_BBCOO,...
    Max_percentage_to_increase_y1_BBCOO,Mean_percentage_to_increase_y1_BBCOO,...
    Standard_deviation_percentage_to_increase_y1_BBCOO,...
    interquartile_range_percentage_to_increase_y1_BBCOO];
Stats_percentage_to_increase_y1_BBCOO=Stats_percentage_to_increase_y1_BBCOO.'
 
% statistical measures for percentage_to_increase_y2_BBCOO 
Min_percentage_to_increase_y2_BBCOO=min(percentage_to_increase_y2_BBCOO); 
Max_percentage_to_increase_y2_BBCOO=max(percentage_to_increase_y2_BBCOO);
Mean_percentage_to_increase_y2_BBCOO=mean(percentage_to_increase_y2_BBCOO);
Standard_deviation_percentage_to_increase_y2_BBCOO=std(percentage_to_increase_y2_BBCOO);
interquartile_range_percentage_to_increase_y2_BBCOO=iqr(percentage_to_increase_y2_BBCOO);
 
Stats_percentage_to_increase_y2_BBCOO=[Min_percentage_to_increase_y2_BBCOO,...
    Max_percentage_to_increase_y2_BBCOO,Mean_percentage_to_increase_y2_BBCOO,...
    Standard_deviation_percentage_to_increase_y2_BBCOO,...
    interquartile_range_percentage_to_increase_y2_BBCOO];
 
Stats_percentage_to_increase_y2_BBCOO=Stats_percentage_to_increase_y2_BBCOO.'
    
%table with all the statistical measures
 statistical_measures=["MIN","MAX","MEAN",...
    "STD","INTER_QUARTILE"]
statistical_measures=statistical_measures.'
 
table_statistical_measures_BBCOO=table(statistical_measures,...
    Stats_DEA_Score_BBCOO,...
    Stats_percentage_to_decrease_x1_BBCOO,...
    Stats_percentage_to_decrease_x2_BBCOO,...
    Stats_percentage_to_increase_y1_BBCOO,...
    Stats_percentage_to_increase_y2_BBCOO,...
    Stats_slack_x1_BBCOO,...
    Stats_slack_x2_BBCOO,...
    Stats_slack_y1_BBCOO,...
    Stats_slack_y2_BBCOO)
table_statistical_measures_BBCOO.Properties.RowNames={'MIN','MAX','MEAN',...
    'STD','INTER_QUARTILE'}
 
filename = 'BBCOO_table_Statistical_Measures.xlsx';
writetable(table_statistical_measures_BBCOO,filename);
    
 
% visualization table: in order to visualize all the important parts in one
% table 
table_visualization_table_BBCOO=table(nameDMU,x1_input,x2_input,y1_input,...
    y2_input,DEA_score_vector_BBCOO,...
    slack_x1_vector_BBCOO,... 
    slack_x2_vector_BBCOO,slack_y1_vector_BBCOO,slack_y2_vector_BBCOO,...
    sum_of_x_slacks_BBCOO,sum_of_y_slacks_BBCOO,projections_x1_BBCOO,...
    projections_x2_BBCOO,projections_y1_BBCOO,projections_y2_BBCOO,...
    percentage_to_decrease_x1_BBCOO,percentage_to_decrease_x2_BBCOO,...
    percentage_to_increase_y1_BBCOO,percentage_to_increase_y2_BBCOO);
 
filename = 'BBCOO_Visualization_Table.xlsx';
writetable(table_visualization_table_BBCOO,filename)


%% Generates the output file "DEA_T3_OO_Output_Data_File.table";

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

%% Generates the file 'DEA_T3_OO_Output_Data_File.table' with the results; 

fid = fopen('DEA_T3_OO_Output_Data_File.table','wt');
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