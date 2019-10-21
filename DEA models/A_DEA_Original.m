% File DEA;
%clear all;

%% Model selection;
model='CCR';               % Choice of model    model= 'add' for additive, 'BCC' or 'CCR' ;
orientation='oo';          % orientation = 'io' or 'oo' (for input or output oriented)    ;

epsilon=.001;              % epsilon (non archimedian number for the BCC & CCR models)    ;   ;

%% Enters the inputs and outputs of all DMUs ;
%  Alternatively, use the "putmatrix" function in Excel to generate X and Y;

% Input matrix (one row per DMU);
% X = [
% 2 3 4;
% 3 2 6;
% 6 6 3;
% 9 1 2;
% 5 9 5;
% 4 3 3;
% 5 2 6;
% 7 6 1;
% 9 2 4;
% 6 9 9;
% 8 3 2;
% 10 3 3]

% Output matrix (one row per DMU);
% Y = [
% 3 2;
% 6 6;
% 9 1;
% 5 9;
% 4 3;
% 2 2;
% 5 1;
% 7 3;
% 8 2;
% 3 2;
% 1 1;
% 7 4]

X1 = xlsread('DEA_Input_Data_File_2018.xls', 'B:B');

X = X1(:, 2:3);
Y = X1(:, 4:5);

% extracts the number of DMUs, inputs and outputs;
[n,m] = size(X);
[n,s] = size(Y);

%% Computes the results from the selected model;
%  The results from the selected model are in matrix Z;

switch model;

    % Additive model;
%     case ('add')
%         Z = zeros(n,n+m+s);
%         
%         f = [zeros(1,n) -ones(1,s+m)]             % Objective function of the additive model: min(-1*s-1*s-);
%         
%         lblambda = zeros(n,1);                    % Lower bounds for (n) lambdas;
%         lboutput = zeros(s,1);                    % Lower bounds for (s) outputs;
%         lbinput  = zeros(m,1);                    % Lower bounds for (m) inputs ;
%         lb = [lblambda; lboutput; lbinput];       % Lower bounds for lambdas, outputs (s+) and inputs (s-);
%         for j=1:n
%             Aeq = [Y', -eye(s,s), zeros(s,m);
%                   -X', zeros(m,s), -eye(m,m);
%                   ones(1,n), zeros(1,s+m)];
%             beq = [Y(j,:)';-X(j,:)';1];
%             z = linprog(f,[],[],Aeq,beq,lb);
%             Z(j,:) = z;
%         end
%         Z
    
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
    % CCR model;
    case ('CCR')
    switch orientation;
        % Input oriented;
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

%% Generates the output file  "DEAresults.table";

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
temp1 = '%12.2f';
formatlambda = [];
for ind = 1:n
    formatlambda = [formatlambda,temp1];
end
clear temp1;

% Generates the file 'DEAresults.table' with the results; 

fid = fopen('DEAresults.table','wt');
fprintf(fid,'DEA results for the ');

switch model
%     case 'add';
%     fprintf(fid,'Additive model\n\n');
%     fmtW = [formatX, formatY, formatlambda, formatX, formatY, '\n'];
%     W = [X, Y, Z ];
%     fprintf(fid,['DMU name ', labelX, labelY, labellambda, labelx_slack, labely_slack '\n']);

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