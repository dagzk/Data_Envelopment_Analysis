%% DEA Task 2 Dual Form CCR-IO  

clear  

X1 = xlsread('DEA_Input_Data_File_2018.xls', 'B:B');
X=X1(:,2:3);    % Input Matrix X
Y=X1(:,4:5);    % Output Matrix Y 
n=size(X,1); 
m=size(X',1); 
s=size(Y',1); 
epsilon=1E-6;  % Epsilon Value being changed
% CCR_IO
LB=epsilon*ones(m+s,1); UB=[]; 
for i=1:n 
f= [zeros(1,m) -Y(i,:)];    
A=[-X  Y]; 
b=zeros(n, 1); 
Aeq=[X(i,:) zeros(1,s)];
beq=1;
w=linprog(f,A,b,Aeq,beq,LB,UB); % Get the best weight(w) for DMUi 
W(i,:)=w;
E(i, i)=Y(i,:)*W(i,m+1:m+s)';   % Get Efficiency(Eii) for DMUi
end  
W    % Get the best weight for all DMUs
E    % Efficiency for all DMUs
Omega=W(:,1:m)    % Input Weights
mu=W(:,m+1:m+s)   % Output Weights

%% Generate the output file "DEA_T2_IO_Output_Data_File.table";


temp1 = 'DMU_';
temp2 = ones(n,1)*temp1;
temp3 = char(temp2);
temp4 = num2str((1:n)');
nameDMU = [temp3, temp4];
clear temp1 temp2 temp3 temp4

DEA_score_CCRIO=diag(E)
table_names=["DMUs","Input_1","Input_2","Output_1","Output_2", "Weight_v1",...
    "Weight_v2","Weight_u1","Weight_u2","DEA_Score"]
DEA_T2_IO_Output_Data_File_table=table(nameDMU,X(:,1),X(:,2),Y(:,1),Y(:,2),...
    W(:,1),W(:,2),W(:,3),W(:,4),DEA_score_CCRIO)
DEA_T2_IO_Output_Data_File_table.Properties.VariableNames = {'DMUs','Input_1','Input_2','Output_1','Output_2', 'Weight_v1',...
    'Weight_v2','Weight_u1','Weight_u2','DEA_Score'}
filedata ='DEA_T2_IO_Output_Data.xlsx'
writetable(DEA_T2_IO_Output_Data_File_table,filedata)

% statistical measures for weight_v1
Min_weight_v1=min(W(:,1)); 
Max_weight_v1=max(W(:,1));
Mean_weight_v1=mean(W(:,1));
Standard_deviation_weight_v1=std(W(:,1));
interquartile_range_weight_v1=iqr(W(:,1));
 
Stats_weight_v1=[Min_weight_v1,...
Max_weight_v1,Mean_weight_v1,...
    Standard_deviation_weight_v1,...
    interquartile_range_weight_v1];
Stats_weight_v1=Stats_weight_v1.'

% statistical measures for weight_v2
Min_weight_v2=min(W(:,2)); 
Max_weight_v2=max(W(:,2));
Mean_weight_v2=mean(W(:,2));
Standard_deviation_weight_v2=std(W(:,2));
interquartile_range_weight_v2=iqr(W(:,2));
 
Stats_weight_v2=[Min_weight_v2,...
Max_weight_v2,Mean_weight_v2,...
    Standard_deviation_weight_v2,...
    interquartile_range_weight_v2];
Stats_weight_v2=Stats_weight_v2.'

% statistical measures for weight_u1
Min_weight_u1=min(W(:,3)); 
Max_weight_u1=max(W(:,3));
Mean_weight_u1=mean(W(:,3));
Standard_deviation_weight_u1=std(W(:,3));
interquartile_range_weight_u1=iqr(W(:,3));
 
Stats_weight_u1=[Min_weight_u1,...
Max_weight_u1,Mean_weight_u1,...
    Standard_deviation_weight_u1,...
    interquartile_range_weight_u1];
Stats_weight_u1=Stats_weight_u1.'

% statistical measures for weight_u2
Min_weight_u2=min(W(:,4)); 
Max_weight_u2=max(W(:,4));
Mean_weight_u2=mean(W(:,4));
Standard_deviation_weight_u2=std(W(:,4));
interquartile_range_weight_u2=iqr(W(:,4));
 
Stats_weight_u2=[Min_weight_u2,...
Max_weight_u2,Mean_weight_u2,...
    Standard_deviation_weight_u2,...
    interquartile_range_weight_u2];
Stats_weight_u2=Stats_weight_u2.'

 statistical_measures=["MIN","MAX","MEAN",...
    "STD","INTER_QUARTILE"]
statistical_measures=statistical_measures.'

table_statistical_measures_Dual_Form=table(statistical_measures,...
    Stats_weight_v1,Stats_weight_v2, Stats_weight_u1,Stats_weight_u2)

filename = 'DEA_T2_IO_Statistics_Weights.xlsx';
writetable(table_statistical_measures_Dual_Form,filename);

