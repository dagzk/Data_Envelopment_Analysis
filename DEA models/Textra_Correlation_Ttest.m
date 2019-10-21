clear all
clear

% ttest formulation: In addition to the previous statistics, we calculate a
% ttest analysis to imporve the scope of our analysis. 
%the following codes are made to calculate the T-tes of our DEA analysis.

% Here we get the DEA Scores from the previous analysis. 
Dataset_CCROO=xlsread("CCROO_Slacks_and_DEA_scores.xlsx");
[~,n]=size(Dataset_CCROO)
DEA_SCORES_CCROO=Dataset_CCROO(:,n)

Dataset_CCRIO=xlsread("CCRIO_Slacks_and_DEA_scores.xlsx");
[~,n]=size(Dataset_CCRIO)
DEA_SCORES_CCRIO=Dataset_CCRIO(:,n)

Dataset_BBCOO=xlsread("BBCOO_Slacks_and_DEA_scores.xlsx");
[~,n]=size(Dataset_BBCOO)
DEA_SCORES_BBCOO=Dataset_BBCOO(:,n)

Dataset_BBCIO=xlsread("BBCIO_Slacks_and_DEA_scores.xlsx");
[~,n]=size(Dataset_BBCIO)
DEA_SCORES_BBCIO=Dataset_BBCIO(:,n)

% Proceed with  calculation of correlation of the DEA Scores between the 
% IO and OO version of each model.
% Calculate the t-test 

% Correlation code
Correlation_CCROO_CCRIO = corrcoef(DEA_SCORES_CCROO,DEA_SCORES_CCRIO);
Correlation_BBCOO_BBCIO = corrcoef(DEA_SCORES_BBCOO,DEA_SCORES_BBCIO);
Correlation_CCROO_CCRIO=Correlation_CCROO_CCRIO(1,2);
Correlation_BBCOO_BBCIO=Correlation_BBCOO_BBCIO(1,2);

% T test Statistical measure code 
[N,~]=size(DEA_SCORES_BBCIO);
vector_ones=ones(N,1);
d1=vector_ones-DEA_SCORES_CCRIO;
d2=DEA_SCORES_CCROO-vector_ones; 
e1=vector_ones-DEA_SCORES_BBCIO;
e2=DEA_SCORES_BBCOO-vector_ones;

[h,p,ci,stats]=ttest(d1,d2) 
[h1,p1,ci1,stats1]=ttest(e1,e2)





