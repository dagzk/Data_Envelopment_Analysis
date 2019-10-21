X1=DEAInputDataFile2018(:,2);
X2=DEAInputDataFile2018(:,3);
Y1=DEAInputDataFile2018(:,4);
Y2=DEAInputDataFile2018(:,5);
% scatter(X1,X2,'k')
hold on
Px1=Analysistask1and3(:,1);
Px2=Analysistask1and3(:,2);
Py1=Analysistask1and3(:,3);
Py2=Analysistask1and3(:,4);
% scatter(Px1,Px2,'r')

Input=[X1,X2];
Output=[Y1,Y2];
%CMD for Input
D1=pdist(Input);
[A1,eigvals1] = cmdscale(D1);

%CMD for Output
D2=pdist(Output);
[A2,eigvals2] = cmdscale(D2);

%IO projections
Px=Analysistask1and3(:,1:2);
Py=Analysistask1and3(:,3:4);
Dp1=pdist(Px);
[Ap1,eigvals_p1] = cmdscale(Dp1);

Dp2=pdist(Py);
[Ap2,eigvals_p2] = cmdscale(Dp2);

NewX=A1(:,1);
NewY=A2(:,1);
plot(NewX,NewY,'k')
hold on

% NewPx=Ap1(:,1);
% NewPy=Ap2(:,1);
% scatter(NewPx,NewPy,'r')


