 

load('COIL100_Obj.mat');	
 nClass = length(unique(gnd));
 dataset='COIL100_Obj'; 
 [a,b]=sort(gnd,'ascend');
 gnd=a;
 fea=fea(b,:);
 fea = NormalizeFea(fea); 
 SelectClasses=nClass;
 LabelsRatio=0.1;
 [X,Xgnd,LabeledNum]=Creat_SampleDatasets(fea,SelectClasses,gnd,nClass,LabelsRatio);    
 % YaleB  UNIST  COIL100  MNIST
 %p  2       3      2       6
 %mu 0.5     0.3    0.2     0.5
 Options.maxIter=20;
 Options.mu=0.2;
 Options.p =2;
 Options.LabeledNum=LabeledNum;
 Options.Xgnd=Xgnd;
 Options.SelectClasses=SelectClasses;
 Options.nClass=nClass;
 [~,V,~]=HLCF(X',Options);
 [~, PreLabels] = max(V');
 result =Clustering8Measure(Xgnd,PreLabels);