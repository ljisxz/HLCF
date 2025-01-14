 function  [SamplsFea,XSmpgnd,gndSmpNum]=Creat_SampleDatasets(fea,KClass,gnd,nClass,LabelsRatio)
Tags=tabulate(gnd);
SmpNumPerClass=Tags(:,2);
Per=cumsum(SmpNumPerClass);
%Per=All/nClass;         % 每个类中含有的样本数
Kmix=randperm(nClass);  % 打乱标签的顺序
Klei=Kmix(1:KClass);    % 取随机乱序后的 Klei 类(所要的类数Klei=2,3，……，10)
% 取出这 Klei 类的特征和标签组成原始分解矩阵 X
% 例如：
% 如果是 KClass=2,则取出特征矩阵如下：
% TempFea=[fea(1+(Klei(1)-1)*Per:Klei(1)*Per,:);fea(1+(Klei(2)-1)*Per:Klei(2)*Per,:)];
TempFea=[];
TempOrders=[];
TempLabelsOrders=[];
for i=1:KClass
    if Klei(i)==1
        SmpNumPerClassOrders=1:Per(Klei(i));
       % TempFea=cat(1,TempFea,fea(1:Klei(i),:));
        TempOrders=cat(2,TempOrders,SmpNumPerClassOrders);
         NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio);
%        if strcmp(dataset,'myORL32')        % 判断两个字符串相等，如果是数据集ORL则每个类中取两个样本
%          NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio);
%        else
%          NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio);     % 每类中提取 10% 的样本作为每个类中所含的标签的样本
%        end
        %SmpNumPerClassOrders=Per(1:Klei(i));
        Orders=randperm(length(SmpNumPerClassOrders));
        TempLabelsOrders=[TempLabelsOrders,SmpNumPerClassOrders(Orders(1:NumP))];
    else
        %TempFea=cat(1,TempFea,fea(Per(Klei(i)-1)+1:Per(Klei(i)),:));
        SmpNumPerClassOrders=Per(Klei(i)-1)+1:Per(Klei(i));
        TempOrders=cat(2,TempOrders,SmpNumPerClassOrders);
%          if strcmp(dataset,'myORL32')        % 判断两个字符串相等，如果是数据集ORL则每个类中取两个样本
%              NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio);
%          else
%              NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio);     % 每类中提取 10% 的样本作为每个类中所含的标签的样本
%          end
         NumP=ceil(SmpNumPerClass(Klei(i))*LabelsRatio); 
         Orders=randperm(length(SmpNumPerClassOrders));
         TempLabelsOrders=[TempLabelsOrders,SmpNumPerClassOrders(Orders(1:NumP))];
       
    end
  
end
gndSmpNum=length(TempLabelsOrders);
SampleOrders=[TempLabelsOrders,setdiff(TempOrders,TempLabelsOrders)];
SamplsFea=fea(SampleOrders,:);% 获得新的特征矩阵
%SamplsFea = NormalizeFea(SamplsFea);
SamplsGnd=gnd(SampleOrders);% 获得新的特征矩阵

% 重新规定 gnd 中的标签
Temp=SamplsGnd;
for i=1:KClass
    SamplsGnd(Temp==Klei(i))=i;
end
XSmpgnd=SamplsGnd;
%X=SamplsFea;
 end
 
