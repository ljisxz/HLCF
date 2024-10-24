function  [W,V,obj]=HLCF(X,Options)
  
 %..............................................
        [~, N]= size(X);
%........................................................
        Y=zeros(Options.LabeledNum, Options.SelectClasses);
        for i =1:Options.LabeledNum
              Y(i,Options.Xgnd(i)) = 1;
        end
      
        W = rand(N, Options.SelectClasses);
        V =  rand(N, Options.SelectClasses);
        V(1:Options.LabeledNum,:)=V(1:Options.LabeledNum,:).*Y;
        obj=[];

        opts.maxIter=20;
        opts.p=Options.p;
        opts.mu= Options.mu;
        opts.gndSmpNum=Options.LabeledNum;
        opts.labels= Options.Xgnd;
        opts.K=Options.SelectClasses;
        [S, ~, ~]= createHyperS(X',opts);
        I=eye(N);
        obj(1)=trace((I-W*V')'*S*(I-W*V'));
        
        for iter = 1:Options.maxIter
          W =W.*(S*V)./(S*W*V'*V+eps);
          V =V.*((S*W)./(V*W'*S*W+eps));
          V(1:Options.LabeledNum,:)=V(1:Options.LabeledNum,:).*Y;
          obj(iter+1)=trace((I-W*V')'*S*(I-W*V'));
          
       end
     
end
