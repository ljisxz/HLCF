function   [ revS, Dv, L]= createHyperS(X,opts)
n=size(X,1);
Z = zeros(n);
for i = 1:opts.gndSmpNum

    for j = 1:opts.gndSmpNum

        if opts.labels(i)==opts.labels(j)
            Z(i,j)=1;
         else
           Z(i,j)=-1; 
        end
    
    end
end
[revS, Dv, L] = PCP(X,Z,opts);
end


function [revS, Dv, L] = PCP(X,Z,opts)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

p = opts.p;
mu = opts.mu;
%view = length(X);

%revS = cell(1,view);
%Dv = cell(1,view);
%L = cell(1,view);


%for i=1:view

    [H, W] = constructHW(X,p);

    invDe = diag(1./sum(H,2));

    S = H'*W*invDe*H;
    S=(S+S')/2;

    revS = reviseHyperW(S,Z,mu);
    
    diag_revS = diag(revS);
    diag_revS_half = diag(diag_revS.^(-0.5));
    tmp_resS = diag_revS_half*revS*diag_revS_half;
    revS = (tmp_resS + tmp_resS') / 2;    
    
    Dv = diag(sum(revS,1));
    L = Dv- revS;
    
end
function [H, W] = constructHW(fea, p)

% Each row of H consists of a hyper edge.
% Each column of H means the sample is used to consist of the corresponding
% hyper edge.

[nSmp, ~] = size(fea);
G = zeros(nSmp*(p+1),3);

dist = EuDist2(fea);
D = dist.^2;
sigma = mean(mean(dist));
%sigma=sqrt(2);
%sigma=2;
A = exp(D/(-sigma.^2));

dump = zeros(nSmp, p+1);
idx = dump;

for j = 1:p+1
    [dump(:,j),idx(:,j)] = min(dist,[],2);
    temp = (idx(:,j)-1)*nSmp+[1:nSmp]';
    dist(temp) = 1e100;
end

G(1:nSmp*(p+1),1) = repmat([1:nSmp]',[p+1,1]);
G(1:nSmp*(p+1),2) = idx(:);
G(1:nSmp*(p+1),3) = 1;

H = sparse(G(:,1),G(:,2),G(:,3),nSmp,nSmp);
W = diag(sum(A.*H'));

end

function reW = reviseHyperW(W, Z, alpha)

cp_num = 30;

n = size(W,1);
gd = full(sum(W,2));
d_mhalf = spdiags(gd.^-.5,0,n,n);
L = d_mhalf*W*d_mhalf;
Fv = 2*rand(size(W))-1;
Fh = 2*rand(size(W))-1;

for i=1:cp_num
    Fv = alpha*L*Fv + (1-alpha)*Z;  
end

for i=1:cp_num
    Fh = alpha*Fh*L + (1-alpha)*Fv;  
end
F = Fh;
idx = (F>=0);
diagW = diag(W);
W_ii = repmat(diagW,1,n);
reW = ((W_ii-(1-F).*(W_ii-W)) .* idx) +  (((1+F).*W) .* (~idx));
end
