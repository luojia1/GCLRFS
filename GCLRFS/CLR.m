% min_{S>=0, S*1=1, F'*F=I}  ||S - A||^2 + r*||S||^2 + 2*lambda*trace(F'*L*F)
% or
% min_{S>=0, S*1=1, F'*F=I}  ||S - A||_1 + r*||S||^2 + 2*lambda*trace(F'*L*F)
% function [y, S, evs, cs] = CLR(A0, c, isrobust, islocal)
function [V, S, summ] = CLR(X, A0, AU, c, alpha, beta, theta);
% A0: the given affinity matrix
% c: cluster number
% isrobust: solving the second (L1 based) problem if isrobust=1
% islocal: only update the similarities of neighbors if islocal=1
% y: the final clustering result, cluster indicator vector
% S: learned symmetric similarity matrix
% evs: eigenvalues of learned graph Laplacian in the iterations
% cs: suggested cluster numbers, effective only when the cluster structure is clear

% Ref:
% Feiping Nie, Xiaoqian Wang, Michael I. Jordan, Heng Huang.
% The Constrained Laplacian Rank Algorithm for Graph-Based Clustering.
% The 30th Conference on Artificial Intelligence (\textbf{AAAI}), Phoenix, USA, 2016.


NITER = 100;  %迭代次数
zr = 1e-20;  %阈值 
% lambda = 10;  %参数
% lambda1 = lambda;
% beta=0.001;
% mu=1;

[dim,num] = size(X);    %1024*1440 dim是维数，num是图片个数
summ = zeros(NITER,1); 
S = zeros(size(A0));  % 假设A0代表了S的合适初始尺寸
SU = zeros(size(AU));
[d,n]=size(X);
% X=X';
% c=10;
% T=abs(rand(c,d));
T=1/d*eye(d);
W=abs(rand(d,c));
V=abs(rand(c,n));
%数据图
% A0 = A0-diag(diag(A0));   %diag(diag(A0))hi一个全为0的矩阵
num = size(A0,1);   %
% A0 = lambda*A0;
A10 = (A0+A0')/2;
D10 = diag(sum(A10));
L0 = D10 - A10;  % 拉普拉斯矩阵
% automatically determine the cluster number
[F0, ~, evs]=eig1(L0, num, 0);  % eigs or eig1 ?? %evs为从小大大进行排列的特征值 F0为对应的特征向量
V = F0(:,1:c);%注意，这里的F就是矩阵V，通常NMF分解是X-UV',这里的UV是通过随机生成或者SVD分解得到的，
%但是我这个代码中的U和V是基于数据空间和特征空间中最优相似度矩阵中非零特征值对应的特征向量学习到的U和V
if sum(evs(1:c+1)) < zr%如果从第1个元素到第c+1个元素的累计和小于某个阈值zr
    [clusternum, y]=graphconncomp(sparse(A0)); y = y';%使用graphconncomp函数找出稀疏矩阵A0的连通分量数量clusternum和每个节点所属的连通分量标识y
end;
if sum(evs(1:c)) < zr
    [clusternum, y]=graphconncomp(sparse(A10)); y = y';
    S = A0;
    return;
end;

%特征图
AU = AU-diag(diag(AU));   %diag(diag(A0))hi一个全为0的矩阵
dim = size(AU,1);   %
% A1U = lambda*AU;
A1U = (AU+AU')/2;
D1U = diag(sum(A1U));
L0U = D1U - A1U;  % 拉普拉斯矩阵
% automatically determine the cluster number
[F0U, ~, evsU]=eig1(L0U, dim, 0);  % eigs or eig1 ?? %evs为从小大大进行排列的特征值 F0为对应的特征向量
W = F0U(:,1:c);%注意，我这里的FU就是矩阵U
T=eye(dim)/dim;
% R=abs(rand(c,c));
% V=rand(c,n);

if sum(evsU(1:c+1)) < zr
    [clusternum, y]=graphconncomp(sparse(A10)); y = y';
end;
if sum(evsU(1:c)) < zr
    [clusternum, y]=graphconncomp(sparse(A1U)); y = y';
    SU = AU;
    return;
end;


% U = abs(rand(dim,c));
% F=abs(F);  % V
% F = F.*((mu*X'*U+A10*F)./(mu*F*U'*U+D10*F+eps));
% U = FU.*((X*F+A1U*FU)./(FU*F'*F+D1U*FU+eps));

% 先通过初始的相似图对矩阵F,U，T进行初始的更新
% F = F.*((X'*T'*FU+lambda*A10*F)./(F*FU'*FU+lambda*D10*F+beta*F+eps));
% U = FU.*((T*X*F+lambda1*A1U*FU)./(FU*F'*F+lambda1*D1U*FU+theta*V*FU+eps));
% Ui=sqrt(sum(U.*U,2)+eps);
% d=0.25./(Ui.^3/2);
    Hi=sqrt(sum(V.*V,2)+eps);h=0.5./Hi;H=diag(h);
    Qi=sqrt(sum(W.*W,2)+eps);q=0.5./Qi;Q=diag(q);
    Mi=sqrt(sum(T.*T,2)+eps);m=0.5./Mi;M=diag(m);
    W=W.*(T*X*V'*H+theta*W)./(T*X*X'*T'*W*H+theta*Q*W);
    T=T.*(W*H*V*X'+theta*T)./(W*H*W'*T*X*X'+theta*M*T);
    V=V.*(H'*W'*T*X)./(H*V);
  
   
% 更新R,T
% R=R.*((U'*T*X*F)./(U'*U*R*F'*F+eps));
% T=T.*((U*F'*X')./(T*X*X'+eps));
% %     F = F.*((mu*X'*FU+lambda*A10*F)./(mu*F*FU'*FU+lambda*D10*F+beta*F+eps));
% %     
% %     U = FU.*((X*F+lambda1*A1U*FU)./(FU*F'*F+lambda1*D1U*FU+beta*FU+eps));
%     F=abs(F);
%     U=abs(FU);
for iter = 1:NITER
    
    dist = L2_distance_1(F',F');
    S = zeros(num);
    for i=1:num
        a0 = A0(i,:);
        idxa0 = find(a0>0);
        ai = a0(idxa0);
        di = dist(i,idxa0);
        ad = ai-0.5*lambda*di; 
        S(i,idxa0) = EProjSimplex_new(ad);
    end;
%     S = lambda*S;
    A = S;

    A = (A+A')/2;
    D = diag(sum(A));
    
    L = D-A;
    F_old = V;
    [V, ~, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;  
        F = F_old;
    else
        break;
    end;
    
%     F = F.*((mu*X'*FU+lambda*A*F+beta*F)./(mu*F*FU'*FU+lambda*D*F+beta*F*F'*F+eps));
%     F=abs(F);
    
%特征图更新
distU = L2_distance_1(W',W');
    SU = zeros(dim);
    for i=1:dim
        a0U = AU(i,:);
        idxa0U = find(a0U>0);
        aiU = a0U(idxa0U);
        diU = distU(i,idxa0U);
        adU = aiU-0.5*lambda1*diU; 
        SU(i,idxa0U) = EProjSimplex_new(adU);
    end;
%     SU = lambda*SU;
    AU = SU;
%     AU = lambda*AU;
    AU = (AU+AU')/2;
    DU = diag(sum(AU));
    LU = DU-AU;
    F_oldU = FU;
    [W, ~, evU]=eig1(LU, c, 0);
    evsU(:,iter+1) = evU;

    fn1U = sum(evU(1:c));
    fn2U = sum(evU(1:c+1));
    if fn1U > zr
        lambda1 = 2*lambda1;
    elseif fn2U < zr
        lambda1 = lambda1/2;  
        W = F_oldU;
    else
        break;
    end;
%     F = F.*((X'*T'*FU*R+lambda*A10*F)./(F*R'*FU'*FU*R+lambda*D10*F+beta*F+eps));
%     U = FU.*((T*X*F*R'+lambda1*A1U*FU)./(FU*R*F'*F*R'+lambda1*D1U*FU+beta*FU+theta*V*FU+eps));
% Ui=sqrt(sum(U.*U,2)+eps);
% d=0.25./(Ui.^3/2);
   %%%根据学习到的最优图对F,U,T，R进行更新 这里F就是矩阵V
% F = F.*((X'*T'*U+lambda*A*F)./(F*U'*U+lambda*D*F+beta*F+eps));
% U = FU.*((T*X*F+lambda1*AU*FU)./(FU*F'*F+lambda1*DU*FU+theta*V*U+eps));
% Ui=sqrt(sum(U.*U,2)+eps);
% d=0.25./(Ui.^3/2);
    Hi=sqrt(sum(V.*V,2)+eps);h=0.5./Hi;H=diag(h);
    Qi=sqrt(sum(W.*W,2)+eps);q=0.5./Qi;Q=diag(q);
    Mi=sqrt(sum(T.*T,2)+eps);m=0.5./Mi;M=diag(m);
    W=W.*(T*X*V'*H+theta*W)./(T*X*X'*T'*W*H+theta*Q*W);
    T=T.*(W*H*V*X'+theta*T)./(W*H*W'*T*X*X'+theta*M*T);
    V=V.*(H'*W'*T*X)./(H*V);

% R=R.*((U'*T*X*F)./(U'*U*R*F'*F+eps));
% T=T.*((U*F'*X')./(T*X*X'+eps));
%     F=abs(F);
%     U=abs(U);
%     FU=abs(FU);
    
%     U = U.*((X*F)./(U*F'*F+eps));
%     F = F.*((mu*X'*U+lambda*A*F)./(mu*F*U'*U+lambda*D*F+eps));
    sum1=alpha*(norm(A-S,'fro')^2+trace(W'*LS*W));
    sum2=beta*((norm(AU-SU,'fro')^2)+trace(W'*X'*LW*W'*X));
    sum3=trace(X'*T'*W*H*W'*T*X-2*X'*T'*W*H*V+V'*H*V);
    sum4=theta(trace(W'*Q*W-W'*W)+trace(T'*M*T-T'*T));
%     sum1=(norm(T*X-U*F','fro')^2);
%     sum2=lambda*(F'*L*F);
%     sum3=alpha*norm(A-S,'fro')^2;
%     sum4=lambda1*trace(U'*LU*U);
%     sum5=norm(AU-SU,'fro')^2;
%     sum6=theta*sum(sqrt(Ui));
%     sum7=beta*(trace(eye(c)-F'*F));
%     summ(iter) = sum1 + sum2 + sum3 + sum4 + sum5 + sum6+ sum7; 
      summ(iter) = sum1 + sum2 + sum3 + sum4; 
        if iter>2
        if abs(summ(iter)-summ(iter-1))<0.001
            break
        end
    end
   iter = iter + 1;
end;
end

%[labv, tem, y] = unique(round(0.1*round(1000*F)),'rows');
% [clusternum, y]=graphconncomp(sparse(A)); y = y';
% if clusternum ~= c
%     sprintf('Can not find the correct cluster number: %d', c)
% end;


