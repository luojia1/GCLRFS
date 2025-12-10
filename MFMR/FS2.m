function [X_new,W,V,iter,obj,D] = FS2(X,Y,D,alpha,lamda,p,maxIter)

C=length(unique(Y));
[m,n]=size(X);
eta=0.1;

[~,V]=litekmeans(X,C);V=V'; 

U=eye(m);

W=(X*X'+eta*eye(m))\(X*V);
W=max(W,1e-8);
X=X';
for iter=1:maxIter
%     G=sqrt(diag(sum(W'*W,2)));
%     W=W*pinv(G);
%     H=G*H;
 
    W=W.*((X'*D'*V)./(X'*D'*D*X*W+alpha*U*W+eps));
    V=V.*((D*X*W+ 2*lamda*D*X*X'*D'*V)./(V+ 2*V*V'*V+eps));
%     D=D.*((V*W'*X'+ 2*lamda*V*V'*D*X*X')./(D*X*W*W'*X'+ 2*lamda*lamda*D*X*X'*D'*D*X*X'+eps));
    
    Wi = sqrt(sum(W.*W,2)+eps);
    u = 0.5./Wi;
    U = diag(u);

    
    obj(iter)=trace((D*X*W-V)*(D*X*W-V)')+trace((V*V'-lamda*D*X*X'*D')*(V*V'-lamda*D*X*X'*D')')+alpha*sum(Wi);

    if iter>2
        if abs(obj(iter)-obj(iter-1))/obj(iter-1)<1e-6
            break
        end
    end
end

score= sqrt(sum(W.*W,2));
[res, idx] = sort(score,'descend');
X_new = X (:,idx(1:p));

end