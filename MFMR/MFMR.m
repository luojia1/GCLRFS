function [X_new] = DRFSMFMR(X, alpha,beta, lamda, p, maxIter,c)

% rng('default');
[n,m]=size(X);
I1=ones(m,m);
I2=ones(c,c);
W=abs(rand(m,c));
H=abs(rand(c,m));
for iter=1:maxIter
 
    W=W.*sqrt((X'*X*H'+beta*W)./(X'*X*W*H*H'+alpha*X'*X*W*I2*I2'+beta*I1*W));
    H=H.*sqrt((W'*X'*X+ lamda*H)./(W'*X'*X*W*H+lamda*H*I1));
    
    obj(iter)=trace((X*W*H-X)*(X*W*H-X)')+0.5*alpha*trace(I2'*W'*X'*X*W*I2)+0.5*beta*(trace(I1*W*W')-trace(W*W'))+0.5*lamda*(trace(I1*H'*H)-trace(H'*H));

    if iter>2
        if abs(obj(iter)-obj(iter-1))/obj(iter-1)<1e-6
            break
        end
    end
end

score=sqrt(sum(W.*W,2));
[~, idx]=sort(score,'descend');
X_new=X (:,idx(1:p));

end