clc;
clear;
load('Isolet.mat');
X=fea;
gnd=gnd;
nClusts = length(unique(gnd));
NIter=30;
resultsMatrix = [];
%%%%%%%%%%%%%%%%%%%%%%%%实验开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AAAA=[];
pp=[ 40 30 20];
% pp=[20 30 40 50 60 70 80 90 100];
% pp=[100 90 80 70 60 50 40 30 20];
alp=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
be=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
the=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
lam=[0];
for ppp=1:9
    p=pp(:,ppp);
    for aaa=1:9
        alpha=alp(:,aaa);
        for bbb=1:9
            beta=be(:,bbb);
            for ccc=1:1
                lamda=lam(:,ccc);
                for ddd=1:9
                    theta=the(:,ddd);
                    rng('default');
                    tic
                    [X_new] = DLDSR(X',alpha,beta,theta,p);
                    toc
                    for i=1:20   
                        label=litekmeans(X_new',nClusts,'MaxIter',100,'Replicates',10);
                        newres = bestMap(gnd,label);
                        AC = length(find(gnd == newres))/length(gnd);
                        MIhat=MutualInfo(gnd,label);
                        result(i,:)=[AC,MIhat];
                    end
                    
                    for j=1:2
                        a=result(:,j);
                        ll=length(a);
                        temp=[];
                        for i=1:ll
                            if i<ll-18
                                b=sum(a(i:i+19));
                                temp=[temp;b];
                            end
                        end
                        [e,f]=max(temp);
                        e=e./20;
                        MEAN(j,:)=[e,f];
                        STD(j,:)=std(result(f:f+19,j));
                        rr(:,j)=sort(result(:,j));
                        BEST(j,:)=rr(end,j);
                    end
                    AAAA=MEAN(1,1);
                    BBBB=MEAN(2,1);
                    AAAA1=STD(1,1);
                    BBBB2=STD(2,1);
                    disp('STD = ')
                    disp(STD)
                    disp('BEST = ')
                    disp(BEST)
                    disp('MEAN = ')
                    disp(MEAN)
                    resultsMatrix = [resultsMatrix; alpha, beta, theta,lamda,p,AAAA, BBBB,AAAA1,BBBB2];
                    xlswrite('C:\Users\Administrator\Desktop\DLDSRisolet70以后.xlsx', resultsMatrix, 'Sheet2');
                end   
            end
        end
    end
end