clc;
clear;
load('lung_discrete.mat');
X=X;
gnd=Y;
% X=fea;
% gnd=gnd;
% % X=data;
% gnd=labels;
% nClusts = length(unique(gnd));
c = length(unique(gnd));
NIter=30;
resultsMatrix = [];
m=5;%邻居数
A0 = InitialGraphCLR(X', m, 1);
AU = InitialGraphCLR(X, m, 1);
%%%%%%%%%%%%%%%%%%%%%%%%实验开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AAAA=[];
%lung_discrete:60,-2,-1,+2; ecoli_uni:100,-1,-4,+1; %UMIST_fac:70,+0,+2,+1
%imm40:50,-4,-3,-3; Yale_32:70,-3,+4,-4; ORL:100,-2,+3,-3; YaleB_32x32:20,-4,+3,-2
%k1a:30,-3,0,-3,30; yeast:100.-2,+4,+2; jaffe:90,-4,-3,+2; Isolet:100,+1,-1,-1; COIL20:90,-2,+4,+4
p=60;
alpha=1e-2;
beta=1e-1;
theta=1e+2;
% pp=[100 90 80 70 60 50 40 30 20];
% alp=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
% be=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
% the=[1e-4 1e-3 1e-2 1e-1 1e+0 1e+1 1e+2 1e+3 1e4];
% for ppp=1:9
%     p=pp(:,ppp);
%     for aaa=1:9
%         alpha=alp(:,aaa);
%         for bbb=1:9
%             beta=be(:,bbb);
%                 for ccc=1:9
%                     theta=the(:,ccc);
                    rng('default');
                    tic
                   [X_new,obj]= GCLRFS(X' ,A0, AU, c, alpha, beta, theta,p);
                    toc
                    for i=1:40   
                        label=litekmeans(X_new',c,'MaxIter',100,'Replicates',10);
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
                    T=toc;
                    resultsMatrix = [resultsMatrix; alpha, beta, theta, p, AAAA, BBBB, AAAA1, BBBB2, T];
                    % xlswrite('C:\Users\Administrator\Desktop\实验结果\GCLRFS\lung_discrete.xlsx', resultsMatrix, 'Sheet1');
%                 end   
%             end
%         end
% end