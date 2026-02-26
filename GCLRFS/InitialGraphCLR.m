function A = InitialGraphCLR(data, m, issymmetric)
%% 
% by Lance, May 2016
% Learning An Initial Graph
% ref: Nie, Feiping, et al. "The Constrained Laplacian Rank Algorithm for Graph-Based Clustering." (2016).
%%

n = size(data,2);   % n是data的列数  data=1024*1440
E = pdist2(data', data').^2;   % 计算data'中每队数据点之间的距离 data'=1440*1024 E=1440*1440
% E的第一行为第一个数据点遍历其他数据点的距离
[SortedE, Ind] = sort(E); % sort函数使得矩阵E的每一列都按照升序排列，其中SortedE排列后的矩阵，Ind是排列后的索引
SortedE = SortedE';   % 转置
Ind = Ind';   % 转置
SortedE = [SortedE(:,2:end),SortedE(:,1)];   % 把全为0的列放在最后一列
Ind = [Ind(:,2:end),Ind(:,1)];   % 对应的索引也要改变
Eim1 = repmat(SortedE(:,m+1),[1,n]);   
Eim1Ind = Ind(:,m+1);

%*****
%make a Indicator Mask to record j<=m or j>m
IndMask = zeros(n); 
for i = 1:n
    IndMask(i, Ind(i,1:m)) = 1;
end
% 对A进行归一化
E_numerator = E.*IndMask;  %取出来前k行的元素
Eim1_numerator = Eim1.*IndMask;  %取出来第k+1个元素

denominator = m*Eim1 - repmat(sum(SortedE(:,1:m),2),[1,n]) + eps;%计算分母
A = (Eim1_numerator - E_numerator)./denominator;

if issymmetric == 1
    A = (A+A')/2;
end;
    