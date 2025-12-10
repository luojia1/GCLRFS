% function [fea, gnd] = select_class(k, fea, gnd)
%     nClass = length(unique(gnd));
%     a = randperm(nClass);
%     a = a(1:k)
%     b = [];
%     for i = 1:k
%         temp = find(gnd == a(i));
%         b = [b; temp];
%         clear temp;
%     end
%     fea = fea(b, :);
%     gnd = gnd(b);
%     a = unique(gnd);
%     b = length(a);
%     for i = 1:b
%         idx = gnd==a(i);
%         gnd(idx) = i;
%     end
% end
% 定义函数select_class，输入参数为类别数量k、特征矩阵fea和对应的标签向量gnd
function [fea, gnd] = select_class(k, fea, gnd)
    % 计算gnd中唯一元素的数量，即总类别数
    nClass = length(unique(gnd));
    % 生成一个随机排列的nClass长度的向量
    a = randperm(nClass);
    % 从这个随机向量中选择前k个元素，代表将要选取的类别索引
    a = a(1:k);
    % 初始化一个空数组b，用于存储最终要选取的所有样本的索引
    b = [];
    % 遍历每一个选定的类别
    for i = 1:k
        % 找出gnd中等于当前选定类别a(i)的所有元素的索引
        temp = find(gnd == a(i));
        % 将找到的索引添加到数组b中
        b = [b; temp];
        % 清除临时变量temp
        clear temp;
    end
    % 使用索引b来重新组织特征矩阵fea，只保留选定类别的样本
    fea = fea(b, :);
    % 同样，使用索引b来重新组织标签向量gnd
    gnd = gnd(b);
    % 获取gnd中所有唯一元素（即当前所含的类别标签）
    a = unique(gnd);
    % 计算现在有多少个不同的类别
    b = length(a);
    % 为了简化类别标签，将它们映射为1到b的连续整数
    for i = 1:b
        % 找出所有属于当前类别a(i)的样本
        idx = gnd == a(i);
        % 将这些样本的标签设置为连续的整数i
        gnd(idx) = i;
    end
end
