function output = smaller_dist(input, neib, dist)
%SMALLER_DIST 
%   输入输出都是索引
%   neib元胞
%   dist距离
    [mdist, rank] = min(dist(neib{input}));
    if mdist<=dist(input)
        output = neib{input}(rank);
    else
        output = input;
    end
end

