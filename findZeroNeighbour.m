function set = findZeroNeighbour(i, set, G, dist, threshold)
%findZeroNeighbour 寻找最大可能0连通
%   递归算法
%   i for source node
%   set for result
%   G for graph from mesh
%   dist for distance check
    nei = neighbors(G,i);

%%
%wuzihan
    for j = 1:length(nei)
        if ~ismember(nei(j),set) & dist(nei(j)) == 0
           set = [set,nei(j)];
           set = findZeroNeighbour(nei(j), set, G, dist, threshold);
        end
    end

end