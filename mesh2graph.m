function G = mesh2graph(k)
%mesh2graph transform mesh to graph
%   easy to find neighbours
    E=[];
    for i = 1 : size(k,1)
        E=[E;k(i,1),k(i,2)];
        E=[E;k(i,1),k(i,3)];
        E=[E;k(i,2),k(i,3)];
    end
    E = sort(E')'
    E = unique(E,"rows");
    G = graph(E(:,1),E(:,2));
end