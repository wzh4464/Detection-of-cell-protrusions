function G=meshs2graph(TR)

e=TR.ConnectivityList;
temp=[e(:,1),e(:,2)];
temp=[temp; e(:,1),e(:,3)];
temp=[temp; e(:,2),e(:,3)];
E = unique(temp,'rows');

G = graph(E(:,1), E(:,2));
end