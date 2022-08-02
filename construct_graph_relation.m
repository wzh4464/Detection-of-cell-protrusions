function face_graph = construct_graph_relation(neib,face)
%CONSTRUCT_GRAPH_RELATION construct the face graph
%   neib, face
face_num = size(face,1);
face_graph = zeros(face_num,3);
parfor i=1:face_num
% for i=1:face_num
    for j=1:3
        newpoint = intersect(neib{face(i,mod3(j+1))},neib{face(i,mod3(j+2))});
        newpoint(newpoint==face(i,j))=[];
        temp = face==sort([face(i,mod3(j+1)),face(i,mod3(j+2)),newpoint]);
%         find(temp[1,])
        face_graph(i,j)=find(sum(temp,2)==3);
    end
    
%     face_graph(i)=[]
end
end

function re = mod3(x)
re = mod(x,3);
if re
    return
else
    re = re+3;
end
end