function face_extend_index=extends_protrusion_point(TR,face_neib,point_index,num,depth)
%%
%     face_extend_index=[];
%     tri_id = vertexAttachments(TR,point_index);%点的所有相邻面
%     nei_1 = unique(face_neib(tri_id{1,1},:));%相邻面的所有相邻面
%     face_extend_index=union(tri_id{1,1},nei_1);
% %     nei_1 = setxor(tri_id{1,1},nei_1);%相邻面的所有相邻面
%     face_extend_index=[face_extend_index nei_1'];
%     for i=1:num-1
%         nei_2 = unique(face_neib(nei_1,:));
%         nei_2 = setxor(nei_1,nei_2);%相邻面的所有相邻面
%         nei_1=nei_2;
%         face_extend_index=[face_extend_index nei_2'];
%     end

%%
%Roger没带自适应的代码
    tri_id = vertexAttachments(TR,point_index);%点的所有相邻面
    nei_1 = unique(face_neib(tri_id{1,1},:));%相邻面的所有相邻面
    searched_index = tri_id{1,1};%查找过的面索引
    nei_1 = setdiff(nei_1,searched_index);%还没查找过的面索引
    searched_index = [searched_index nei_1'];
    for i=1:num-1
        nei_2 = unique(face_neib(nei_1,:));
        nei_2 = setdiff(nei_2,searched_index);
        searched_index = [searched_index nei_2'];
        nei_1=nei_2;
    end
    face_extend_index=searched_index;
 

%%
%吴子涵带自适应的代码（失败）
% tri_id = vertexAttachments(TR,point_index);%点的所有相邻面
%     nei_1 = unique(face_neib(tri_id{1,1},:));%相邻面的所有相邻面
%     searched_index = tri_id{1,1};%查找过的面索引
%     nei_1 = setdiff(nei_1,searched_index);%还没查找过的面索引
%     searched_index = [searched_index nei_1'];
%     searched_box=zeros(size(face_neib,1),1);
%     searched_box(searched_index)=1;
%     for i=1:num-1
%         nei_2=zeros(size(face_neib,1),1);
%         nei_2(nei_1)=1;
%         for j=1:length(nei_1)
%             for k=1:3
%                 if depth(face_neib(nei_1(j),k))>depth(nei_1(j)) && searched_box(face_neib(nei_1(j),k))==0
%                     nei_2(face_neib(nei_1(j),k))=1;
%                     searched_box(face_neib(nei_1(j),k))=1;
%                 end
%             end
%         end
%         nei_1 = find(nei_2==1);
% 
%     end
%     face_extend_index=find(searched_box==1);
end