function face_extend_index=extend_protrusion_point(ori_faces,point_index,num,neib)
    if num==0
        face_extend_index=[];
        return;
    end
    face_extend_index=[];
    index=neib{point_index};%寻找包含邻居的点
    index(index==point_index)=[];
    for k=1:length(index)
        [face_index_temp,~]=find(ori_faces==index(k));%包含某点的全部面
        for l=1:length(face_index_temp)
            if ~ismember(face_index_temp(l),face_extend_index)%find(face_extend_index~=face_index_temp(l))
                face_extend_index=[face_extend_index face_index_temp(l)];
            end
        end
        face_extend_index1=extend_protrusion_point(ori_faces,index(k),num-1,neib);
        face_extend_index=[face_extend_index face_extend_index1];
    end
end