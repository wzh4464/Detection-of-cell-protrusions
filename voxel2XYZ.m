function [X,Y,Z]=voxel2XYZ(V)
%三维图像转XYZ坐标
    idx=find(V);
    [X,Y,Z]=ind2sub(size(V),idx);
end