clear
addpath(genpath(pwd))

%%
% file = gunzip('datasets\Sample07_088_segCell.nii.gz');%change
% V = niftiread('datasets\data\Sample04_004_segCell.nii');
% index = unique(V);
% 
% %3 1608 2796 3931 5543 6591
% % for i=2:length(index)
%     V1=V;
% %     V1(find(V~=index(i)))=0;
%     V1(find(V~=3))=0;
%     [X,Y,Z]=voxel2XYZ(V1);
% %     V2 = imresize3(V1, 1.5);
% %     B=edge3(V1);
%     [hh,TT,X1,Y1,Z1,CC,AA] = voxelSurf(V1,false);
%     figure
%     pcshow([X(:),Y(:),Z(:)]);
%     k = boundary(X(:),Y(:),Z(:));
%     figure
%     trisurf(k,X(:),Y(:),Z(:),'Facecolor','red','FaceAlpha',0.1)
% %     figure
% %     g=volshow(V1);
% % end

%%
%导入三维点，开始计算凸包，计算距离
load('points.mat')
X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[k1,av1] = convhull(X,Y,Z);
figure
trisurf(k1,X,Y,Z,'FaceColor','cyan')
axis equal

solve=[];
for i=1:size(k1,1)
    A=[X(k1(i,1)) Y(k1(i,1)) Z(k1(i,1))
    X(k1(i,2)) Y(k1(i,2)) Z(k1(i,2))
    X(k1(i,3)) Y(k1(i,3)) Z(k1(i,3))];
    B=[1;1;1];
    solve(i,1:3)=(A\B)';
end

for i=1:length(X)
    for j=1:size(k1,1)
        temp(i,j)=abs(solve(j,1:3)*[X(i) Y(i) Z(i)]'-1)/norm(solve(j,1:3));
    end
    [dist(i),plane(i)]=min(temp(i,:));
    t(i)=(solve(plane(i),1:3)*[X(i) Y(i) Z(i)]'-1)/(norm(solve(plane(i),1:3))^2);
    proj(i,1:3)=[X(i) Y(i) Z(i)]-solve(plane(i),1:3)*t(i);
end

% convhullIdx=unique(k1);
% Cpoints=points(convhullIdx,:);

%k1叫不到原来points里的冗余点

%球 theta phi
center=sum(points)/size(points,1);
for i=1:size(points,1)
    r(i)=norm(points(i,:)-center);
    theta(i)=asin(Z(i)/r(i));
    phi(i)=real(acos(X(i)/r(i)/cos(theta(i))));
end

figure
pcshow([theta',phi',dist'],"MarkerSize",40)
% figure
% surf(theta',phi',dist')

%三维点生成mesh
tri=delaunay(theta,phi);
trimesh(tri,theta,phi,dist);
TO = triangulation(tri,theta',phi',dist');

%找谷
%周期性延拓
distM=[theta' phi' dist'];
% distMpro = repmat(distM,9,1);
distMpro=[];
for i=-1:1
    for j=-1:1
        tmp=distM;
        tmp=[tmp(:,1)+i*pi,tmp(:,2)+j*pi,tmp(:,3)];
        distMpro=[distMpro;tmp];
    end
end
tri=delaunay(distMpro);
trimesh(tri,distMpro);
% fx = gradient(dist)




% map=dist;
% colormap(map)
% figure
% plot3(X,Y,Z)

figure;
pcshow([X(:),Y(:),Z(:)],dist,"MarkerSize",40);

% plot(1:length(dist),dist)


%三维图像转XYZ坐标
function [X,Y,Z]=voxel2XYZ(V)
    idx=find(V);
    [X,Y,Z]=ind2sub(size(V),idx);
end



