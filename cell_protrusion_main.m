clear
close all
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
load('Jun15cell.mat')
X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[k1,av1] = convhull(X,Y,Z);
figure
trisurf(k1,X,Y,Z,'FaceColor','cyan')
axis equal
title('convexhull')

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
% for i=1:size(points,1)
%     r(i)=norm(points(i,:)-center);
%     theta(i)=real(asin(Z(i)/r(i)));
%     phi(i)=real(acos(X(i)/r(i)/cos(theta(i))));
% end
tep=points-repmat(center,size(X,1),1);
[phi,theta,~] = cart2sph(tep(:,1),tep(:,2),tep(:,3));

figure
pcshow([theta,phi,dist'],"MarkerSize",40)
title('dist_theta_phi')
% figure
% surf(theta',phi',dist')

%三维点生成mesh
% tri=delaunay(theta,phi);
% figure
% trimesh(tri,theta,phi,dist);
% TO = triangulation(tri,theta,phi,dist');

%找谷
%周期性延拓
distM=[theta phi dist'];
% distMpro = repmat(distM,9,1);
distMpro=[];
for i=-1:1
    for j=-1:1
        tmp=distM;
        tmp=[tmp(:,1)+i*pi,tmp(:,2)+j*2*pi,tmp(:,3)];
        distMpro=[distMpro;tmp];
    end
end
% tri=delaunay(distMpro);
% trimesh(tri,distMpro);
% fx = gradient(dist)
tri=delaunay(distMpro(:,1),distMpro(:,2));
% trimesh(tri,distMpro(:,1),distMpro(:,2),distMpro(:,3));
[zx,zy] = trigradient(distMpro(:,1),distMpro(:,2),distMpro(:,3),tri); 
z = sqrt(zx.^2+zy.^2);

figure;
pcshow(distMpro,z,"MarkerSize",100);
title('distMpro')



z_center=z(size(X,1)*4+1:size(X,1)*5,:);% 原距离的标量梯度
figure
sorted=sort(z_center);
plot(1:length(z_center),sorted);

% figure;
% pcshow(distM,z_center,"MarkerSize",100);
% title('distMpro')

figure;
scatter3(distM(:,1),distM(:,2),z_center);
title('distM')

% threshold = 0.05;
% I=find(z_center<threshold & z_center~=0);
% 
% figure;
% scatter3(distM(I,1),distM(I,2),distM(I,3));

threshold = 0.05;
I=find(z<threshold & z~=0);
J=find(distMpro(:,3)~=0);
I=intersect(I,J);

for i=1:length(I)
    find(k1==I(i))
    if I(i)
end

figure;
scatter3(distMpro(I,1),distMpro(I,2),distMpro(I,3));
% find(tri=)
% hold on
% quiver3(distMpro(:,1),distMpro(:,2),distMpro(:,3),zx,zy,zeros(size(zx,1),1))
% hold off



% map=dist;
% colormap(map)
% figure
% plot3(X,Y,Z)

figure;
pcshow([X(:),Y(:),Z(:)],dist,"MarkerSize",40);
title('xyz_dist')

% plot(1:length(dist),dist)


%三维图像转XYZ坐标
function [X,Y,Z]=voxel2XYZ(V)
    idx=find(V);
    [X,Y,Z]=ind2sub(size(V),idx);
end



