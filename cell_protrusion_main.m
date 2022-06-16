% @Author: Huang Roger
% @Date:   2022-06-15 15:26:47
% @Last Modified by:   WU Zihan
% @Last Modified time: 2022-06-15 19:04:14
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
%cell name ABalapa
% V = niftiread('datasets\data\191108plc1p1_053_segCell.nii');
% index = unique(V);
% V1=V;
% V1(find(V~=65))=0;
% % figure
% % h=volshow(V1);
% [X,Y,Z]=voxel2XYZ(V1);
% %     figure
% %     pcshow([X(:),Y(:),Z(:)]);
% k = boundary(X(:),Y(:),Z(:));%cell surface mesh
% bind=unique(k);
% % figure
% % pcshow([X(bind),Y(bind),Z(bind)],"MarkerSize",100);
% figure
% trisurf(k,X(:),Y(:),Z(:),'Facecolor','cyan','FaceAlpha',0.1)
% axis equal
% points=[X(bind),Y(bind),Z(bind)];

%%
%导入三维点，开始计算凸包，计算距离
load('k.mat')
load('ABalapa.mat')
X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[k1,av1] = convhull(X,Y,Z);
figure
trisurf(k1,X,Y,Z,'FaceColor','cyan')
axis equal
title('convexhull')


% 求解Convex Hull平面的方程 Ax + By + Cz = 1
solve=[];
for i=1:size(k1,1)
    A=[X(k1(i,1)) Y(k1(i,1)) Z(k1(i,1))
    X(k1(i,2)) Y(k1(i,2)) Z(k1(i,2))
    X(k1(i,3)) Y(k1(i,3)) Z(k1(i,3))];
    B=[1;1;1];
    solve(i,1:3)=(A\B)';
end

% 寻找射影点
for i=1:length(X)
    for j=1:size(k1,1)
        temp(i,j)=abs(solve(j,1:3)*[X(i) Y(i) Z(i)]'-1)/norm(solve(j,1:3));
    end
    [dist(i),plane(i)]=min(temp(i,:));
    t(i)=(solve(plane(i),1:3)*[X(i) Y(i) Z(i)]'-1)/(norm(solve(plane(i),1:3))^2);
    proj(i,1:3)=[X(i) Y(i) Z(i)]-solve(plane(i),1:3)*t(i);
end

% %把老的mesh索引转到new_mesh上面
% for i=1:size(k,1)
%     for j=1:3
%         new_mesh(i,j)=find(bind==k(i,j));
%     end
% end
% k=new_mesh;
% tri=delaunay(points(:,1),points(:,2));
% trimesh(tri,points(:,1),points(:,2),points(:,3));
%合并点，找protrusions
coinPointsInd=find(dist==0);%原图与convexhull重合点的索引
protrusions=[];
for i=1:length(coinPointsInd)
    [row,col]=find(k==coinPointsInd(i));%找dist=0的点所在的同一个mesh中的点
    temp_k=k(row,:);%找第i点所在的同一个mesh中的点
    pointsInd_temp=unique(temp_k(find(temp_k~=coinPointsInd(i))));%与第i个点同在一个mesh的所有点
    zeroNum_temp=length(find(dist(pointsInd_temp)==0));%距离为0的点的数量
    nonzeroNum_temp=length(find(dist(pointsInd_temp)~=0));
    if zeroNum_temp<3%||nonzeroNum_temp/(zeroNum_temp+nonzeroNum_temp)>0.1 %nonzeroNum_temp>=5
        protrusions=[protrusions;coinPointsInd(i)];
    end
end

mycolor=zeros(1,size(dist,2));
mycolor(protrusions)=200;
figure;
pcshow(points,mycolor,"MarkerSize",100);
% title('xyz_protrusions')


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


% 把投影点转化为球坐标
tep=proj-repmat(center,size(X,1),1);
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

% 延拓后的 mesh
tri=delaunay(distMpro(:,1),distMpro(:,2));
% trimesh(tri,distMpro(:,1),distMpro(:,2),distMpro(:,3));
[zx,zy] = trigradient(distMpro(:,1),distMpro(:,2),distMpro(:,3),tri); % 梯度
z = sqrt(zx.^2+zy.^2); % 标量梯度

figure;
pcshow(distMpro,z,"MarkerSize",100);
title('梯度温度')



z_center=z(size(X,1)*4+1:size(X,1)*5,:);% 原距离的标量梯度
% figure
% sorted=sort(z_center);
% plot(1:length(z_center),sorted);

% figure;
% pcshow(distM,z_center,"MarkerSize",100);
% title('distMpro')

figure;
scatter3(distM(:,1),distM(:,2),z_center);
title('scalar gradient')

% threshold = 0.05;
% I=find(z_center<threshold & z_center~=0);
% 
% figure;
% scatter3(distM(I,1),distM(I,2),distM(I,3));

threshold = 0.05;
I=find(z<threshold & z~=0); % index of 梯度小 -> 去掉未连接的点
J=find(distMpro(:,3)~=0); % index of 高度小 -> 凸包上的点
I=intersect(I,J);

% for i=1:length(I)
%     find(k1==I(i))
%     if I(i)
% end

figure;
scatter3(distMpro(I,1),distMpro(I,2),distMpro(I,3));
title('zero points height')
% find(tri=)
% hold on
% quiver3(distMpro(:,1),distMpro(:,2),distMpro(:,3),zx,zy,zeros(size(zx,1),1))
% hold off



% map=dist;
% colormap(map)
% figure
% plot3(X,Y,Z)

figure;
pcshow(points,dist,"MarkerSize",40);
title('xyz_dist')

% plot(1:length(dist),dist)






