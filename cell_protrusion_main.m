% @Author: Huang Roger
% @Date:   2022-06-15 15:26:47
% @Last Modified by:   WU Zihan
% @Last Modified time: 2022-06-15 19:04:14
clear
close all
addpath(genpath(pwd))
threshold = 0.05;

%% mat from nature
load('surface_1_1.mat');
vertices = surface.vertices;
ori_faces = surface.faces;
points = vertices;
% for i = 1:size(raw_points,1)
%     for j = 1:size(raw_points,2)
%         raw_points(i,j) = raw_points(i,j)+0.01*rand();
%     end
% end
% k = boundary(raw_points);
% bind=unique(k);
% raw_points=[X(bind),Y(bind),Z(bind)];
figure
pcshow(points,"MarkerSize",100);
figure
trisurf(ori_faces,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)

%% load from tif
% close all
% V = tiffreadVolume('1_CAM01_000000.tif');
% figure
% slice(double(V),size(V,2)/2,size(V,1)/2,size(V,3)/2)
% colormap gray 
% shading interp
% thr = [0.01]
% for i = -5:1:0
%     J = imbinarize(V,'adaptive','ForegroundPolarity','bright','Sensitivity',10^(i));
%     % J = imbinarize(V,'global')*255;
%     % J = 255-V;
%     % figure
%     % slice(double(J),size(J,2)/2,size(J,1)/2,size(J,3)/2)
%     % colormap gray
%     % shading interp
%     figure
%     volshow(J);
% end
% BW = imbinarize(V)*255;
% h=volshow(BW);
% V = tiffreadVolume('1_CAM01_000000.tif','PixelRegion',{rows,columns,slices});
% v1 = V(:,:,1);
% b1 = BW(:,:,75);
% [x,y,z]=find(BW == 0);
%% nii to XYZ(boundary)
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

%% whole XYZ to boundary
% % points read todo
% load Sample06_024_ABalp.txt
% points = Sample06_024_ABalp;
% for i = 1:size(points,1)
%     for j = 1:size(points,2)
%         points(i,j) = points(i,j)+0.01*rand();
%     end
% end
% k = boundary(points);
% bind=unique(k);
% X=points(:,1);
% Y=points(:,2);
% Z=points(:,3);
% % points=[X(bind),Y(bind),Z(bind)];
% figure
% pcshow([X(bind),Y(bind),Z(bind)],"MarkerSize",100);
% figure
% trisurf(k,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)
% points=[X(bind),Y(bind),Z(bind)];

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
% load('k.mat')
% load('ABalapa.mat')
X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[cv_mesh,~] = convhull(points(:,1),points(:,2),points(:,3));
cv_face_num = size(cv_mesh,1);
figure
trisurf(cv_mesh,X,Y,Z,'FaceColor','cyan')
axis equal
title('convexhull')

%%
% 求解Convex Hull平面的方程 Ax + By + Cz = 1
solve=zeros(cv_face_num,3);
for i=1:cv_face_num
    A=[X(cv_mesh(i,1)) Y(cv_mesh(i,1)) Z(cv_mesh(i,1))
        X(cv_mesh(i,2)) Y(cv_mesh(i,2)) Z(cv_mesh(i,2))
        X(cv_mesh(i,3)) Y(cv_mesh(i,3)) Z(cv_mesh(i,3))];
    B=[1;1;1];
    solve(i,:)=(A\B)';
end

%求距离和投影
point_num = size(points,1);
tmp = abs(solve*points'-1);
nrm = (sum(abs(solve).^2,2)).^(1/2);%一行一行（逐行）norm
% temp = tmp./repmat(nrm,1,point_num);%
temp = tmp./nrm;
[dist,plane]=min(temp);%细胞表面的点到最近的一个凸包面的距离
% nrm1 = (sum(abs(solve(plane,:)).^2,2)).^(1/2);%一行一行（逐行）norm
t=(sum(solve(plane,:).*points,2)-1)./nrm(plane).^2;%除了两次
proj=points-solve(plane,:).*t;

%%
neib=cell(point_num,1);
p = parpool(8);
parfor i=1:point_num
    [temp1,~] = find(ori_faces==i);
    neib{i} = unique(ori_faces(temp1,:));%找每个点的相同面上的点的索引
    neib{i}(neib{i}==i) = [];
end
delete(p);
%% 
best = zeros(1,5000);
j=0;
for i=1:point_num
    [temp3,~] = min(dist(neib{i}));%每个点相同面上的距离最小的点
%     diff(i)=dist(i)-dist(temp4);
    if dist(i)<=temp3
        j = j+1;
        best(j)=i;
    end
end
best = best(1:j);
figure
pcshow(points(best,:))

% %% 把老的mesh索引转到new_mesh上面
% for i=1:size(ori_faces,1)
%     for j=1:3
%         new_mesh(i,j)=find(bind==ori_faces(i,j));
%     end
% end
% ori_faces=new_mesh;

%%
% tri=delaunay(points(:,1),points(:,2));
% trimesh(tri,points(:,1),points(:,2),points(:,3));

% %合并点，找protrusions
% coinPointsInd=find(dist==0);%原图与convexhull重合点的索引
% protrusions=[];
% 
% G = mesh2graph(ori_faces);
% zeroNeibour={};
% for i = 1:length(coinPointsInd)
%     zeroNeibour{coinPointsInd(i)}= [];
%     zeroNeibour{coinPointsInd(i)}= findZeroNeighbour(coinPointsInd(i),zeroNeibour{coinPointsInd(i)},G,dist,threshold);
% end

% for i=1:length(coinPointsInd)
%     [row,col]=find(k==coinPointsInd(i));%找dist=0的点所在的同一个mesh中的点
%     temp_k=k(row,:);%找第i点所在的所有mesh面
%     pointsInd_temp=unique(temp_k(temp_k~=coinPointsInd(i)));%与第i个点同在一个mesh的所有点
%     zeroNum_temp=length(find(dist(pointsInd_temp)==0));%距离为0的点的数量
%     nonzeroNum_temp=length(find(dist(pointsInd_temp)~=0));
%     if zeroNum_temp<3%||nonzeroNum_temp/(zeroNum_temp+nonzeroNum_temp)>0.1 %nonzeroNum_temp>=5
%         protrusions=[protrusions;coinPointsInd(i)];
%     end
% end
% %% drawing
% mycolor=zeros(1,size(dist,2));
% for i=1:size(zeroNeibour,2)
%     mycolor(zeroNeibour{i})=i;
% end
% 
% 
% % mycolor=zeros(1,size(dist,2));
% % mycolor(protrusions)=200;
% figure;
% pcshow(points,mycolor,"MarkerSize",100);
% % title('xyz_protrusions')
% 
% figure
% trisurf(ori_faces,X,Y,Z,mycolor)
% axis equal
% % title('convexhull')
% 
% 
% % convhullIdx=unique(k1);
% % Cpoints=points(convhullIdx,:);
% 
% %k1叫不到原来points里的冗余点
% 
% %球 theta phi
% center=sum(points)/size(points,1);
% % for i=1:size(points,1)
% %     r(i)=norm(points(i,:)-center);
% %     theta(i)=real(asin(Z(i)/r(i)));
% %     phi(i)=real(acos(X(i)/r(i)/cos(theta(i))));
% % end
% 
% 
% % 把投影点转化为球坐标
% tep=proj-repmat(center,size(X,1),1);
% [phi,theta,~] = cart2sph(tep(:,1),tep(:,2),tep(:,3));
% 
% figure
% pcshow([theta,phi,dist'],"MarkerSize",40)
% title('dist_theta_phi')
% % figure
% % surf(theta',phi',dist')
% 
% %三维点生成mesh
% % tri=delaunay(theta,phi);
% % figure
% % trimesh(tri,theta,phi,dist);
% % TO = triangulation(tri,theta,phi,dist');
% 
% %找谷
% %周期性延拓
% distM=[theta phi dist'];
% % distMpro = repmat(distM,9,1);
% distMpro=[];
% for i=-1:1
%     for j=-1:1
%         tmp=distM;
%         tmp=[tmp(:,1)+i*pi,tmp(:,2)+j*2*pi,tmp(:,3)];
%         distMpro=[distMpro;tmp];
%     end
% end
% % tri=delaunay(distMpro);
% % trimesh(tri,distMpro);
% % fx = gradient(dist)
% 
% % 延拓后的 mesh
% tri=delaunay(distMpro(:,1),distMpro(:,2));
% % trimesh(tri,distMpro(:,1),distMpro(:,2),distMpro(:,3));
% [zx,zy] = trigradient(distMpro(:,1),distMpro(:,2),distMpro(:,3),tri); % 梯度
% z = sqrt(zx.^2+zy.^2); % 标量梯度
% 
% figure;
% pcshow(distMpro,z,"MarkerSize",100);
% title('梯度温度')
% 
% 
% 
% z_center=z(size(X,1)*4+1:size(X,1)*5,:);% 原距离的标量梯度
% % figure
% % sorted=sort(z_center);
% % plot(1:length(z_center),sorted);
% 
% % figure;
% % pcshow(distM,z_center,"MarkerSize",100);
% % title('distMpro')
% 
% figure;
% scatter3(distM(:,1),distM(:,2),z_center);
% title('scalar gradient')
% 
% % threshold = 0.05;
% % I=find(z_center<threshold & z_center~=0);
% %
% % figure;
% % scatter3(distM(I,1),distM(I,2),distM(I,3));
% 
% threshold = 0.05;
% I=find(z<threshold & z~=0); % index of 梯度小 -> 去掉未连接的点
% J=find(distMpro(:,3)~=0); % index of 高度小 -> 凸包上的点
% I=intersect(I,J);
% 
% % for i=1:length(I)
% %     find(k1==I(i))
% %     if I(i)
% % end
% 
% figure;
% scatter3(distMpro(I,1),distMpro(I,2),distMpro(I,3));
% title('zero points height')
% % find(tri=)
% % hold on
% % quiver3(distMpro(:,1),distMpro(:,2),distMpro(:,3),zx,zy,zeros(size(zx,1),1))
% % hold off
% 
% 
% 
% % map=dist;
% % colormap(map)
% % figure
% % plot3(X,Y,Z)
% 
% figure;
% pcshow(points,dist,"MarkerSize",40);
% title('xyz_dist')
% 
% % plot(1:length(dist),dist)






