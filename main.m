clear
close all
addpath(genpath(pwd))
threshold = 0.05;

%% mat from nature
load('surface_1_1.mat');
points = surface.vertices;
ori_faces = surface.faces;
ori_faces = sort(ori_faces,2);
figure
pcshow(points,"MarkerSize",100);
figure
trisurf(ori_faces,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)

%%
%导入三维点，开始计算凸包，计算距离
X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[cv_mesh,~] = convhull(points(:,1),points(:,2),points(:,3));
cv_face_num = size(cv_mesh,1);
% figure
% trisurf(cv_mesh,X,Y,Z,'FaceColor','cyan')
% axis equal
% title('convexhull')

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


%% Union the connected parts
[disjoint_set,size_DS]=union_zero_fast(best,neib);
protrusion_rep = unique(disjoint_set);
protrusion_rep(protrusion_rep==0)=[];

% for i=1:length(protrusion_rep)
%     construct_graph_relation

%%
% face_neib = construct_graph_relation(neib,ori_faces); %% face_num * 3
TR = triangulation(ori_faces,points);
face_neib = neighbors(TR);