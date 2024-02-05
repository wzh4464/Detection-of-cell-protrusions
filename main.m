%% init and set path
clear
close all
addpath(genpath(pwd))
threshold = 0.05;

%% mat from nature

% surface is a mesh struct
% surface is a struct with fields: vertices, faces
% vertices is a n*3 matrix, n is the number of points
% each row is a point, the first column is x, the second column is y, the third column is z
% faces is a m*3 matrix, m is the number of faces
% each row is a face, each column is the index of the point in vertices

load('surface_1_1.mat');

points = surface.vertices;
ori_faces = surface.faces;
ori_faces = sort(ori_faces,2);

% figure
% use pcshow to show the point cloud
% pcshow(points,"MarkerSize",100);

% figure
% use trisurf to show the surface as a mesh
% trisurf(ori_faces,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)

%% calculate the convex hull

X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[cv_mesh,~] = convhull(X,Y,Z); % cv_mesh is the index of the points in the convex hull
cv_face_num = size(cv_mesh,1); % number of faces of the convex hull

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

%% calculate the distance and projection
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
num_cpu = feature('numCores');
p = parpool(num_cpu);
parfor i=1:point_num
    [temp1,~] = find(ori_faces==i);
    neib{i} = unique(ori_faces(temp1,:));%找每个点的相同面上的点的索引
    neib{i}(neib{i}==i) = [];
end
delete(p);

%%Roger best

%寻找邻接点，都比自己高的就是最小值

TR = triangulation(ori_faces,points);
face_neib = neighbors(TR);%所有面的相邻面
G=meshs2graph(TR);
neib=cell(1,point_num);
best=[];
for i=1:point_num
    temp2 = neighbors(G,i);
    neib{i}=temp2;
    temp3 = min(dist(temp2));%每个点相同面上的距离最小的点
    if dist(i)<=temp3
        best=[best i];
    end
end

%%
best = zeros(1,size(points,1));
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
f = figure('visible','off');
pcshow(points(best,:))
saveas(f, 'best_points.png')


%% Union the connected parts
[disjoint_set,size_DS]=union_zero_fast(best,neib);
protrusion_rep = unique(disjoint_set);
protrusion_rep(protrusion_rep==0)=[];

% for i=1:length(protrusion_rep)
%     construct_graph_relation

%%
% face_neib = construct_graph_relation(neib,ori_faces); %% face_num * 3
% TR = triangulation(ori_faces,points);
% face_neib = neighbors(TR);
%以一个点为中心，扩展n轮面

%% Roger 模范
%合并连在一起的best点
[disjoint_set,size_DS]=union_zero_fast(best,neib);

blebSegment=zeros(size(ori_faces,1),1);
points_color=unique(disjoint_set);%总共有x个点需要上色
points_color(1)=[];
face_extend_index=cell(1,length(points_color));

for i=1:size(ori_faces,1)
    depth(i)=sum(dist(ori_faces(i,:)))/3;
end

for i=1:length(points_color)%0不要
    face_extend_index{i}=extends_protrusion_point(TR,face_neib,points_color(i),10,depth);
end

cancel_j=[];
for i=1:length(points_color)-1
    for j=i+1:length(points_color)
        if ~isempty(intersect(face_extend_index{i}, face_extend_index{j}))
            %             face_extend_index{i}=union(face_extend_index{i},face_extend_index{j});
            face_extend_index{j}=[];
            cancel_j=[cancel_j j];
            %             points_color(j)=points_color(i);
        end
    end
end
points_color(cancel_j)=0;

for i=1:length(points_color)%0不要
    blebSegment(face_extend_index{i}')=points_color(i);
end
