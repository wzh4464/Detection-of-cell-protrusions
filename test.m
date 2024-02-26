%% init and set path
clear
close all
addpath(genpath(pwd))
threshold = 0.05;

% if parrallel pool is not open, open a new one
% if isempty(gcp('nocreate'))
%     num_cores = feature('numCores');
%     parpool(num_cores);
% end

%% mat from nature

% surface is a mesh struct
% surface is a struct with fields: vertices, faces
% vertices is a n*3 matrix, n is the number of points
% each row is a point, the first column is x, the second column is y, the third column is z
% faces is a m*3 matrix, m is the number of faces
% each row is a face, each column is the index of the point in vertices

% load('surface_18.mat');
% load('torus.mat');

% load from self_defined_mesh/surface.mat
load('self_defined_mesh/surface.mat');

points = surface.vertices;
faces = surface.faces;
faces = sort(faces,2); % face: f * 3, f is the number of faces, each row is the index of the points in the face
edges = get_uni_edges(faces); % edges: e * 2, e is the number of edges, each row is the index of the points in the edge

%% convex hull
[con_faces,~] = convhull(points);
con_edges = get_uni_edges(con_faces);
con_points_ind = unique(con_faces(:));
cell_only_points = setdiff(1:size(points, 1), con_points_ind); % points that are not on the convex hull boundary but on the cell surface

save_mesh_figure(con_faces, points, 'convex_hull.png');
save_mesh_figure(faces, points, 'cell_surface.png');

%% more init
vertexToEdges = construct_vertex_to_edges(edges); % vertexToEdges{i} is a list of edges that are connected to vertex i
% cell for neighbors (unique([vertexToEdges{edges(i, 1)}; vertexToEdges{edges(i, 2)}]);)
edge_neighbors = cellfun(@(x) unique([vertexToEdges{edges(x, 1)}; vertexToEdges{edges(x, 2)}]), num2cell(1:size(edges, 1)), 'UniformOutput', false);

% color the edges of the cell
color_cell_edges = color_edges(edges, con_edges, edge_neighbors);
color_cell_faces = max(color_cell_edges(faces), [], 2);

% color the edges of the convex hull
color_con_edges = color_edges(con_edges, edges, edge_neighbors);

A = color_con_edges(con_faces(:,1));
B = color_con_edges(con_faces(:,2));
C = color_con_edges(con_faces(:,3));
color_con_faces = max([A,B,C],[],2);

% get the number of regions
num_cell_regions = max(color_cell_edges);
num_con_regions = max(color_con_edges);

% get the distance matrix, dist_matrix(i, j) is the number of different regions between i-th cell region and j-th convex hull region
dist_matrix = zeros(num_cell_regions, num_con_regions);

cell_face_center_set = zeros(num_cell_regions, 3);
con_face_center_set = zeros(num_con_regions, 3);

% get the center of each region
for i = 1:num_cell_regions
    cell_face_center_set(i, :) = mean(points(color_cell_faces == i, :));
end

for i = 1:num_con_regions
    con_face_center_set(i, :) = mean(points(color_con_faces == i, :));
end


for i = 1:num_cell_regions
    for j = 1:num_con_regions
        % get the distance between the centers of the two regions
        dist_matrix(i, j) = hausdorffDistance(points, cell_face_center_set(i, :), con_face_center_set(j, :));
    end
end

% find the nearest region for each cell region
[~, nearest_con_region] = min(dist_matrix, [], 2);

% here, we have the nearest convex hull region for each cell region

%% projection
% get the projection of the points on the cell surface to the convex hull boundary

for i = 1:length(cell_only_points)
    point = points(cell_only_points(i), :);
    % find the nearest face of the convex hull
    [~, nearest_face] = min(pdist2(point, points(con_faces, :)));
    nearest_edge = con_edges(any(con_faces == nearest_face, 2), :);
    nearest_edge = nearest_edge(1, :);
    nearest_point = points(nearest_edge, :);
    % project the point to the line
    projected_point = project_point_to_line(point, nearest_point(1, :), nearest_point(2, :));
    points(cell_only_points(i), :) = projected_point;
end



%% functions

function uni_edges = get_uni_edges(faces)
    edges1 = faces(:, [1, 2]);
    edges2 = faces(:, [2, 3]);
    edges3 = faces(:, [3, 1]);
    edges = [edges1; edges2; edges3];
    edges = sort(edges, 2);
    [uni_edges, ~] = unique(edges, 'rows');
end

% function euler_char = get_euler_char(points, faces)
%     uniqueEdges = get_uni_edges(faces);
%     euler_char = size(points, 1) - size(uniqueEdges, 1) + size(faces, 1);
% end

% Step 1: Initialize a mapping from vertices to edges
function vertexToEdges = construct_vertex_to_edges(edges)
    % vertexToEdges{i} is a list of edges that are connected to vertex i
    % edges(vertexToEdges{i}, :) are two vertices that are connected to vertex i
    vertexToEdges = cell(max(edges(:)), 1);

    % profile this
    
    begin_time = tic;

    for i = 1:size(edges, 1)
        vertexToEdges{edges(i, 1)} = [vertexToEdges{edges(i, 1)}; i];
        vertexToEdges{edges(i, 2)} = [vertexToEdges{edges(i, 2)}; i];
    end

    end_time = toc(begin_time);
    fprintf('Time to construct vertex to edges: %f\n', end_time);
end

function color_cell_edges = color_edges(edges, con_edges, edge_neighbors)

    % initialize the color of the edges with -1
    color_cell_edges = zeros(size(edges, 1), 1);
    color_cell_edges = color_cell_edges - 1;
    
    % find the boundary edges
    for i = 1:size(edges, 1)
        if ismember(i, con_edges)
            color_cell_edges(i) = 0;
        end
    end
    
    % for non-boundary edges, find neighbor edges and make them the same color
    currentColor = 1;

    % Iterate over all edges to color non-boundary edges
    for i = 1:size(edges, 1)
        if color_cell_edges(i) == -1
            % find edges from the same region, which means:
            % 1. the edge is not a boundary edge
            % 2. the edge is a neighbor of i
            color_cell_edges = color_the_region(i, currentColor, color_cell_edges, edge_neighbors);
            currentColor = currentColor + 1;

        end
    end

end

function color_cell_edges = color_the_region(edge, color, color_cell_edges, edge_neighbors)
    % 预估可能需要的栈空间大小，这里假设为边的总数
    maxEdges = numel(color_cell_edges);
    stack = zeros(1, maxEdges); % 预分配栈空间
    stackTop = 1;
    stack(stackTop) = edge;

    while stackTop > 0
        currentEdge = stack(stackTop);
        stackTop = stackTop - 1;

        if color_cell_edges(currentEdge) ~= -1
            continue;
        end

        color_cell_edges(currentEdge) = color;

        neighbors = edge_neighbors{currentEdge};
        % 逻辑索引找出所有未上色的邻边
        uncoloredNeighbors = neighbors(color_cell_edges(neighbors) == -1);

        % 批量将未上色的邻边加入栈
        numUncolored = numel(uncoloredNeighbors);
        stack(stackTop + 1:stackTop + numUncolored) = uncoloredNeighbors;
        stackTop = stackTop + numUncolored;
    end
end


% function color_points = color_points(faces_this, faces_other, 
%     color_points = zeros(size(faces, 1), 1);

function save_mesh_figure(faces, vertices, filename)
    f = figure;
    patch('Faces',faces,'Vertices',vertices,'FaceColor','red','EdgeColor','black', 'FaceAlpha', 0.1);
    axis equal;
    view(3);
    saveas(f, filename);
end