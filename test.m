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

load('surface_18.mat');
% load('torus.mat');

points = surface.vertices;
faces = surface.faces;
faces = sort(faces,2);

edges = get_uni_edges(faces);

X=points(:,1);
Y=points(:,2);
Z=points(:,3);
[con_faces,~] = convhull(X,Y,Z); % con_faces is the index of the points in the convex hull

con_edges = get_uni_edges(con_faces);
con_points = unique(con_faces(:));
cell_only_points = setdiff(1:size(points, 1), con_points);

vertexToEdges = construct_vertex_to_edges(edges);
% cell for neighbors (unique([vertexToEdges{edges(i, 1)}; vertexToEdges{edges(i, 2)}]);)
edge_neighbors = cellfun(@(x) unique([vertexToEdges{edges(x, 1)}; vertexToEdges{edges(x, 2)}]), num2cell(1:size(edges, 1)), 'UniformOutput', false);

% color the edges of the cell
color_cell_edges = color_edges(edges, con_edges, edge_neighbors);
color_cell_faces = max(color_cell_edges(faces), [], 2);

% color the edges of the convex hull
color_con_edges = color_edges(con_edges, edges, edge_neighbors);
color_con_faces = max(color_con_edges(con_faces), [], 2);

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
    % initialize the stack
    stack = [edge];

    % iterate over the stack
    while ~isempty(stack)
        % pop the last element 
        currentEdge = stack(end);
        stack(end) = [];

        % if the edge is already colored, continue
        if color_cell_edges(currentEdge) ~= -1
            continue;
        end

        % color the edge 
        color_cell_edges(currentEdge) = color;

       % get the neighbors of the edge
        neighbors = edge_neighbors{currentEdge};

        % push the neighbors to the stack, if they are not colored
        for i = 1:length(neighbors)
            neighborEdge = neighbors(i);
            % 
            if color_cell_edges(neighborEdge) == -1
                stack(end + 1) = neighborEdge;
            end
        end
    end
end

% function color_points = color_points(faces_this, faces_other, 
%     color_points = zeros(size(faces, 1), 1);
    