% Parameters
R = 5; % Major radius
r = 2; % Minor radius
nTheta = 36; % Number of divisions in theta direction
nPhi = 18; % Number of divisions in phi direction

% Mesh grid for theta and phi
[theta, phi] = meshgrid(linspace(0, 2*pi, nTheta), linspace(0, 2*pi, nPhi));

% Parametric equations for the torus
x = (R + r*cos(phi)).*cos(theta);
y = (R + r*cos(phi)).*sin(theta);
z = r*sin(phi);

% Plotting the torus
% figure;
% surf(x, y, z, 'FaceColor', 'cyan', 'EdgeColor', 'none');
% camlight left; lighting phong
% axis equal
% title('Torus Mesh');

% Flatten the matrices to create a vertex list
vertices = [x(:), y(:), z(:)];

% Initialize the faces matrix
faces = [];

% Number of vertices in each direction
nVerticesTheta = nTheta;
nVerticesPhi = nPhi;

% Create faces
for i = 1:nVerticesTheta-1
    for j = 1:nVerticesPhi-1
        % Calculate vertex indices for the corners of the square
        v1 = (i-1)*nVerticesPhi + j;
        v2 = (i-1)*nVerticesPhi + j + 1;
        v3 = i*nVerticesPhi + j + 1;
        v4 = i*nVerticesPhi + j;
        
        % Define two triangles for the current square
        faces = [faces; 
                 v1, v2, v3; % Triangle 1
                 v1, v3, v4]; % Triangle 2
    end
end

% Adjust for periodic boundary conditions
for i = 1:nVerticesTheta-1
    faces = [faces;
             i*nVerticesPhi, (i-1)*nVerticesPhi + 1, i*nVerticesPhi + 1;
             i*nVerticesPhi, (i+1)*nVerticesPhi, (i-1)*nVerticesPhi + 1];
end
for j = 1:nVerticesPhi-1
    faces = [faces;
             j, (nVerticesTheta-1)*nVerticesPhi + j + 1, j+1;
             j, (nVerticesTheta-1)*nVerticesPhi + j, (nVerticesTheta-1)*nVerticesPhi + j + 1];
end
faces = [faces;
         nVerticesPhi, nVerticesTheta*nVerticesPhi, 1;
         nVerticesPhi, (nVerticesTheta-1)*nVerticesPhi, nVerticesTheta*nVerticesPhi];

% Plot the mesh with face indexing
% figure;
% trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'black');
% axis equal
% title('Torus Mesh with Face Indexing');

% Save the mesh to a file

% create a structure named surface
% surface.vertices is a matrix of size n x 3
% surface.faces is a matrix of size m x 3

surface.vertices = vertices;
surface.faces = faces;

save('torus.mat', 'surface');
