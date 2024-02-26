% 大立方体的顶点
verticesLarge = [-1 -1 -1; -1 -1 1; -1 1 1; -1 1 -1; 1 -1 -1; 1 -1 1; 1 1 1; 1 1 -1];

% 小立方体的顶点
verticesSmall = [-0.5 -0.5 -0.5; -0.5 -0.5 0.5; -0.5 0.5 0.5; -0.5 0.5 -0.5; 0.5 -0.5 -0.5; 0.5 -0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 -0.5];

vertices = [verticesLarge; verticesSmall];

% add little noise to vertices
vertices = vertices + 0.01 * randn(size(vertices));

small_faces = [9 10 11 12; 9 10 14 13; 10 11 15 14; 11 12 16 15; 12 9 13 16; 13 14 15 16];

big_faces = [1 9 10 2; 2 10 11 3; 3 11 12 4; 4 12 9 1; 1 2 3 4; 5 13 14 6; 6 14 15 7; 7 15 16 8; 8 16 13 5; 5 6 7 8; 1 5 13 9; 2 6 14 10; 3 7 15 11; 4 8 16 12];

faces = [big_faces; small_faces];

% figure;

% % 绘制立方体
% patch('Vertices',vertices,'Faces',faces,'FaceColor','r', 'FaceAlpha', 0.1);

% axis equal; % 保证各轴等比缩放
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% view(3); % 设置3D视图

% % 在每个顶点位置添加顶点编号
% for i = 1:size(vertices, 1)
%     text(vertices(i,1), vertices(i,2), vertices(i,3), num2str(i), ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','bottom');
% end

surface.vertices = vertices;
surface.faces = faces;

save('surface.mat', 'surface'); % 保存立方体的顶点和面信息
