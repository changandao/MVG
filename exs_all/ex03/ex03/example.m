%% load & display model
close all;
clear all;
[V,F] = openOFF('model.off');

%% close figure
T = [0.5,0.2,0.1]';
[V1, V2] = rotation(45/180*pi,0,120/180*pi,T,V);

%% display model again (possibly with changed vertex positions V)
figure,
P = patch('Vertices', V1, 'Faces', F, 'FaceVertexCData',0.3*ones(size(V,1),3));
axis equal;
shading interp;
camlight right;
camlight left;

figure,
P = patch('Vertices', V2, 'Faces', F, 'FaceVertexCData',0.3*ones(size(V,1),3));
axis equal;
shading interp;
camlight right;
camlight left;