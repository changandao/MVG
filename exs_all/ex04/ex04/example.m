%% load & display model
close all;
clear all;
[V,F] = openOFF('model.off');

%% close figure


%% display model again (possibly with changed vertex positions V)
figure,
P = patch('Vertices', V, 'Faces', F, 'FaceVertexCData',0.3*ones(size(V,1),3));
axis equal;
shading interp;
camlight right;
camlight left;

%% projection

K = [540,0,320;
    0,540,240;
    0,0,1];
T = [-0.5,-0.5,1]';
V1 = transformation(0,0,0,T,V);
V_project = (K * V1')';
V_uv = V_project(:,1:2)./V_project(:,3);

figure,
P = patch('Vertices', V_uv, 'Faces', F, 'FaceVertexCData',0.3*ones(size(V,1),3));
axis equal, axis([0 640 0 480]-0.5);


%% ortogonal projection
V_para_project = (K * [V1(:,1:2), ones(size(V,1),1)]')';
V_para_uv = V_para_project(:,1:2)./V_para_project(:,3);
figure,
P = patch('Vertices', V_para_uv, 'Faces', F, 'FaceVertexCData',0.3*ones(size(V,1),3));
axis equal, axis([0 640 0 480]-0.5);

%%distorted image
Kd1 = [388.6, 0, 343.7; 
    0, 389.4, 234.6;
    0, 0, 1];
w = 0.92646;
I1 = imreadbw('img1.jpg');
imagesc(I1), axis image, colormap gray;

%% b 
[u,v]= meshgrid(0:1023,0:767);
Knew = [250,0,512;
    0,250,384;
    0,0,1];
N_img = 1024*768;
uv_hom = [u(:), v(:), ones(N_img,1)];

X_hom = Knew\uv_hom';

r = sqrt(X_hom(1,:).^2 + X_hom(2,:).^2);
project_new1 = (1./(w*r)) .* atan(2*r*tan(w/2)) ;

Xd1_hom = [project_new1 .* X_hom(1:2,:); X_hom(3,:)];

uvd1_hom = Kd1 * Xd1_hom; 

[md, nd] = size(I1);
[xd, yd] = meshgrid(0:nd-1, 0:md-1);
Inew = interp2(xd,yd,I1, uvd1_hom(1,:), uvd1_hom(2,:), 'linear', 0);
Inew = reshape(Inew, size(u));

figure,
imagesc(Inew), axis image, colormap gray
title('Undistorted image')
imwrite(Inew,'img1_undist.jpg')

%% c

Kd2 = [279.7, 0, 347.3; 
    0, 279.7, 235.0;
    0, 0, 1];



project_new2 = 1 - 0.3407 * r + 0.057*r.^2 - 0.0046*r.^3 + 0.00014*r.^4 ;
Xd2_hom = [project_new2 .* X_hom(1:2,:); X_hom(3,:)];
uvd2_hom = Kd2 * Xd2_hom; 

Inew = interp2(xd,yd,I1, uvd2_hom(1,:), uvd2_hom(2,:), 'linear', 0);
Inew = reshape(Inew, size(u));

figure,
imagesc(Inew), axis image, colormap gray
title('Undistorted image')
imwrite(Inew,'img2_undist.jpg')

