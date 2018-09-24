clear all
I1 = imreadbw('img2.jpg');
imagesc(I1), axis image, colormap gray;

%% b 
w = 0.92646;
[u,v]= meshgrid(0:1023,0:767);
Knew = [250,0,512;
    0,250,384;
    0,0,1];
N_img = 1024*768;
uv_hom = [u(:), v(:), ones(N_img,1)];

X_hom = Knew\uv_hom';

r = sqrt(X_hom(1,:).^2 + X_hom(2,:).^2);
%project_new1 = (1./(w*r)) .* atan(2*r*tan(w/2));
project_new1 = 1 - 0.3407 * r + 0.057*r.^2 - 0.0046*r.^3 + 0.00014*r.^4 ;

Xd1_hom = [project_new1 .* X_hom(1:2,:); X_hom(3,:)];

Kd1 = [279.7, 0, 347.3; 
    0, 279.7, 235.0;
    0, 0, 1]

uvd1_hom = Kd1 * Xd1_hom; 

[md, nd] = size(I1);
[xd, yd] = meshgrid(0:nd-1, 0:md-1);
Inew = interp2(xd,yd,I1, uvd1_hom(1,:), uvd1_hom(2,:), 'linear', 0);
Inew = reshape(Inew, size(u));

figure,
imagesc(Inew), axis image, colormap gray
title('Undistorted image')
imwrite(Inew,'img1_undist.jpg')