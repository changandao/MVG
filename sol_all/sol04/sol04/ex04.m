%% Multiple View Geometry 2018, Exercise Sheet 4
% Prof. Daniel Cremers, Nikolaus Demmel, Marvin Eisenberger

%% Exercise 1: Image Formation

% Load data
[V, F] = openOFF('model.off','');
N_vert = size(V,1);

% (a)

% Transformation for shift
T = [eye(3) [-0.5 -0.5 1]'];

% convert to homogenious coordinates (and transpose for easier vector-wise
% multiplication)
V_hom = [V'; ones(1, N_vert)];

% shift and do perspective projection (cut 4th coordinate)
V_shifted = T * V_hom;

% (b)

% Kamera intrinsics
K = [540 0 320; ...
     0 540 240; ...
     0   0   1];

% apply camera matrix and convert from homogenious to pixel coordinates
V_proj_hom = K * V_shifted;
V_proj_uv = V_proj_hom(1:2,:) ./ V_proj_hom(3, :);

% visualize (patch also works with 2d coordinates)
figure(1)
subplot(121)
patch('Vertices', V_proj_uv', 'Faces', F)
% NOTE: It is only a minor detail, but as discussed in the tutorial, 
% to get an image that corresponds to our convention of the center of the 
% top-left pixel having coordinate (0,0), we set the viewport (axis) to a 
% range that starts at -0.5.
axis equal, axis([0 640 0 480]-0.5)
title('Perspective projection')

% (c)

% NOTE: In the case of parallel projection K(1,1) and K(2,2) are not
% interpreted as a focal length * scale factor, but instead just as the 
% scale factors that map camera coordinates to pixel coordinates.

% do parallel projection (set z to 1), apply camera intrinsics and convert
% to pixel coordinates
V_proj_parallel = [V_shifted(1:2,:); ones(1, N_vert)];
V_proj_parallel_hom = K * V_proj_parallel;
V_proj_parallel_uv = V_proj_parallel_hom(1:2,:) ./ V_proj_parallel_hom(3, :);


% visualize result
subplot(122)
patch('Vertices', V_proj_parallel_uv', 'Faces', F)
axis equal, axis([0 640 0 480]-0.5)
title('Parallel projection')


%% Exercise 2: Radial Distortion and Image Rectification

% (a)

% NOTE: The given camera intrinsics assume pixel coordinates that start at
% (0,0) in the top-left corner. If you want to use indexing starting at
% 1, you need to add +1 to each of the center-offsets o_x and o_y.

% load image and given camera intrinsics (camera matrix and FOV distortion
% parameter) and visualize image

Id1 = imreadbw('img1.jpg');

Kd1 = [388.6 0 343.7; ...
       0 389.4 234.6; ...
       0     0     1];
w1 = 0.92646;

figure(2)
subplot(121)
% NOTE: imagesc(Id1) in MATLAB will plot the image in the range 
%       (0.5, width+0.5) x (0.5, height+0.5) such that the center of the
%       top-left pixel is at (1,1). For visualization this is not a
%       problem, but later, when we lookup/interpolate pixel values, we
%       need to account for this convetion vs (0,0) for the top-left pixel.
imagesc(Id1), axis image, colormap gray
title('Distorted image')



% (b)

% distortion function for FOV/ATAN model with parameter w1; note that this
% can be applied to a whole vector of r-values.
g_ATAN_1 = @(r) (1./(w1*r) .* atan(2*tan(w1/2)*r));

% desired camera intrinsics for rectified image
K_new = [250 0 512; ...
         0 250 384; ...
         0   0   1];

% tic / toc for checking the runtime of rectification
tic
     
% meshgrid, starting with 0,0 in the top-left. This generates
% pixel-coordinates for all pixels in the image we want to create. For each
% of then we need to later lookup an (interpolated) intensity value in the
% distorted image. We then put all pixels as homogenous coordinates in one 
% long vector uv_hom for easier manipulation.
[u,v] = meshgrid(0:1023,0:767);
N_img = 1024 * 768;
uv_hom = [u(:) v(:) ones(N_img, 1)]';

% unproject image coordinates of ideal pinhole camera to generic image 
% plane (at Z=1). NOTE: inv(K_new) * uv_hom also works.
X_generic = K_new\uv_hom;

% compute the norm of the undistorted image coordinates
r = sqrt(X_generic(1,:).^2 + X_generic(2,:).^2);

% apply distortion, we can ignore the z coordinate of X_generic, since we
% know it is 1 for all points
X_d1 = [g_ATAN_1(r) .* X_generic(1:2,:); ones(1, N_img)];

% project the distorted coordinates to the actual image
uv_d1_hom = Kd1 * X_d1;

% Now find the pixel values for each point in uv_d1_hom by linear
% interpolation. Again we need to take care to ensure that the top-left
% corner has coordinates (0, 0). Also, we ignore the z coordinate of the
% homogenious vectors in uv_d1_hom, since we know they will be 1. Pixels
% outside the original image are set to black (0). Finally we reshape the
% vector of pixel values to a rectangular image.
[Hd1, Wd1] = size(Id1);
[grid_u_d1, grid_v_d1] = meshgrid(0:Wd1-1, 0:Hd1-1);
Inew = interp2(grid_u_d1, grid_v_d1, Id1, uv_d1_hom(1,:), uv_d1_hom(2,:), 'linear', 0);
Inew = reshape(Inew, size(u));

% prints the time since last tic
toc

% visualize the rectified image and save to disc
subplot(122)
imagesc(Inew), axis image, colormap gray
title('Undistorted image')
imwrite(Inew,'img1_undist.jpg')

%%

%(c)

% load image and given camera intrinsics (camera matrix and distortion 
% function) and visualize image

Id2 = imreadbw('img2.jpg');

Kd2 = [279.7 0 347.3; ...
       0 279.7 235.0; ...
       0     0     1];

% polinomial distortion function (which works on vector of r values)
g_pol_2 = @(r) 1 - 0.3407*r + 0.057*r.^2 - 0.0046*r.^3 + 0.00014*r.^4;
   
figure(3)
subplot(121)
imagesc(Id2), axis image, colormap gray
title('Distorted image')

tic

% apply distortion, we can ignore the z coordinate of X_generic, since we
% know it is 1 for all points
X_d2 = [g_pol_2(r) .* X_generic(1:2,:); ones(1, N_img)];

% project the distorted coordinates to the actual image
uv_d2_hom = Kd2 * X_d2;

% Now find the pixel values for each point in uv_d1_hom by linear
% interpolation. Again we need to take care to ensure that the top-left
% corner has coordinates (0, 0). Also, we ignore the z coordinate of the
% homogenious vectors in uv_d1_hom, since we know they will be 1. Pixels
% outside the original image are set to black (0). Finally we reshape the
% vector of pixel values to a rectangular image.
[Hd2, Wd2] = size(Id2);
[grid_u_d2, grid_v_d2] = meshgrid(0:Wd2-1, 0:Hd2-1);
Inew2 = interp2(grid_u_d2, grid_v_d2, Id2, uv_d2_hom(1,:), uv_d2_hom(2,:), 'linear', 0);
Inew2 = reshape(Inew2, size(u));

toc

% visualize the rectified image and save to disc
subplot(122)
imagesc(Inew2), axis image, colormap gray
title('Undistorted image')
imwrite(Inew2,'img2_undist.jpg')


%%

% (f)

% Inverse FOV distortion function for parameter w1
f_ATAN_1 = @(r) tan(w1*r) ./ (2*tan(w1/2)*r);

tic
     
% meshgrid as in (b), but now the dimensions are as for img1, or target
% camera
[H_img1, W_img1] = size(Id1);
[u,v] = meshgrid(0:W_img1-1,0:H_img1-1);
N_img1 = H_img1 * W_img1;
uv_hom = [u(:) v(:) ones(N_img1, 1)]';

% unproject image coordinates of target camera to distorted coordinates in
% the generic (Z=1) image plane
X_generic_d = Kd1\uv_hom;

% compute the norm of the distorted image coordinates
r_d = sqrt(X_generic_d(1,:).^2 + X_generic_d(2,:).^2);

% undistort the points using the inverse distorion function; z remains 1
X_generic = [f_ATAN_1(r_d) .* X_generic_d(1:2,:); ones(1, N_img1)];

% compute the norm of the undistorted image coordinates
r = f_ATAN_1(r_d) .* r_d;

% apply distortion of the second camera; as before z remains 1
X_d2 = [g_pol_2(r) .* X_generic(1:2,:); ones(1, N_img1)];

% project the distorted coordinates to image coordinates in the second
% camera
uv_d2_hom = Kd2 * X_d2;

% Now linear interpolation as before
Inew3 = interp2(grid_u_d2, grid_v_d2, Id2, uv_d2_hom(1,:), uv_d2_hom(2,:), 'linear', 0);
Inew3 = reshape(Inew3, size(u));

% prints the time since last tic
toc

% visualize the rectified image and save to disc
figure(4)
subplot(121)
imagesc(Id2), axis image, colormap gray
title('Original distorted image from cam2')
subplot(122)
imagesc(Inew3), axis image, colormap gray
title('New virtual image as if taking with cam1')
imwrite(Inew,'img2_cam1.jpg')
