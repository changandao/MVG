clear all;
close all;

I1 = imreadbw('img1.png');
I2 = imreadbw('img2.png');
[M,N] = size(I1);

[vx,vy] = getFlow(I1, I2, 2);

figure,
% subplot(3,1,1)
% imagesc(vx);
% subplot(3,1,2)
% imagesc(vy);
% subplot(3,1,3)
hold off
imagesc(I2);
colormap gray
hold on
quiver(vx,vy,0);
