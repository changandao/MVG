clear all;
close all;

I = imreadbw('img1.png');

[score, points] = getHarrisCorners(I, 1, 0.05, 10e-7);
figure,
drawPts(I, points);