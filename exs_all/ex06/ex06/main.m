clear 
close all;

img1 = imread('batinria0.tif');
img2 = imread('batinria1.tif');

points1 = load('points1.mat');
x1 = [points1.x1, points1.y1, ones(8,1)];

points2 = load('points2.mat');
x2 = [points2.x2, points2.y2, ones(8,1)];

recmat = zeros(8,9);


