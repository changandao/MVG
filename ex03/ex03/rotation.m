function [V1,V2] = rotation(alpha, belta, gamma ,translate, V)
rotatex = [1,0,0;0,cos(alpha),-sin(alpha);0, sin(alpha), cos(alpha)];
rotatey = [cos(belta), 0, sin(belta);
            0, 1, 0;
            -sin(belta),0, cos(belta)];
rotatez = [cos(gamma), -sin(gamma), 0;
    sin(gamma), cos(gamma), 0;
    0, 0, 1];

rotate1 = rotatex*rotatey*rotatez;
rotate2 = rotatez* rotatey* rotatex;
g1 = [rotate1, translate;
    zeros(1,3),1];
g2 = [rotate2, translate;
    zeros(1,3),1];

V_mean = mean(V);
V = V - ones(size(V,1),1) * V_mean;

V = [V,ones(size(V,1),1)];
V1 = (g1 * V')';
V2 = (g2 * V')';

V1 = V1(:, 1:3)./V1(:, 4);
V2 = V2(:, 1:3)./V2(:, 4);
end