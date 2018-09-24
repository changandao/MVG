function V = transformation(alpha, belta, gamma ,translate, V)
rotatex = [1,0,0;0,cos(alpha),-sin(alpha);0, sin(alpha), cos(alpha)];
rotatey = [cos(belta), 0, sin(belta);
            0, 1, 0;
            -sin(belta),0, cos(belta)];
rotatez = [cos(gamma), -sin(gamma), 0;
    sin(gamma), cos(gamma), 0;
    0, 0, 1];

rotate = rotatex * rotatey * rotatez;

V_mean = mean(V);

g_submean = [eye(3), -V_mean';
    zeros(1,3),1];
g1 = [rotate, translate;
    zeros(1,3),1];
g_addmean = [eye(3), V_mean';
    zeros(1,3),1];

V = [V,ones(size(V,1),1)];
g = g_addmean * g1 * g_submean;
V = (g*V')';

V = V(:, 1:3)./V(:, 4);
end