function [M11, M12, M22, q1, q2] = getMq(I1, I2, sigma)
differ = [1,0,-1];

sobely = differ' * ones(1,3);
sobelx = sobely';

Ix = conv2(I1, sobelx, 'same');
Iy = conv2(I1, sobely, 'same');

k = 4*sigma+1;
Gaus_kernel = zeros(k,k);
for i = 1:2*k+1
    for j = 1:2*k+1
    Gaus_kernel(i,j) = exp(-((i+j)^2)/(2*sigma));
    end
end
Gaus_sum = sum(sum(Gaus_kernel));
Gaus_kernel = Gaus_kernel/Gaus_sum;

M11 = conv2(Ix.*Ix, Gaus_kernel, 'same');
M22 = conv2(Iy.*Iy, Gaus_kernel, 'same');
M12 = conv2(Iy.*Ix, Gaus_kernel, 'same');

%% q
It = I2 - I1;

q1 = conv2(Ix.*It, Gaus_kernel, 'same');
q2 = conv2(Iy.*It, Gaus_kernel, 'same');


end