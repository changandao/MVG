function [M11, M12, M22] = getM(I, sigma)
differ = [1,0,-1];

sobely = differ' * ones(1,3);
sobelx = sobely';

Ix = conv2(I, sobelx, 'same');
Iy = conv2(I, sobely, 'same');

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

end