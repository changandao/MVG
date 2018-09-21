function [score, points] = getHarrisCorners(I, sigma, kappa, theta)
[M11, M12, M22] = getM(I, sigma);

% TODO: compute score using det and trace
score = M11.*M22 - M12.*M12 - kappa * (M11 + M22).^2;

% TODO: display score
figure,
imagesc(sign(score));
colormap gray;
% TODO: find corners (variable points)
[row,colum] = find((score - theta)>0);
points = [colum, row];

end