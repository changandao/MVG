function w = twist(R)
norm_w = acos((trace(R)-1)/2),

w = (norm_w / (2* sin(norm_w)))*[R(3,2)-R(3,2), R(1,3)- R(3,1), R(2,1)- R(1,2)];

end