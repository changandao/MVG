function R = retwist(w)

wbar = trans(w);
R = eye(3) + wbar/norm(w)*sin(norm(w)) + (w'*w/(norm(w)*norm(w))-eye(3))*(1-cos(norm(w))); 

end

%%test


