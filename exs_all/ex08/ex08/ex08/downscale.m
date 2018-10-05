function [ Id, Dd, Kd ] = downscale( I, D, K, level )
    % downscale intensity image, depth map and camera intrinsics (recursively)

    if(level <= 1)
        % coarsest pyramid level
        Id = I;
        Dd = D;
        Kd = K;
        return;
    end

    % downscale camera intrinsics
    % TODO
    %Kd = ...
    Kd = [K(1,1)/2, 0 ,K(1,3)/2 - 1/4; 0, K(2,1)/2, K(2,3)/2 - 1/4; 0,0,1];

    % downscale intensity image
    % TODO
    %Id = ...
    Id = (I(0+1:2:end,0+1:2:end)+...
        I(1+1:2:end,0+1:2:end)+...
        I(0+1:2:end,1+1:2:end)+...
        I(1+1:2:end,1+1:2:end))/4;
    

    % downscale depth map
    % TODO
    %Dd = ...
    Dvalid = (sign(D(0+1:2:end,0+1:2:end))+...
        sign(D(1+1:2:end,0+1:2:end))+...
        sign(D(0+1:2:end,1+1:2:end))+...
        sign(D(1+1:2:end,1+1:2:end)))/4;
    
    Dd = (D(0+1:2:end,0+1:2:end)+...
        D(1+1:2:end,0+1:2:end)+...
        D(0+1:2:end,1+1:2:end)+...
        D(1+1:2:end,1+1:2:end))./Dvalid; 
    Dd(isnan(Dd)) = 0;
    % recursively downscale on next pyramid level
    [Id, Dd, Kd] = downscale( Id, Dd, Kd, level - 1);    
end