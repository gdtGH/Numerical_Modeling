function [x, u, u_ex] = PointwiseSol(femregion, uh, u_ex, ~)
    x  = femregion.coord(:,1);
    u  = full(uh);
end