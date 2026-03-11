function [x, u, u_ex] = Snapshot(femregion, uh, ~, ~, u_ex)
    if nargin < 5, u_ex = []; end
    x  = femregion.coord(:,1);
    u  = full(uh);
end